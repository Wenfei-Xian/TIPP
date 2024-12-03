#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <cmath>
#include <mutex>
#include <algorithm>
#include <set>
#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"

//v1.1 May 29, 2024
//v1.2 May 31, 2024 ouput the high frequency IDs
//v1.3 Dec 03, 2024 accept the parameters for lkc and hkc

std::mutex median_mtx;
std::mutex validSeqIDs_mtx;
std::mutex tmpFileMtx;
std::vector<double> medians;
std::set<std::string> validSeqIDs;
std::ofstream tmpFile;
std::ofstream highPercentFile;

CKMCFile kmcDb;

std::string reverseComplement(const std::string& kmer) {
    std::string rc(kmer.size(), ' ');
    for (size_t i = 0; i < kmer.size(); ++i) {
        switch (kmer[i]) {
            case 'A': rc[kmer.size() - 1 - i] = 'T'; break;
            case 'T': rc[kmer.size() - 1 - i] = 'A'; break;
            case 'G': rc[kmer.size() - 1 - i] = 'C'; break;
            case 'C': rc[kmer.size() - 1 - i] = 'G'; break;
        }
    }
    return rc;
}

void processRead(const std::string& read, const std::string& seqID) {
    uint32_t kmerLength = kmcDb.KmerLength();
    CKmerAPI kmer(kmerLength);

    std::vector<uint32_t> counts;
    for (size_t i = 0; i <= read.length() - kmerLength; ++i) {
        std::string kmerStr = read.substr(i, kmerLength);
        kmer.from_string(kmerStr.c_str());
        uint32_t count;
        if (kmcDb.CheckKmer(kmer, count)) { //canonical
            counts.push_back(count);
        }
	else {
	    std::string rcStr = reverseComplement(kmerStr);
	    kmer.from_string(rcStr.c_str());
	    if (kmcDb.CheckKmer(kmer, count)) {
		    counts.push_back(count);
            }
        }
    }

    double median = 0.0;
    if (!counts.empty()) {
        size_t n = counts.size() / 2;
        nth_element(counts.begin(), counts.begin() + n, counts.end());
        median = counts[n];
        if (counts.size() % 2 == 0) {
            nth_element(counts.begin(), counts.begin() + n - 1, counts.end());
            median = (median + counts[n - 1]) / 2.0;
        }
    }

    std::lock_guard<std::mutex> median_lock(median_mtx);
    medians.push_back(median);

    std::lock_guard<std::mutex> tmpFile_lock(tmpFileMtx);
    tmpFile << seqID;
    for (auto count : counts) {
        tmpFile << " " << count;
    }
    tmpFile << std::endl;
}

std::string getFileExtension(const std::string& filename) {
    size_t dotPos = filename.find_last_of('.');
    if (dotPos != std::string::npos)
        return filename.substr(dotPos);
    return "";
}

void analyzeTmpFile (const std::string& tmpFilename, double fileMedian, float lkc, float hkc ) {
    std::ifstream tmpFile(tmpFilename);
    if (!tmpFile.is_open()) {
        std::cerr << "Could not open temporary file for reading: " << tmpFilename << std::endl;
        return;
    }

    std::string highPercentFilename = tmpFilename + ".highfrq.ID";
    highPercentFile.open(highPercentFilename);
    if (!highPercentFile.is_open()) {
	    std::cerr << "Could not open file for writing high percent IDs: " << highPercentFilename << std::endl;
	    return;
    }

    std::string line;
    while (std::getline(tmpFile, line)) {
        std::istringstream iss(line);
        std::string seqID;
        uint32_t count;
        int totalCount = 0;
        int lowCount = 0;
	int highCount = 0;

        iss >> seqID;
        while (iss >> count) {
            totalCount++;
            if (count < fileMedian * lkc ) { // organelle genome 3 times larger than nuclear
                lowCount++;
            }
	    else if( count > fileMedian * hkc ) { // high repetative region
		highCount++;
	    }
        }

        double lowPercent = 100.0 * lowCount / totalCount;
	double highPercent = 100.0 * highCount/ totalCount;

        if ( lowPercent < 20.0 && highPercent < 30.0 ) { // 1, remove low quality reads. 2, remove low frequency reads (origin from nuclear)
            std::lock_guard<std::mutex> validSeqIDs_lock(validSeqIDs_mtx);
            validSeqIDs.insert(seqID);
        }
	else if( highPercent >= 30.0 ){
	    highPercentFile << seqID << std::endl;
	}
    }

    tmpFile.close();
    highPercentFile.close();
}

void outputValidSequences(const std::string& inputFile, const std::string& outputFilename) {
    std::ifstream file(inputFile);
    std::ofstream outputFile(outputFilename);
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << inputFile << std::endl;
        return;
    }
    if (!outputFile.is_open()) {
        std::cerr << "Could not open output file: " << outputFilename << std::endl;
        return;
    }

    std::string line, read, seqID;
    bool output = false;

    while (std::getline(file, line)) {
        if (line.front() == '>' || line.front() == '@') {
            if (output) {
                outputFile << read << std::endl;
            }
            seqID = line.substr(1);
            read.clear();
            output = false;

            std::lock_guard<std::mutex> validSeqIDs_lock(validSeqIDs_mtx);
            if (validSeqIDs.find(seqID) != validSeqIDs.end()) {
                outputFile << line << std::endl;
                output = true;
            }
        } else if (!line.empty() && output) {
            read += line;
        }
    }

    if (output) {
        outputFile << read << std::endl;
    }

    file.close();
    outputFile.close();
}

void countKmersInReadsMultiThreaded(const std::string& inputFile, int maxThreads) {
    std::ifstream file(inputFile);
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << inputFile << std::endl;
        return;
    }

    std::string extension = getFileExtension(inputFile);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

    std::string line, read, seqID;
    int lineCounter = 0;
    std::vector<std::thread> threads;

    while (std::getline(file, line)) {
        if (extension == ".fq" || extension == ".fastq") {
            lineCounter++;

            if (lineCounter == 1) {
                size_t firstSpace = line.find(' ');
                seqID = (firstSpace != std::string::npos) ? line.substr(1, firstSpace - 1) : line.substr(1);
            } else if (lineCounter == 2) {
                read = line;
            } else if (lineCounter == 4) {
                lineCounter = 0;
                threads.push_back(std::thread(processRead, read, seqID));
            }
        } else if (extension == ".fa" || extension == ".fasta") {
            if (line.front() == '>') {
                if (!read.empty()) {
                    threads.push_back(std::thread(processRead, read, seqID));
                    read.clear();
                }
                seqID = line.substr(1);
            } else {
                read += line;
            }
        }

        if (threads.size() >= maxThreads) {
            for (auto& t : threads) {
                if (t.joinable()) {
                    t.join();
                }
            }
            threads.clear();
        }
    }

    if (!read.empty() && (extension == ".fa" || extension == ".fasta")) {
        threads.push_back(std::thread(processRead, read, seqID));
    }

    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }

    file.close();
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <kmc_database> <input_file> <max_threads> <lkc> <hkc>" << std::endl;
        return 1;
    }

    const std::string kmcDatabase = argv[1];
    const std::string inputFile = argv[2];
    //int peak = std::stoi(argv[3]);
    int maxThreads = std::stoi(argv[3]);
    float lkc = std::stof(argv[4]);
    float hkc= std::stof(argv[5]);

    if (!kmcDb.OpenForRA(kmcDatabase)) {
        std::cerr << "Could not open KMC database: " << kmcDatabase << std::endl;
        return 1;
    }

    std::string tmpFilename = inputFile + ".kmercount.txt";
    tmpFile.open(tmpFilename);
    if (!tmpFile.is_open()) {
        std::cerr << "Could not open temporary file for writing: " << tmpFilename << std::endl;
        return 1;
    }

    countKmersInReadsMultiThreaded(inputFile, maxThreads);

    tmpFile.close();

    double fileMedian = 0.0;
    
    //if( peak == 0 ){
        if (!medians.empty()) {
            size_t m = medians.size() / 2;
            nth_element(medians.begin(), medians.begin() + m, medians.end());
            fileMedian = medians[m];
            if (medians.size() % 2 == 0) {
                nth_element(medians.begin(), medians.begin() + m - 1, medians.end());
                fileMedian = (fileMedian + medians[m - 1]) / 2.0;
            }
	}
    //}
    //else{
    //	    fileMedian = peak;
    //}

    //for (size_t i = 0; i < medians.size(); ++i) {
    //    std::cout << medians[i] << " ";
    //}
    //std::cout << std::endl;

    std::cout << "File Median Count: " << fileMedian << std::endl;

    analyzeTmpFile(tmpFilename, fileMedian, lkc, hkc);

    std::string outputFilename = inputFile + ".filter.fasta";
    outputValidSequences(inputFile, outputFilename);

    kmcDb.Close();

    return 0;
}

