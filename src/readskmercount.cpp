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

std::mutex median_mtx;
std::mutex validSeqIDs_mtx;
std::mutex tmpFileMtx;
std::vector<double> medians;
std::set<std::string> validSeqIDs;
std::ofstream tmpFile;

void processRead(const std::string& kmcDatabase, const std::string& read, const std::string& seqID) {
    CKMCFile kmcDb;
    if (!kmcDb.OpenForRA(kmcDatabase)) {
        std::cerr << "Could not open KMC database: " << kmcDatabase << std::endl;
        return;
    }

    uint32_t kmerLength = kmcDb.KmerLength();
    CKmerAPI kmer(kmerLength);

    std::vector<uint32_t> counts;
    for (size_t i = 0; i <= read.length() - kmerLength; ++i) {
        std::string kmerStr = read.substr(i, kmerLength);
        kmer.from_string(kmerStr.c_str());
        uint32_t count;
        if (kmcDb.CheckKmer(kmer, count)) {
            counts.push_back(count);
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

    median_mtx.lock();
    medians.push_back(median);
    median_mtx.unlock();

    tmpFileMtx.lock();
    tmpFile << seqID;
    for (auto count : counts) {
        tmpFile << " " << count;
    }
    tmpFile << std::endl;
    tmpFileMtx.unlock();

    kmcDb.Close();
}

std::string getFileExtension(const std::string& filename) {
    size_t dotPos = filename.find_last_of('.');
    if (dotPos != std::string::npos)
        return filename.substr(dotPos);
    return "";
}

void analyzeTmpFile(const std::string& tmpFilename, double fileMedian) {
    std::ifstream tmpFile(tmpFilename);
    if (!tmpFile.is_open()) {
        std::cerr << "Could not open temporary file for reading: " << tmpFilename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(tmpFile, line)) {
        std::istringstream iss(line);
        std::string seqID;
        uint32_t count;
        int totalCount = 0;
        int lowCount = 0;

        iss >> seqID;
        while (iss >> count) {
            totalCount++;
            if (count < fileMedian * 0.5) {
                lowCount++;
            }
        }

        double lowPercent = 100.0 * lowCount / totalCount;
        //if (lowPercent < 10.0) { //potiential sequencing error
        if (lowPercent < 20.0 && lowPercent >1.0 ) { //potiential sequencing error, sequencing error should below 1%
            validSeqIDs_mtx.lock();
            validSeqIDs.insert(seqID);
            validSeqIDs_mtx.unlock();
        }
    }

    tmpFile.close();
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

            validSeqIDs_mtx.lock();
            if (validSeqIDs.find(seqID) != validSeqIDs.end()) {
                outputFile << line << std::endl;
                output = true;
            }
            validSeqIDs_mtx.unlock();
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

void countKmersInReadsMultiThreaded(const std::string& kmcDatabase, const std::string& inputFile, int maxThreads) {
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
                threads.push_back(std::thread(processRead, kmcDatabase, read, seqID));
            }
        } else if (extension == ".fa" || extension == ".fasta") {
            if (line.front() == '>') {
                if (!read.empty()) {
                    threads.push_back(std::thread(processRead, kmcDatabase, read, seqID));
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
        threads.push_back(std::thread(processRead, kmcDatabase, read, seqID));
    }

    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }

    file.close();
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <kmc_database> <input_file> <max_threads>" << std::endl;
        return 1;
    }

    const std::string kmcDatabase = argv[1];
    const std::string inputFile = argv[2];
    int maxThreads = std::stoi(argv[3]);

    std::string tmpFilename = inputFile + ".tmp";
    tmpFile.open(tmpFilename);
    if (!tmpFile.is_open()) {
        std::cerr << "Could not open temporary file for writing: " << tmpFilename << std::endl;
        return 1;
    }

    countKmersInReadsMultiThreaded(kmcDatabase, inputFile, maxThreads);

    tmpFile.close();

    double fileMedian = 0.0;
    if (!medians.empty()) {
        size_t m = medians.size() / 2;
        nth_element(medians.begin(), medians.begin() + m, medians.end());
        fileMedian = medians[m];
        if (medians.size() % 2 == 0) {
            nth_element(medians.begin(), medians.begin() + m - 1, medians.end());
            fileMedian = (fileMedian + medians[m - 1]) / 2.0;
        }
    }

    std::cout << "File Median Count: " << fileMedian << std::endl;

    analyzeTmpFile(tmpFilename, fileMedian);

    std::string outputFilename = inputFile + ".filter.fa";
    outputValidSequences(inputFile, outputFilename);

    return 0;
}
