## Name of TIPP
When I heard that the Tübingen International PhD Program (TIPP) had stopped recruiting new PhD students, I was in the process of naming three small tools I was working on. Coincidentally, I realized that the three small tools I was developing could be combined to form Telomere local assembly, Improved whole genome polish, and Plastid assembly (TIPP), which I decided to use to commemorate the program that allowed me to continue my PhD studies.   
Although the TIPP program has ceased, the IMPRS PhD program continues! Those interested can read more at the following link: https://www.phd.tuebingen.mpg.de/imprs.   

## Dependencies for TIPP_telomere
BCFtools https://github.com/samtools/bcftools   
SPOA https://github.com/rvaser/spoa    
MCL https://github.com/micans/mcl   
Samtools https://github.com/samtools/samtools   
Minimap2 https://github.com/lh3/minimap2   
seqtk https://github.com/Wenfei-Xian/seqtk (forked from lh3/seqtk)  
TRF https://github.com/Benson-Genomics-Lab/TRF   

## Dependency for TIPP_plastid
KMC3 https://github.com/refresh-bio/KMC   
Flye https://github.com/fenderglass/Flye   
Graphaligner https://github.com/maickrau/GraphAligner   
TIARA https://github.com/ibe-uw/tiara   

## Installation
```
conda create -n TIPP
conda activate TIPP python=3.7
pip install tiara (TIPP_plastid)
conda install graphaligner (TIPP_plastid)
conda install flye (TIPP_plastid)
conda install -c bioconda minimap2 (TIPP_telomere)
conda install -c bioconda spoa (TIPP_telomere)
conda install -c bioconda mcl (TIPP_telomere)
conda install -c bioconda r-pheatmap (TIPP_telomere)
conda install -c conda-forge r-igraph (TIPP_telomere)
conda install -c bioconda bioconductor-biostrings (TIPP_telomere)
conda install -c r r-stringdist (TIPP_telomere)
conda install -c conda-forge r-ggplot2 (TIPP_telomere)
conda install -c bioconda trf (TIPP_telomere)   

git clone https://github.com/Wenfei-Xian/TIPP.git
cd TIPP
bash install.sh (Compile KMC3, seqtk, and readskmercount)
CURRENT_PATH=$(pwd)/src
echo "export PATH=$CURRENT_PATH:\$PATH" >> ~/.bashrc
echo "export PATH=$CURRENT_PATH/seqtk:\$PATH" >> ~/.bashrc
echo "export PATH=$CURRENT_PATH/kmc3/bin:\$PATH" >> ~/.bashrc
```

## Chloroplast and Mitochondrial Database for TIPP_plastid  (only for v1)  
```
wget https://figshare.com/ndownloader/files/44102183 -O Plant_chlo_mito.fa.gz
```
or visit
https://figshare.com/articles/dataset/Plant_Chloroplast_and_Mitochondrial_Genomes/25018469   


## Usage   
### TIPP_plastid   
```
Usage: TIPP_plastid.v2.1.pl [options]
  -h: Show this help message.
  -f: HiFi reads (required).
  -g: chloroplast or organelle (default: organelle).
  -t: Threads for tiara, flye, KMC3 and readskmercount.
  -n: Number of reads in each downsample for chloroplast.
  -r: Number of random downsamplings (default: 5).
  -p: Sequencing platform - either 'pacbio' or 'ont'. Only Q20 reads are accepted (default: pacbio).
  -i: Assume the presence of the inverted repeats (default: 1).
  -v: version.

```
### TIPP_telomere   
```
Usage: TIPP_telomere.pl
-h: show this help message.
-u: telomere unit.
-f: hifi reads.
-e: extend the contigs with new assembled telomere sequences.(default=0, no extend)
-c: contigs. If the extension is not specified, it will be used as the output name.
-t: threads for minimap2.
-m: minimum length of uniq sequence (without telomere)

```

## Errors I met
1) /tmp is full   
2) lower version of singularity
3) when coverage is high, sopa will be killed because of high memery usage, downsample could be a good option to solve it.
4) storage is full.
5) g++: error: unrecognized command line option '-std=c++14', please update the gcc - conda install gcc_linux-64 (KMC3)
6) plase use seqtk https://github.com/Wenfei-Xian/seqtk (forked from lh3/seqtk) 

## Citation
Citation for TIPP_plastid: https://www.biorxiv.org/content/10.1101/2024.01.29.577798v1

### Please cite the dependencies if you use TIPP_plastid:   
Flye: https://www.nature.com/articles/s41587-019-0072-8   
KMC: https://academic.oup.com/bioinformatics/article/33/17/2759/3796399   
Minimap2: https://academic.oup.com/bioinformatics/article/34/18/3094/4994778    
Diamond: https://www.nature.com/articles/s41592-021-01101-x
Graphaligner: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02157-2
