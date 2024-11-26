# TIPPo: A User-Friendly Tool for De Novo Assembly of Organellar Genomes with HiFi Data
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/tipp) [![Downloads]([https://anaconda.org/bioconda/tipp/badges/downloads.svg](https://anaconda.org/bioconda/tipp/badges/downloads.svg))](https://anaconda.org/bioconda/tipp)
### Note:- Conda currently includes version 1.1.0, which does not have TIPPo.v2.2.pl. We will update it soon

## Contents ##

* [General Introduction](#Name-of-TIPP)
* [Dependencies for TIPP_telomere](#Dependencies-for-TIPP_telomere)
* [Dependency for TIPP_plastid_v2.1](#Dependency-for-TIPP_plastid_v2.1)
* [Installation and Run](#Installation)
    * [General](#Installation)
    * [Note for Windows User](#Note-for-Windows-User)
    * [Chloroplast and Mitochondrial Database for TIPP_plastid  (only for v1)](#Chloroplast-and-Mitochondrial-Database-for-TIPP_plastid  (only for v1))
    * [Docker and Singularity (recommended)](#Docker-and-Singularity (recommended))
      * [Run TIPP with Docker](#Run-TIPP-with-Docker)
      * [Run TIPP with Apptainer/Singularity](#Run-TIPP-with-Apptainer/Singularity)
* [Usage](#Usage)
   * [TIPP_plastid](#TIPP_plastid)
   * [TIPP_telomere](#TIPP_telomere)
* [Errors I met](#Errors-I-met)
* [Citation](#Citation)
* [Please cite the dependencies if you use TIPP_plastid](#Please-cite-the-dependencies-if-you-use-TIPP_plastid) 


## Name of TIPP ##
<div align="justify">

When I heard that the Tübingen International PhD Program (TIPP) had stopped recruiting new PhD students, I was in the process of naming three small tools I was working on. Coincidentally, I realized that the three small tools I was developing could be combined to form Telomere local assembly, Improved whole genome polish, and Plastid assembly (TIPP), which I decided to use to commemorate the program that allowed me to continue my PhD studies.   
Although the TIPP program has ceased, the IMPRS PhD program continues! Those interested can read more at the following link: https://www.phd.tuebingen.mpg.de/imprs.
</div>

## Dependencies for TIPP_telomere
BCFtools https://github.com/samtools/bcftools   
SPOA https://github.com/rvaser/spoa    
MCL https://github.com/micans/mcl   
Samtools https://github.com/samtools/samtools   
Minimap2 https://github.com/lh3/minimap2   
seqtk https://github.com/Wenfei-Xian/seqtk (forked from lh3/seqtk)  
TRF https://github.com/Benson-Genomics-Lab/TRF   

## Dependency for TIPP_plastid_v2.1
KMC3 https://github.com/refresh-bio/KMC   
Flye https://github.com/fenderglass/Flye   
Graphaligner https://github.com/maickrau/GraphAligner   
TIARA https://github.com/ibe-uw/tiara   

## Reference free approach TIPP_plastid_v2.1, reference rely TIPP_plastid_v1


 

## Installation

## Via Conda [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/tipp)

### Note:- Conda currently includes version 1.1.0, which does not have TIPPo.v2.2.pl. We will update it soon

```
conda create -n TIPP
conda activate TIPP
conda install bioconda::tipp
```
## Manual

```
conda create -n TIPP python=3.7
conda activate TIPP
pip install tiara (TIPP_plastid)
conda install graphaligner (TIPP_plastid)
conda install bioconda::trf (TIPP_plastid)
conda install bioconda::flye (TIPP_plastid)
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

## Note for Windows User

<div align="justify">
  
If you are utilizing a Windows-based system or laptop, it is recommended to run a Linux emulator for optimal compatibility. One approach is to set up a Virtual Machine with a Linux distribution such as Ubuntu. To accomplish this, you can explore the usage of [VirtualBox](https://www.virtualbox.org/), a widely-used virtualization software. VirtualBox allows you to create and manage virtual machines on your Windows operating system. For more information and to download VirtualBox, please visit their official website at [VirtualBox](https://www.virtualbox.org/). 

Moreover, it's worth mentioning that Ubuntu is also available as [Windows Subsystem for Linux (WSL)](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support#1-overview). WSL allows you to install a complete Ubuntu terminal environment on your Windows machine, enabling you to seamlessly develop and run cross-platform applications without leaving the Windows environment. Installing WSL is a breeze, whether you opt for PowerShell or the Windows Store app. As of November 2022, WSL is now available as a Windows Store app for both Windows 10 and Windows 11. This means you have multiple convenient options to set up Ubuntu on your Windows system for an enhanced and hassle-free experience.

Please note that the availability of Ubuntu through WSL offers flexibility and ease of use, allowing you to run TIPP_plastid even as a Windows user.

So, whether you choose VirtualBox or WSL, integrating Linux into your Windows environment will enable you to seamlessly carry out the TIPP_plastid and ensure a smooth and efficient analysis process. </div>

## Chloroplast and Mitochondrial Database for TIPP_plastid  (only for v1)  
```
wget https://figshare.com/ndownloader/files/44102183 -O Plant_chlo_mito.fa.gz
```
or visit
https://figshare.com/articles/dataset/Plant_Chloroplast_and_Mitochondrial_Genomes/25018469   

## Docker and Singularity (recommended)
<div align="justify">
   
The easiest way to install TIPP is with a Docker container. A image is available at [Dockerhub](https://hub.docker.com/r/weigelworld/tipp) and can be run via Docker or Singularity/Apptainer. A Dockerfile is included in the Github repository, if you prefer to build the image yourself.
</div>

### Run TIPP with Docker
```
# Pull the image from Dockerhub
docker pull weigelworld/tipp:latest

# Start an interactive shell in the container
docker run -it weigelworld/tipp:latest
```

Inside the container, the tools described below are directly accessible via command line.
An example run of TIPP_plastid would work as follows:
```
(base) root@84037bc1f36c:~# cd
# Download an example dataset
wget -O Arabidopsis_thaliana.4X.fastq.gz https://figshare.com/ndownloader/files/47427487
# Run TIPP_plastid command
TIPP_plastid.v2.1.pl -f Arabidopsis_thaliana.4X.fastq.gz
```
The above example does not pass data to or from the host system -- please check the [Docker documentation](https://docs.docker.com/storage/bind-mounts/) on how to exchange data between the host system and the container.

As an example, the following command would run an interactive shell in the container and mount the directory `/data/experiments` on the host system and make it available under `/mnt` inside the container:
```
docker run -it -v /data/experiments:/mnt weigelworld/tipp:latest
```

### Run TIPP with Apptainer/Singularity
Apptainer, or formerly Singularity can be used to easily run TIPP when Docker is not available. If it is not installed on your system, it can be installed without super user rights, e.g. via conda from conda-forge.

The procedure to run TIPP via Singularity/Apptainer is similar to the above. Here, we show an example with apptainer, but the same parameters work with singularity.

```
# Pull the image from Dockerhub
apptainer pull docker://weigelworld/tipp:latest

# Start an interactive shell in the container
apptainer run tipp_latest.sif
```
Inside the container, the same example commands shown above can be run:
```
# Download an example dataset
wget -O Arabidopsis_thaliana.4X.fastq.gz https://figshare.com/ndownloader/files/47427487
# Run TIPP_plastid command
TIPP_plastid.v2.1.pl -f Arabidopsis_thaliana.4X.fastq.gz
```

By default, apptainer mounts your home and the current working directory inside the container. Thus, your data will be stored on the host file system. You can mount additional volumes as required (see [Apptainer documentation](https://apptainer.org/docs/user/main/bind_paths_and_mounts.html) ).


## Usage   
### TIPPo   
```
Usage: ./TIPPo.v2.2.pl [options]
  -h: Show this help message.
  -f: HiFi reads (required).
  -g: chloroplast or organelle (default: organelle).
  -t: Threads for tiara, flye, KMC3 and readskmercount.
  -n: Number of reads in each downsample for chloroplast.
  -r: Number of random downsamplings (default: 5).
  -p: Sequencing platform - either 'pacbio' or 'ont'. Only Q20 reads are accepted (default: pacbio).
  -i: Assume the presence of the inverted repeats (default: 1).
  -l: lower kmer count - lkc (default: 0.3).
  -c: high kmer count - hkc (default: 5).
  -m: minimum overlap in repeat graph construction (default:800)
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
3) when coverage is high, sopa will be killed because of high memory usage, downsample could be a good option to solve it.
4) storage is full.
5) g++: error: unrecognized command line option '-std=c++14', please update the gcc - conda install gcc_linux-64 (KMC3)
6) plase use seqtk https://github.com/Wenfei-Xian/seqtk (forked from lh3/seqtk) 

## Citation

Citation for TIPP_plastid: [https://www.biorxiv.org/content/10.1101/2024.01.29.577798v2](https://www.biorxiv.org/content/10.1101/2024.01.29.577798v2)

### Please cite the dependencies if you use TIPP_plastid:
Flye: https://www.nature.com/articles/s41587-019-0072-8   
KMC: https://academic.oup.com/bioinformatics/article/33/17/2759/3796399   
Graphaligner: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02157-2   
Tiara: https://academic.oup.com/bioinformatics/article/38/2/344/6375939
