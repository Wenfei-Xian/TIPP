# TIPPo: A User-Friendly Tool for De Novo Assembly of Organellar Genomes with HiFi Data
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/tipp) [![image](https://anaconda.org/bioconda/tipp/badges/downloads.svg)](https://anaconda.org/bioconda/tipp)

## Contents ##

* [General Introduction](#Name-of-TIPP)
* [Dependencies for TIPP_telomere](#Dependencies-for-TIPP_telomere)
* [Dependency for TIPPo](#Dependency-for-TIPPo)
* [Installation and Run](#Installation)
    * [General](#Installation)
    * [Note for Windows User](#Note-for-Windows-User)
    * [Chloroplast and Mitochondrial Database for TIPP_plastid  (only for v1)](#Chloroplast-and-Mitochondrial-Database-for-TIPP_plastid  (only for v1))
    * [Docker and Singularity (recommended)](#Docker-and-Singularity (recommended))
      * [Run TIPP with Docker](#Run-TIPP-with-Docker)
      * [Run TIPP with Apptainer/Singularity](#Run-TIPP-with-Apptainer/Singularity)
* [Usage](#Usage)
   * [TIPPo](#TIPPo-usage)
   * [TIPP_telomere](#TIPP_telomere)
* [Output](#Output)
   * [TIPPo](#TIPPo-output)
* [Errors I met](#Errors-I-met)
* [Citation](#Citation)
* [Please cite the dependencies if you use TIPPo](#Please-cite-the-dependencies-if-you-use-TIPPo) 


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

## Dependency for TIPPo
KMC3 https://github.com/refresh-bio/KMC   
Flye https://github.com/fenderglass/Flye   
Graphaligner https://github.com/maickrau/GraphAligner   
TIARA https://github.com/ibe-uw/tiara   
TRF https://github.com/Benson-Genomics-Lab/TRF  
Minimap2 https://github.com/lh3/minimap2  

## Installation

## Via Conda [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/tipp)

```
conda create -n TIPP python=3.8 #please specific the python version to 3.8 :)
conda activate TIPP
conda install bioconda::tipp
```
## Manual

```
conda create -n TIPP python=3.8
conda activate TIPP
pip install tiara (TIPPo)
conda install bioconda::graphaligner (TIPPo)
conda install bioconda::trf (TIPPo & TIPP_telomere)
conda install bioconda::flye (TIPPo)
conda install bioconda::minimap2 (TIPPo & TIPP_telomere)
conda install -c bioconda spoa (TIPP_telomere)
conda install -c bioconda mcl (TIPP_telomere)
conda install -c bioconda r-pheatmap (TIPP_telomere)
conda install -c conda-forge r-igraph (TIPP_telomere)
conda install -c bioconda bioconductor-biostrings (TIPP_telomere)
conda install -c r r-stringdist (TIPP_telomere)
conda install -c conda-forge r-ggplot2 (TIPP_telomere)

git clone https://github.com/Wenfei-Xian/TIPP.git
cd TIPP
bash install.sh (Compile KMC3, seqtk, and readskmercount)
CURRENT_PATH=$(pwd)/src
echo "export PATH=$CURRENT_PATH:\$PATH" >> ~/.bashrc
echo "export PATH=$CURRENT_PATH/seqtk:\$PATH" >> ~/.bashrc
echo "export PATH=$CURRENT_PATH/kmc3/bin:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

## Note for Windows User

<div align="justify">
  
If you are utilizing a Windows-based system or laptop, it is recommended to run a Linux emulator for optimal compatibility. One approach is to set up a Virtual Machine with a Linux distribution such as Ubuntu. To accomplish this, you can explore the usage of [VirtualBox](https://www.virtualbox.org/), a widely-used virtualization software. VirtualBox allows you to create and manage virtual machines on your Windows operating system. For more information and to download VirtualBox, please visit their official website at [VirtualBox](https://www.virtualbox.org/). 

Moreover, it's worth mentioning that Ubuntu is also available as [Windows Subsystem for Linux (WSL)](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support#1-overview). WSL allows you to install a complete Ubuntu terminal environment on your Windows machine, enabling you to seamlessly develop and run cross-platform applications without leaving the Windows environment. Installing WSL is a breeze, whether you opt for PowerShell or the Windows Store app. As of November 2022, WSL is now available as a Windows Store app for both Windows 10 and Windows 11. This means you have multiple convenient options to set up Ubuntu on your Windows system for an enhanced and hassle-free experience.

Please note that the availability of Ubuntu through WSL offers flexibility and ease of use, allowing you to run TIPPo even as a Windows user.

So, whether you choose VirtualBox or WSL, integrating Linux into your Windows environment will enable you to seamlessly carry out the TIPPo and ensure a smooth and efficient analysis process. </div>

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
An example run of TIPPo would work as follows:
```
(base) root@84037bc1f36c:~# cd
# Download an example dataset
wget -O Arabidopsis_thaliana.4X.fastq.gz https://figshare.com/ndownloader/files/47427487
# Run TIPPo command
TIPPo.v2.3.pl -f Arabidopsis_thaliana.4X.fastq.gz
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
# Run TIPPo command
TIPPo.v2.3.pl -f Arabidopsis_thaliana.4X.fastq.gz
```

By default, apptainer mounts your home and the current working directory inside the container. Thus, your data will be stored on the host file system. You can mount additional volumes as required (see [Apptainer documentation](https://apptainer.org/docs/user/main/bind_paths_and_mounts.html) ).


## Usage   
### TIPPo-usage   
```
TIPPo.v2.4.pl 
Usage: /tmp/global2/wxian/conda/envs/tipp13/bin/TIPPo.v2.4.pl [options]
  -h: Show this help message.
  -f: Long reads (required).
  -d: reference sequence (default: None).
  -g: chloroplast or organelle (default: organelle).
  -t: Threads for tiara, flye, KMC3 and readskmercount.
  -n: Number of reads in each downsample for chloroplast.
  -r: Number of random downsamplings (default: 5).
  -p: Sequence technology - 'hifi','clr','ont', 'onthq' (default: hifi)
  -i: Assume the presence of the inverted repeats in the chloroplast genome (default: 1).
  -l: lower kmer count - lkc (default: 0.3).
  -c: high kmer count - hkc (default: 5).
  -m: minimum overlap in repeat graph construction (default:800)
  --trf: remove the reads are tandem repeats, only avaliable for reference-free and hifi/onthq reads
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

## Output
### TIPPo-output
For the **chloroplast genome**, you can find the file at:   
```
your_input_file.chloroplast.fasta.filter.800.round1.edge_*.edge_*.edge_*.organelle.chloroplast.fasta
```
For the **mitochondrial genome**, the file is located at:   
```
your_input_file.mitochondrial.fasta.filter.fasta.flye/assembly_graph.gfa
```
If you are interested in exploring the **complexity of the mitochondrial genome**, you can find the file at:
```
your_input_file.mitochondrial.fasta.filter.fasta.flye/50.repeat-graph/graph_before_rr.gfa
```
### GFA → Linearized FASTA
#### Detects circular paths in a GFA assembly graph and outputs a linearized FASTA sequence starting at the longest segment.
```
```


## Errors I met

1) /tmp is full   
2) lower version of singularity
3) when coverage is high, sopa will be killed because of high memory usage, downsample could be a good option to solve it.
4) storage is full.
5) g++: error: unrecognized command line option '-std=c++14', please update the gcc - conda install gcc_linux-64 (KMC3)
6) plase use seqtk https://github.com/Wenfei-Xian/seqtk (forked from lh3/seqtk) 

## Citation

Citation for TIPPo: Wenfei Xian, Ilja Bezrukov, Zhigui Bao, Sebastian Vorbrugg, Anupam Gautam, Detlef Weigel, TIPPo: A User-Friendly Tool for De Novo Assembly of Organellar Genomes with High-Fidelity Data, ***Molecular Biology and Evolution***, Volume 42, Issue 1, January 2025, msae247,[https://doi.org/10.1093/molbev/msae247](https://doi.org/10.1093/molbev/msae247)

### Please cite the dependencies if you use TIPPo:
Flye: https://www.nature.com/articles/s41587-019-0072-8   
Graphaligner: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02157-2   

Reference free:  
Tiara: https://academic.oup.com/bioinformatics/article/38/2/344/6375939  
KMC: https://academic.oup.com/bioinformatics/article/33/17/2759/3796399  
Reference based:  
Minimap2: https://academic.oup.com/bioinformatics/article/34/18/3094/4994778  
