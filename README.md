# TIPP
Telomere local assembly, Improved whole genome polish, and Plastid assembly

## Dependencies
Diamond https://github.com/bbuchfink/diamond   
Minimap2 https://github.com/lh3/minimap2   
BCFtools https://github.com/samtools/bcftools   
SPOA https://github.com/rvaser/spoa    
MCL https://github.com/micans/mcl   
Samtools https://github.com/samtools/samtools   
Hifiasm https://github.com/chhylp123/hifiasm  
seqtk https://github.com/Wenfei-Xian/seqtk (forked from lh3/seqtk)  
Deepvariant https://github.com/google/deepvariant   
Singularity https://github.com/sylabs/singularity   

## Installation
```
conda create -n TIPP
conda activate TIPP
conda clean --all # optional
conda install -c bioconda diamond
conda install -c bioconda minimap2
conda install -c bioconda spoa
conda install -c bioconda mcl
conda install -c bioconda hifiasm
conda install -c bioconda r-pheatmap

#conda install -c conda-forge singularity # Until October 24, 2023, the available version of Singularity is 3.8.7. Ensure the version you downloaded is not lower than 3.8.7.
#conda install -c bioconda samtools# not recommond via conda
#conda install -c bioconda bcftools# not recommond via conda

git clone https://github.com/Wenfei-Xian/seqtk.git
cd seqtk
make
pwd #(you will get the current path)
export PATH=current_path:$PATH #or you can paste it to your .bashrc, but it will replace the origin seqtk if it's already in your PATH.
cd ..

git clone https://github.com/Wenfei-Xian/TIPP.git
cd TIPP
pwd #(you will get the current path)
export PATH=current_path:$PATH #paste it to your .bashrc
BIN_VERSION="1.6.0"
singularity cache clean # before pull down, I usually clean the cache
singularity pull docker://google/deepvariant:"${BIN_VERSION}"   
```
While there are quite a few dependencies, they are all easy to install, and I believe some of the software in there should be familiar to you :)

## Usage   
conda activate TIPP


## Potiential error I meet
1) /tmp is full   
2) lower version of singularity
