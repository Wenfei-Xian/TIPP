# TIPP
Telomere local assembly, Improved whole genome polish, and Plastid assembly (TIPP)

## Dependencies
Minimap2 https://github.com/lh3/minimap2   
BCFtools https://github.com/samtools/bcftools   
SPOA https://github.com/rvaser/spoa    
MCL https://github.com/micans/mcl   
Samtools https://github.com/samtools/samtools   
Hifiasm https://github.com/chhylp123/hifiasm  
seqtk https://github.com/Wenfei-Xian/seqtk (forked from lh3/seqtk)  
TRF https://github.com/Benson-Genomics-Lab/TRF

## Installation
```
conda create -n TIPP
conda activate TIPP
conda install -c bioconda minimap2
conda install -c bioconda spoa
conda install -c bioconda mcl
conda install -c bioconda r-pheatmap
conda install -c conda-forge r-igraph
conda install -c bioconda bioconductor-biostrings
conda install -c r r-stringdist
conda install -c conda-forge r-ggplot2
conda install -c bioconda trf

git clone https://github.com/Wenfei-Xian/seqtk.git
cd seqtk
make
pwd #(you will get the current path)
export PATH=current_path:$PATH #or you can paste it to your .bashrc, but it will replace the origin seqtk if it's already in your PATH.
cd ..

git clone https://github.com/Wenfei-Xian/TIPP.git
cd TIPP
pwd #(you will get the current path)
export PATH=current_path/src:$PATH #paste it to your .bashrc
```

## Usage   
conda activate TIPP


## Errors I met
1) /tmp is full   
2) lower version of singularity
3) when coverage is high, sopa will be killed because of high memery usage, downsample could be a good option to solve it.
