## Name of TIPP
Upon learning that the Tübingen International PhD Program (TIPP) was no longer enrolling new PhD students, I realized that the three tools I was developing could be combined into a toolkit named Telomere local assembly, Improved whole genome polish, and Plastid assembly (TIPP). This was a fitting tribute to the program. I am profoundly grateful to the TIPP PhD program for providing me the opportunity to embark on my doctoral journey at the Max Planck Institute, and I extend my thanks to my mentor, Detlef Weigel. Although TIPP has concluded, the IMPRS PhD program remains active. Interested individuals can learn more by following this link: https://www.phd.tuebingen.mpg.de/imprs.   

## Dependency for all
Minimap2 https://github.com/lh3/minimap2

## Dependencies for TIPP_telo
BCFtools https://github.com/samtools/bcftools   
SPOA https://github.com/rvaser/spoa    
MCL https://github.com/micans/mcl   
Samtools https://github.com/samtools/samtools   
seqtk https://github.com/Wenfei-Xian/seqtk (forked from lh3/seqtk)  
TRF https://github.com/Benson-Genomics-Lab/TRF   

## Dependency for TIPP_plastid
KMC3 https://github.com/refresh-bio/KMC

## Installation
```
conda create -n TIPP
conda activate TIPP
conda install -c bioconda minimap2 (TIPP_telo,TIPP_plastid)
conda install -c bioconda spoa (TIPP_telo)
conda install -c bioconda mcl (TIPP_telo)
conda install -c bioconda r-pheatmap (TIPP_telo)
conda install -c conda-forge r-igraph (TIPP_telo)
conda install -c bioconda bioconductor-biostrings (TIPP_telo)
conda install -c r r-stringdist (TIPP_telo)
conda install -c conda-forge r-ggplot2 (TIPP_telo)
conda install -c bioconda trf (TIPP_telo)   

git clone https://github.com/Wenfei-Xian/TIPP.git
cd TIPP
bash install.sh
pwd #(you will get the current path)
export PATH=current_path/src:$PATH #paste it to your .bashrc
```

## Chloroplast Database for TIPP_plastid   
https://figshare.com/articles/dataset/Chloroplast_genome_Database_for_TIPP_plastid/24600084   

## Usage   
conda activate TIPP


## Errors I met
1) /tmp is full   
2) lower version of singularity
3) when coverage is high, sopa will be killed because of high memery usage, downsample could be a good option to solve it.
