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
conda install -c bioconda diamond
conda install -c bioconda minimap2
conda install -c bioconda bcftools
conda install -c bioconda spoa
conda install -c bioconda mcl
conda install -c bioconda samtools
conda install -c bioconda hifiasm
conda install -c conda-forge singularity
git clone https://github.com/Wenfei-Xian/seqtk.git
cd seqtk
make
export PATH=path/seqtk:$PATH
git clone https://github.com/Wenfei-Xian/TIPP.git
cd TIPP
singularity pull docker://google/deepvariant:1.5.0
```

