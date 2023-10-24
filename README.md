# TIPP
Telomere local assembly, Improved whole genome polish, and Plastid assembly

# Installation
## Dependencies
Diamond https://github.com/bbuchfink/diamond   
Minimap2 https://github.com/lh3/minimap2   
BCFtools https://github.com/samtools/bcftools   
SPOA https://github.com/rvaser/spoa    
MCL https://github.com/micans/mcl   
Samtools https://github.com/samtools/samtools   
Hifiasm https://github.com/chhylp123/hifiasm  
Above dependencies can be easily installed using conda, and I believe some of them have already been added to your PATH environment variable :)   

seqtk https://github.com/Wenfei-Xian/seqtk (forked from lh3/seqtk)  
```
git clone https://github.com/Wenfei-Xian/seqtk.git
cd seqtk
make
export PATH=path/seqtk:$PATH
```
