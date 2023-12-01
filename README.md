## Name of TIPP
When I heard that the Tübingen International PhD Program (TIPP) had stopped recruiting new PhD students, I was in the process of naming three small tools I was working on. Coincidentally, I realized that the three small tools I was developing could be combined to form Telomere local assembly, Improved whole genome polish, and Plastid assembly (TIPP), which I decided to use to commemorate the program that allowed me to continue my PhD studies.   
Although the TIPP program has ceased, the IMPRS PhD program continues! Those interested can read more at the following link: https://www.phd.tuebingen.mpg.de/imprs.   

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
Flye https://github.com/fenderglass/Flye

## Installation
```
conda create -n TIPP
conda activate TIPP
conda install -c bioconda minimap2 (TIPP_telomere,TIPP_plastid)
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
pwd #(you will get the current path)
export PATH=current_path/src:$PATH #paste it to your .bashrc
```

## Chloroplast Database for TIPP_plastid   
https://figshare.com/articles/dataset/Chloroplast_genome_Database_for_TIPP_plastid/24600084   

## Usage   
### TIPP_plastid   
```
Usage: TIPP_plastid.pl [options]
  -h: Show this help message.
  -d: Chloroplast database (required).
  -f: HiFi reads (required).
  -g: Chloroplast or Mitochondrion (default: Chloroplast).
  -t: Threads for Minimap2, Flye, KMC3 and readskmercount.
  -n: Number of reads in each downsample, if the read length is short (<=15kb) or the genome size is extremely large, please increase this value (default: 2000).
  -r: Number of random downsamplings (default: 5).
  -m: Sequencing platform - either 'pacbio' or 'ont'. Only Q20 reads are accepted (default: pacbio).
  -i: Assume the presence of the inverted repeats (default: 1).
  -c: Maximum number of candidate reads used (default: 60000).
  -p: The proportion of total number of M in cigar / the length of reads, greater than this value is considered a match (default: 0.3).

```
### TIPP_telo   
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
