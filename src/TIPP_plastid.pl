#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Long;

my$help;
my$vcf_file;
my$unit;
my$fastq;
my$fastq_name;
my$contigs;
my$contigs_name;
my$bam;
my$vcf_file_name;
my$threads=40;
my$current_directory = cwd();
my$reads=2000;

GetOptions(
    'h'   => \$help,
    'd=s' => \$vcf_file,
    #'u=s' => \$unit,
    'f=s' => \$fastq,
    'c=s' => \$contigs,
    'b=s' => \$bam,
    't=s' => \$threads,
    'n=s' => \$reads,
);

if ($help or not defined $vcf_file and not defined $unit and not defined $fastq and not defined $contigs and not defined $bam) {
    print "Usage: $0\n";
    print "-h: show this help message.\n";
    #print "-d: specify the VCF file from deepvariant.\n";
    #print "-b: bam file\n";
    #print "-u: telomere unit.\n";
    print "-f: hifi reads.\n";
    print "-c: contigs.\n";
    print "-t: threads for minimap2,hifiasm,diamond.\n";
    exit;
}

for my $tool ("minimap2","samtools","diamond","hifiasm","pigz") {
	my $return_val = system("which $tool > /dev/null 2>&1");
	if ($return_val != 0) {
		die "Error: $tool does not seem to be installed or is not in the PATH. Please check.\n";
	}
}

my$script_dir = ($0 =~ m{(.*/)?})[0];

$contigs_name=(split /\//,$contigs)[-1];
#$fastq_name=(split /\//,$fastq)[-1];

unless (-d "$contigs_name.chloroplast") {
        system("mkdir $contigs_name.chloroplast");
}

if( defined $contigs ){
	$fastq_name=(split /\//,$fastq)[-1];
	system("minimap2 -x map-hifi -I 20G -t $threads -c --eqx $contigs $fastq -o $contigs_name.chloroplast/$contigs_name.$fastq_name.paf");
	my$bam_file = "$contigs_name.chloroplast/$contigs_name.$fastq_name.paf";
	my$full_map = "$contigs_name.chloroplast/$contigs_name.$fastq_name.fullmap.ID";
	system("awk '{if(\$3<100 && \$4+100>\$2)print \$1 }' $bam_file | sort | uniq -u > $full_map");
}
else{
	my$full_map = "$contigs_name.chloroplast/$contigs_name.$fastq_name.fullmap.ID";
	system("touch $full_map");
}

if (-e "$script_dir/coregene.fa" && -e "$script_dir/coregene.fa.dmnd" ) {
	print "Chloroplast seed found!\n";
	print "Diamond search starting!\n";
	system("diamond blastx --threads $threads --db $script_dir/coregene.fa --query $fastq --max-target-seqs 50 --outfmt 6 qseqid full_qseq --out $contigs_name.chloroplast/$fastq_name.chloroplast.diamond.out");
}
else{
	print "Seed does not exist!";
	exit;
}

system("pigz $contigs_name.chloroplast/$fastq_name.chloroplast.diamond.out");
system("zcat $contigs_name.chloroplast/$fastq_name.chloroplast.diamond.out | awk '{print \$1}' | sort | uniq -c |awk '{if(\$1 >=3) print \$2}' > $contigs_name.chloroplast/$fastq_name.chloroplast.candidate.ID");

open my$diamond_canID,"$contigs_name.chloroplast/$fastq_name.chloroplast.candidate.ID" or die "Can't open $contigs_name.chloroplast/$fastq_name.chloroplast.candidate.ID";
my%diamond_canID_hash;
while(<$diamond_canID>){
	chomp;
	$diamond_canID_hash{$_}=$_;
}
close $diamond_canID;

open my$full_map_ID,"$contigs_name.chloroplast/$contigs_name.$fastq_name.fullmap.ID" or die "Can't open $contigs_name.chloroplast/$contigs_name.$fastq_name.fullmap.ID";
my%full_map_hash;
while(<$full_map_ID>){
	chomp;
	$full_map_hash{$_}=$_;
}
close $full_map_ID;

open my $diamond_out, "gunzip -dc $contigs_name.chloroplast/$fastq_name.chloroplast.diamond.out.gz|" or die "Can't open $contigs_name.chloroplast/$fastq_name.chloroplast.diamond.out.gz";
open my $chloroplast_out, ">$contigs_name.chloroplast/$fastq_name.chloroplast.reads.fasta";
my%hash_single;
while(<$diamond_out>){
	chomp;
	if (/(\S+)\s+(\S+)/) {
		my ($diamond_out_id, $diamond_out_seq) = ($1, $2);
		if (exists $diamond_canID_hash{$diamond_out_id} && !exists $full_map_hash{$diamond_out_id}) {
			if( !exists $hash_single{$diamond_out_id} ){
				print $chloroplast_out ">$diamond_out_id\n$diamond_out_seq\n";
				$hash_single{$diamond_out_id}=$diamond_out_id;
			}
		}
	}
}
close $diamond_out;
close $chloroplast_out;

open my $chloroplast, ">$contigs_name.chloroplast/Chloroplast.fa" or die "Cannot open $contigs_name.chloroplast/Chloroplast.fa: $!\n";

for (my $round=1; $round<50; $round++){
    print "Starting downsampling of chloroplast reads, round$round...\n";
    system("seqtk sample -s$round $contigs_name.chloroplast/$fastq_name.chloroplast.reads.fasta $reads > $contigs_name.chloroplast/$fastq_name.chloroplast.reads.subsample$reads.round$round.fasta");

    system("hifiasm -D 1000 -o $contigs_name.chloroplast/$fastq_name.chloroplast.reads.subsample$reads.round$round -t 40 $contigs_name.chloroplast/$fastq_name.chloroplast.reads.subsample$reads.round$round.fasta");

    my $result = `grep '^S' $contigs_name.chloroplast/$fastq_name.chloroplast.reads.subsample$reads.round$round.bp.p_ctg.gfa | grep 'c\\s'`;

    if ($result) {
        print "Circular contig found in round$round.\n";
        my ($ID, $seq) = (split /\t/, $result)[1,2];
        print $chloroplast ">$ID\n$seq\n";
        last;
    }
}
close $chloroplast;
