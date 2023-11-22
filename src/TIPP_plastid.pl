#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Long;

my$help;
my$fastq;
my$fastq_name;
my$threads=40;
my$current_directory = cwd();
my$reads_num=2000;
my$db;
my$round_cutoff=5;
my$no_IR=0;

GetOptions(
    'h'   => \$help,
    'f=s' => \$fastq,
    't=s' => \$threads,
    'n=s' => \$reads_num,
    'd=s' => \$db,
    'i=s' => \$no_IR,
);

if ($help or not defined $fastq or not defined $db) {
    print "Usage: $0\n";
    print "-h: Show this help message.\n";
    print "-d: Chloroplast database (required).\n";
    print "-f: HiFi reads (required).\n";
    print "-t: Threads for Minimap2, Flye, KMC3 and readskmercount.\n";
    print "-n: Number of reads in each downsample (default:2000).\n";
    print "-i: Assume the presence of the inverted repeats (default: 1).\n ";
    exit;
}

for my $tool ("minimap2","samtools","flye") {
	my $return_val = system("which $tool > /dev/null 2>&1");
	if ($return_val != 0) {
		die "Error: $tool does not seem to be installed or is not in the PATH. Please check.\n";
	}
}

my$script_dir = ($0 =~ m{(.*/)?})[0];

$fastq_name=(split /\//,$fastq)[-1];
my$fastq_name_reads="$fastq_name.Chloroplast.reads.fa";

unless (-d "$fastq_name.Chloroplast") {
        system("mkdir $fastq_name.Chloroplast");
}

my%hash_uniq_seq;
open IN1,"minimap2 --secondary no --sam-hit-only -ax map-hifi -t $threads $db $fastq_name |";
open OUT1,">$fastq_name.Chloroplast/$fastq_name_reads";
while(<IN1>){
        chomp;
        next if(/^@/);
        my@sam=split /\t/,$_;
        next if( exists $hash_uniq_seq{$sam[0]} );
        $hash_uniq_seq{$sam[0]}=0;
        my$total_matches = calculate_total_length_of_matches($sam[5]);
        my$query_len=length$sam[5];
        if($total_matches >= $query_len*0.5){
                print OUT1 ">$sam[0]\n$sam[9]\n";
        }
}

sub calculate_total_length_of_matches {
    my $cigar = shift;
    my $total_length = 0;

    while ($cigar =~ /(\d+)M/g) {
        $total_length += $1;
    }

    return $total_length;
}

system("mkdir $fastq_name.Chloroplast/$fastq_name_reads.tmp");
system("$script_dir/kmc3/bin/kmc -k31 -cs999999 -ci3 -fa $fastq_name.Chloroplast/$fastq_name_reads $fastq_name.Chloroplast/$fastq_name_reads $fastq_name.Chloroplast/$fastq_name_reads.tmp");
system("rm -rf $fastq_name.Chloroplast/$fastq_name_reads.tmp");

system("$script_dir/readskmercount $fastq_name.Chloroplast/$fastq_name_reads $fastq_name.Chloroplast/$fastq_name_reads $threads");

my%hash_assembly_round;
my$found=0;
for (my $round=1; $round<=$round_cutoff; $round++){

	system("seqtk sample -s$round $fastq_name.Chloroplast/$fastq_name_reads.filter.fa $reads_num > $fastq_name.Chloroplast/$fastq_name_reads.filter.$reads_num.round$round.fa");
	
	system("flye --pacbio-hifi $fastq_name.Chloroplast/$fastq_name_reads.filter.$reads_num.round$round.fa --threads 40 -o $fastq_name.Chloroplast/$fastq_name_reads.filter.$reads_num.round$round.fa.flye");
	
	if( $no_IR == 0 ){ # assumed that the chloroplast carries inverted repeats
	
		open my$chloro_assem,"$fastq_name.Chloroplast/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa" or die "Can't open $fastq_name.Chloroplast/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa";
		my%hash_gfa;
	
		while(<$chloro_assem>){
			chomp;
			if($_=~m/^L/){
				my@linkage=split /\t/,$_;
				$hash_gfa{$linkage[1]}+=1;
				$hash_gfa{$linkage[3]}+=1;
			}
		}

		my $number_edge = keys %hash_gfa;

		if ($number_edge == 3) {
			my %value_counts;
			foreach my $value (values %hash_gfa) {
				$value_counts{$value}++;
			}

			if ($value_counts{2} == 2 && $value_counts{4} == 1) {
 				print "Chloroplast genome found\n";
				$found=1;
				last;
			}
		}
	}
	elsif( $no_IR == 1 ){
		open my$chloro_assem,"$fastq_name.Chloroplast/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa" or die "Can't open $fastq_name.Chloroplast/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa";
		my%hash_gfa;
		while(<$chloro_assem>){
			chomp;
			if($_=~m/^L/){
				my@linkage=split /\t/,$_;
				$hash_gfa{$linkage[1]}+=1;
				$hash_gfa{$linkage[3]}+=1;
			}
		}

		my $number_edge = keys %hash_gfa;
		if ($number_edge == 1) {
			if($no_IR == 1){
				print "Chloroplast genome found\n";
				$found=1;
				last;
			}
		}
	}
	else{
		open my$chloro_assem,"$fastq_name.Chloroplast/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa" or die "Can't open $fastq_name.Chloroplast/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa";
		my%hash_gfa;
		while(<$chloro_assem>){
			chomp;
			if($_=~m/^L/){
				my@linkage=split /\t/,$_;
				$hash_gfa{$linkage[1]}+=1;
				$hash_gfa{$linkage[3]}+=1;
			}
		}

		my $number_edge = keys %hash_gfa;
		if ($number_edge == 1) {
			$hash_assembly_round{$round}=$round;
		}
	}
}

if ($found == 0 && keys %hash_assembly_round > 0) {
	print "After $round_cutoff rounds of downsampling and assembly, we couldn't find an assembly with an inverted repeat, but we found an assembly without an inverted repeat.\n";
	print "You can find the assembly in round ";
	foreach my $round_assem (keys %hash_assembly_round) {
		print "$round_assem ";
	}
	print "\n";
}
elsif ($found == 0 && keys %hash_assembly_round == 0) {
	print "After $round_cutoff rounds of downsampling and assembly, we couldn't find any assembly that meets our expectations.\n";
	print "You can increase the value of -n and try again :)\n";
}
