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
my$platform='pacbio';
my$organelle='Chloroplast';
my$cap=60000;

GetOptions(
    'h'   => \$help,
    'f=s' => \$fastq,
    't=s' => \$threads,
    'n=s' => \$reads_num,
    'd=s' => \$db,
    'i=s' => \$no_IR,
    'm=s' => \$platform,
    'c=s' => \$cap,
    'r=s' => \$round_cutoff,
    'g=s' => \$organelle,
);

if( $platform ne 'pacbio' && $platform ne 'ont' ){
        print "Error, only pacbio or ont can be accecpted for -m\n";
        exit;
}

if( $organelle ne 'Chloroplast' && $organelle ne 'Mitochondrion' ){
        print "Error, only Chloroplast or Mitochondrion can be accecpted for -g\n";
        exit;
}

if ($help or not defined $fastq or not defined $db) {
    print "Usage: $0 [options]\n";
    print "  -h: Show this help message.\n";
    print "  -d: Chloroplast database (required).\n";
    print "  -f: HiFi reads (required).\n";
    print "  -g: Chloroplast or Mitochondrion (default: Chloroplast).\n";
    print "  -t: Threads for Minimap2, Flye, KMC3 and readskmercount.\n";
    print "  -n: Number of reads in each downsample, if the read length is short (<=15kb) or the genome size is extremely large, please increase this value (default: 2000).\n";
    print "  -r: Number of random downsamplings (default: 5).\n";
    print "  -m: Sequencing platform - either 'pacbio' or 'ont'. Only Q20 reads are accepted (default: pacbio).\n";
    print "  -i: Assume the presence of the inverted repeats (default: 1).\n";
    print "  -c: Maximum number of candidate reads used (default: 60000).\n";
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
my$fastq_name_reads="$fastq_name.$organelle.reads.fa";

unless (-d "$fastq_name.$organelle") {
        system("mkdir $fastq_name.$organelle");
}

my%hash_uniq_seq;
my$number_of_reads=0;

=header
if( $platform eq 'pacbio' ){
	open IN1,"minimap2 --secondary no --sam-hit-only -ax map-hifi -t $threads $db $fastq_name |";
}
else{
	open IN1,"minimap2 --secondary no --sam-hit-only -ax map-ont -t $threads $db $fastq_name |";
}
open OUT1,">$fastq_name.$organelle/$fastq_name_reads";
while(<IN1>){
        chomp;
        next if(/^@/);
        my@sam=split /\t/,$_;
        next if exists $hash_uniq_seq{$sam[0]};
        $hash_uniq_seq{$sam[0]} = undef;
	#my$total_matches = calculate_total_length_of_matches($sam[5]);
        my$query_len=length$sam[9];
	#if( $total_matches >= $query_len*0.1 && $sam[2]=~m/$organelle/ ){
	if( $sam[2]=~m/$organelle/ ){
		$number_of_reads++;
                print OUT1 ">$sam[0]\n$sam[9]\n";
		if( $number_of_reads > $cap ){
			last;
		}
	}
}
=cut

sub calculate_total_length_of_matches {
    my $cigar = shift;
    my $total_length = 0;

    while ($cigar =~ /(\d+)M/g) {
        $total_length += $1;
    }

    return $total_length;
}

#system("mkdir $fastq_name.$organelle/$fastq_name_reads.tmp");
#system("$script_dir/kmc3/bin/kmc -k31 -cs999999 -ci3 -fa $fastq_name.$organelle/$fastq_name_reads $fastq_name.$organelle/$fastq_name_reads $fastq_name.$organelle/$fastq_name_reads.tmp");
#system("rm -rf $fastq_name.$organelle/$fastq_name_reads.tmp");


#if($threads >24 ){
#	system("$script_dir/readskmercount $fastq_name.$organelle/$fastq_name_reads $fastq_name.$organelle/$fastq_name_reads 24");
#}
#else{
#	system("$script_dir/readskmercount $fastq_name.$organelle/$fastq_name_reads $fastq_name.$organelle/$fastq_name_reads $threads");
#}

my%hash_assembly_round;
my$found=0;
for (my $round=1; $round<=$round_cutoff; $round++){

	#system("seqtk sample -s$round $fastq_name.$organelle/$fastq_name_reads.filter.fa $reads_num > $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa");
	
	if( $platform eq 'pacbio' ){
		#system("flye --pacbio-hifi $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa --threads $threads -o $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa.flye");
	}
	else{
		#system("flye --nano-hq $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa --threads $threads -o $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa.flye");
	}

	if( $organelle eq 'Chloroplast' ){
		if( $no_IR == 0 ){ # assumed that the chloroplast carries inverted repeats
	
			open my$chloro_assem,"$fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa" or die "Can't open $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa";
			my%hash_gfa;
			my%hash_edge_seq;
			my%hash_edge;
			while(<$chloro_assem>){
				chomp;
				if($_=~m/^L/){
					my@linkage=split /\t/,$_;
					$hash_gfa{$linkage[1]}+=1;
					$hash_gfa{$linkage[3]}+=1;
					push @{$hash_edge{$linkage[1]}},$linkage[3];
					push @{$hash_edge{$linkage[3]}},$linkage[1];
				}
				elsif($_=~m/^S/){
					my@S_edge=split /\t/,$_;
					$hash_edge_seq{$S_edge[1]}=$S_edge[2];
				}
			}

			my $number_edge = keys %hash_gfa;

			###Looking for IR edge
			my$IR;
			my$LSC;
			my$SSC;
			my$ring=0;
			foreach my$Edge (sort keys %hash_edge){
				if( scalar @{$hash_edge{$Edge}} == 4 ){ # IR edge
					my%hash;
					foreach my$item (@{$hash_edge{$Edge}}){
						$hash{$item} = 1;
					}
					my@unique_array = keys %hash;
					if( scalar @unique_array == 2 ){
						if( scalar @{$hash_edge{$unique_array[0]}} == 2 && scalar @{$hash_edge{$unique_array[1]}} == 2 ){
							print "Chloroplast genome with inverted repeats found\n";
							$found=1;
							$IR=$hash_edge_seq{$Edge};
							if( length$hash_edge_seq{$unique_array[0]} > length$hash_edge_seq{$unique_array[1]} ){
								$LSC=$hash_edge_seq{$unique_array[0]};
								$SSC=$hash_edge_seq{$unique_array[1]};
							}
							else{
								$LSC=$hash_edge_seq{$unique_array[1]};
								$SSC=$hash_edge_seq{$unique_array[0]};
							}
							my$SSC_RC=reverse_complement($SSC);
							my$IR_R=reverse_complement($IR);
							open my$heteroplasmy,">$fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.$Edge.$unique_array[0].$unique_array[1].$organelle.fa" or die "Can't open $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.$Edge.$unique_array[0].$unique_array[1].$organelle.fa";
							print $heteroplasmy ">Heteroplasmy1\n$LSC$IR$SSC$IR_R\n";
							print $heteroplasmy ">Heteroplasmy2\n$LSC$IR$SSC_RC$IR_R\n";
							system("ln -s $fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa $fastq_name.Chloroplast/$fastq_name_reads.filter.$reads_num.round$round.$organelle.gfa");
						}
					}
					elsif( scalar @unique_array == 1 ){
						if( keys %hash == 1 && $unique_array[0] eq $Edge ){
							$ring=1;
						}
					}

				}	
				elsif( scalar @{$hash_edge{$Edge}} == 2 ){
					my%hash;
					foreach my$item (@{$hash_edge{$Edge}}){
						$hash{$item} = 1;
					}
					my @unique_keys = keys %hash;
					if( keys %hash == 1 && $unique_keys[0] eq $Edge ){
						$ring=1;
					}	
				}
			}
		
			if( $found == 1 ){
				last;
			}
			elsif( $found == 0 ){ # No inverted repeats genome found
				if ($ring == 1) {
					$hash_assembly_round{$round}=$round;
					print "Chloroplast genome without inverted repeats found in round$round\n";
					system("ln -s $fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly.fasta $fastq_name.Chloroplast/$fastq_name_reads.filter.$reads_num.round$round.$organelle.fa");
				}
			}
		}
		elsif( $no_IR == 1 ){
			open my$chloro_assem,"$fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa" or die "Can't open $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa";
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
					system("ln -s $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly.fasta $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.$organelle.fa");
					$found=1;
					last;
				}
			}
		}
	}
	else{
		open my$chloro_assem,"$fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa" or die "Can't open $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly_graph.gfa";
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
			print "Mitochondrion genome found\n";
			system("ln -s $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fa.flye/assembly.fasta $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.$organelle.fa");
			$found=1;
			last;
		}
	}
}

if( $organelle eq 'Chloroplast' ){
	if ($found == 0 && keys %hash_assembly_round > 0) {
		print "After $round_cutoff rounds of downsampling and assembly, we couldn't find an assembly with an inverted repeat, but we found an assembly without an inverted repeat.\n";
		print "You can find the assembly in round ";
		foreach my $round_assem (sort keys %hash_assembly_round) {
			print "$round_assem ";
		}
		print "\n";
	}
	elsif ($found == 0 && keys %hash_assembly_round == 0) {
		print "After $round_cutoff rounds of downsampling and assembly, we couldn't find any assembly that meets our expectations.\n";
		print "You can increase the value of -n and try again :)\n";
	}
}

sub reverse_complement {
    my ($dna) = @_;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}
