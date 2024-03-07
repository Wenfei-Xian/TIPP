#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Long;

my$help;
my$fastq;
my$fastq_name;
my$threads=128;
my$current_directory = cwd();
my$reads_num=2000;
my$db;
my$round_cutoff=5;
my$no_IR=0;
my$platform='pacbio';
my$organelle='both';
my$cap=60000;
my$percentage=0.05;
my$mito_core;
my$version="v1.0.0";
#v1.0.0: March 07, 2024
my$show_version = 0;

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
    'p=s' => \$percentage,
    'b=s' => \$mito_core,
    'v' => \$show_version,
);

if ($show_version) {
    print "Version: $version\n";
    exit;
}

if( $platform ne 'pacbio' && $platform ne 'ont' ){
        print "Error, only pacbio or ont can be accecpted for -m\n";
        exit;
}

if( $organelle ne 'chloroplast' && $organelle ne 'both' ){
        print "Error, only Chloroplast or Both can be accecpted for -g\n";
        exit;
}

if ($help or not defined $fastq or not defined $db) {
    print "Usage: $0 [options]\n";
    print "  -h: Show this help message.\n";
    print "  -d: Chloroplast database (required).\n";
    print "  -f: HiFi reads (required).\n";
    print "  -g: chloroplast or both (default: both).\n";
    print "  -t: Threads for Minimap2, Flye, KMC3, Diamond and readskmercount.\n";
    print "  -n: Number of reads in each downsample, if the read length is short (<=15kb) or the genome size is extremely large, please increase this value (default: 2000 for chlo; mito will be 2*n).\n";
    print "  -r: Number of random downsamplings (default: 5).\n";
    print "  -m: Sequencing platform - either 'pacbio' or 'ont'. Only Q20 reads are accepted (default: pacbio).\n";
    print "  -i: Assume the presence of the inverted repeats (default: 1).\n";
    print "  -c: Maximum number of candidate reads used (default: 60000).\n";
    print "  -p: The proportion of total number of M in cigar / the length of reads, greater than this value is considered a match (default: 0.05, mito will be 2*p).\n";
    print "  -v: version.\n";
    exit;
}

for my $tool ("minimap2","samtools","flye","diamond") {
	my $return_val = system("which $tool > /dev/null 2>&1");
	if ($return_val != 0) {
		die "Error: $tool does not seem to be installed or is not in the PATH. Please check.\n";
	}
}

sub reverse_complement {
    my ($dna) = @_;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

my$script_dir = ($0 =~ m{(.*/)?})[0];
$mito_core="$script_dir/mito.485.functional.proteins.fasta.gz";

if( $organelle eq 'both' ){
	if ( -e "$mito_core" ) {
	}
	else{
        	print "Mito core genes does not exist!\n";
		print "Please check the path:$mito_core:)\n";
        	exit;
	}
}

$fastq_name=(split /\//,$fastq)[-1];
my$fastq_name_reads="$fastq_name.chloroplast.reads.fasta";

unless (-d "$fastq_name.$organelle") {
        system("mkdir $fastq_name.$organelle");
}

my%hash_uniq_seq;
my$number_of_chlo_reads=0;
my$number_of_mito_reads=0;

print "###TIPP_plastid for plastid start:\n";

if( $platform eq 'pacbio' ){
	open IN1,"minimap2 --secondary no --sam-hit-only -ax map-hifi -t $threads $db $fastq |";
}
else{
	open IN1,"minimap2 --secondary no --sam-hit-only -ax map-ont -t $threads $db $fastq |";
}

my%chloroplast_map_hash;
my%mito_minimap_hash;
open OUT1,">$fastq_name.$organelle/$fastq_name_reads";
open OUT2,">$fastq_name.$organelle/$fastq_name_reads.sam";
open OUT3,">$fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta";
open OUT4,">$fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta.sam";
while(<IN1>){
	chomp;
	if(/^@/){
		print OUT2 "$_\n";
	}
	else{
		my@sam=split /\t/,$_;
		next if exists $hash_uniq_seq{$sam[0]};
		$hash_uniq_seq{$sam[0]} = undef;
		my$total_matches = calculate_total_length_of_matches($sam[5]);
		my$query_len=length$sam[9];
		if( $sam[2]=~m/Chloroplast/ && $total_matches >= $query_len*$percentage ){
			$number_of_chlo_reads++;
			if( $number_of_chlo_reads < $cap ){
				print OUT1 ">$sam[0]\n$sam[9]\n";
				print OUT2 "$_\n";
				$chloroplast_map_hash{$sam[0]}=$sam[0];
			}
			#if( $number_of_reads > $cap ){
			#	last;
			#}
		}
		elsif( $sam[2]=~m/Mitochondrion/ && $total_matches >= $query_len*$percentage*2 && $organelle eq 'both' ){
			$number_of_mito_reads++;
			if( $number_of_mito_reads < $cap ){
				print OUT3 ">$sam[0]\n$sam[9]\n";
				print OUT4 "$_\n";
				$mito_minimap_hash{$sam[0]}=$sam[0];
			}
		}
		if($number_of_chlo_reads > $cap && $number_of_mito_reads > $cap){
			last;
		}
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

system("mkdir $fastq_name.$organelle/$fastq_name_reads.tmp");
system("$script_dir/kmc3/bin/kmc -t$threads -k31 -cs999999 -ci3 -fa $fastq_name.$organelle/$fastq_name_reads $fastq_name.$organelle/$fastq_name_reads $fastq_name.$organelle/$fastq_name_reads.tmp");
system("rm -rf $fastq_name.$organelle/$fastq_name_reads.tmp");


if($threads >12 ){
	system("$script_dir/readskmercount $fastq_name.$organelle/$fastq_name_reads $fastq_name.$organelle/$fastq_name_reads 12");
}
else{
	system("$script_dir/readskmercount $fastq_name.$organelle/$fastq_name_reads $fastq_name.$organelle/$fastq_name_reads $threads");
}

my%hash_assembly_round;
my$found=0;
for (my $round=1; $round<=$round_cutoff; $round++){

	system("$script_dir/seqtk/seqtk sample -s$round $fastq_name.$organelle/$fastq_name_reads.filter.fasta $reads_num > $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta");
	
	if( $platform eq 'pacbio' ){
		system("flye --pacbio-hifi $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta --threads $threads -o $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta.chloroplast.flye");
	}
	elsif( $platform eq 'ont' ){
		system("flye --nano-hq $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta --threads $threads -o $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta.chloroplast.flye");
	}

	if( $no_IR == 0 ){ # assumed that the chloroplast carries inverted repeats
	
		open my$chloro_assem,"$fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta.chloroplast.flye/assembly_graph.gfa" or die "Can't open $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta.chloroplast.flye/assembly_graph.gfa";
		my%hash_gfa;
		my%hash_edge_seq;
		my%hash_edge;
		my$S=0;
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
				$S++;
			}
		}

		my $number_edge = keys %hash_gfa;
		next unless ($S == 3 || $S == 1);
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
						open my$heteroplasmy,">$fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.$Edge.$unique_array[0].$unique_array[1].$organelle.chloroplast.fasta" or die "Can't open $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.$Edge.$unique_array[0].$unique_array[1].$organelle.chloroplast.fasta";
						print $heteroplasmy ">Heteroplasmy1\n$LSC$IR$SSC$IR_R\n";
						print $heteroplasmy ">Heteroplasmy2\n$LSC$IR$SSC_RC$IR_R\n";
						system("ln -s $fastq_name_reads.filter.$reads_num.round$round.fasta.chloroplast.flye/assembly_graph.gfa $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.$organelle.chloroplast.gfa");
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
				system("ln -s $fastq_name_reads.filter.$reads_num.round$round.fasta.chloroplast.flye/assembly.fasta $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.$organelle.chloroplast.fasta");
			}
		}
	}
	elsif( $no_IR == 1 ){
		open my$chloro_assem,"$fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta.chloroplast.flye/assembly_graph.gfa" or die "Can't open $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta.chloroplast.flye/assembly_graph.gfa";
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
				system("ln -s $fastq_name_reads.filter.$reads_num.round$round.fasta.chloroplast.flye/assembly.fasta $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.$organelle.chloroplast.fasta");
				$found=1;
				last;
			}
		}
	}
}

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

print "###TIPP_plastid for plastid end.\n\n\n";

if( $organelle eq 'both' ){
	
	my $command = "grep \">\" $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta | wc -l";
	my $count = `$command`;
	chomp $count;
	
	if( $count > 0 ){
		print "###TIPP_plastid for Mitochondrion start:\n";
		print "###Starting to process the reads retrieved from minimap2\n";

		system("mkdir $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta.tmp");
		system("$script_dir/kmc3/bin/kmc -t$threads -k31 -cs999999 -ci3 -fa $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta $fastq_name.$organelle/$fastq.mitochondrial.reads.minimap2.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta.tmp");
		system("rm -rf $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta.tmp");

		if($threads >12 ){
			system("$script_dir/readskmercount $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta 12");
		}
		else{
			system("$script_dir/readskmercount $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta $threads");
		}

		system("ln -s $fastq_name.mitochondrial.reads.minimap2.fasta.filter.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta");
		$reads_num=$reads_num*3;
		for (my $round=1; $round<=3; $round++){

			system("$script_dir/seqtk/seqtk sample -s$round $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta $reads_num > $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta.$reads_num.round$round.fasta");

			if( $platform eq 'pacbio' ){
				system("flye --pacbio-hifi $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta.$reads_num.round$round.fasta --threads $threads -o $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta.$reads_num.round$round.fasta.mitochondrial.flye");
			}
			elsif( $platform eq 'ont' ){
				system("flye --nano-hq $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta.$reads_num.round$round.fasta --threads $threads -o $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta.$reads_num.round$round.fasta.mitochondrial.flye");
			}
		}
	}
	else{

		my%hash_mito_minimap2_ID;
		#open my$mito_minimap2,"$fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta.filter.fasta" or die "Can't open $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta.filter.fasta";
		#$/=">";<$mito_minimap2>;
		#while(<$mito_minimap2>){
		#	my($mito_minimap2_ID,$mito_minimap2_seq)=(split /\n/,$_,2)[0,1];
		#	$hash_mito_minimap2_ID{$mito_minimap2_ID}=$mito_minimap2_ID;
		#}
		#$/="\n";

		my%hash_mito_minimap_false_positive_ID;
		#open my$mito_minimap2_raw,"$fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta" or die "Can't open $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta";
		#$/=">";<$mito_minimap2_raw>;
		#while(<$mito_minimap2_raw>){
		#	my($mito_minimap2_raw_ID,$mito_minimap2_raw_seq)=(split /\n/,$_,2)[0,1];
		#	next if exists $hash_mito_minimap2_ID{$mito_minimap2_raw_ID};
		#	$hash_mito_minimap_false_positive_ID{$mito_minimap2_raw_ID}=$mito_minimap2_raw_ID;
		#}
		#$/="\n";

		###
		#$fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta.filter.fasta
		###
	
		print "###Since mitochondrial genome are diverse, but genes are conserved\n";
		print "###Following stepwill use protein sequences to recruit the mitochondrial reads\n";
	
		open my $diamond_out,"diamond blastx --threads $threads --evalue 1e-5 --db $mito_core --query $fastq --max-target-seqs 10 --outfmt 6 qseqid full_qseq|";
		open my $mito_out, ">$fastq_name.$organelle/$fastq_name.mitochondrial.reads.diamond.fasta";
		my%hash_diamond_ID;
		my%hash_diamond_hits;
		my$number_of_reads_mito=0;
		while(<$diamond_out>){
			chomp;
			if (/(\S+)\s+(\S+)/) {
				my ($diamond_out_id, $diamond_out_seq) = ($1, $2);
				$hash_diamond_ID{$diamond_out_id}=$diamond_out_id;
				$hash_diamond_hits{$diamond_out_id}++;
				next if exists $hash_mito_minimap2_ID{$diamond_out_id};#de-redundancy
				next if exists $hash_mito_minimap_false_positive_ID{$diamond_out_id};
				next if exists $chloroplast_map_hash{$diamond_out_id};#remove the reads origin from chloroplast
				if( $hash_diamond_hits{$diamond_out_id} == 10 ){
					print $mito_out ">$diamond_out_id\n$diamond_out_seq\n";
				}
				$number_of_reads_mito++;
			}
		}
		close $diamond_out;

		system("mkdir $fastq_name.$organelle/$fastq_name.mitochondrial.reads.diamond.tmp");
		system("$script_dir/kmc3/bin/kmc -t$threads -k31 -cs999999 -ci3 -fa $fastq_name.$organelle/$fastq_name.mitochondrial.reads.diamond.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.reads.diamond.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.reads.diamond.tmp");
		system("rm -rf $fastq_name.$organelle/$fastq_name.mitochondrial.reads.diamond.tmp");

		if($threads >12 ){
			system("$script_dir/readskmercount $fastq_name.$organelle/$fastq_name.mitochondrial.reads.diamond.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.reads.diamond.fasta 12");
		}
		else{
			system("$script_dir/readskmercount $fastq_name.$organelle/$fastq_name.mitochondrial.reads.diamond.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.reads.diamond.fasta $threads");
		}

		#merge the reads using minimap2 and diamond
		system("cat $fastq_name.$organelle/$fastq_name.mitochondrial.reads.diamond.fasta.filter.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2.fasta > $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2diamond.fasta");

		my%hash_multiple_extension;
		for(my$extend=0;$extend<=5;$extend++){#extend, strict thresold
		
			my$new_ID=0;
			if( $platform eq 'pacbio' ){
				open IN1_seed,"minimap2 --secondary no --sam-hit-only -ax map-hifi --MD --eqx -t $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2diamond.fasta $fastq |";
			}
			else{
				open IN1_seed,"minimap2 --secondary no --sam-hit-only -ax map-ont --MD --eqx -t $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2diamond.fasta $fastq |";
			}

			while(<IN1_seed>){
				chomp;
				if(/^@/){
					next;
				}
				else{
					my@sam=split /\t/,$_;
					next if exists $chloroplast_map_hash{$sam[0]};#skip the reads belong to chloroplast;
					next if exists $hash_uniq_seq{$sam[0]}; #only process the best hit
					next if exists $hash_diamond_ID{$sam[0]}; #skip the reads already present in the diamond step;
					next if exists $hash_mito_minimap2_ID{$sam[0]}; #skip the reads already present the minimap2 step;
					next if exists $hash_mito_minimap_false_positive_ID{$sam[0]};#skip the reads already marked as false positive
					$hash_uniq_seq{$sam[0]} = undef;
					my ($matches, $mismatches, $insertions, $deletions, $right_clip, $left_clip) = analyze_cigar($sam[5]);
					next if( $right_clip > 100 && $left_clip > 100 ); #no double clip
					next if( $mismatches >= $matches*0.0001 ); #no large snps
					if( $matches >=5000 ){
						$hash_multiple_extension{$sam[0]}=$sam[9];
						$new_ID++;
					}
				}
			}
			print "Extension $extend: there are $new_ID added\n";
			if($new_ID<100){
				print "In extension $extend, new reads are smaller than 100, extension finished\n";
				last;
			}
		}

		open OUT_EXTENSION,">$fastq_name.$organelle/$fastq_name.mitochondrial.reads.extension.fasta";
		foreach my$reads_extension (keys %hash_multiple_extension){
			print OUT_EXTENSION  ">$reads_extension\n$hash_multiple_extension{$reads_extension}\n";
		}
	
		system("cat $fastq_name.$organelle/$fastq_name.mitochondrial.reads.minimap2diamond.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.reads.filter.fasta.extension.fasta > $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta");
		system("ln -s $fastq_name.mitochondrial.reads.minimap2.fasta.filter.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta");
		$reads_num=$reads_num*3;
		for (my $round=1; $round<=3; $round++){
		
			system("$script_dir/seqtk/seqtk sample -s$round $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta $reads_num > $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta.$reads_num.round$round.fasta");

			if( $platform eq 'pacbio' ){
				system("flye --pacbio-hifi $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta.$reads_num.round$round.fasta --threads $threads -o $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta.$reads_num.round$round.fasta.mitochondrial.flye");
			}
        		elsif( $platform eq 'ont' ){
				system("flye --nano-hq $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta.$reads_num.round$round.fasta --threads $threads -o $fastq_name.$organelle/$fastq_name.mitochondrial.reads.merged.fasta.$reads_num.round$round.fasta.mitochondrial.flye");
			}
		}
	}
}

sub analyze_cigar {
	my ($cigar) = @_;
	my $matches = 0;
	my $mismatches = 0;
	my $insertions = 0;
	my $deletions = 0;
	my $right_clip = 0;
	my $left_clip = 0;

	while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
		my $len = $1;
		my $op = $2;

		if ($op eq 'M') {
			$matches += $len;
		}
		elsif ($op eq '=') {
			$matches += $len;
		}
		elsif ($op eq 'X') {
			$mismatches += $len;
		}
		elsif ($op eq 'I') {
			$insertions += $len;
		}
		elsif ($op eq 'D') {
			$deletions += $len;
		}
		elsif ($op eq 'S') {
			if ($cigar =~ /^$len$op/) {
				$left_clip += $len;
			}
			else {
				$right_clip += $len;
			}
		}
	}
	return ($matches, $mismatches, $insertions, $deletions, $right_clip, $left_clip);
}
