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
my$reads_num=800;
my$round_cutoff=5;
my$no_IR=0;
#my$platform='pacbio';
my$organelle='organelle';
my$version="v2.3";
my$lkc=0.3;
my$hkc=5;
my$flyep="--pacbio-hifi";
my$minimap2p="map-hifi";
my$minoverlap=800;
my$reference;
my$trf;
#v1.0: March 07, 2024
#v1.1: May 09, 2024; 
#v2.0: May 20, 2024; reference free
#v2.1: May 26, 2024; exclude ultrl high frequency reads
#v2.2: Nov 11, 2024; accept the file with the absoult path, reconstructe the repeat graph with smaller overlap length;
#v2.3: Dec 03, 2024; intergrated the reference based approach 
#v2.4: Dec 04, 2024; excluded the tandem repeats reads (HiFi only)
my$show_version = 0;

GetOptions(
    'h'   => \$help,
    'f=s' => \$fastq,
    't=s' => \$threads,
    'n=s' => \$reads_num,
    'i=s' => \$no_IR,
    #'p=s' => \$platform,
    'r=s' => \$round_cutoff,
    'g=s' => \$organelle,
    'v' => \$show_version,
    'l=s' => \$lkc,
    'c=s' => \$hkc,
    'y=s' => \$flyep,
    'a=s' => \$minimap2p,
    'm=s' => \$minoverlap,
    'b=s' => \$reference,
    'trf' => \$trf,
);

if ($show_version) {
    print "Version: $version\n";
    exit;
}

#if( $platform ne 'pacbio' && $platform ne 'ont' ){
#        print "Error, only pacbio or ont can be accecpted for -m\n";
#        exit;
#}

if( $organelle ne 'chloroplast' && $organelle ne 'organelle' ){
        print "Error, only chloroplast or organelle can be accecpted for -g\n";
        exit;
}

if ($help or not defined $fastq ) {
	print "Usage: $0 [options]\n";
	print "  -h: Show this help message.\n";
	print "  -f: HiFi reads (required).\n";
	print "  -g: chloroplast or organelle (default: organelle).\n";
	print "  -t: Threads for tiara, flye, KMC3 and readskmercount.\n";
	print "  -n: Number of reads in each downsample for chloroplast.\n";
	print "  -r: Number of random downsamplings (default: 5).\n";
	#print "  -p: Sequencing platform - either 'pacbio' or 'ont'. Only Q20 reads are accepted (default: pacbio).\n";
	print "  -i: Assume the presence of the inverted repeats (default: 1).\n";
	print "  -l: lower kmer count - lkc (default: 0.3).\n";
	print "  -c: high kmer count - hkc (default: 5).\n";
	print "  -y: parameter for flye (default: --pacbio-hifi).\n";
	print "  -a: parameter for minimap2 (default: map-hifi).\n";
	print "  -m: minimum overlap in repeat graph construction (default:800)\n";
	print "  -b: reference sequence (default: No).\n";
	print "  --trf: remove the reads are tandem repeats (default: used for HiFi reads)\n";
	print "  -v: version.\n";
	exit;
}

for my $tool ("flye","tiara","GraphAligner","trf") {
	my $return_val = system("which $tool > /dev/null 2>&1");
	if ($return_val != 0) {
		die "Error: $tool does not seem to be installed or is not in the PATH. Please check.\n";
	}
}

#my$flyep=join " ",@flyep_array;
#my$minimap2p=join " ",@minimap2_array;

my$script_dir = ($0 =~ m{(.*/)?})[0];

$fastq_name=(split /\//,$fastq)[-1];
my$fastq_name_reads="$fastq_name.chloroplast.fasta";

unless (-d "$fastq_name.$organelle") {
        system("mkdir $fastq_name.$organelle");
}

print "###TIPPo for plastid start:\n";

my$input=$fastq;
if( $fastq=~m/\.fastq.gz$/ || $fastq=~m/\.fq.gz$/ || $fastq=~m/\.fastq$/ || $fastq=~m/\.fq$/){
	system("$script_dir/seqtk/seqtk seq -A $fastq | awk '{print \$1}' > $fastq.fasta");
	$input="$fastq.fasta";
}
elsif( $fastq=~m/\.fasta.gz$/ || $fastq=~m/\.fa.gz$/ || $fastq=~m/\.fasta$/ || $fastq=~m/\.fa$/ ){
	system("$script_dir/seqtk/seqtk seq -A $fastq | awk '{print \$1}' > $fastq.fasta");
	$input="$fastq.fasta";
}

my$input_name=(split /\//,$input)[-1];

if( defined $reference ){
	
	my%hash_uniq_seq;
	my%chloroplast_map_hash;
	my%mito_minimap_hash;

	print "minimap2 --secondary no --sam-hit-only -ax $minimap2p -t $threads $reference $input |";
	open ALIGN,"minimap2 --secondary no --sam-hit-only -ax $minimap2p -t $threads $reference $input |";
	open OUT1_CHLO,">$fastq_name.$organelle/$fastq_name_reads.filter.fasta";
	open OUT2_CHLO,">$fastq_name.$organelle/$fastq_name_reads.filter.fasta.sam";
	open OUT3_MITO,">$fastq_name.$organelle/$fastq_name.mitochondrial_tag.fasta"; #$fastq_name.mitochondrial_tag.fasta
	open OUT4_MITO,">$fastq_name.$organelle/$fastq_name.mitochondrial_tag.fasta.sam";
	open OUT5_MITO,">$fastq_name.$organelle/$fastq_name.mitochondrial_tag.ID";
	while(<ALIGN>){
		chomp;
		if(/^@/){
			#print OUT2 "$_\n";
		}
		else{
			my@sam=split /\t/,$_;
			next if exists $hash_uniq_seq{$sam[0]};
			$hash_uniq_seq{$sam[0]} = undef;
			my$total_matches = calculate_total_length_of_matches($sam[5]);
			my$query_len=length$sam[9];
			if( $sam[2]=~m/Chloroplast/ && $total_matches >= $query_len*0.3 ){
				next if( exists $chloroplast_map_hash{$sam[0]} );
				print OUT1_CHLO ">$sam[0]\n$sam[9]\n";
				print OUT2_CHLO "$_\n";
				$chloroplast_map_hash{$sam[0]}=$sam[0];
			}
			elsif( $sam[2]=~m/Mitochondrion/ && $total_matches >= $query_len*0.3 && $organelle eq 'organelle' ){
				next if( exists $mito_minimap_hash{$sam[0]} );
				print OUT3_MITO ">$sam[0]\n$sam[9]\n";
				print OUT4_MITO "$_\n";
				print OUT5_MITO "$sam[0]\n";
				$mito_minimap_hash{$sam[0]}=$sam[0];
			}
		}
	}
}
else{

	system("mkdir $fastq_name.$organelle/$fastq_name.tmp");
	#system("$script_dir/kmc3/bin/kmc -t$threads -k31 -cs999999 -ci2 -fa $input $fastq_name.$organelle/$fastq_name $fastq_name.$organelle/$fastq_name.tmp");
	system("kmc -t$threads -k31 -cs999999 -ci2 -fa $input $fastq_name.$organelle/$fastq_name $fastq_name.$organelle/$fastq_name.tmp");
 	system("rm -rf $fastq_name.$organelle/$fastq_name.tmp");

	system("tiara -i $input -o $fastq_name.$organelle/$fastq_name.tiara.out --to_fasta mit pla  -t $threads");
	#system("cat $fastq_name.$organelle/mitochondrion_$input $fastq_name.$organelle/plastid_$input > $fastq_name.$organelle/$fastq_name.organelle.fasta");
	system("grep mitochondrion $fastq_name.$organelle/$fastq_name.tiara.out | awk '{print \$1}' > $fastq_name.$organelle/$fastq_name.mitochondrial_tag.ID");
	system("mv $fastq_name.$organelle/plastid_$input_name $fastq_name.$organelle/$fastq_name_reads");
	system("mv $fastq_name.$organelle/mitochondrion_$input_name $fastq_name.$organelle/$fastq_name.mitochondrial_tag.fasta");

	system("$script_dir/readskmercount $fastq_name.$organelle/$fastq_name $fastq_name.$organelle/$fastq_name_reads $threads $lkc $hkc");

}

my%hash_assembly_round;
my$found=0;
my$found_round=0;
for (my $round=1; $round<=$round_cutoff; $round++){

	system("$script_dir/seqtk/seqtk sample -s$round $fastq_name.$organelle/$fastq_name_reads.filter.fasta $reads_num > $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta");
	
	system("flye $flyep $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta --threads $threads -o $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round.fasta.chloroplast.flye");

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
			$found_round=$round;
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

if ( $found == 1 ){
	system( "GraphAligner -g $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$found_round.$organelle.chloroplast.gfa -f $fastq_name.$organelle/$fastq_name.mitochondrial_tag.fasta -x vg -a $fastq_name.$organelle/$fastq_name_reads.graphaligner.chloroplast.gfa.gaf --precise-clipping 0.9 -t $threads");
	
	open GA,"$fastq_name.$organelle/$fastq_name_reads.graphaligner.chloroplast.gfa.gaf" or die "Can't open $fastq_name.$organelle/$fastq_name_reads.graphaligner.chloroplast.gfa.gaf\n";
        open GAOUT,">$fastq_name.$organelle/$fastq_name_reads.mito_chloroplast.ID";
        my%hash_GA_paf;
        while(<GA>){
                chomp;
                my@GA_paf=split /\t/,$_;
                next if( exists $hash_GA_paf{$GA_paf[0]} );
                if( $GA_paf[2] <= 100 && $GA_paf[3]+100 >= $GA_paf[1] && $GA_paf[9] > $GA_paf[10]*0.95 ){
                        $hash_GA_paf{$GA_paf[0]}=0;
                        print GAOUT "$GA_paf[0]\n";
                }
        }
        print "\n";

}
elsif ($found == 0 && keys %hash_assembly_round > 0) {
	print "After $round_cutoff rounds of downsampling and assembly, we couldn't find an assembly with an inverted repeat, but we found an assembly without an inverted repeat.\n";
	print "You can find the assembly in round ";
	my$round_found0;
	foreach my $round_assem (sort keys %hash_assembly_round) {
		print "$round_assem ";
		$round_found0=$round_assem;
	}
	system( "GraphAligner -g $fastq_name.$organelle/$fastq_name_reads.filter.$reads_num.round$round_found0.fasta.chloroplast.flye/assembly_graph.gfa -f $fastq_name.$organelle/$fastq_name.mitochondrial_tag.fasta -x vg -a $fastq_name.$organelle/$fastq_name_reads.graphaligner.chloroplast.gfa.gaf --precise-clipping 0.9 -t $threads");
        
	print "###Start to filter out the GraphAligner result\n";
	
	open GA,"$fastq_name.$organelle/$fastq_name_reads.graphaligner.chloroplast.gfa.gaf" or die "Can't open $fastq_name.$organelle/$fastq_name_reads.graphaligner.chloroplast.gfa.gaf\n";
	open GAOUT,">$fastq_name.$organelle/$fastq_name_reads.mito_chloroplast.ID";
	my%hash_GA_paf;
	while(<GA>){
		chomp;
		my@GA_paf=split /\t/,$_;
		next if( exists $hash_GA_paf{$GA_paf[0]} );
		if( $GA_paf[2] <= 100 && $GA_paf[3]+100 >= $GA_paf[1] && $GA_paf[9] > $GA_paf[10]*0.95 ){
			$hash_GA_paf{$GA_paf[0]}=0;
			print GAOUT "$GA_paf[0]\n";
		}
	}
	print "\n";
}
elsif ($found == 0 && keys %hash_assembly_round == 0) {
	print "After $round_cutoff rounds of downsampling and assembly, we couldn't find any assembly that meets our expectations.\n";
	print "You can increase the value of -n and try again :)\n";
	system("touch $fastq_name.$organelle/$fastq_name_reads.mito_chloroplast.ID");
}

print "###TIPPo for plastid end.\n\n\n";

if( $organelle eq 'organelle' ){
	
	print "###TIPPo for Mitochondrion start:\n";

	system("sort $fastq_name.$organelle/$fastq_name.mitochondrial_tag.ID $fastq_name.$organelle/$fastq_name_reads.mito_chloroplast.ID $fastq_name.$organelle/$fastq_name_reads.mito_chloroplast.ID | uniq -u > $fastq_name.$organelle/$fastq_name.mitochondrial.ID");
	system("$script_dir/seqtk/seqtk subseq $fastq_name.$organelle/$fastq_name.mitochondrial_tag.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.ID > $fastq_name.$organelle/$fastq_name.mitochondrial.fasta");

	if( defined $reference ){
		system("cp $fastq_name.$organelle/$fastq_name.mitochondrial.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.filter.fasta");
	}	
	else{
		if($trf) {
			system("sort $fastq_name.$organelle/$fastq_name.mitochondrial_tag.ID $fastq_name.$organelle/$fastq_name_reads.mito_chloroplast.ID $fastq_name.$organelle/$fastq_name_reads.mito_chloroplast.ID | uniq -u > $fastq_name.$organelle/$fastq_name.mitochondrial.ID");
			system("$script_dir/seqtk/seqtk subseq $fastq_name.$organelle/$fastq_name.mitochondrial_tag.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.ID > $fastq_name.$organelle/$fastq_name.mitochondrial.fasta");

			system("seqtk split -n $threads $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.split $fastq_name.$organelle/$fastq_name.mitochondrial.fasta");
			my @split_files = glob("$fastq_name.$organelle/$fastq_name.mitochondrial.fasta.split*.fa");
			my $trf_cmd = join("\n", map { "trf $_ 2 7 7 80 10 100 500 -h -l 6 -ngs > $_.trf.out" } @split_files);
			system("echo \"$trf_cmd\" | parallel -j $threads");
			system("cat $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.split*.trf.out > $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.trf.out");
			system("rm $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.split*.fa");
			system("rm $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.split*.trf.out");
			#system("trf $fastq_name.$organelle/$fastq_name.mitochondrial.fasta 2 7 7 80 10 100 500 -h -l 6 -ngs > $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.trf.out");
	
			open MITO,"$fastq_name.$organelle/$fastq_name.mitochondrial.fasta";
			$/=">";<MITO>;
			my%hash_mito_len;
			my%hash_mito_seq;
			while(<MITO>){
				chomp;
				my($mito_ID,$mito_seq)=(split /\n/,$_,2)[0,1];
				$mito_seq=~s/>//;
				$mito_seq=~s/\n//g;
				my$mito_len=length$mito_seq;
				$hash_mito_len{$mito_ID}=$mito_len;
				$hash_mito_seq{$mito_ID}=$mito_seq;
			}
			$/="\n";

			open TRF,"$fastq_name.$organelle/$fastq_name.mitochondrial.fasta.trf.out";
			my$chr="";
			my%hash_mito_satellite;
			while(<TRF>){
				chomp;
				if(/@/){
					$chr=$_;
					$chr=~s/@//;
				}
				else{
					my@trf_out=split /\s/,$_;
					if( $trf_out[1]-$trf_out[0] > $hash_mito_len{$chr}*0.5 && $trf_out[2] > 50 && $trf_out[2] < 300 ){
						$hash_mito_satellite{$chr}="";
					}
				}
			}

			open MITOOUT,">$fastq_name.$organelle/$fastq_name.mitochondrial.fasta.trf.fasta";
			foreach my$mito_reads (keys %hash_mito_seq){
				if( exists $hash_mito_satellite{$mito_reads} ){
				}
				else{
					print MITOOUT ">$mito_reads\n$hash_mito_seq{$mito_reads}\n";
				}
			}
			
			system("$script_dir/readskmercount $fastq_name.$organelle/$fastq_name $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.trf.fasta $threads $lkc $hkc");
			system("mv $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.trf.fasta $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.filter.fasta");

		}else{
	
			system("$script_dir/readskmercount $fastq_name.$organelle/$fastq_name $fastq_name.$organelle/$fastq_name.mitochondrial.fasta $threads $lkc $hkc");
		}
	}
	#for (my $round=1; $round<=$round_cutoff_mito; $round++){
	#system("$script_dir/seqtk/seqtk sample -s$round $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.filter.fasta $reads_num_m > $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.filter.$reads_num_m.round$round.fasta");

	system("flye $flyep $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.filter.fasta --threads $threads -o $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.filter.fasta.flye");
	
	system("mkdir $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.filter.fasta.flye/50.repeat-graph");
	system("flye-modules repeat --disjointigs $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.filter.fasta.flye/10-consensus/consensus.fasta --reads $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.filter.fasta --out-dir $fastq_name.$organelle/$fastq_name.mitochondrial.fasta.filter.fasta.flye/50.repeat-graph --config $script_dir/asm_hifi.cfg --min-ovlp $minoverlap --threads $threads");
}

sub reverse_complement {
    my ($dna) = @_;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

sub calculate_total_length_of_matches {
        my $cigar = shift;
        my $total_length = 0;

        while ($cigar =~ /(\d+)M/g) {
                $total_length += $1;
        }

        return $total_length;
}
