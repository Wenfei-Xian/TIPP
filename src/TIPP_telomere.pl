#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Long;

##########
#Although the Perl language has been abandoned by many bioinformaticians, it remains a good language for pipeline construction.
#Below, I have added detailed comments, hoping that users unfamiliar with Perl can also understand what each step is doing :)
#Author: Wenfei Xian
##########

my$help;
my$vcf_file;
my$unit;
my$fastq;
my$contigs;
my$bam;
my$threads=20;
my$vcf_file_name;
my$current_directory = cwd();
my$contigs_name;
my$fastq_name;
my$extend=0;
my$minimum_uniq=2000;
my$threads_spoa=2;
my$two_telomere=1;

GetOptions(
    'h'   => \$help,
    'd=s' => \$vcf_file,
    'u=s' => \$unit,
    'f=s' => \$fastq,
    'c=s' => \$contigs,
    'b=s' => \$bam,
    't=s' => \$threads,
    'e=s' => \$extend,
    'm=s' => \$minimum_uniq,
    'p=s' => \$threads_spoa,
);

###Help message
if ($help or not defined $vcf_file and not defined $unit and not defined $fastq and not defined $contigs and not defined $bam) {
    print "Usage: $0\n";
    print "-h: show this help message.\n";
    print "-u: telomere unit.\n";
    print "-f: hifi reads.\n";
    print "-e: extend the contigs with new assembled telomere sequences.(default=0, no extend)\n";
    print "-c: contigs. If the extension is not specified, it will be used as the output name.(default=TIPP_telo)\n";
    print "-t: threads for minimap2.\n";
    print "-m: minimum length of uniq sequence (without telomere)\n";
    exit;
}

###Check dependencies
for my $tool ("spoa", "seqtk", "minimap2", "mcl", "samtools","trf") {
	my $return_val = system("which $tool > /dev/null 2>&1");
	if ($return_val != 0) {
		die "Error: $tool does not seem to be installed or is not in the PATH. Please check.\n";
	}
}

###Extract the name of the contig file, ensuring that the generated file path is not confused by the directory structure.
if( defined $contigs ){
	$contigs_name=(split /\//,$contigs)[-1];
}
else{
	$contigs_name="TIPP_telo";
}
$fastq_name=(split /\//,$fastq)[-1];

###Determine the path where the script resides to ensure related scripts are found correctly (like the Rscripts). 
my$script_dir = ($0 =~ m{(.*/)?})[0];

###Start
print "Extracting telomere reads, classifying into clusters, and generating a consensus for each cluster\n";

###Create the output directory.
unless (-d "$contigs_name.telomere") {
        system("mkdir $contigs_name.telomere");
}

###Use seqtk to extract HiFi reads with telomeres unit.
system("seqtk telo -m $unit $fastq  > $contigs_name.telomere/$fastq_name.telo.tmp.fa");

if( $two_telomere == 1 ){
	system("grep '>' $contigs_name.telomere/$fastq_name.telo.tmp.fa | sed 's/>//' | awk '{print \$1}'| sort | uniq -c |awk '{if(\$1 == 1)print \$2}' > $contigs_name.telomere/$fastq_name.telo.fa.single.ID");
	system("seqtk subseq $contigs_name.telomere/$fastq_name.telo.tmp.fa $contigs_name.telomere/$fastq_name.telo.fa.single.ID > $contigs_name.telomere/$fastq_name.telo.fa");
}
else{
	system("cp $contigs_name.telomere/$fastq_name.telo.tmp.fa $contigs_name.telomere/$fastq_name.telo.fa");
}

###Process the seqtk output to extract and output sequences with telomeres removed.
system("grep -A 1 'remove' --no-group-separator $contigs_name.telomere/$fastq_name.telo.fa > $contigs_name.telomere/$fastq_name.remove.telo.fa");

###Employ minimap2 to conduct an all-vs-all alignment of sequences with telomeres removed.
system("minimap2 --eqx -c -x asm20 -X -t $threads $contigs_name.telomere/$fastq_name.remove.telo.fa $contigs_name.telomere/$fastq_name.remove.telo.fa -o $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf");

###Extract potential read pairs stemming from the same telomere based on alignment results, treating each read as a node and each pair as an edge, for use as input to the MCL algorithm.
open my$paf,"$contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf" or die "Can't open $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf\n";
open my$abc,">$contigs_name.telomere/${fastq_name}.remove.telo.fa.self.paf.abc" or die "Can't open $contigs_name.telomere/${fastq_name}.remove.telo.fa.self.paf.abc\n";
###In theory, if the two hifi reads origin from the same locus, we should detect more Indels, instead of more mismatches.
my%paf_unit;
my%paf_pair;
while(<$paf>){
	chomp;
	my@array=split /\t/,$_;
	my$num_I=0;
	my$num_D=0;
	my$num_X=0;

	foreach my$col (@array){
		if($col=~m/cg:/){
			while ($col =~ /(\d+)([IDX])/g) {
				my $count = $1;
				my $type = $2;
				if ($type eq 'I') {
					$num_I+=$count;
				}
				elsif ($type eq 'D') {
					$num_D+=$count;
				}
				elsif ($type eq 'X') {
					$num_X+=$count;
				}
			}
		}
	}

	if( !exists $paf_pair{"$array[0]\t$array[5]"} ){
		$paf_pair{"$array[0]\t$array[5]"}=0;
		if( $array[9]/$array[10] > 0.98  && $array[1] >= $minimum_uniq && $array[6] >= $minimum_uniq && $array[0] ne $array[5] && $num_X < ($num_I+$num_D) ){
			#MMMMMMMXMMMMMMMMMMMMMMMMM Q20
			#MMMMMMMMMMMMXMMMMMMMMMMMM Q20
			#0.99 * 0.99 = 0.98
			if( $array[1] > $array[6] ){
				#if( $array[2] <= 500 || $array[3] >= $array[1]-500 ){
					if( $array[8]-$array[7] >= $array[6]*0.9 ){
						if( ($array[2] == 0 && $array[7] == 0 && $array[8]+200 < $array[6]) || ( $array[1] == $array[3] && $array[6] == $array[8] && $array[7] > 200 ) || ( ($array[7] > 200) && ($array[8]+200 > $array[6] )) ){
							#if the length of X, unmap sequence longer than the cutoff, two sequences will be considered to be origined from two telomeres
							#
							# MMMMMMMMMMMMMMMMMMMMMXXXXXX
							# MMMMMMMMMMMMMMMMMMMMMXXXXXXX
							#
							# XXXXMMMMMMMMMMMMMMMMMMM
							# XXXXMMMMMMMMMMMMMMMMMMM
							#
							# XXXXXMMMMMMMMMMMMXXXXXXXXXX
							# XXXXXMMMMMMMMMMMXXXXXXX
						}
						else{
							if( !exists $paf_unit{"$array[0]\t$array[5]"} ){
								print $abc "$array[0]\t$array[5]\n";
								$paf_unit{"$array[0]\t$array[5]"}=0;
							}
						}
					}
				#}
			}
			else{
				#if( $array[7] <=500 || $array[8] >= $array[6]-500 ){
					if( $array[3]-$array[2] >= $array[1]*0.9 ){
						if( ($array[2] == 0 && $array[7] == 0 && $array[3]+200 < $array[1]) || ( $array[1] == $array[3] && $array[6] == $array[8] && $array[2] > 200 ) || ( ($array[2]>200) && ($array[3]+200 > $array[1] )) ){
						}
						else{
							if( !exists $paf_unit{"$array[0]\t$array[5]"} ){
								print $abc "$array[0]\t$array[5]\n";
								$paf_unit{"$array[0]\t$array[5]"}=0;
							}
						}
					}
				#}
			}
		}
	}
}

###Perform clustering analysis using the MCL algorithm.
system("mcl $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc --abc -o $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc.mcl");

###Visualize the clustering results with an R, allowing for manual inspection of the clusters by the user.
system("Rscript $script_dir/graph.plot.r $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc.mcl $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc.mcl.graph.pdf");

###Use a hash to store the output results from seqtk.
open my$telo,"$contigs_name.telomere/$fastq_name.telo.fa" or die "Can't open $contigs_name.telomere/$fastq_name.telo.fa";
$/=">";<$telo>;
my%hash_seq;
while(<$telo>){
        chomp;
        my($ID,$seq)=(split /\n/,$_,2)[0,1];
	if( $ID=~m/_remove/){
		next;
	}
	else{
		my($id,$start_telo,$end_telo,$length_telo)=($1,$2,$3,$4) if($ID=~m/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/);
		$ID=(split / /,$ID)[0];
		$ID=(split /\t/,$ID)[0];
		$seq=~s/\n//g;
		$seq=~s/>//;
		#print "$id\n";
		#print "$start_telo\n";
		#print "$end_telo\n";
		#print "$length_telo\n";
		if( $length_telo-$end_telo < $start_telo ){
			my$revcomp = reverse $seq;
			$revcomp =~ tr/ACGTacgt/TGCAtgca/;
			$seq=$revcomp;
		}
        	$hash_seq{$ID}=$seq;
	}
}
$/="\n";
close $telo;

my$telomere=0;
###Open the MCL clustering results.
open my$mcl,"$contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc.mcl" or die "Can't open $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc.mcl";
###Create spoa command lines
open my$spoa,">$contigs_name.telomere/$fastq_name.remove.telo.fa.spoa.sh";
while(<$mcl>){
	chomp;
	open my$cluster,">$contigs_name.telomere/$fastq_name.telomere$telomere.fa";
	my@reads=split /\t/,$_;
	foreach my$read (@reads){
		$read=~s/_remove//;###obtain the sequences with telomere
		print $cluster ">$read\n$hash_seq{$read}\n";
	}
	print $spoa "spoa -r 2 $contigs_name.telomere/$fastq_name.telomere$telomere.fa > $contigs_name.telomere/$fastq_name.telomere$telomere.cons.MSA.fa\n";
	$telomere++;
}

###Because SPOA can consume a high amount of memory, therefore only run two instances concurrently:
system("cat $contigs_name.telomere/$fastq_name.remove.telo.fa.spoa.sh | xargs -I {} -P $threads_spoa bash -c {}");

sub calculate_overlap {
    my ($start1, $end1, $start2, $end2) = @_;

    # Sort the coordinates to make sure start is always less than end
    ($start1, $end1) = ($start1 < $end1) ? ($start1, $end1) : ($end1, $start1);
    ($start2, $end2) = ($start2 < $end2) ? ($start2, $end2) : ($end2, $start2);

    # Find the maximum of the two starts and the minimum of the two ends
    my $overlap_start = ($start1 > $start2) ? $start1 : $start2;
    my $overlap_end = ($end1 < $end2) ? $end1 : $end2;

    # Calculate overlap
    my $overlap = $overlap_end - $overlap_start;

    # If the overlap is negative, there is no overlap
    return $overlap > 0 ? $overlap : 0;
}

###Using a loop to process each telomere sequences.
for (my $j=0; $j<$telomere; $j++) {

	###Change the name of consensus sequence.
	system("grep -A 1 'Consensus' $contigs_name.telomere/$fastq_name.telomere$j.cons.MSA.fa | sed 's/>Consensus/>telomere$j/' | sed 's/-//g' > $contigs_name.telomere/$fastq_name.telomere$j.cons.fa");
	
	###For the convenience of later visualization, adjust the telomere sequences to the start of the sequence.
	system("seqtk telo -m $unit $contigs_name.telomere/$fastq_name.telomere$j.cons.fa 1> $contigs_name.telomere/$fastq_name.telomere$j.cons.fa.stat 2>$contigs_name.telomere/$fastq_name.telomere$j.cons.fa.stat2");
	system("sed -i '1!d' \"$contigs_name.telomere/$fastq_name.telomere$j.cons.fa.stat\" ");
	
	open my$cons_seq, "$contigs_name.telomere/$fastq_name.telomere$j.cons.fa.stat" or die "Can't open $contigs_name.telomere/$fastq_name.telomere$j.cons.fa.stat";
	$/ = "\n";
	while (my$line=<$cons_seq>) {
		chomp$line;
		my@fields = split /\t/,$line;
		my$half=$fields[3]/2;
		my$overlap_start = calculate_overlap($fields[1], $fields[2], 0, $half);
		my$overlap_end = calculate_overlap($fields[1], $fields[2], $half, $fields[3]);		
		if ( $overlap_start >= $overlap_end ) { ###No change needed since the telomere is already at the beginning of the sequence.
			system("cp $contigs_name.telomere/$fastq_name.telomere$j.cons.fa $contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa");
			system("trf $contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa 2 5 7 80 10 50 2000 -l 10");
			system("mv $fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.*.txt.html $contigs_name.telomere");
			system("mv $fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.*.html $contigs_name.telomere");
			system("cat $contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.*.txt.html > $contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.merged.txt.html");
			system("perl ${script_dir}html2repeatbed.pl $contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.merged.txt.html > $contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.merged.txt.html.bed");
		}
		else{ ###Convert the sequence to its reverse complement.
			my$filename_cons = "$contigs_name.telomere/$fastq_name.telomere$j.cons.fa";
			my$out_filename_telomere_start = "$contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa";
			open my$in_cons, $filename_cons or die "Cannot open $filename_cons\n";
			open my$out_cons_telomere_start, ">$out_filename_telomere_start" or die "Cannot open $out_filename_telomere_start";
			$/=">";<$in_cons>;
			while (my$cons_line = <$in_cons>) {
				my($ID_cons,$seq_cons)=(split /\n/,$cons_line,2)[0,1];
				$seq_cons=~s/\n//g;
				$seq_cons=~s/>//g;
				my$revcomp = reverse $seq_cons;
				$revcomp =~ tr/ACGTacgt/TGCAtgca/;
				print $out_cons_telomere_start ">$ID_cons\n$revcomp\n";
			}
			$/="\n";
			close($in_cons);
			close($out_cons_telomere_start);
			system("trf $contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa 2 5 7 80 10 50 2000 -l 10");
			system("mv $fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.*.txt.html $contigs_name.telomere");
			system("mv $fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.*.html $contigs_name.telomere");
			system("cat $contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.*.txt.html > $contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.merged.txt.html");
			system("perl ${script_dir}html2repeatbed.pl $contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.merged.txt.html > $contigs_name.telomere/$fastq_name.telomere$j.cons.telomere_start.fa.2.5.7.80.10.50.2000.merged.txt.html.bed");
		}
	}
	close $cons_seq;

	open my$MSA,"$contigs_name.telomere/$fastq_name.telomere$j.cons.MSA.fa" or die "Can't open $contigs_name.telomere/$fastq_name.telomere$j.cons.MSA.fa";
	open my$MSA01,">$contigs_name.telomere/$fastq_name.telomere$j.cons.MSA.fa.matrix" or die "Can't open $contigs_name.telomere/$fastq_name.telomere$j.cons.MSA.fa.matrix";
	$/=">";<$MSA>;
	while(<$MSA>){
		chomp;
		my($ID,$seq)=(split /\n/,$_,2)[0,1];
		$ID=(split / /,$ID)[0];
		$ID=(split /\t/,$ID)[0];
		if($ID eq "Consensus"){
			$ID="telomere$j.cons.fa";
		}
		$seq=~s/>//;
		$seq=~s/\n//g;
		$seq=~tr/atcgATCG-/123412340/;
		my@split=(split //,$seq);
		unshift @split,$ID;
		print $MSA01 join(" ", @split);
		print $MSA01 "\n";
	}
	close $MSA;
	close $MSA01;

	system("Rscript ${script_dir}MSA.plot.r $contigs_name.telomere/$fastq_name.telomere$j.cons.MSA.fa.matrix $contigs_name.telomere/$fastq_name.telomere$j.cons.MSA.fa.matrix.pdf");
}

system("cat $contigs_name.telomere/*bed > $contigs_name.telomere/$fastq_name.telomere.cons.telomere_start.bed");
system("cat $contigs_name.telomere/*cons.telomere_start.fa > $contigs_name.telomere/$fastq_name.telomere.cons.telomere_start.fa");
system("Rscript ${script_dir}telomeres.visulization.r $contigs_name.telomere/$fastq_name.telomere.cons.telomere_start.fa $contigs_name.telomere/$fastq_name.telomere.cons.telomere_start.bed $contigs_name.telomere/$fastq_name.telomere.cons.telomere_start.pdf");

if( $telomere == 0 ){
	print "   Sorry, no telomeres have been assembled. :(\n";
}
else{
	print "   $telomere telomere regions have been assembled. :)\n";
}

if( $extend == 0 ){
	exit;
}

my$alltelomere="$contigs_name.all.telomere.cons.fa";
my$alltelomerepaf="$contigs_name.all.telomere.cons.fa.paf";
system("rm $contigs_name.telomere/$contigs_name.all.telomere.cons.fa");
system("cat $contigs_name.telomere/*cons.fa > $contigs_name.telomere/$alltelomere");
system("minimap2 -I 20G -x asm5 -t $threads -c --eqx $contigs $contigs_name.telomere/$alltelomere -o $contigs_name.telomere/$alltelomerepaf");

sub reverse_complement {
        my $sequence = shift;
        $sequence = reverse($sequence);
        $sequence =~ tr/ATCGatcg/TAGCtagc/;
        return $sequence;
}

open my$query_file,"$contigs_name.telomere/$alltelomere" or die "Can't open $contigs_name.telomere/$alltelomere"; #query.fa
my%hash_query;
$/=">";<$query_file>;
while(<$query_file>){
        chomp;
        my($ID,$seq)=(split /\n/,$_,2)[0,1];
        $ID=(split / /,$ID)[0];
        $ID=(split /\t/,$ID)[0];
        $seq=~s/>//;
        $seq=~s/\n//g;
        $hash_query{$ID}=$seq;
}
close $query_file;
$/="\n";

open my$target_file,"$contigs" or die "Can't open $contigs"; #target.fa
my%hash_target;
my%hash_targetlen;
$/=">";<$target_file>;
while(<$target_file>){
        chomp;
        my($ID,$seq)=(split /\n/,$_,2)[0,1];
        $ID=(split / /,$ID)[0];
        $ID=(split /\t/,$ID)[0];
        $seq=~s/>//;
        $seq=~s/\n//g;
        $hash_target{$ID}=$seq;
        $hash_targetlen{$ID}=length$seq;
}
close $target_file;
$/="\n";

open my$paf_file,"$contigs_name.telomere/$alltelomerepaf" or die "Can't open $contigs_name.telomere/$alltelomerepaf"; #paf
my%hash_paf;
while(<$paf_file>){
        chomp;
        my@fields=split /\t/,$_;
        my ($q_id, $t_id, $t_start, $t_end, $strand, $mapq,$q_len, $t_len, $matches) = ($fields[0], $fields[5], $fields[7], $fields[8], $fields[4], $fields[11],$fields[1], $fields[6], $fields[9]);
        next unless $mapq == 60; #uniq mapping
        next unless $matches >= $q_len*0.8; #true telomere
        next unless ($t_start < 50000 || $t_end + 50000 > $t_len);
	push @{$hash_paf{$t_id}},"$_";
}
close $paf_file;

my%hash_target_change;
foreach my$target (sort keys %hash_paf){
	my$num=@{$hash_paf{$target}};
	if($num == 1){
		my@fields=split /\t/,${$hash_paf{$target}}[0];
		my($q_id, $t_id, $t_start, $t_end, $strand, $mapq,$t_len) = ($fields[0], $fields[5], $fields[7], $fields[8], $fields[4], $fields[11], $fields[6]);
		my$q_seq=$hash_query{$q_id};
		if( $strand eq '-' ){
			$q_seq=reverse_complement($q_seq);
		}
		if( $t_end < $t_len/2 ){#start
			my$seq1=$q_seq;
			my$seq2=substr($hash_target{$t_id},$t_end+1);
			$hash_target_change{$t_id}="$seq1$seq2";
		}
		else{ #end
			my$seq1=substr($hash_target{$t_id},0,$t_start);
			my$seq2=$q_seq;
			$hash_target_change{$t_id}="$seq1$seq2";
		}
	}
	elsif($num == 2){
		my@fields1=split /\t/,${$hash_paf{$target}}[0];
		my@fields2=split /\t/,${$hash_paf{$target}}[1];
		my ($q_id1, $t_id1, $t_start1, $t_end1, $strand1, $mapq1, $t_len1) = ($fields1[0], $fields1[5], $fields1[7], $fields1[8], $fields1[4], $fields1[11], $fields1[6]);
		my ($q_id2, $t_id2, $t_start2, $t_end2, $strand2, $mapq2, $t_len2) = ($fields2[0], $fields2[5], $fields2[7], $fields2[8], $fields2[4], $fields2[11], $fields1[6]);
		my$q_seq1=$hash_query{$q_id1};
		my$q_seq2=$hash_query{$q_id2};
		if( $strand1 eq '-' ){
			$q_seq1=reverse_complement($q_seq1);
		}
		if( $strand2 eq '-' ){
			$q_seq2=reverse_complement($q_seq2);
		}
		if( $t_end1 < $t_len1/2 && $t_end2 > $t_len2/2 ){#q_id1 in start; q_id2 in end;
			my$seq1=$q_seq1;
			my$seq2=substr($hash_target{$t_id1},$t_start1, $t_len1-$t_start2+1 );
			my$seq3=$q_seq2;
			$hash_target_change{$t_id1}="$seq1$seq2$seq3";
		}
		elsif( $t_end1 > $t_len1/2 && $t_end2 < $t_len2/2 ){#q_id1 in end; q_id2 in start;
			my$seq1=$q_seq2;
			my$seq2=substr($hash_target{$t_id1},$t_start2,$t_len1-$t_start1+1 );
			my$seq3=$q_seq1;
			$hash_target_change{$t_id1}="$seq1$seq2$seq3";
		}
	}
}

my$tla="$contigs_name.TLA.fa";
open my$out,">$contigs_name.telomere/$tla" or die "Can't open $contigs_name.telomere/$tla\n";
foreach my$ID (sort keys %hash_target){
	if(exists $hash_target_change{$ID} ){
		#print "Changed:$ID\n";
		$hash_target_change{$ID}=~s/\n//g;
		print $out ">$ID\n$hash_target_change{$ID}\n";
	}
	else{
		#print "No Changed:$ID\n";
		$hash_target{$ID}=~s/\n//g;
		print $out ">$ID\n$hash_target{$ID}\n";
	}
}
