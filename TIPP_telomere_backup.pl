#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Long;

my$help;
my$vcf_file;
my$unit;
my$fastq;
my$contigs;
my$bam;
my$threads=40;
my$vcf_file_name;
my$current_directory = cwd();
my$contigs_name;
my$fastq_name;
my$extend=0;
my$telomere_length=1000;

GetOptions(
    'h'   => \$help,
    'd=s' => \$vcf_file,
    'u=s' => \$unit,
    'f=s' => \$fastq,
    'c=s' => \$contigs,
    'b=s' => \$bam,
    't=s' => \$threads,
    'e=s' => \$extend,
    'l=s' => \$telomere_length,
);

if ($help or not defined $vcf_file and not defined $unit and not defined $fastq and not defined $contigs and not defined $bam) {
    print "Usage: $0\n";
    print "-h: show this help message.\n";
    #print "-d: specify the VCF file from deepvariant.\n";
    #print "-b: bam file\n";
    print "-u: telomere unit.\n";
    print "-f: hifi reads.\n";
    print "-e: extend the contigs with new assembled telomere sequences.(default=0, no extend)\n";
    print "-c: contigs. If the extension is not specified, it will be used as the output name.\n";
    print "-t: threads for minimap2.\n";
    print "-l: telomere length used in seqtk -d.(default=1000)\n";
    exit;
}

for my $tool ("spoa", "seqtk", "minimap2", "mcl", "bgzip", "bcftools","samtools") {
	my $return_val = system("which $tool > /dev/null 2>&1");
	if ($return_val != 0) {
		die "Error: $tool does not seem to be installed or is not in the PATH. Please check.\n";
	}
}

$contigs_name=(split /\//,$contigs)[-1];
$fastq_name=(split /\//,$fastq)[-1];

my$script_dir = ($0 =~ m{(.*/)?})[0];
my$BIN_VERSION = "1.6.0"; 

print "Extracting telomere reads, classifying into clusters, and generating a consensus for each cluster\n";

unless (-d "$contigs_name.telomere") {
        system("mkdir $contigs_name.telomere");
}

system("seqtk telo -m $unit $fastq -d $telomere_length > $contigs_name.telomere/$fastq_name.telo.fa");
system("grep -A 1 'remove' --no-group-separator $contigs_name.telomere/$fastq_name.telo.fa > $contigs_name.telomere/$fastq_name.remove.telo.fa");
system("minimap2 --eqx -c -x ava-pb -X -t $threads $contigs_name.telomere/$fastq_name.remove.telo.fa $contigs_name.telomere/$fastq_name.remove.telo.fa -o $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf");

#system("awk '{if( \$10>\$11*0.98 && ( (\$4-\$3)>\$2*0.9 && \$2>=1000 ) || (\$9-\$8)>\$7*0.9 && \$7>=1000 )print}' 01.telomere_local_assembly/${fastq_name}.remove.telo.fa.self.paf | awk '{print \$1\"\\t\"\$6}' > 01.telomere_local_assembly/${fastq_name}.remove.telo.fa.self.paf.abc");
system("awk '{if( ((\$10/\$11)>0.98) && (( (\$4-\$3)>\$2*0.9 && \$2>=1000 ) || ((\$9-\$8)>\$7*0.9 && \$7>=1000)) && (\$1 != \$6) ) print \$1\"\\t\"\$6}' $contigs_name.telomere/${fastq_name}.remove.telo.fa.self.paf | sort | uniq > $contigs_name.telomere/${fastq_name}.remove.telo.fa.self.paf.abc");

system("mcl $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc --abc -o $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc.mcl");
system("Rscript $script_dir/graph.plot.r $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc.mcl $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc.mcl.graph.pdf");

open my$telo,"$contigs_name.telomere/$fastq_name.telo.fa" or die "Can't open $contigs_name.telomere/$fastq_name.telo.fa";
$/=">";<$telo>;
my%hash_seq;
while(<$telo>){
        chomp;
        my($ID,$seq)=(split /\n/,$_,2)[0,1];
        $ID=(split / /,$ID)[0];
        $ID=(split /\t/,$ID)[0];
        $seq=~s/\n//g;
        $seq=~s/>//;
        $hash_seq{$ID}=$seq;
}
$/="\n";
close $telo;

my$telomere=0;
open my$mcl,"$contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc.mcl" or die "Can't open $contigs_name.telomere/$fastq_name.remove.telo.fa.self.paf.abc.mcl";
open my$spoa,">$contigs_name.telomere/$fastq_name.remove.telo.fa.spoa.sh";
while(<$mcl>){
        chomp;
        open my$cluster,">$contigs_name.telomere/$fastq_name.telomere$telomere.fa";
        my@reads=split /\t/,$_;
        foreach my$read (@reads){
                $read=~s/_remove//;
                print $cluster ">$read\n$hash_seq{$read}\n";
        }
        print $spoa "spoa -r 2 -s $contigs_name.telomere/$fastq_name.telomere$telomere.fa > $contigs_name.telomere/$fastq_name.telomere$telomere.cons.MSA.fa\n";
	$telomere++;
}

system("cat $contigs_name.telomere/$fastq_name.remove.telo.fa.spoa.sh | xargs -I {} -P 2 bash -c {}");

for (my $j=0; $j<$telomere; $j++) {
	system("grep -A 1 'Consensus' $contigs_name.telomere/$fastq_name.telomere$j.cons.MSA.fa | sed 's/>Consensus/>telomere$j/' | sed 's/-//g' > $contigs_name.telomere/$fastq_name.telomere$j.cons.fa");
	
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
