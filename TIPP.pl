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

GetOptions(
    'h'   => \$help,
    'd=s' => \$vcf_file,
    'u=s' => \$unit,
    'f=s' => \$fastq,
    'c=s' => \$contigs,
    'b=s' => \$bam,
    't=s' => \$threads,
);

if ($help or not defined $vcf_file and not defined $unit and not defined $fastq and not defined $contigs and not defined $bam) {
    print "Usage: $0\n";
    print "-h: show this help message.\n";
    #print "-d: specify the VCF file from deepvariant.\n";
    #print "-b: bam file\n";
    print "-u: telomere unit.\n";
    print "-f: hifi reads.\n";
    print "-c: contigs.\n";
    print "-t: threads for minimap2,hifiasm,diamond.\n";
    exit;
}

for my $tool ("spoa", "seqtk", "minimap2", "mcl", "parallel", "bgzip", "bcftools","samtools","diamond","hifiasm") {
	my $return_val = system("which $tool > /dev/null 2>&1");
	if ($return_val != 0) {
		die "Error: $tool does not seem to be installed or is not in the PATH. Please check.\n";
	}
}

$contigs_name=(split /\//,$contigs)[-1];
$fastq_name=(split /\//,$fastq)[-1];

my$script_dir = ($0 =~ m{(.*/)?})[0];
my$BIN_VERSION = "1.6.0"; 

print "Step 1: Extracting telomere reads, classifying into clusters, and generating a consensus for each cluster\n";

unless (-d "01.telomere_local_assembly") {
        system("mkdir 01.telomere_local_assembly");
}

system("seqtk telo -m $unit $fastq -d 1000 > 01.telomere_local_assembly/$fastq_name.telo.fa");
system("grep -A 1 'remove' --no-group-separator 01.telomere_local_assembly/$fastq_name.telo.fa > 01.telomere_local_assembly/$fastq_name.remove.telo.fa");
system("minimap2 --eqx -c -x ava-pb -t $threads 01.telomere_local_assembly/$fastq_name.remove.telo.fa 01.telomere_local_assembly/$fastq_name.remove.telo.fa -o 01.telomere_local_assembly/$fastq_name.remove.telo.fa.self.paf");

system("awk '{if( \$10>\$11*0.98 && ( (\$4-\$3)>\$2*0.9 && \$2>=1000 ) || (\$9-\$8)>\$7*0.9 && \$7>=1000 )print}' 01.telomere_local_assembly/${fastq_name}.remove.telo.fa.self.paf | awk '{print \$1\"\\t\"\$6}' > 01.telomere_local_assembly/${fastq_name}.remove.telo.fa.self.paf.abc");

system("mcl 01.telomere_local_assembly/$fastq_name.remove.telo.fa.self.paf.abc --abc -o 01.telomere_local_assembly/$fastq_name.remove.telo.fa.self.paf.abc.mcl");

open my$telo,"01.telomere_local_assembly/$fastq_name.telo.fa" or die "Can't open 01.telomere_local_assembly/$fastq_name.telo.fa";
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
open my$mcl,"01.telomere_local_assembly/$fastq_name.remove.telo.fa.self.paf.abc.mcl" or die "Can't open 01.telomere_local_assembly/$fastq_name.remove.telo.fa.self.paf.abc.mcl";
open my$spoa,">01.telomere_local_assembly/$fastq_name.remove.telo.fa.spoa.sh";
while(<$mcl>){
        chomp;
        open my$cluster,">01.telomere_local_assembly/$fastq_name.telomere$telomere.fa";
        my@reads=split /\t/,$_;
        foreach my$read (@reads){
                $read=~s/_remove//;
                print $cluster ">$read\n$hash_seq{$read}\n";
        }
        print $spoa "spoa -r 2 -s 01.telomere_local_assembly/$fastq_name.telomere$telomere.fa > 01.telomere_local_assembly/$fastq_name.telomere$telomere.cons.MSA.fa\n";
	$telomere++;
}

system("cat 01.telomere_local_assembly/$fastq_name.remove.telo.fa.spoa.sh | xargs -I {} -P 2 bash -c {}");

for (my $j=0; $j<$telomere; $j++) {
	system("grep -A 1 'Consensus' 01.telomere_local_assembly/$fastq_name.telomere$j.cons.MSA.fa | sed 's/>Consensus/>telomere$j/' | sed 's/-//g' > 01.telomere_local_assembly/$fastq_name.telomere$j.cons.fa");
	
	open my$MSA,"01.telomere_local_assembly/$fastq_name.telomere$j.cons.MSA.fa" or die "Can't open 01.telomere_local_assembly/$fastq_name.telomere$j.cons.MSA.fa";
	open my$MSA01,">01.telomere_local_assembly/$fastq_name.telomere$j.cons.MSA.fa.matrix" or die "Can't open 01.telomere_local_assembly/$fastq_name.telomere$j.cons.MSA.fa.matrix";
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

	system("Rscript ${script_dir}MSA.plot.r 01.telomere_local_assembly/$fastq_name.telomere$j.cons.MSA.fa.matrix 01.telomere_local_assembly/$fastq_name.telomere$j.cons.MSA.fa.matrix.pdf");
}

if( $telomere == 0 ){
	print "   Sorry, no telomeres have been assembled. :(\n";
}
else{
	print "   $telomere telomere regions have been assembled. :)\n";
}

my$alltelomere="$contigs_name.all.telomere.cons.fa";
my$alltelomerepaf="$contigs_name.all.telomere.cons.fa.paf";
system("rm 01.telomere_local_assembly/$contigs_name.all.telomere.cons.fa");
system("cat 01.telomere_local_assembly/*cons.fa > 01.telomere_local_assembly/$alltelomere");
system("minimap2 -x asm5 -I 20G -t $threads -c --eqx $contigs 01.telomere_local_assembly/$alltelomere -o 01.telomere_local_assembly/$alltelomerepaf");

sub reverse_complement {
        my $sequence = shift;
        $sequence = reverse($sequence);
        $sequence =~ tr/ATCGatcg/TAGCtagc/;
        return $sequence;
}

open my$query_file,"01.telomere_local_assembly/$alltelomere" or die "Can't open 01.telomere_local_assembly/$alltelomere"; #query.fa
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

open my$paf_file,"01.telomere_local_assembly/$alltelomerepaf" or die "Can't open 01.telomere_local_assembly/$alltelomerepaf"; #paf
my%hash_paf;
while(<$paf_file>){
        chomp;
        my@fields=split /\t/,$_;
        my ($q_id, $t_id, $t_start, $t_end, $strand, $mapq,$q_len, $t_len, $matches) = ($fields[0], $fields[5], $fields[7], $fields[8], $fields[4], $fields[11],$fields[1], $fields[6], $fields[9]);
        next unless $mapq == 60; #uniq mapping
        next unless $matches >= $q_len*0.8; #true telomere
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
open my$out,">01.telomere_local_assembly/$tla" or die "Can't open 01.telomere_local_assembly/$tla\n";
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


unless (-d "02.polish") {
	system("mkdir 02.polish");
}

print "Step2. Polishing\n";
print "Starting alignment\n";
system("minimap2 -ax map-hifi -I 20G -t $threads 01.telomere_local_assembly/$tla $fastq | samtools sort -\@12 -T 02.polish/$tla.$fastq_name.samtoolstmp -o 02.polish/$tla.$fastq_name.sorted.bam -");
system("samtools index -\@12 02.polish/$tla.$fastq_name.sorted.bam");

if (-e "$script_dir/deepvariant_$BIN_VERSION.sif") {
	print "Required tools found\n";
	print "Starting Deepvariant\n";
	$vcf_file="$contigs_name.deepvariant.vcf.gz";
	$vcf_file_name="$contigs_name.deepvariant.vcf.gz";
	system("samtools faidx 01.telomere_local_assembly/$tla");
	system("singularity exec -B /usr/lib/locale/:/usr/lib/locale/ $script_dir/deepvariant_$BIN_VERSION.sif /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref 01.telomere_local_assembly/$tla --reads 02.polish/$tla.$fastq_name.sorted.bam --output_vcf 02.polish/$vcf_file --num_shards $threads");
}
else{
	print "The deepvariant singularity container is not found.\n";
	print "Please navigate to $script_dir\n";
	print "Use the following command to download it: singularity pull docker://google/deepvariant:\"${BIN_VERSION}\"\n";
	exit;
}

print "Filtering the VCF file, and polishing the contigs\n";

if ($vcf_file =~ /\.gz$/) {
	open IN1, "gunzip -dc 02.polish/$vcf_file_name |" or die "Can't open 02.polish/$vcf_file_name";
}
else{
	open IN1, "02.polish/$vcf_file_name" or die "Can't open 02.polish/$vcf_file_name";
}

my$variant_retained=0;

my $output_file = "02.polish/$vcf_file_name.filter.vcf.gz";
open my $OUT, "| bgzip > $output_file" or die "Can't open output file: $output_file";

while (<IN1>) {
	chomp;
	if(/^#/) {
		print $OUT "$_\n";
	}
	else{
		my@fields = split(/\t/, $_);
		next unless $fields[6] eq 'PASS';

		my@alts = split(/,/, $fields[4]);
		next unless scalar @alts == 1;

		next unless length($alts[0]) <= 50;

		my@formats = split(/:/, $fields[8]);
		my@sample_vals = split(/:/, $fields[9]);
		my ($GQ_idx, $DP_idx, $AD_idx, $VAF_idx) = (-1, -1, -1, -1);

		for (my $i = 0; $i < scalar @formats; $i++) {
			if($formats[$i] eq "GQ") {
				$GQ_idx = $i;
			}
			elsif($formats[$i] eq "DP") {
                		$DP_idx = $i;
            		}
			elsif ($formats[$i] eq "AD") {
				$AD_idx = $i;
			}
			elsif ($formats[$i] eq "VAF") {
				$VAF_idx = $i;
			}
		}
		if ($GQ_idx == -1 || $DP_idx == -1 || $AD_idx == -1 || $VAF_idx == -1) {
			die "One or more required fields (GQ, DP, AD, VAF) were not found in the VCF format column.\n";
		}

		my($ref_depth, $alt_depth) = (split /,/, $sample_vals[$AD_idx])[0, 1];
		next unless $ref_depth <= 3;
		next unless $sample_vals[$GQ_idx] >= 3;
		next unless $sample_vals[$VAF_idx] >= 0.5;
		print $OUT "$_\n";
		$variant_retained++;
	}
}
close $OUT;

if ($variant_retained == 0) {
    print "   After filtering, no variants were retained :(\n";
}
elsif ($variant_retained > 0 && -e $output_file) {
    print "   After filtering, $variant_retained variants were retained :)\n";
}

my$tlapolish_fa="$contigs_name.tla.polish.fa";

system("bcftools index $output_file");
system("bcftools consensus -o 02.polish/$tlapolish_fa -f 01.telomere_local_assembly/$tla $output_file");

#chloroplast
unless (-d "03.chloroplast") {
        system("mkdir 03.chloroplast");
}

my$bam_file = "02.polish/$tla.$fastq_name.sorted.bam";

open my $fh, "samtools view -H $bam_file |" or die "Failed to open pipe from samtools: $!";
my @large_refs;
while (<$fh>) {
    if (/^\@SQ\tSN:(\S+)\tLN:(\d+)/) {
        push @large_refs, $1 if $2 > 1000000;
    }
}
close $fh;


open my$full_map,">03.chloroplast/$fastq_name.fullmap.ID";
foreach my $ref (@large_refs) {
    open my $fh, "samtools view $bam_file $ref |" or die "Failed to open pipe from samtools: $!";
    while (<$fh>) {
        my @fields = split /\t/;
        my $read_id = $fields[0];
        my $cigar = $fields[5];

        my $clip_condition = 1;
        while ($cigar =~ /(\d+)([SH])/g) {
            if ($1 >= 100) {
                $clip_condition = 0;
                last;
            }
        }
        print $full_map "$read_id\n" if $clip_condition;
    }
    close $fh;
}
close $full_map;

if (-e "$script_dir/coregene.fa" && -e "$script_dir/coregene.fa.dmnd" ) {
	print "Chloroplast seed found!\n";
	print "Diamond search starting!\n";
	system("diamond blastx --threads $threads --db $script_dir/coregene.fa --query $fastq --max-target-seqs 50 --outfmt 6 qseqid full_qseq --out 03.chloroplast/$fastq_name.chloroplast.diamond.out");
}
else{
	print "Seed does not exist!";
	exit;
}

system("awk '{print \$1}' 03.chloroplast/$fastq_name.chloroplast.diamond.out | sort | uniq -c |awk '{if(\$1 >=3) print \$2}' > 03.chloroplast/$fastq_name.chloroplast.candidate.ID");

open my$diamond_canID,"03.chloroplast/$fastq_name.chloroplast.candidate.ID" or die "Can't open 03.chloroplast/$fastq_name.chloroplast.candidate.ID";
my%diamond_canID_hash;
while(<$diamond_canID>){
	chomp;
	$diamond_canID_hash{$_}=$_;
}
close $diamond_canID;

open my$full_map_ID,"03.chloroplast/$fastq_name.fullmap.ID" or die "Can't open 03.chloroplast/$fastq_name.fullmap.ID";
my%full_map_hash;
while(<$full_map_ID>){
	chomp;
	$full_map_hash{$_}=$_;
}
close $full_map_ID;

open my $diamond_out, "03.chloroplast/$fastq_name.chloroplast.diamond.out" or die "Can't open 03.chloroplast/$fastq_name.chloroplast.diamond.out";
open my $chloroplast_out, ">03.chloroplast/$fastq_name.chloroplast.reads.fasta";
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

open my $chloroplast, ">03.chloroplast/Chloroplast.fa" or die "Cannot open 03.chloroplast/Chloroplast.fa: $!\n";

for (my $round=1; $round<50; $round++){
    print "Starting downsampling of chloroplast reads, round$round...\n";

    system("seqtk sample -s$round 03.chloroplast/$fastq_name.chloroplast.reads.fasta 2000 > 03.chloroplast/$fastq_name.chloroplast.reads.subsample2000.round$round.fasta");

    system("hifiasm -D 1000 -o 03.chloroplast/$fastq_name.chloroplast.reads.subsample2000.round$round -t $threads 03.chloroplast/$fastq_name.chloroplast.reads.subsample2000.round$round.fasta");

    my $result = `grep '^S' 03.chloroplast/$fastq_name.chloroplast.reads.subsample2000.round$round.bp.p_ctg.gfa | grep 'c\\s'`;

    if ($result) {
        print "Circular contig found in round$round.\n";
        my ($ID, $seq) = (split /\t/, $result)[1,2];
        print $chloroplast ">$ID\n$seq\n";
        last;
    }
}
close $chloroplast;

