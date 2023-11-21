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
my$vcf_file_name;
my$current_directory = cwd();
my$contigs_name;
my$fastq_name;
my$bam_name;
my$threads=48;

GetOptions(
    'h'   => \$help,
    'd=s' => \$vcf_file,
    'u=s' => \$unit,
    'f=s' => \$fastq,
    'c=s' => \$contigs,
    'b=s' => \$bam,
    't=s' => \$threads,
);

if ($help or not defined $fastq and not defined $contigs) {
    print "Usage: $0\n";
    print "-h: show this help message.\n";
    #print "-d: specify the VCF file from deepvariant.\n";
    print "-b: bam file\n";
    #print "-u: telomere unit.\n";
    print "-f: hifi reads.\n";
    print "-c: contigs.\n";
    print "-t: threads.\n";
    exit;
}

$contigs_name=(split /\//,$contigs)[-1];

for my $tool ("minimap2","bgzip", "bcftools","samtools") {
	my $return_val = system("which $tool > /dev/null 2>&1");
	if ($return_val != 0) {
		die "Error: $tool does not seem to be installed or is not in the PATH. Please check.\n";
	}
}

my$script_dir = ($0 =~ m{(.*/)?})[0];
$script_dir =~ s/\/$//;
my$BIN_VERSION = "1.6.0"; 

unless (-d "$contigs_name.polish") {
        system("mkdir $contigs_name.polish");
}

if( not defined $bam ){
	$fastq_name=(split /\//,$fastq)[-1];
	$bam_name="$contigs_name.$fastq_name.sorted.bam";
	system("minimap2 -I 20G -ax map-hifi -t $threads $contigs $fastq | samtools sort -\@12 -T $contigs_name.polish/$contigs_name.$fastq_name.samtoolstmp -o $contigs_name.polish/$bam_name -");
	system("samtools index -\@12 $contigs_name.polish/$bam_name");
}
else{
	$bam_name=(split /\//,$bam)[-1];
	system("ln -s $current_directory/$bam $contigs_name.polish/$bam_name");
	system("ln -s $current_directory/$bam.bai $contigs_name.polish/$bam_name.bai");
}

if (-e "$script_dir/deepvariant_$BIN_VERSION.sif") {
	print "Required tools found\n";
	print "Starting Deepvariant\n";
	#$vcf_file="$contigs_name.deepvariant.vcf.gz";
	$vcf_file_name="$contigs_name.deepvariant.vcf.gz";
	system("samtools faidx $contigs");
	system("singularity exec -B /usr/lib/locale/:/usr/lib/locale/ $script_dir/deepvariant_$BIN_VERSION.sif /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref $contigs --reads $contigs_name.polish/$bam_name --output_vcf $contigs_name.polish/$vcf_file_name --num_shards $threads");
}
else{
	print "The deepvariant singularity container is not found.\n";
	print "Please navigate to $script_dir\n";
	print "Use the following command to download it: singularity pull docker://google/deepvariant:\"${BIN_VERSION}\"\n";
	exit;
}

print "Filtering the VCF file, and polishing the contigs\n";

open IN1, "gunzip -dc $contigs_name.polish/$vcf_file_name |" or die "Can't open 02.polish/$vcf_file_name";
my$variant_retained=0;
my $output_file = "$contigs_name.polish/$vcf_file_name.filter.vcf.gz";
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

my$tlapolish_fa="$contigs_name.polished.fa";
print "bcftools index\n";
system("bcftools index $output_file");
print "bcftools consensus\n";
system("bcftools consensus -o $contigs_name.polish/$tlapolish_fa -f $contigs $output_file");
