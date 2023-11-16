use strict;
use warnings;

##############
#overlapped records (>50%) will be removed
##############

open IN0, "$ARGV[0]" or die "Can't open $ARGV[0]";
$/ = "Alignment explanation";
my %hash_overlap_long;

while (<IN0>) {
    my @unit_f = split /\n/, $_;
    my $all_f = @unit_f;
    my $start_f;
    my $end_f;
    my $range_f;

    for (my $j = 0; $j < $all_f; $j++) {
        if ($unit_f[$j] =~ m/Indices/) {
            ($range_f) = $unit_f[$j] =~ /Indices:\s+(\S+?)\s+/;
            ($start_f, $end_f) = split /--/, $range_f;
	    my$skip=0;

            foreach my $bed (keys %hash_overlap_long) {
                my ($start1, $end1) = split /\t/, $bed;

                # Calculate overlap
                my $overlap = 0;
                my $start_overlap = $start_f > $start1 ? $start_f : $start1;
                my $end_overlap = $end_f < $end1 ? $end_f : $end1;
                $overlap = $end_overlap - $start_overlap if $end_overlap > $start_overlap;

                # Calculate lengths
                my $length_f = $end_f - $start_f;
                my $length1 = $end1 - $start1;

                # Check conditions and replace or skip as necessary
                if ($overlap / $length1 > 0.5 && $length_f > $length1) {
                    # Replace the key in %hash_overlap_long with start_f and end_f
                    delete $hash_overlap_long{$bed};
                    $hash_overlap_long{"$start_f\t$end_f"} = 1;
                } elsif ($overlap / $length_f > 0.5 && $length_f < $length1) {
                    # Skip to the next $bed
		    $skip=1;
		    last;
                }
            }

            # If start_f and end_f did not match any $bed, add it to the hash
            if ($skip==0) {
                $hash_overlap_long{"$start_f\t$end_f"} = 1;
            }
        }
    }
}

open IN1,"$ARGV[0]" or die "Can't open $ARGV[0]";
$/="Alignment explanation";
my$first=0;
my$seq_name;
while(<IN1>){
	if($first == 0){
		$first=1;
		my@unit=split /\n/,$_;
		foreach my$line (@unit){
			if($line=~m/Sequence:/){
				($seq_name) = $line =~ m/Sequence:\s+(\S+)/;
			}
		}
	}
	else{
		my@unit=split /\n/,$_;
		my$all=@unit;
		my$start=0;
		my$end=0;
		my$range;
		my$First_one=0;
		my$total_len=0;
		for(my$i=0; $i< $all; $i++){
			if( $unit[$i]=~m/Indices/ ){
				($range) = $unit[$i] =~/Indices:\s+(\S+?)\s+/;
				($start,$end) = (split /--/,$range)[0,1];
				if( !exists $hash_overlap_long{"$start\t$end"} ){
					last;
				}
				$total_len=$end-$start;
				if( $total_len < 100){
					last;
				}
			}
			else{
				#my@columns = split /\s+/, $unit[$i];
				#my$num=0;
				#$num=@columns;
				#if( $unit[$i]=~m/^\s+?1\s+?(\S+?)$/ ){
				if( $unit[$i]=~m/^\s+?1\s+?([A-Za-z\s-]+)$/ ){
					if( $First_one == 0 ){
						if( $unit[$i-1] =~ /[ATCGatcg]/i ){
							my($Start,$seq)=($1,$2) if( $unit[$i-1] =~ /\s+(\d+)\s+(\S+)/ );
							print "$seq_name\t$Start\t";
							$First_one = 1;
						}
					}
					else{
						if( $unit[$i-1] =~ /[ATCGatcg]/i ){
							my($Start,$seq)=($1,$2) if( $unit[$i-1] =~ /\s+(\d+)\s+(\S+)/ );
							my$End=$Start-1;
							print "$End\t$range\n";
							print "$seq_name\t$Start\t";
						}
					}
				}
			}
		}
		if($total_len >= 100 && exists $hash_overlap_long{"$start\t$end"} ){
			print "$end\t$range\n";
		}
	}
}
