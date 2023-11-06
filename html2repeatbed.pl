use strict;
use warnings;
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
		my$total_len;
		for(my$i=0; $i< $all; $i++){
			if( $unit[$i]=~m/Indices/ ){
				($range) = $unit[$i] =~/Indices:\s+(\S+?)\s+/;
				($start,$end) = (split /--/,$range)[0,1];
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
		if($total_len >= 100 ){
			print "$end\t$range\n";
		}
	}
}
