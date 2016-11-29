#!/usr/bin/env perl
use strict;
use warnings;

my @start;
my @end;
my @name;
my $ifh;
my @arrow_start;
my @arrow_end;
my $sed_start;
my $sed_end;
my @orientation;
my $l;
my @enrichment;
my @average;
my$x;	

open($ifh,"targets.tsv	") || die "Failed: $!\n";
	while (my $l = <$ifh> )
	{
  	chomp( $l );    
  	my @s = split( /\t/, $l );	
  	push(@name, $s[6]);
  	push(@start, $s[0]);
  	push(@end, $s[1]);
  	push(@arrow_start, $s[3]);
  	push(@arrow_end, $s[4]);  	
  #	push(@orientation, $s[6])
  	}
close ($ifh);
my $len = @start; 
for (my $i=0; $i < $len; $i++){
				$sed_start = $start[$i] - 200;
				$sed_end = $end[$i] + 200;				
				# open (my $xfh,"sed -n $sed_start,$sed_end"."p wg-results.txt | ");
# 				while (my $line = <$xfh> ){
#   					chomp( $line );    
#   					my @s = split( /\t/, $line );
#   					push(@enrichment, $s[5]);
# 				}
# 				$l = @enrichment;
# 				if ($orientation[$i] == 1){
# 							for (my $j=0; $j < $l; $j++){
# 							$average[$j] += shift(@enrichment);
# 							}
# 				}
# 				else {		for (my $j=0; $j < $l; $j++){
# 							$average[$j] += pop(@enrichment);
# 							}
# 				}
				system("gnuplot -e \"as=$arrow_start[$i];ae=$arrow_end[$i];s=$start[$i];e=$end[$i];of=\'$name[$i].pdf\';\" plot.gpl" );
				print ("$i \t as=$arrow_start[$i];ae=$arrow_end[$i];s=$start[$i];e=$end[$i];of=$name[$i].pdf; \n");
}

open (my $ofh,">","average.txt");
foreach my $k (@average){
	$x = $k / $len;
	print ($ofh "$x \n");
}
close ($ofh);