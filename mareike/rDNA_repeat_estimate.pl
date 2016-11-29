#!/usr/bin/env perl
# 
# Author:       mh23
# Maintainer:   mh23
# Created:      16.08.2014
# Name:			rDNA_repeat_estimate.pl
# Test:			perl rDNA_repeat_estimate.pl -i SC_MFY5784126.B.wig.XII.rDNA
#				for x in *.rDNA; do echo $x; perl ../rDNA_repeat_estimate.pl $x >> Rate_140816.txt ; done

# Description:
# This program takes a .rDNA coverage file and sums the coverage from 3rd to the 400th entry (440051 to 449976)
# and the 460th to  1152nd entry (451701 to 469010). It then divides the former by 398 and the latter 693.
# The second result is divided by the first and multiplied by two. This is our rDNA repeat estimate.
# It will print:
# Sample_name \t Repeat_number_estimate

#standard use
use strict;
use warnings;
use Getopt::Long;


#Get input from the command line
my ($input);

GetOptions
(
'i|input=s'         => \$input,
);

( $input && -f $input ) or die qq[Usage: $0 -i <input vcf>\n];



my $ifh;
if( $input =~ /\.gz/ ){open($ifh, qq[gunzip -c $input|]);}else{open($ifh, $input ) or die $!;}

#declare hash to store coverage data
my %cov;

#build the coverage into a hash to draw from later
#Go through the coverage files
while( my $l = <$ifh> ){
	#read the coverage into a hash
	chomp( $l );
    my @s = split( / /, $l );
	$cov{$s[0]} = $s[1];
}
close( $ifh );

#This is what we want
#mean: 3 bis 400		440051 to 449976	divide by 398
#rDNA: 460 bis 1152		451701 to 469010	divide by 693

#declare necessary variables
my $mean_cov;
my $mean_sum = 0;
my $rDNA_cov;
my $rDNA_sum = 0;
my $cont1 = 0;
my $cont2 = 0;
#Calculate the mean coverage
for (my $i=440045; $i <= 449980; $i++) {
    #check teh value exists; only every 25th should
    if (exists $cov{$i}){
    #sum them up
    $mean_sum = $mean_sum + $cov{$i}; 
    $cont1 = $cont1 +1 ; 
}}
#get the average coverage
$mean_cov = $mean_sum/$cont1;

#Calculate the rDNA coverage
#repeat as above
for (my $i=452000; $i <= 459000; $i++) { #FABIO: I edited the region used to calculate the coverage. will explain
    if (exists $cov{$i}){
    $rDNA_sum = $rDNA_sum + $cov{$i};   
    $cont2 = $cont2 +1 ; 
}}
$rDNA_cov = $rDNA_sum/$cont2;

#claculate the repeat number
my $repeats;
$repeats = $rDNA_cov/$mean_cov*2;

#print("mean=$mean_cov\t rDNA=$rDNA_cov\n");

#split the file name on the dot to get just the SC number
my @name = split( /\./, $input );

#print sample name and repeat
print("$name[0]\t$repeats\t");

		
