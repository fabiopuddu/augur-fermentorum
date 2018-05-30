#!/usr/bin/env perl

# Author:       	Fabio Puddu  
# Maintainer:   	Fabio Puddu
# Created:      	Oct 2016
# Description:	This script takes yeast fastq files in input and returns the number of telomeric reads per million

use strict;
use warnings;

#my $i=1;
my $tresh1= shift @ARGV; #get first threshold: minimum number of GT/AC repeats
#my $tresh2= shift@ARGV; #get second threshold: minimum GT% in the read
my $infastq= shift @ARGV; #get path to fastq file; can be omitted, it will read from STDIN
#my %base_numbers; #counter to get total number of base in each read

my $fh;
my $is_stdin = 0; #boolean value to determine wheter STDIN should be used
if (defined $infastq){
  open my $fh, "<", $infastq or die $!;  #open provided fastq file
} 
else {
  $fh = *STDIN; #open STDIN 
  $is_stdin++; #STDIN is being used
}

my $totcount=0;
my $telocount=0;
while(my $row=<$fh>){ #read file or STDIN line by line
	chomp $row;
	next unless $row =~ /^[TAGCN]+$/; #skip line if it does not contain DNA
	next if (length $row < 100); #skip reads that are shorter than expected
	$totcount++;
	#next unless $row =~ /TG{1,3}|C{1,3}A/; #skip the line if the DNA does not contain any telomeric sequence
	next unless $row =~ /(TG){1,3}TG{2,3}|C{2,3}A(CA){1,3}/; #skip the line if the DNA does not contain any telomeric sequence
#	my %base_numbers= (
#  			  	"A" => "0",	#
#    				"C" => "0",	#
#    				"G" => "0",	# initialise counters#
#				"T" => "0",	#
#				"N" => "0",	#
#	);
#	my @sequence = split '', $row;
#	foreach my $base (@sequence){		
#		   $base_numbers{"$base"}++; #count each base		
#	}
#	my $gt_perc = ($base_numbers{'G'}+$base_numbers{'T'})/@sequence*100; #calculate GT percent
	my $count=0; #initialise counter of repeat instances	
	#while ($row =~ /TG{1,3}|C{1,3}A/g){
	while ($row =~ /(TG){1,3}TG{2,3}|C{2,3}A(CA){1,3}/g){
	$count++	#count repeat instances
	}
#        printf "$row\n" if ($count>=$tresh1 && ($gt_perc > $tresh2 || $gt_perc < (100-$tresh2))); #output the reads if the number of matches is higher than thresh1 and the GT or AC percent is higher than tresh2
	#printf "$row\n" if $count>=$tresh1;
	$telocount++ if $count>=$tresh1;
}

my $telo_estimate = $telocount/$totcount*1000000; #telomeric reads per million reads
printf "$telo_estimate"

