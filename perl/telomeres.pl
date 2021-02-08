#!/usr/bin/env perl

# Author:       	Fabio Puddu
# Maintainer:   	Fabio Puddu
# Created:      	Oct 2016
# Description:	This script takes yeast fastq files in input and returns the number of telomeric reads per million
use strict;
use warnings;
my $tresh1= shift @ARGV;                        #get first threshold: minimum number of GT/AC repeats
my $tresh2= shift@ARGV;                         #get read threshold: minimum read length
my $infastq= shift @ARGV;                       #get path to fastq file; can be omitted, it will read from STDIN
my $fh;
my $is_stdin = 0;                               #boolean value to determine wheter STDIN should be used
if (defined $infastq){                          #if the user defined a fastq file
        open my $fh, "<", $infastq or die $!;   #open provided fastq file
}
else {                                          #otherwise
  $fh = *STDIN;                                 #open STDIN
  $is_stdin++;                                  #STDIN is being used
}

my $totcount=0;
my $telocount=0;
while(my $row=<$fh>){                           #read file or STDIN line by line
	chomp $row;
	next unless $row =~ /^[TAGCN]+$/;      #skip line if it does not contain DNA
	next if (length $row < $tresh2);           #skip reads that are significantly shorter than expected
	$totcount++;
	next unless $row =~ /(TG){1,3}TG{2,3}|C{2,3}A(CA){1,3}/; #skip the line if the DNA does not contain any telomeric sequence
	my $count=0;                                           #initialise counter of repeat instances
	while ($row =~ /(TG){1,3}TG{2,3}|C{2,3}A(CA){1,3}/g){
	               $count++	                               #count repeat instances
	}
	$telocount++ if $count>=$tresh1;
}
my $telo_estimate = $telocount/$totcount*1000000; #telomeric reads per million reads
printf "$telo_estimate"
