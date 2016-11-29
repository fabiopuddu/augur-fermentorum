#!/usr/bin/env perl
# 
# Author:       mh23	
# Maintainer:   mh23
# Created: 		23.02.2015

#Test: perl /nfs/users/nfs_m/mh23/Scripts/bam_stats_table.pl -i test.stats
#description: A script turn the output from  samtools stats in a table of
#sample_name total_bases mapped_bases mapped_bases-(duplicated_reads*read_length) (mapped_bases-(duplicated_reads*read_length))/12.1Mb
 
         
use Carp;
use strict;
use warnings;
use Getopt::Long;

my ($input);

GetOptions
(
'i|input=s'         => \$input,
);

( $input && -f $input ) or die qq[Usage: $0 -i <input vcf>\n];



my $ifh;
open(F, $input ) or die ("Unable to open file $input: $!\n" );
my @s = split( /\./, $input );
my $sample_name = $s[0];
my $total_bases = 0;
my $mapped_bases = 0;
my $dup_reads = 0;
my $read_legth = 0;

 #go through file line by line
while ( my $l = <F>) {
	chomp $l;
	next if($l !~ /^SN/);
	if ($l =~ /total length/) {my @k = split( /\t/, $l ); $total_bases = $k[2]; }
	if ($l =~ /cigar/ && $l =~ /bases mapped/ && $l !~ /mismatches/) {my @k = split( /\t/, $l );  $mapped_bases = $k[2]; }
	if ($l =~ /reads duplicated/) {my @k = split( /\t/, $l );$dup_reads = $k[2];}
	if ($l =~ /average length/) {my @k = split( /\t/, $l );  $read_legth = $k[2];}
	if ($l =~ /First Fragment Qualitites/) { last; }
}
close( F );

my $four = $total_bases-($dup_reads*$read_legth);
my $five = ($total_bases-($dup_reads*$read_legth))/12000000; 
print("$sample_name\n$total_bases\n$mapped_bases\n$four\n$five\n\n");
