#!/usr/bin/env perl
# 
# Author:       mh23
# Maintainer:   mh23
# Created: 		11.02.2015

#Test: perl /nfs/users/nfs_m/mh23/Scripts/filter_exons.pl - i<vcf>
    

use Carp;
use strict;
use warnings;
use Getopt::Long;

#define the storage arrays for the exon start and end nucleotide locations for each chromosome
my @chr1_s; my @chr1_e;
my @chr2_s; my @chr2_e;
my @chr3_s; my @chr3_e;
my @chr4_s; my @chr4_e;
my @chr5_s; my @chr5_e;
my @chr6_s; my @chr6_e;
my @chr7_s; my @chr7_e;
my @chr8_s; my @chr8_e;
my @chr9_s; my @chr9_e;
my @chr10_s; my @chr10_e;
my @chr11_s; my @chr11_e;
my @chr12_s; my @chr12_e;
my @chr13_s; my @chr13_e;
my @chr14_s; my @chr14_e;
my @chr15_s; my @chr15_e;
my @chr16_s; my @chr16_e;
my @chr17_s; my @chr17_e;
my @chr18_s; my @chr18_e;
my @chr19_s; my @chr19_e;
my @chrX_s; my @chrX_e;

#give location of the file that contains all the info about exon location: 3073253 3074322
my $conversion_file='/nfs/users/nfs_m/mh23/ReferenceGenomes/Mus_musculus.GRCm38.78.gtf.gz';

my ($input);

GetOptions
(
'i|input=s'         => \$input,
);

( $input && -f $input ) or die qq[Usage: $0 -i <input vcf>\n];


my $ifi;
# Open the Conversion file and make arrays for the exonic locations:
if( $conversion_file =~ /\.gz/ ){open($ifi, qq[gunzip -c $conversion_file|]);}else{open($ifi, $conversion_file ) or die $!;}
        #go through file line by line
        while( my $line = <$ifi> ) {
                next if ($line =~ /^#/); #header lines skipped
                #split the line into its columns
                #Example:
                #1       havana  gene    3073253 3074322 .       +       .       gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "RP23-271O17.1"; gene_source "havana"; gene_biotype "TEC";
                my $ln = $line;
                my($f1, $f2, $f3, $f4, $f5, $f6, $f7, $f8) = split '\t', $line;
        		# make specfic arrays for each chromosome
        		#f1 is the chromosome, f4 the start and f5 the end of an exon
        		if ($f1 eq '1') {push @chr1_s, $f4 ; push @chr1_e,$f5; $f1=0; $f4=0; $f5=0;}
        		elsif ($f1 eq '2') {push @chr2_s, $f4 ; push @chr2_e,$f5; $f1=0; $f4=0; $f5=0;}
        		elsif ($f1 eq 3) {push @chr3_s, $f4 ; push @chr3_e, $f5;}
        		elsif ($f1 eq 4) {push @chr4_s, $f4 ; push @chr4_e, $f5;}
        		elsif  ($f1 eq 5) {push @chr5_s, $f4 ; push @chr5_e, $f5;}
        		elsif  ($f1 eq 6) {push @chr6_s, $f4 ; push @chr6_e, $f5;}
        		elsif ($f1 eq 7) {push @chr7_s, $f4 ; push @chr7_e, $f5;}
        		elsif ($f1 eq 8) {push @chr8_s, $f4 ; push @chr8_e, $f5;}
        		elsif ($f1 eq 9) {push @chr9_s, $f4 ; push @chr9_e, $f5;}
        		elsif ($f1 eq 10) {push @chr10_s, $f4 ; push @chr10_e, $f5;}
        		elsif ($f1 eq 11) {push @chr11_s, $f4 ; push @chr11_e, $f5;}
        		elsif ($f1 eq 12) {push @chr12_s, $f4 ; push @chr12_e, $f5;}
        		elsif ($f1 eq 13) {push @chr13_s, $f4 ; push @chr13_e, $f5;}
        		elsif ($f1 eq 14) {push @chr14_s, $f4 ; push @chr14_e, $f5;}
        		elsif ($f1 eq 15) {push @chr15_s, $f4 ; push @chr15_e, $f5;}
        		elsif ($f1 eq 16) {push @chr16_s, $f4 ; push @chr16_e, $f5;}
        		elsif ($f1 eq 17) {push @chr17_s, $f4 ; push @chr17_e, $f5;}
        		elsif ($f1 eq 18) {push @chr18_s, $f4 ; push @chr18_e, $f5;}
        		elsif ($f1 eq 19) {push @chr19_s, $f4 ; push @chr19_e, $f5;}
        		elsif ($f1 eq 'X') {push @chrX_s, $f4 ; push @chrX_e, $f5;}
        		#make sure that the variables are emptied at the end of each cycle
        		$f1=0; $f4=0; $f5=0;
        }# close while
close( $ifi );

my $ifh;
if( $input =~ /\.gz/ ){open($ifh, qq[gunzip -c $input|]);}else{open($ifh, $input ) or die $!;}

while( my $l = <$ifh> )
{
    
    chomp( $l );
    my $f = $l;
    if ( $l =~ /^#/ ) {print("$f\n");}
    next if( $l =~ /^#/ || $l =~ /^#CHROM/);
    #k1 is the chromosome and k2 is the position
    my($k1, $k2) = split '\t', $l;
    if ($k1 eq '1') {for(my $i=0;$i<@chr1_s;$i++){if ($k2 >= $chr1_s[$i] && $k2 <= $chr1_e[$i]) {print("$f\n");last;}   }}
    if ($k1 eq '2') {for(my $i=0;$i<@chr2_s;$i++){if ($k2 >= $chr2_s[$i] && $k2 <= $chr2_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '3') {for(my $i=0;$i<@chr3_s;$i++){if ($k2 >= $chr3_s[$i] && $k2 <= $chr3_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '4') {for(my $i=0;$i<@chr4_s;$i++){if ($k2 >= $chr4_s[$i] && $k2 <= $chr4_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '5') {for(my $i=0;$i<@chr5_s;$i++){if ($k2 >= $chr5_s[$i] && $k2 <= $chr5_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '6') {for(my $i=0;$i<@chr6_s;$i++){if ($k2 >= $chr6_s[$i] && $k2 <= $chr6_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '7') {for(my $i=0;$i<@chr7_s;$i++){if ($k2 >= $chr7_s[$i] && $k2 <= $chr7_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '8') {for(my $i=0;$i<@chr8_s;$i++){if ($k2 >= $chr8_s[$i] && $k2 <= $chr8_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '9') {for(my $i=0;$i<@chr9_s;$i++){if ($k2 >= $chr9_s[$i] && $k2 <= $chr9_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '10') {for(my $i=0;$i<@chr10_s;$i++){if ($k2 >= $chr10_s[$i] && $k2 <= $chr10_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '11') {for(my $i=0;$i<@chr11_s;$i++){if ($k2 >= $chr11_s[$i] && $k2 <= $chr11_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '12') {for(my $i=0;$i<@chr12_s;$i++){if ($k2 >= $chr12_s[$i] && $k2 <= $chr12_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '13') {for(my $i=0;$i<@chr13_s;$i++){if ($k2 >= $chr13_s[$i] && $k2 <= $chr13_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '14') {for(my $i=0;$i<@chr14_s;$i++){if ($k2 >= $chr14_s[$i] && $k2 <= $chr14_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '15') {for(my $i=0;$i<@chr15_s;$i++){if ($k2 >= $chr15_s[$i] && $k2 <= $chr15_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '16') {for(my $i=0;$i<@chr16_s;$i++){if ($k2 >= $chr16_s[$i] && $k2 <= $chr16_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '17') {for(my $i=0;$i<@chr17_s;$i++){if ($k2 >= $chr17_s[$i] && $k2 <= $chr17_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '18') {for(my $i=0;$i<@chr18_s;$i++){if ($k2 >= $chr18_s[$i] && $k2 <= $chr18_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq '19') {for(my $i=0;$i<@chr19_s;$i++){if ($k2 >= $chr19_s[$i] && $k2 <= $chr19_e[$i]) {print("$f\n"); last;}   }}
    if ($k1 eq 'X') {for(my $i=0;$i<@chrX_s;$i++){if ($k2 >= $chrX_s[$i] && $k2 <= $chrX_e[$i]) {print("$f\n"); last;}   }}
    
    
    }
close( $ifh );

