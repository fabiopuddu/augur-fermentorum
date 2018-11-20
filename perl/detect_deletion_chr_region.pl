#!/usr/bin/env perl

# Author:		Mareike Herzog
# Maintainer:	Fabio Puddu
# Created: 	July 2014
# Description:
# Usage: 		perl detect_deletion_chr_region.pl <gene> <bam_file>
	
use Carp;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Data::Dumper;

#Get the bam file as input
my $gene_name = shift;
my $input = shift;

#Get the path of the current script to localise gene_name list file
my $path = __FILE__;
my @path_components = split( /\//, $path );
my $rem = pop @path_components; 
$rem = pop @path_components;
$path = (join "/",  @path_components);
my $gene_file="$path".'/defaults/all_yeast_genes.txt';

#Exit if the input file does not exist
if ( ! -e $input) {print "Could not find input file\n";  exit};

#This file contains information like this:
#Ensembl Gene ID	Chromosome Name	Gene Start (bp)	Gene End (bp)	Associated Gene Name
#YHR055C		VIII		214533		214718		CUP1-2

my $gene_oi; my $chrom; my $st; my $en;

my $ifi;
# Open the Conversion file and make an array for the location of the gene:
if( $gene_file =~ /\.gz/ ){open($ifi, qq[gunzip -c $gene_file|]);}else{open($ifi, $gene_file ) or die $!;}

#go through file line by line
while( my $line = <$ifi> ) {
#        	chomp $line;	
          next if $line =~ /^#/ ; #header lines skipped
          if ($line =~ /$gene_name\s/){
               	my($f1, $f2, $f3, $f4, $f5) = split '\t', $line; #split the line into its columns
                	my $f6 = $f4 - $f3 +1;
                	$gene_oi = $f1; $chrom = $f2; $st = $f3; $en = $f4;
          }          
}
close( $ifi );

#Get the chromosome-wide median
my @mp_out = `samtools view -b $input -F 0x0400 \'$chrom\' | genomeCoverageBed -dz -ibam stdin -g  `; 
my @cwm_coverage;
for my $line(@mp_out){
	chomp $line;
	my($s0, $s1, $s2) = split '\t', $line;
	push @cwm_coverage, $s2;
	
}
my $CWM=median(@cwm_coverage);

my $threshold=$CWM * 0.05;

#Get the coverage of the gene of interest:

my $range_goi = "$chrom".':'."$st".'-'."$en";
#print("$range_goi\n");
@mp_out = `samtools view -@ 8 -b $input -F 0x0400 \'$range_goi\'| genomeCoverageBed -dz -ibam stdin -g  `; 
#Get the coverage of the gene of interest:
my %cov;
for my $line(@mp_out){
	chomp $line;
	my($s0, $s1, $s2) = split '\t', $line;
	$cov{$s1}=$s2;
}
my @coverage;
my $cons_bas=0;
my $longest_cons_bas=0;
my $prev_cov=1000;

print Dumper \%cov;
for (my $i=$st; $i<=$en; $i++){
		my $s2=$cov{$i} || 0; # This is required because genomeCoverageBed -dz does not output positions with zero coverage
		print "$s2 ";
                push @coverage, $s2;	
		#calculate if there are gaps in the coverage that are smaller than the whole gene
		if ($s2 < $threshold and $prev_cov < $threshold){
			$cons_bas++;
			$longest_cons_bas=$cons_bas if $cons_bas>$longest_cons_bas;
		}
		else{
			$cons_bas=0;
		}		
		$prev_cov=$s2;
}
#Calculate statistics
my $cov_of_interest = median(@coverage);
my $perc_cov = sprintf("%.1f", $cov_of_interest/$CWM*100);
$cov_of_interest = sprintf("%.1f", $cov_of_interest);

#Output the result
print ("\n$gene_oi\t$gene_name\t$cov_of_interest\t$perc_cov %\tDeleted:");
if ($perc_cov < 15) {print "Y\n";} elsif ($longest_cons_bas > 9) {print "Yp (p=partial)\n";} else {print "N\n";}


sub median{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
