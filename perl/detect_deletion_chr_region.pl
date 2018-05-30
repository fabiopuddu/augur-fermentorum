#!/usr/bin/env perl

# Author:		Mareike Herzog
# Maintainer:	Mareike Herzog
# Created: 	July 2014
# Description:
# Test: 		perl /nfs/users/nfs_m/mh23/Scripts/detect_deletion.pl <gene> <bam_file>
# Test: 		perl /lustre/scratch109/sanger/mh23/Scripts/detect_deletion_chr_region.pl YML061C SD0735b.bam
	
use Carp;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);



# define a subroutine to calculate the mean
sub mean {
    return sum(@_)/@_;
}

#Get the bam file as input
my $gene_name = shift;
my $input = shift;
#Get the path to say that the gene_name list is in the same folder
my $path = __FILE__;
my @path_components = split( /\//, $path );
my $rem = pop @path_components; 
$rem = pop @path_components;
$path = (join "/",  @path_components);
my $gene_file="$path".'/defaults/all_yeast_genes.txt';
#my $gene_file = '/lustre/scratch109/sanger/mh23/Scripts/all_yeast_genes.txt';


#Output if the file doe snot exist
if ( ! -e $input) {print "File missing\n";  exit};



my %results_hash;
my $results_file;


#This file contains information like this:
#Ensembl Gene ID	Chromosome Name	Gene Start (bp)	Gene End (bp)	Associated Gene Name
#YHR055C	VIII	214533	214718	CUP1-2

#gene of interest:
my $gene_oi; my $chrom; my $st; my $en;
my $risposta = 0;
my $linecounter = 0;
my $check=0; my $check_count=0;
my $length_of_gene=0;
my $ifi;
# Open the Conversion file and make arrays for the exonic location of the gene:
if( $gene_file =~ /\.gz/ ){open($ifi, qq[gunzip -c $gene_file|]);}else{open($ifi, $gene_file ) or die $!;}
        #go through file line by line
        while( my $line = <$ifi> ) {
                $linecounter++;
                next if ($linecounter == 1); #header lines skipped
                #we are only checking one gene
                if ($line =~ /$gene_name\s/){
                #split the line into its columns
                my $ln = $line;
                my($f1, $f2, $f3, $f4, $f5) = split '\t', $line;
                my $f6 = $f4 - $f3 +1;
                $gene_oi = $f1; $chrom = $f2; $st = $f3; $en = $f4;
                $length_of_gene = $en - $st +1;
                $check = 1;
                }#close if          
        }# close while
close( $ifi );

#Get the coverage of the gene of interest:
my $range_goi = "$chrom".':'."$st".'-'."$en";
#print("$range_goi\n");
my @mp_out = `samtools mpileup  -A -Q 0 -r $range_goi -f $path/mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa $input 2>/dev/null`; 
#Get the coverage of the gene of interest:
my $coverage=0; my $checker = 0; my $break = 0;
for (my $i =0; $i < @mp_out; $i++) {
	my $l = $mp_out[$i];
	my($s0, $s1, $s2, $s3, $s4) = split '\t', $l;	
	$coverage = $coverage + $s3;
	#we also want to check whether there are continous stretches without coverage
	#for that in the first line we will just read in the position
	if ($i == 0){$checker = $s1; next;}
	#in the next line we will calculate the difference between the previous position ($checker) and the current position
	my $diff = $s1 - $checker;
	#if the difference between the two is bigger than 9 (which means we have at least 9nt with 0 coverage)
	#we will set a variable from 0 to 1
	#at the end of the programme this will be used to decide whether the gene is deleted
	if ($diff > 9) {
	$risposta = 1;
	}
	#at the end of each line the position is updated
	$checker = $s1; 
}
my $cov_of_interest = $coverage / $length_of_gene;


#Get the genomic coverage of several 10kb bits of chromosome
my @ranges = (

'I:126501-136,500',
'II:371950-381,950',
'III:180329-190329',
'IV:673307-683307',
'V:224934-234934',
'VI:87738-97738',
'VII:344556-354556',
'VIII:304608-314608',
'IX:203999-213998',
'X:213998-223997',
'XI:223997-233996',
'XII:542742-552742',
'XIII:552742-562742',
'XIV:305008-315008',
'XV:558831-568831',
'XVI:391338-401338'

);

my @coverages;
# 1- loop through the array
foreach my $range (@ranges) {
#alternative: do it on a chunk of chromosome
	
	
	# 7- empty array withs samtools results
	# define the array to store the samtools mpileup output
	my @mpileup_out;

	# 3- use samtools mpileup to get the coverage for that position
	# execute samtools mpileup as a systems command and store the output in an array
	# 4- capture results in array
	@mpileup_out =  `samtools mpileup  -A -Q 0 -r $range -f $path/mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa $input 2>/dev/null`; 

	# 5- calculate average coverage of gene
	#The output comes in the format: 
	#chr is $s[0] & the position is $s[1] & coverage is $s[3]
	
	my $coverage=0; my $length_of_gene=0;
	foreach my $l (@mpileup_out) {
		my($s0, $s1, $s2, $s3, $s4) = split '\t', $l;	
		$length_of_gene++;
		$coverage = $coverage + $s3;
	}
	next if ($length_of_gene == 0);
	my $average_cov = $coverage/$length_of_gene;
	push @coverages, $average_cov;
	
	#print("$range\t$average_cov\n");
	
}


#11 - Determine relative coverage #
#adding up average coverages and dividing
my $num_o_cov=0;
my $control_sum_cov=0;
my $removed_cov_from_gene_of_interest = shift @coverages;
foreach my $cov (@coverages) {
	$num_o_cov++;
	$control_sum_cov = $control_sum_cov + $cov;
}

my $control_cov = $control_sum_cov / $num_o_cov; 

my $final_cov = sprintf("%.1f", $cov_of_interest/$control_cov*100);
$cov_of_interest = sprintf("%.1f", $cov_of_interest);
print ("$gene_oi\t$gene_name\t$cov_of_interest\t$final_cov %\tDeleted:");
if ($final_cov < 15) {print "Y\n";} elsif ($risposta == 1) {print "Yp (p=partial)\n";} else {print "N\n";}



