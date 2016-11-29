#!/usr/bin/env perl

# Author:       mh23
# Maintainer:   mh23
# Created: 		July 2014
# Name:			rDNA_cov_extract.pl
# Test: 		perl /nfs/users/nfs_m/mh23/Scripts/rDNA_attempt.pl <bam_file>
#
# Two control files:
#   /lustre/scratch105/projects/SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants_2/REL-140723/ERP001366/SC_MFY5784126/10332485/13416_1#33/13416_1#33.bam
#   /lustre/scratch105/projects/SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants_2/REL-140723/ERP001366/SC_MFY5784126/10332485/13416_2#33/13416_2#33.bam
#   /lustre/scratch105/projects/SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants_2/REL-140723/ERP001366/SC_MFY5784127/10332497/13416_1#34/13416_1#34.bam
#   /lustre/scratch105/projects/SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants_2/REL-140723/ERP001366/SC_MFY5784127/10332497/13416_2#34/13416_2#34.bam
#
# Command to generate files
# REL-140723
# > for x in /lustre/scratch105/projects/SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants_2/REL-140723/ERP*/*/*/*/*.bam; do n=$(echo $x | sed 's/^.*SC/SC/'); m=$(echo $n | sed 's/\(SC_MFY[0-9]*\).*/\1/g'); echo $m; bsub-1gb -o sub.o -e sub.e "perl /nfs/users/nfs_m/mh23/Scripts/rDNA_attempt.pl -i $x > $m.wig.XII.rDNA"; done
# REL-140401
# > for x in /lustre/scratch105/projects/SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants_2/REL-140401/ERP001366/SC*/*.bam; do n=$(echo $x | sed 's/^.*SC/SC/'); m=$(echo $n | sed 's/\(SC_MFY[0-9]*\).*/\1/g'); echo $m; bsub-1gb -o sub.o -e sub.e "perl /nfs/users/nfs_m/mh23/Scripts/rDNA_attempt.pl -i $x > $m.wig.XII.rDNA"; done
#
# Description:
# This script reads in a S. cerevisiae Bam file and the employes samtools mpileup to make a vcf file
# for the region "XII     440000  495947". Nucleotide position and coverage are kept in an array. 
# The coverage is averaged, 25 at a time and position and average coverage is given as output. 


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
my ($input);

GetOptions
(
'i|input=s'         => \$input,
);

( $input && -f $input ) or die qq[Usage: $0 -i <input vcf>\n];


# define the array to store the samtools mpileup output
my @mpileup_out;

# execute samtools mpileup as a systems command and store the output in an array
# /nfs/users/nfs_m/mh23/rDNA_region.bed contains only this: XII     440000  495947
@mpileup_out =  `samtools mpileup  -A -Q 0 -r XII:440000-495947 -f ~/sw/bin/PF/mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa $input 2>/dev/null`; 

# declare variables
my %mpileup_hash;
my @positions;
my @values;
 
#go through the mpileup array and extract the position and coverage 
foreach my $l (@mpileup_out){
        chomp( $l );    
        my @s = split( /\t/, $l );
        #the position is $s[1] & coverage is $s[3]
        #make it into a hash and arrays
        $mpileup_hash{$s[1]} = $s[3];
        push(@positions, $s[1]); 
        push(@values, $s[3]); 
 }
 

#go through the hash and take the mean of each 25
#get length of array to determine how long the for loop should be
my $length_ar = scalar @positions;
my $length = $length_ar/25;

#declare variable 
my %output;

#calculate the average every 25 
for (my $i=0; $i < $length; $i++) {
        #get a subset of 25 values at a time
        my @mean_values=();
        @mean_values = splice (@values, 0, 25);
        my @mean_positions=();
        @mean_positions = splice (@positions, 0, 25);
        
        #get the first position from the list of 25
        my $pos = shift @mean_positions;
        
        #get the mean of the 25 coverage values
        my $mean = eval(join("+", @mean_values)) / @mean_values;
        
        #put the output into a hash
        $output{$pos}=$mean;

 } 
 
 #print output of the hash
 #position and coverage average
print "$_ $output{$_}\n" for (sort { $a cmp $b } keys %output);
