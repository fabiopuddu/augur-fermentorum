#!/usr/bin/perl -w

#Author: Mareike Herzog	
#Date: 20.4.2017
#Name: rep_median_adjust_wt-plate_median.pl
#Description: A script written for the normalisation of repetitive DNA measurements by the sequencing lane they were on.
			# Based on rep_median_adjust_telomere_only.pl.
			# Recent changes include: it determines the plate containing the wild-type samples, calculate that plate's median and sets that as a target median for every plate.
			
			# Input files needed: rep.txt
			#SD1862b 128     12      10      39      0       0       0       0       81      878.978969256493        ERS001999       Del1000_YNL120C
			#SD1863b2        135     12      7       39      12      2       3       1       94      767.640538213103        ERS002002       Del1001_TOM70
			#SD1863b 127     12      8       39      11      2       3       1       85      834.366835889768        ERS002001       Del1001_TOM70
			
			# Input files needed: name\ conversion.tsv
			#3792STDY6185637 Del1_TDA8       (C001)DN414641I TDA8    SD0863b         ERS000001
			#3792STDY6185638 Del1_TDA8       (C001)DN414641I TDA8    SD0863b2                ERS000002
			#3792STDY6185639 Del2_SCS22      (C001)DN414641I SCS22   SD0864b         ERS000003

			#You will also need to pass a column number to the script to indicate which repeat number you want to adjust
			#Column 1: SD_number (don't pass that one); Column 2: rDNA repeats; Column 3: CUP1 repeats;
			#Column 4: mtDNA copies; Column 5: Ty1 repeats; Column 6: Ty2 repeats;
			#Column 7: Ty3 copies; Column 8: Ty4 repeats; Column 9: Ty5 repeats;
			#Column 10: genome wide median; Column 11: Telomere repeats; Column 12: ERS numbers; Column 13: sample name;
			
			#THE SCRIPT WILL NOT WORK IF THE COLUMNS ARE NOT IN THE CORRECT ORDER


#sample usage: perl rep_median_adjust_wt-plate_median.pl rep.txt name\ conversion.tsv *column_no*

#Test command: perl rep_median_adjust_wt-plate_median.pl rep_example.txt name\ conversion.tsv 11

use strict;
use warnings;

#Load functions needed for median calculations
use POSIX qw(ceil);
use List::Util qw(sum);

######################
######################
###	SCRIPT  LOGIC  ###
######################
######################

#Step A: Find the wild-type plate  

#Step B: Get all the measurements & Compute the target median

#Step C: Get all the measurements plate by plate

#Step D: Compute the adjustment factor for each plate

#Step E: Replace the data line by line


###################
###################
### SCRIPT BODY ###
###################
###################

#############
# Get input #
#############

my $rep_file = shift; #results file
my $plate_file = shift; #name conversion
my $column = shift; # which measure you want to fix

#####################
# Declare variables #
#####################

# Declare other variables 

my $target_median = ''; #The target median all plates should have after adjustments. It will be the median of the plate containing the WT strains
my @plates; #an array containing the name of all plates
my %factors; #a hash to store plate => factor


##########
# Step A #
##########

######
##Find the wild-type plate 

#This command executes an external command to quickly search the plate_file for any line containing the word "WT-" and prints the third column which is the plate. This only works if all wt strains were on the same plate and if they are the only samples with a name containing 'WT-'
my $command = "cat '$plate_file' | grep 'WT-' | awk '{print ".'$3'."}' | sort | uniq | grep 'SD' -v";
my $wt_plate = `$command`;
chomp $wt_plate;

##########
# Step B #
##########

######
## Get all the measurements & Compute the target median

#For the $wt_plate get all the sample names into an array e.g. SD followed by number and b followed by a number e.g. SD5584b2
my @samples = @{get_SD_numbers($wt_plate)}; 

#Get all measurements for the column in question for those particular SD samples into an array
my @wt_measurements; #array to store all the wild type measurements
#the subroutine to extract the numbers needs to be passed the sample names and the column in question 
@wt_measurements = @{mes_per_plate(\@samples, $column)}; 
#remove any empty lines
@wt_measurements = grep { $_ ne '' } @wt_measurements;
#Compute the median of the plate containing the wild types
$target_median = median (@wt_measurements);


##########
# Step C #
##########

######
# Get a list of all plates
my $command2 = "cat '$plate_file' | awk '{print ".'$3'."}' | sort | uniq | grep 'SD' -v";
@plates = `$command2`;
chomp @plates;

######
# Get all the measurements plate by plate to determine the factor

foreach my $plate (@plates){
	#Skip the mock plates
	next if $plate =~ /mock/;
	#For each plate get the SD numbers
	my @samples_per_plate = @{get_SD_numbers($plate)}; 
	#Go sample by sample to get the measurement in question for each sample
	#the mes_per_plate function needs to be passed references to two arrays (input and output) as well as the column of the repeats file we are dealing with
	my @measurements = @{mes_per_plate(\@samples_per_plate, $column)}; 
	@measurements = grep { $_ ne '' } @measurements;
		
	##########
	# Step D #
	##########
	
	######
	# Compute the adjustment factor for each plate
	
	#Skip the plate if it has not been sequenced yet or is in the results
	my $length = scalar @measurements;
	next if $length == 0;
	#print "Plate: $plate; Length: $length; Samples:@samples_per_plate; Measures: @measurements\n\n";
	
	#Compute the median of all measurements on the plate
	my $plate_median = median(@measurements);
	#Compute the plate factor
	#print "Plate: $plate\tMedian:$plate_median\n";
	my $factor = $target_median/$plate_median;
	#store the factors determined for each plate in a hash
	$factors{$plate}=$factor;
	
}



##########
# Step E #
##########

######
# Replace the data line by line

#Step C: Replace the data line by line
my $fh;
open $fh, "<", $rep_file or die $!;  #open provided fastq file

while(my $row=<$fh>){ #read file or STDIN line by line
	chomp $row;
	
	#Split the line on the tab 
	my @tags = split "\t", $row;
	
	#Get the sample name and the value to be normalised
	my $index = $column - 1;
	my $sample = $tags[0];
	my $old_measure = $tags[$index];

	#Get the plate number
	my $plate = get_plate_number($sample);
	
	#Get the plate factor
	my $factor = $factors{$plate};
	
	#Get the adjusted measurement
	my $adj_measure = $old_measure * $factor;
	
	#Replace value and print the new line 
	$tags[$index] = $adj_measure;
	#print "$row\n"; #if one wants to compare the row as it was
	print join("\t",@tags),"\n";
	
	
	}
close $fh;



######################
######################
###	SUB  ROUTINES  ###
######################
######################

sub get_SD_numbers {
	my $pl = shift; #This subroutines first input variable is the plate number in question
	my $command3 = "cat '$plate_file' | grep '$pl' | awk '{print ".'$5'."}'";
	my @SDs = `$command3`;
	chomp @SDs;
	return \@SDs; #just return the reference to the array
}

sub mes_per_plate {
	#This subroutine needs to be passed an array reference followed by the column of interest 
	my @out_array = '';
	my $array_ref = shift;
	my $column_number = shift; #The column number is just a number like 11
	my $column = '$'.$column_number; #In the awk command that is to follow it needs to be preceeded by a $
	for my $samp (@{$array_ref}){
		my $command4 = "cat '$rep_file' | grep -w '$samp' | awk '{print ".$column."}'";
		my @ts_per_sample = `$command4`;
		chomp @ts_per_sample;
		if (defined $ts_per_sample[0] and $ts_per_sample[0] =~ /[0-9]+/ ) {push @out_array, $ts_per_sample[0]};
	}
	return \@out_array;
}

sub median {
  my $med = sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
  return $med;
}

sub get_plate_number {
	my $sample = shift;
	my $command2 = "cat '$plate_file' | grep -w $sample | awk '{print ".'$3'."}'";
    my @results=`$command2`;
    next if !defined $results[0];
    #my $plate = substr $results[0], 1, 4;
    my $plate = $results[0];
    chomp $plate;
	return $plate;
}

