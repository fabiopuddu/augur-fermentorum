#!/usr/bin/env perl
#Author: 	Mareike Herzog
#Maintainer: 	Fabio Puddu
#Created: 	Apr 2017
#Description: 	A script written for the normalisation of repetitive DNA measurements by the sequencing lane they were on.
# 		it determines the plate containing the wild-type samples, calculate that plate's median and sets that as a target median for every plate.
# 		Input files needed: rep.txt
# 		Input files needed: name\ conversion.tsv
#
#	THE SCRIPT WILL NOT WORK IF THE COLUMNS ARE NOT IN THE CORRECT ORDER
#sample usage: 	perl rep_median_adjust_wt-plate_median.pl rep.txt name_conversion.tsv *column_no*
#Test command: 	perl rep_median_adjust_wt-plate_median.pl rep_example.txt name_conversion.tsv 11
use strict;
use warnings;

#Load functions needed for median calculations
use POSIX qw(ceil);
use List::Util qw(sum);
use Data::Dumper qw(Dumper);

##########################
##########################
###	SCRIPT  LOGIC  ###
##########################
##########################
#Step A: Read in data into a hash of hashes
#Step B: Find the wild-type plate , Get all the measurements & Compute the target median
#Step C: Get all the measurements plate by plate
#Step D: Compute the adjustment factor for each plate
#Step E: Replace the data line by line

#############
# Get input #
#############

my $rep_file = shift; #results file
my $plate_file = shift; #name conversion

#which measure you want to fix (columns in the rep.txt file)
my @columns = (2,3,4,5,6,7,8,9,10,13);
#print "@columns\n";

#####################
# Declare variables #
#####################
my $target_median = ''; 	#The target median all plates should have after adjustments. It will be the median of the plate containing the WT strains
my @plates; 			#an array containing the name of all plates
my %factors; 			#a hash to store plate => factor
my %data;

################
# Read in data #
################
my @header;
my $col_no;

#reading info in repDNA file
print STDERR "Reading repDNA file\n";
my $fh;
open $fh, "<", $rep_file or die $!;
chomp(my @REP_FILE=<$fh>);
close $fh;

#reading info in name conversion file
print STDERR "Reading name conversion file\n";
open $fh, "<", $plate_file or die $!;
chomp(my @PLATE_FILE=<$fh>);
close $fh;
@PLATE_FILE=grep { $_ !~ /mock/ } @PLATE_FILE; #remove odd looking data

#Process  repDNA file
print STDERR "Processing repDNA file:\n ";
my $tot=scalar @REP_FILE;
@header=split "\t",$REP_FILE[0];
$col_no = scalar @header;
my $i=0;
for my $row(@REP_FILE){ 			#read file or STDIN line by lin
	$i++;
	next if $i==1;
	my @file_tags = split "\t", $row;	#Split the line on the tab
	print STDERR "\r $i / $tot";
	for (my $i=0; $i <$col_no; $i++) {	#Set up the hash structure
		#$data{$file_tags[0]}{$header[$i]}=$file_tags[$i]; #this writes the name of teh column into the key
		my $c_n = $i + 1;
		$data{$file_tags[0]}{$c_n}=$file_tags[$i];
	}
}
print STDERR "\n";
#print Dumper \%data;
#exit;

###############################
# Start the correction process #
###############################

##########
# Step A #
##########
##Find the wild-type plate
# This searches the plate_file for any line containing the word "WT-" and prints the third column which is the plate.
#This only works if all wt strains were on the same plate and if they are the only samples with a name containing 'WT-'
my @wt_samples = grep {/WT-/} @PLATE_FILE;
@wt_samples=split ("\t", $wt_samples[0]);
my $wt_plate=$wt_samples[2];

# Get a list of all plates
for my $line(@PLATE_FILE){
        my @riga=split "\t", $line;
	       push @plates, $riga[2] if not $riga[2] ~~ @plates and not $riga[2] eq '';
}

#For each column we want to correct
foreach my $column (@columns){
	print STDERR "Working on column $column\n";
	##########
	# Step B #
	##########
	## Get all the measurements & Compute the target median
	my @samples = @{get_SD_numbers($wt_plate)};			#For the $wt_plate get all the sample names (e.g. SD5584b2) into an array
	my @wt_measurements = @{mes_per_plate(\@samples, $column)};	#Get all measurements for the column in question for those particular wild-type SD numbers into an array
	@wt_measurements = grep { $_ ne '' } @wt_measurements;		#remove any empty lines
	$target_median = median (@wt_measurements);			#Compute the median of the plate containing the wild types
	# Get all the measurements plate by plate to determine the factor
	foreach my $plate (@plates){
		next if $plate =~ /mock/;						#Skip the mock plates
		##########
		# Step C #
		##########
		my @samples_per_plate = @{get_SD_numbers($plate)};			#For each plate get the SD numbers it contains
		my @measurements = @{mes_per_plate(\@samples_per_plate, $column)};	#Get all measurements for the column in question for that list of SD numbers
		@measurements = grep { $_ ne '' } @measurements;			#remove any empty lines

		##########
		# Step D #
		##########
		# Compute the adjustment factor for each plate
		next if scalar @measurements == 0;			#Skip that plate if it contains no data
		#print "Plate: $plate; Length: $length; Samples:@samples_per_plate; Measures: @measurements\n\n";
		my $plate_median = median(@measurements);		#Compute the median of all measurements on the plate
		#print "Plate: $plate\tMedian:$plate_median\n";
		my $factor = $target_median/$plate_median;		#Compute the plate factor
		$factors{$plate}=$factor;				#store the factors determined for each plate in a hash
	}
	##########
	# Step E #
	##########
	# Replace the data line by line
	foreach my $sample (keys %data){
		#print "$sample\n";
		my $old_measure = $data{$sample}{$column};
		#print "Old_measure: $old_measure\n";
		my $plate = get_plate_number($sample);			#Identify which plate that sample was on
		#print "Plate: $plate\t";
		my $factor = $factors{$plate};				#Get the plate factor
		# print STDERR "$old_measure\t$sample\t$column\n";
		my $adj_measure = $old_measure * $factor;		#Get the adjusted measurement
		#print "New measure: $adj_measure\n";
		$data{$sample}{$column}=$adj_measure;
	}
}
#################
# Print results #
#################
print (join "\t", @header);
print "\n";
foreach my $strain (sort keys %data) {
	my @out_line=();
	for (my $i=1; $i <= $col_no; $i++) {
		push @out_line, $data{$strain}{$i};
	}
	print join "\t", @out_line;
	print "\n";
}

###########################
###########################
###	SUBROUTINES     ###
###########################
###########################
#This subroutines gets all the sample names (Sd numbers) present in the plate in question, passed as argument
sub get_SD_numbers {
	my $pl = shift;
	my @SDs;
	#chomp @SDs;
	#print "$pl\n";
	$pl =substr $pl, 6;
	my @samples=grep { $_ =~ /$pl/ } @PLATE_FILE;
	foreach my $sample(@samples){
		my @line=split ("\t", $sample);
		push @SDs, $line[4];
	}
	return \@SDs; 	# return the reference to the array
}
###########################
#This subroutine gets all the measures for a list of strains (passed an array reference) for a single genomic feature passed as argument
sub mes_per_plate {
	my @out_array = '';
	my $array_ref = shift;	#List of strain names (SD numbers)
	my $column = shift; 	#The column number
	for my $samp (@{$array_ref}){
		my $ts_per_sample=$data{$samp}{$column} if exists $data{$samp} and defined $data{$samp}{$column};
		if (defined $ts_per_sample and $ts_per_sample =~ /[0-9]+/ ) {push @out_array, $ts_per_sample};
	}
	return \@out_array;
}
###########################
#This subroutine calculates a median
sub median {
  my $med = sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
  return $med;
}
###########################
#This subroutine returns the plate number corresponding to a sample
sub get_plate_number {
	my $sample = shift;
	my $line=(grep { $_ =~ /\b$sample\b/ } @PLATE_FILE)[0];
	my $results = (split "\t", $line)[2];
#       next if !defined $results;
	#my $plate = substr $results[0], 1, 4;
	my $plate = $results;
	chomp $plate;
	return $plate;
}
###########################
