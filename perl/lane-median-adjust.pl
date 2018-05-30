#!/usr/bin/env perl

#Author: 	Mareike Herzog	
#Maintainer: 	Fabio Puddu
#Created: 		Apr 2017
#Description: 	A script written for the normalisation of repetitive DNA measurements by the sequencing lane they were on.
# 		it determines the plate containing the wild-type samples, calculate that plate's median and sets that as a target median for every plate.
# 		Input files needed: rep.txt/all_combined.tsv
#			SDname rDNA    CUP1    mitochondria    Ty1     Ty2     Ty3     Ty4     Ty5     GWM     Telomeres       ERSno   Deletion        chr01   chr02   chr03   chr04   chr05   chr06   chr07   chr08   chr09   chr10   chr11   chr12   chr13   chr14   chr15   chr16   AneupNumber
#			SD0863b 118.085106382979        11.468085106383 12.8936170212766        39.9787234042553        12.0531914893617        1.90425531914894        3.1063829787234 1.11702127659574        94      614.23391551713 ERS000001       Del1_TDA8       2       2       2       2       2       2       2       2       2       2       2       2       2       2       2       2       0
#			SD0863b2        122.652631578947        12.1052631578947        12.3368421052632        39.7842105263158        12.4210526315789        1.93684210526316        3.23157894736842        0.989473684210526       95      598.304236757884        ERS000002       Del1_TDA8       2       2       2       2       2       2       2       2       2       2       2       2       2       2       2       2       0				
# 		Input files needed: name\ conversion.tsv
#			3792STDY6185637 Del1_TDA8       (C001)DN414641I TDA8    SD0863b         ERS000001
#			3792STDY6185638 Del1_TDA8       (C001)DN414641I TDA8    SD0863b2                ERS000002
#			3792STDY6185639 Del2_SCS22      (C001)DN414641I SCS22   SD0864b         ERS000003
#			3792STDY6185640 Del2_SCS22      (C001)DN414641I SCS22   SD0864b2                ERS000004
#			3792STDY6185641 Del3_SDH8       (C001)DN414641I SDH8    SD0865b         ERS000005
#		You will also need to pass a column number to the script to indicate which repeat number you want to adjust
#			Column 1: SD_number (don't pass that one); Column 2: rDNA repeats; Column 3: CUP1 repeats;
#			Column 4: mtDNA copies; Column 5: Ty1 repeats; Column 6: Ty2 repeats;
#			Column 7: Ty3 copies; Column 8: Ty4 repeats; Column 9: Ty5 repeats;
#			Column 10: genome wide median; Column 11: Telomere repeats; Column 12: ERS numbers; Column 13: sample name/Deletion;
#			Column 14: chr01   Column 15: chr02   Column 16:chr03   Column 17: chr04   Column 18: chr05   
#			Column 19: chr06   Column 20: chr07   Column 21: chr08  Column 22: chr09   Column 23: chr10   
#			Column 24: chr11   Column 25: chr12   Column 26: chr13  Column 27: chr14   Column 28: chr15   
#			Column 29: chr16   Column 30: AneupNumber
#	THE SCRIPT WILL NOT WORK IF THE COLUMNS ARE NOT IN THE CORRECT ORDER
#sample usage: 	perl rep_median_adjust_wt-plate_median.pl rep.txt name\ conversion.tsv *column_no*
#Test command: 	perl rep_median_adjust_wt-plate_median.pl rep_example.txt name\ conversion.tsv 11

use strict;
use warnings;

#Load functions needed for median calculations
use POSIX qw(ceil);
use List::Util qw(sum);
use Data::Dumper qw(Dumper);


######################
######################
###	SCRIPT  LOGIC  ###
######################
######################


#Step A: Read in data into a hash of hashes 

#Step B: Find the wild-type plate , Get all the measurements & Compute the target median

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

#which measure you want to fix

my @columns = (2,3,4,5,6,7,8,9,10,13);


#print "@columns\n";

#exit;

#####################
# Declare variables #
#####################

# Declare other variables 

my $target_median = ''; #The target median all plates should have after adjustments. It will be the median of the plate containing the WT strains
my @plates; #an array containing the name of all plates
my %factors; #a hash to store plate => factor
my %data;

################
# Read in data #
################

my @header;
my $col_no;
print STDERR "Reading repDNA file\n";
my $fh;
open $fh, "<", $rep_file or die $!;  #open provided fastq file
chomp(my @REP_FILE=<$fh>);
close $fh;

print STDERR "Reading name conversion file\n";
open $fh, "<", $plate_file or die $!;  #open provided fastq file
chomp(my @PLATE_FILE=<$fh>);
close $fh;

@PLATE_FILE=grep { $_ !~ /mock/ } @PLATE_FILE;

print STDERR "Processing repDNA file:\n ";
my $tot=scalar @REP_FILE;
@header=split "\t",$REP_FILE[0];
$col_no = scalar @header;
my $i=0;
for my $row(@REP_FILE){ #read file or STDIN line by lin
	$i++;
	next if $i==1;
	#Split the line on the tab 
	my @file_tags = split "\t", $row;
	print STDERR "\r $i / $tot";
	#Set up the hash structure
	for (my $i=0; $i <$col_no; $i++) {
		#$data{$file_tags[0]}{$header[$i]}=$file_tags[$i]; #this writes the name of teh column into the key
		my $c_n = $i + 1;
		$data{$file_tags[0]}{$c_n}=$file_tags[$i];
	}
}

print STDERR "\n";
#print Dumper \%data;
#exit; 

###############################
# Start the iterative process #
###############################

##########
# Step A #
##########

######
##Find the wild-type plate 

#This command executes an external command to quickly search the plate_file for any line containing the word "WT-" and prints the third column which is the plate. This only works if all wt strains were on the same plate and if they are the only samples with a name containing 'WT-'
my @wt_samples=grep { $_ =~ /WT-/ } @PLATE_FILE;
@wt_samples=split ("\t", $wt_samples[0]);
my $wt_plate=$wt_samples[2];



######
# Get a list of all plates
my @plates;
for my $line(@PLATE_FILE){
        my @riga=split "\t", $line;
        push @plates, $riga[2] if not $riga[2] ~~ @plates and not $riga[2] eq '';
}

##########
# Step B #
##########



foreach my $column (@columns){

print STDERR "Working on column $column\n";

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

#Step C: Replace the data in the hash of hashes

	foreach my $sample (keys %data){

		#print "$sample\n";		
		my $old_measure = $data{$sample}{$column};
		#print "Old_measure: $old_measure\n";
		
		my $plate = get_plate_number($sample);
		#print "Plate: $plate\t";
		
		#Get the plate factor
		my $factor = $factors{$plate};
			
	
		#Get the adjusted measurement
#		print STDERR "$old_measure\t$sample\t$column\n";
		my $adj_measure = $old_measure * $factor;
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
	for (my $i=1; $i <= $col_no; $i++) {
	 	print "$data{$strain}{$i}\t";
	}
	print "\n";
}


######################
######################
###	SUB  ROUTINES  ###
######################
######################

sub get_SD_numbers {
	my $pl = shift; #This subroutines first input variable is the plate number in question
	#my $command3 = "cat '$plate_file' | grep '$pl' | awk '{print ".'$5'."}'";
	my @SDs;
	#chomp @SDs;
	#print "$pl\n";
	$pl =substr $pl, 6;
	my @samples=grep { $_ =~ /$pl/ } @PLATE_FILE;
	foreach my $sample(@samples){
		my @line=split ("\t", $sample);
		push @SDs, $line[4];
	}
	return \@SDs; #just return the reference to the array
}

sub mes_per_plate {
	#This subroutine needs to be passed an array reference followed by the column of interest 
	my @out_array = '';
	my $array_ref = shift;
	my $column = shift; #The column number is just a number like 11
	#my $column = '$'.$column_number; #In the awk command that is to follow it needs to be preceeded by a $
	for my $samp (@{$array_ref}){
		#my $command4 = "cat '$rep_file' | grep -w '$samp' | awk '{print ".$column."}'";
		#my @ts_per_sample = `$command4`;
		#chomp @ts_per_sample;
		#my row=(grep { $_ == /$samp/ } @REP_FILE)[0];
		my $ts_per_sample=$data{$samp}{$column} if exists $data{$samp} and defined $data{$samp}{$column};
		if (defined $ts_per_sample and $ts_per_sample =~ /[0-9]+/ ) {push @out_array, $ts_per_sample};

	}
	return \@out_array;
}

sub median {
  my $med = sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
  return $med;
}

sub get_plate_number {
	my $sample = shift;
	my $line=(grep { $_ =~ /$sample/ } @PLATE_FILE)[0];
	#my $command2 = "cat '$plate_file' | grep -w $sample | awk '{print ".'$3'."}'";
    	#my @results=`$command2`;
	my $results = (split "\t", $line)[2];
        next if !defined $results;
    #my $plate = substr $results[0], 1, 4;
    my $plate = $results;
    chomp $plate;
	return $plate;
}

