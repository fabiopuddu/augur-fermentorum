#!/usr/bin/env perl
# 
# Author:       mh23	
# Maintainer:   mh23
# Created: 		03.06.2015

# Test: perl /lustre/scratch109/sanger/mh23/Scripts/merge_bams_samtools.pl bams_for_mplieup
     

use Carp;
use strict;
use warnings;


################################
## Input the bams for mpileup ##
################################

my $input = shift;

################################
## Read the files into  array ##
## and get sample names		  ##
## and build the path back    ##
################################

my @files; my @samples; my @paths;
open(F, $input) or die ("Unable to open input file $input: $!\n");
#go through file line by line
while (my $line = <F>) {	
	chomp($line);
	my $path = ''; 
	# get the file paths into an array
	push(@files,$line);
	# split path on the tab
	my @elements = split '/', $line;
	# cycle through the elements of the path
	foreach my $ele (@elements){
		# if the element starts with SC_MFY or SD
		# put sample names into array
		if ($ele =~ /^SC_MFY/ || $ele =~ /^SD/) {
			push(@samples,$ele); last;
		}
		$path = $path.$ele.'/';
	}
	push(@paths,$path);
}
close F or die "Cannot close $input: $!\n";
#add a last "nonsense" element to the samples array. 
push(@samples, '0');



##############################
# Save the bams_for_mpileup ## 
##############################

system ("mv $input $input".'_old');

#################################
## Define new bams_for_mpileup ##
#################################
my $filename = "$input";

open (my $fh, '>', $filename) or die ("Could not open new file $input: $!");


###########################
## Generate the commands ##
###########################
my $command; my @commands;
#Loop through the sample names and check whether two consecutive samples are the same
my $loop_length = scalar @samples;
my $i=1; my $k; my $flag = 0;

# we loop through the samples names
# we start at the second sample and compare to the previous
while ($i < $loop_length) {
	$k = $i - 1;
	#if the two sample names are not the same we move on 
	if ($samples[$i] ne $samples[$k]) {
		#printing files to bams_for_mpileup
		if ($flag ==1) {
		 print $fh "$paths[$k]"."$samples[$k]".'/'."$samples[$k]".'.merged.bam'."\n"; #print the merged path into bams_for_mpileup
		}
		else {
		 print $fh "$files[$k]\n";
		}
		
		
		#if the samples are not the same, a new command will have to be generated
		#and we save the previously made command in the array
		if ($i != 1 && $flag ==1) { 
		push(@commands, $command);
		$flag = 0; 
		$command = ''; }
		$i++;		
	}
	# if they are we build commands
	else {
		# if a command for these samples has not been started, start generating a command
		if  ($flag == 0) {
			$command = "samtools merge "."$paths[$k]"."$samples[$k]".'/'."$samples[$k]".'.merged.bam '."$files[$k] $files[$i]";	
		}
		# if the command has been started we add the newest sample to the end
		else {
			$command = $command . " $files[$i]"; 	
		}
		
		$flag = 1;
		$i++;
	} 
}
close $fh;	

###############################
## Clean up bams_for_mpileup ##
###############################

system( "cat $input | grep '/' > temp");
system( "mv temp $input");

##########################
## Execute the commands ##
##########################

# we loop array of commands and execute them one by one
foreach my $c (@commands) {
	system("sbatch --wrap=\"$c\" ");	
}
	
