#!/usr/bin/env perl

# Author:       	Mareike Herzog	
# Maintainer:   	Mareike Herzog and Fabio Puddu
# Created: 	Nov 2016
# Rewritten: 	Aug 2017
# Description:

use autodie;
use utf8;
use Carp;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
use Pod::Usage;
use Cwd 'abs_path';
use POSIX qw(ceil);
use POSIX qw(floor);
use IPC::Open3;
use IO::File;
use Math::Complex qw(sqrt);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;


#######################################
##    								 ##
##				USAGE				 ##
##									 ##
#######################################

#Get the input and print usage if help of wrong input

#DEFAULTS
my $input = '0';
my $help = 0;

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions('help|?' => \$help, 
		   'i|input=s' => \$input,
		   ) or pod2usage(2);
pod2usage(1) if $help;
## If no input argument were given, then allow STDIN to be used only
## if it's not connected to a terminal (otherwise print usage)
pod2usage("$0: No input given.")  if (($input eq 0) && (-t STDIN));

#Check that input exists and has a size bigger than 0		
pod2usage("$0: File $input does not exist.")  unless ( -e $input);
pod2usage("$0: File $input is empty.")  if ( -z $input);

##############################################################
## Read the database converting systematic and common names ##
##############################################################

#Get the path of the current script to localise gene_name list file
my $path = __FILE__;
my @path_components = split( /\//, $path );
my $rem = pop @path_components; $rem = pop @path_components;
$path = (join "/",  @path_components);
my $gene_file="$path".'/defaults/all_yeast_genes.txt';

# Open the Conversion file and make an array for the location of the gene:
my $ifi;
if( $gene_file =~ /\.gz/ ){open($ifi, qq[gunzip -c $gene_file|]);}else{open($ifi, $gene_file ) or die $!;}

#go through file line by line
my %gene_db;
while( my $line = <$ifi> ) {
        	chomp $line;	
          next if $line =~ /^#/ ; #header lines skipped
	my($f1, $f2, $f3, $f4, @f5) = split '\t', $line; #split the line into its columns
	if ($f5[0]){
		$gene_db{$f5[0]}=$f1;	
	}
	else {
		$gene_db{$f1}=$f1;
	}
}
close( $ifi );


########################################
## Script logic and parameter setting ##
########################################

### Step1
# Get wild type stats and numbers (to identify the window of significance) 

#Define range of data columns
my $first_col = 2; 
my $last_col = 13; 

#Define the order of the measures
my @columns = qw[rDNA	CUP1	M	2mic	T1		T2		T3		T4		T5	gwm	mat	TEL x x x x x x x x x x x x x x x x x x ANEUP CR]; 

### Step 1a
# Compute the averages of wild types across all measures

### Step 1b
# Compute the standard error of wild types across all measures

### Step 2
# Compute the window of significance using the mean and the standard deviation

#Define the multiplier (the Stdev will be multiplied with this and this will be added/subtracted from the mean to identify the confidence interval) 
my $multiplier = 3;

### Step 3
# Loop through the file and identify which sample measures are outside of the bounds

### Step 3b
# Compile a list of all genes
my $gene_column = 15;

### Step 4
# Loop through all significant results and identify whether more than two measures for that gene are out of bounds

### Step 5 
# Print results to files



#######################
##  STEP1 - WT STATS ##
#######################

# The data columns are set in "Script logic and parameter setting"
# Save results in array and a hash
my @wt_means; my %wt_means_hash; 
my @wt_stdev; my %wt_stdev_hash;
my @wt_stderr; my %wt_stderr_hash;
### Step 1a
# Compute the averages of wild types across all measures
# Loop through all data columns and get the averages for the wild type in all columns
# Loop from the first data column to the last
for my $x (2 .. 15, 32, 33) {
	#Compute the average for the particular column
    my $average = wt_column_average($x);     
    #Push to array
    push @wt_means, $average; 
    #Save in hash with column name as key
    $wt_means_hash{$columns[$x-2]}=$average;   
### Step 1b
# Compute the standard deviation of wild types across all measures
	my $std_dev = wt_column_stdev ($x);
	#save stdev as well
	push @wt_stdev, $std_dev; 
	$wt_stdev_hash{$columns[$x-2]}=$std_dev;
	#save sterr as well
	my $std_err = wt_column_stderr($std_dev);
	push @wt_stderr, $std_err; 
	$wt_stderr_hash{$columns[$x-2]}=$std_err;
	
}


########################
##  STEP2 - Define CI ##
########################

# Compute the window of significance using the mean, the standard deviation and a multiplier
# Save bounds in this hash
my %lower_bounds; my %upper_bounds;
# Loop through the columns and identify the windows
for my $measure (@columns){
	#Multiply StDev by the multiplier
	my $interval = $multiplier * $wt_stdev_hash{$measure};
	#Add and subtract from the mean 
	my $lower_bound = $wt_means_hash{$measure}-$interval;  
	my $upper_bound = $wt_means_hash{$measure}+$interval;
	#Save in hash
	#We round and floor the intervals for distinct units eg everything but Telomeres
	if ($measure eq 'ANEUP') {
		$lower_bound=ceil($lower_bound);
		$upper_bound=floor($upper_bound);
	}
	$lower_bounds{$measure}=$lower_bound; 
	$upper_bounds{$measure}=$upper_bound;
}
####


###############################
##  STEP3 - Find sig samples ##
###############################

#Idea: 
# - get a list of all genes
# - make a hash with gene as key and all SD numbers in a comma separated string or array 
# - determine for each SD number whether it is sig for any of the measures
# - store SD numbers in Hash of arrays - where each column (high and low) has an array containing all SD numbers
# - later we can loop through all genes, extract all SD numbers for that gene and then loop through the hash of arrays to determine whether the gene should be used as output  


#Define a hash for all genes
my %all_genes; 
#Build a significance hash
my %sig_SD_numbers;
#Build a hash containing all significant values
my %values_sig;
#and one that contains all values
my %values;
#Build a complicated structure:
# a hash: measure (key) -> Gene (key) -> SD number (array)
my %results; 

my %counts_strains_failing;

foreach (@columns) {$counts_strains_failing{$_.'_-'}=0; $counts_strains_failing{$_.'_+'}=0;};

#Open rep.txt data file and loop through
my $ifh;
my $tc;
if( $input =~ /\.gz/ ){open($ifh, qq[gunzip -c $input|]);}else{open($ifh, $input ) or die $!;}
while( my $l = <$ifh> ){
	next if $l =~ /^\#/;
	#read the cover into a hash
	chomp( $l );
    # - make a hash with gene as key and all SD numbers in a comma separated string or array 
	my @s = split( /\t/, $l ); #Split the line into its columns and loop through columns
	my $SD_no = $s[0];
	my ($del, $gene) = split('_',$s[$gene_column-1]);
	#add SD number to hash of genes
	$all_genes{$gene}.="$SD_no,"; 
	# - determine for each SD number whether it is sig for any of the measures
	$tc++;
	for my $i (0 .. 14,31,32) {
   		#skip any columns that are not numbers
   		#next if (!looks_like_number($s[$i]));
  		next if $s[$i] =~ /ERS[0-9]{5,}/ or $s[$i] =~ /Del/ or $s[$i] =~ /SD/; #skip ERS column
		#To allow mat field analysis instead of skipping we just remove non numerical characters
		if ($s[$i] =~ /alpha|a/){
			$s[$i] =~ s/[a-z]//g;
			$s[$i] =~ s/\///g;
			$s[$i] =~ s/\(//g;
			$s[$i] =~ s/\)//g;
			$s[$i] =~ s/ //g;
	#		print "$s[$i]\n";				
		}
		#print "$columns[$i-1]\t$s[$i]\n"; #This shows that this is how you can relate column to its number
   		#Check whether the columns are outside the boundaries
		# --> This bit could be coded a bit more elegant. Have a think. Too much repetition.
		#First check the Lower bounds
		my $lower_key = $columns[$i-1].'_-';
		my $upper_key = $columns[$i-1].'_+';
		$values{$upper_key}{$SD_no}=$s[$i];
		$values{$lower_key}{$SD_no}=$s[$i];
		if (  $s[$i] <  $lower_bounds{$columns[$i-1]} ) {
			$counts_strains_failing{$lower_key}++;
			#print "YES:$lower_key  $s[$i] <  $lower_bounds{$columns[$i-1]}\n";
			#Store the significant samples in a hash of arrays
			push(@{$sig_SD_numbers{$lower_key}}, $SD_no); #creates a hash of array like this: push(@{$hash{$key}}, $insert_val); 
			#Store the values in a hash of hashes
			$values_sig{$lower_key}{$SD_no}=$s[$i];
			#Store the results in the complicated hash
			push(@{$results{$lower_key}{$gene}}, $SD_no);
		} 
		#Then check the higher bounds
		elsif (  $s[$i] >  $upper_bounds{$columns[$i-1]} ) {
			 $counts_strains_failing{$upper_key}++;
			#print "YES: $upper_key $s[$i] <  $lower_bounds{$columns[$i-1]}\n";
			#Store the significant samples in a hash of arrays
			push(@{$sig_SD_numbers{$upper_key}}, $SD_no); #creates a hash of array like this: push(@{$hash{$key}}, $insert_val); 
			#Store the values in a hash of hashes
			$values_sig{$upper_key}{$SD_no}=$s[$i];
			#Store the results in the complicated hash
			push(@{$results{$upper_key}{$gene}}, $SD_no);
		}			
    }
}
close( $ifh );


### Step 3b
#Print all genes to a file called ALL.txt
my @genes = sort(keys %all_genes);
open(my $fh1, '>', 'ALL.txt');
for my $gene (@genes){ print $fh1 "$gene\n"; }
close $fh1;

#If you want to check the multi-dimensional data structures:
#print Dumper \%sig_SD_numbers;
#print Dumper \%values;
#print Dumper \%values_sig;
#print Dumper \%results;

#########################################
##  STEP4 - Check all samples per gene ##
#########################################
# Now we want to check whether at least two samples for each gene are significant 
# If there was only one sample per gene

my %counts_gene_failing;

foreach (@columns) {$counts_gene_failing{$_.'_-'}=0; $counts_gene_failing{$_.'_+'}=0;};

my %hits;

#Loop through repeat by repeat (smaller and larger sets)
foreach my $repeat (sort keys %results) {
	#Print the results to file
	my $file_name = $repeat.'.txt';
	open(my $fh2, '>', $file_name);
	#Loop through all the genes that have at least one significant sample for this significance category
	foreach my $gene (sort keys %{$results{$repeat}}){
		my $val_string = '';
		#Check how many significant samples there are
		my $number = scalar (@{$results{$repeat}{$gene}});	
		#if it is two or more we can print it to results immediately
		#if it is only one we need to check how many samples there are for that gene
		my @SDs = split ',', $all_genes{$gene};
		if ($number == 1){
			my $number_of_samples =	scalar @SDs;
			next if ($number_of_samples != 1);
			#Get the actual value associated with the one SD number
			$val_string = $values{$repeat}{$SDs[0]};   		
			 $counts_gene_failing{$repeat}++;
		}
		elsif ($number > 1) {
			$counts_gene_failing{$repeat}++;
			my @estimates='';
			#Get the values
			foreach my $sd (@SDs){
				if (defined $values{$repeat}{$sd}){ 
					my $val = $values{$repeat}{$sd};
					push @estimates, $val;
				}
			} 
			$val_string = join "\t", @estimates;	
		}
		elsif ($number == 0) {
			#There really shouldn't be a gene with 0 samples
			print "Something funny happened here: $repeat\t$gene\t$number! Investigate!\n";
		}
		#Print the gene, followed by the values (all of them, bot just those that are significant) 
		print $fh2 "$gene_db{$gene}\t$gene\t$val_string\n"; 	
		push @{$hits{$gene}}, $repeat;
	}
	#Close the results file
	close $fh2;
}

#print Dumper \%counts_gene_failing;
print "Tot.strains ".scalar @genes."\n\n";
# Print CIs
print "Repeat\tWT_mean\tWT_StDev\tCI(+/-StDev*$multiplier)\tno. genes failing\tno. strains failing\n";
foreach my $rep (@columns){
		next if $rep eq 'x';
		my $tot_fail=$counts_gene_failing{$rep.'_-'}+$counts_gene_failing{$rep.'_+'};
		my $tot_str_fail=$counts_strains_failing{$rep.'_-'}+$counts_strains_failing{$rep.'_+'};
                printf ("%s\t%.2f\t%.2f\t( %.1f - %.1f )\t%d ( %.2f%% )\t%d ( %.2f%% )\n", $rep,$wt_means_hash{$rep},$wt_stdev_hash{$rep},$lower_bounds{$rep},$upper_bounds{$rep}, $tot_fail, ($tot_fail/scalar @genes)*100, $tot_str_fail, ($tot_str_fail /$tc) *100 );

}
print "\n\n";

#Now we want to print a table with all the strains

@columns=qw[rDNA CUP1 M 2mic T1 T2 T3 TEL ANEUP CR];
open(my $out, '>', 'overlaps.tsv');
print $out "Gene\t".join("\t", @columns)."\tSum"."\n";
my @output;
my @keys;
my %stats;
foreach my $gene(@genes){
	 print $out "$gene\t";
	if (exists $hits{$gene} and defined $hits{$gene}){
	 #    print $out "$gene\t"; 
		for my $k (@columns){
			my $out_value=0;
			if ("$k"."_-" ~~ @{$hits{$gene}}){push @output, "-1"; push @keys, $k}
			elsif ("$k"."_+" ~~ @{$hits{$gene}}){push @output,"1";push @keys, $k}
			else {push @output, "0"}
		}
		 print $out join ("\t", @output);
		my $sum=0;
		$sum=$sum+abs foreach @output;
		print $out "\t$sum";
		s/^T[0-9]$/T/g for @keys;	
		s/^rDNA$/Tandem/g for @keys;
		s/^CUP1$/Tandem/g for @keys;
		@keys=uniq(@keys);
		$stats{join ':', @keys}++;
		@output=();
		@keys=();		
	}
	else{
		for my $k (@columns){push @output,"0"}
		print $out join ("\t", @output);
		my $sum=0;
                $sum=$sum+abs foreach @output;
                print $out "\t$sum";
                @output=();

	}
	print $out "\n";
}
close ($out);
#print Dumper \%hits;

open($out, '>', 'stats_hits.tsv');
for my $k (keys %stats){
	print $out "$k\t$stats{$k}\n" if ($k =~ tr/://) >= 2;
}
close ($out);
##################
##	SUBROUTINES ##
##################

sub wt_column_average {
	my $column_no = shift;
	my $command = "cat $input | grep 'WT-' | tr -d \"[a-z]/()\" | awk -F \"\t\" '".'{ sum += $'."$column_no".'; n++ } END { if (n > 0) print'." sum / n; }'";
	my $col_av = `$command`;
	chomp $col_av;
	return $col_av;
}

sub wt_column_stdev {
	my $column_no = shift;
	my $command = "cat $input | grep 'WT-' |tr -d \"[a-z]/()\"| awk -F \"\t\" '". '{sum+=$'."$column_no".'; array[NR]=$'."$column_no".'} END '."{for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'"; 
	my $stdev = `$command`; chomp $stdev;
	return $stdev;
}

sub wt_column_stderr {
	my $stdev = shift;
	my $command = "cat $input | grep 'WT-' | wc -l";
	my $no_samp = `$command`; chomp $no_samp;
	my $root = sqrt($no_samp);
	my $stderr = $stdev / $root;
	return $stderr;
}

##################
##				##
## 	USAGE POD	##
##				##
##################

__END__


=head1 SYNOPSIS

repeat_stats_sig.pl [options] -i <filename>
 
 Options:
   -help	brief help message
   -i		input (repeat file)
   
   

=head1 OPTIONS

=over 4

=item B<-help>

Prints a brief help message and exits.

=item B<-i>

Accepts a path to file with results for DNA repeat estimates.

=back


=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something useful with the contents thereof.

=cut

