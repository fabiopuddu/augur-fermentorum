#!/usr/bin/env perl
# 
# Author:       mh23	
# Maintainer:   mh23
# Created: 		23.02.2015

#Test: perl /nfs/users/nfs_m/mh23/Scripts/vcf_stats_table.pl -i test.mutations.vcf 
#description: A script to turn the output of vcf-stats into sth more useful
     
 
     

use Carp;
use strict;
use warnings;
use Getopt::Long;

#Declare variables
my $input = shift;
my @s; my $sample;
my $indel_count;
my $snp_count; 
$sample='NA';$indel_count='0';$snp_count='0';
my $check=0; 
my $transitions; my $transversions;
my $A_C=0; 
my $A_G=0; 
my $A_T=0; 
my $C_A=0; 
my $C_G=0; 
my $C_T=0; 
my $G_A=0; 
my $G_C=0; 
my $G_T=0; 
my $T_A=0; 
my $T_C=0; 
my $T_G=0; 

#Get command line options
GetOptions
(
    'i|input=s'         => \$input,
);

( $input && -f $input ) or die qq[Usage: $0 -i <input vcf>\n];


my $ifh;
#if( $input =~ /\.gz/ ){open($ifh, qq[gunzip -c $input|]);}else{open($ifh, $input ) or die $!;}
open(F, $input ) or die ("Unable to open file $input: $!\n" );

 #go through file line by line
print("Sample\tINDEL_Count\tSNP_Count\tTransitions\tTransversions\tA>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\n");
while ( my $l = <F>) {
	my $k = $l;
	chomp( $l );
	#The following code will work like this:
	#The line will be matched for the information I am looking for e.g. indel_count
	#The line will then be split on an empty space
	#Because the output of vcf-stats has a lot of tabs in it we will iterate through the split information
	#These the script will try to match to numbers (because we are interested in those)
	#Any number is fine (that's why the match says [0-9]), the + indicates we want one or more digit
	#the () around say that we want to capture the match in the special variable $1
	#So once this match has been found we capture the information $1 in a more permanent variable
	#get the sample name
   	
   	#The other bit is that the logic of the code is a bit backward
   	#We need to find a way to count from sample to sample
   	#Thus everytime a sample is found we set the variable check to 1 
   	#and then when the information is printed we set it to 0
   	#Thus if a sample has already been found and this is a new sample we know any
   	#available information has been captured, we print it and then we reset variables
   	
   	if ($l =~ /E/ && $check==1 || $l =~ /all/) {	
   		$transitions=$A_G+$G_A+$C_T+$T_C;
   		$transversions=$A_T+$A_C+$C_A+$C_G+$G_T+$G_C+$T_A+$T_G;
   		print("$sample\t$indel_count\t$snp_count\t$transitions\t$transversions\t$A_C\t$A_G\t$A_T\t$C_A\t$C_G\t$C_T\t$G_A\t$G_C\t$G_T\t$T_A\t$T_C\t$T_G\n");
   		$sample='NA';$indel_count='0';$snp_count='0'; $A_C=0; $A_G=0; $A_T=0; $C_A=0; 
		$C_G=0; $C_T=0; $G_A=0; $G_C=0; $G_T=0; $T_A=0; $T_C=0; $T_G=0; 
   		$check=0;
   	}
   	
   	#whereas this block captures information:
	if ($l =~ /E/ && $check==0) {$check=1; @s = split( /'/, $l ); foreach my $t (@s){if ($t =~ /E/) {$sample = $t;}}}
	if ($l =~ /indel_count/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$indel_count = $1;} };}
    if ($l =~ /snp_count/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$snp_count = $1;} };}
	if ($l =~ /A>C/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$A_C = $1;} };}
	if ($l =~ /A>G/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$A_G = $1;} };}
	if ($l =~ /A>T/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$A_T = $1;} };}
	if ($l =~ /C>A/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$C_A = $1;} };}
	if ($l =~ /C>G/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$C_G = $1;} };}
	if ($l =~ /C>T/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$C_T = $1;} };}
	if ($l =~ /G>A/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$G_A = $1;} };}
	if ($l =~ /G>C/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$G_C = $1;} };}
	if ($l =~ /G>T/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$G_T = $1;} };}
	if ($l =~ /T>A/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$T_A = $1;} };}
	if ($l =~ /T>C/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$T_C = $1;} };}
	if ($l =~ /T>G/) { @s = split( / /, $l ); foreach my $t (@s){if ($t =~ /([0-9]+)/){$T_G = $1;} };}
      
}# close while
close F or die "Cannot close $input: $!\n";
