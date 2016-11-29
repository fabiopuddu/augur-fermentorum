#!/usr/bin/env perl
#
# Author:       mh23
# Maintainer:   mh23
# Created: 		23.02.2015

#Test: perl /nfs/users/nfs_m/mh23/Scripts/vcf_stats_table_all.pl Mutation_Types_140225.txt
#description: A script to turn the output of vcf-stats into sth more useful




use Carp;
use strict;
use warnings;
use Getopt::Long;
use Cwd;


#Declare variables
my $pwd = cwd(); #current working directory; full path
#if you want only the folder name
my @p = split( /\//, $pwd);
$pwd = $p[-2]; #the directory will be the last element of the array
my @s; my $sample;
my $indel_count;
my $snp_count;
$indel_count='0';$snp_count='0';
my $check=0;
my $transitions; my $transversions;
my @number; my %result;
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


#go through file line by line
#open my $fh, "<", $input or die $!;
my $input = shift;
#my $command =  "cat $input"." | grep ".'"'.'1/1" -v'." | vcf-stats";
my $command = "cat $input"." | grep ".'"'.'1/1\|0/1\|#"'." | vcf-stats"; #edited by FABIO to resolve bug that does eliminates rows where one of the mutation is masked
my @stats_out =  readpipe("$command");
print("INDEL\tSNP\tTs\tTv\tSample Name\tC>T\tA>G\tA>T\tC>G\tG>T\tA>C\t<-3\t-3\t-2\t-1\t1\t2\t3\t>3\tNumber of samples\n");


foreach my $l (@stats_out){
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
   	
   	
   	#whereas this block captures information:
	if ($l =~ /all/) { $indel_count='0';$snp_count='0'; $A_C=0; $A_G=0; $A_T=0; $C_A=0; $C_G=0; $C_T=0; $G_A=0; $G_C=0; $G_T=0; $T_A=0; $T_C=0; $T_G=0; }
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
    
    #get the indel count information
    #the check tracks where the indel information starts and stops
    if ($l =~ /indel/ && $l !~ /indel_count/){$check = 1; }
    #if the check is 1 and the line contains a number it is indel information
    if ($check == 1 && $l =~ /([0-9]+)/) {#get the info
        @s = split( / /, $l ); #split on the empty characters and loop through
        #first try to match a negative indel with the - and push it into arrays
        #if that has happened we don't want it to match again (the one without the hyphen, so next)
        #if no - has been mached we want both numbers to be pushed into the array
        foreach my $k (@s){if ($k =~ /(-[0-9]+)/){push(@number, $1);next;} if ($k =~ /([0-9]+)/){push(@number, $1);}}
        #use the first number as the key and the second as the value
        $result{$number[0]} =$number[1];
        #the array has to be emptied every time so that the first two numbers are always the newest numbers
        my $removed = shift @number; $removed = shift @number;
    }
    #when the closing bracket appears that means the indel bit is over, so check is 0
    if ($l =~ /}/){$check = 0;}
   	
}# close

#transitions
$C_T =	$C_T + $G_A;
$A_G =	$A_G + $T_C;
#transversions
$A_T =	$A_T + $T_A;
$C_G =	$C_G + $G_C;
$G_T =  $G_T + $C_A;
$A_C =	$A_C + $T_G;

$transitions=$A_G+$C_T;
$transversions=$A_T+$A_C+$C_G+$G_T;

#INDELs
my $one=0;
my $two=0;
my $three=0;
my $more=0;
my $minone=0;
my $mintwo=0;
my $minthree=0;
my $less=0;

for (sort { $a cmp $b } keys %result){
	if ($_ == 1){$one = $result{$_};}
	if ($_ == 2){$two = $result{$_};}
	if ($_ == 3){$three = $result{$_};}
	if ($_ == -1){$minone = $result{$_};}
	if ($_ == -2){$mintwo = $result{$_};}
	if ($_ == -3 ){$minthree = $result{$_};}
	if ($_ > 3) {$more = $more + $result{$_};}
	if ($_ < -3) {$less = $less + $result{$_};}
}

print("$indel_count\t$snp_count\t$transitions\t$transversions\t$pwd\t$C_T\t$A_G\t$A_T\t$C_G\t$G_T\t$A_C\t");
print("$less\t$minthree\t$mintwo\t$minone\t$one\t$two\t$three\t$more\t");
