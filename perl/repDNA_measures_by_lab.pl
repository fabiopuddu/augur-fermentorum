#!/usr/bin/env perl

# Author:       	Fabio Puddu
# Maintainer:   	Fabio Puddu
# Created: 	Apr 2017
# Description:	

use strict;
use warnings;
#This scripts runs through the deletion check log file and the name conversion file to determine
#How did undeleted and failed samples distributes across the different studies
#Two known bugs: WT samples and deletion of MF(ALPHA)1 and MF(ALPHA)2 are reported as not deleted but are good 
#at least in the samples before the propagation
open (my $rep, '<', 'rep.txt');
open (my $lab_conversion, '<', 'strain_by_lab.tab');

chomp(my @RD = <$rep>);
chomp(my @LC = <$lab_conversion>);

my @OUT;

my @lab_list=('lab_1', 'lab_2', 'lab_3', 'lab_4', 'lab_5', 'lab_6', 'lab_7', 'lab_8', 'lab_9', 'lab_10','lab_11', 'lab_12', 'lab_13', 'lab_14', 'lab_15', 'lab_16', 'unknown' );
my @SDnumbers;
foreach my $lab (@lab_list){
	@SDnumbers=();
	@OUT=();
	my @matching_samples= grep { $_ =~ $lab } @LC;
	foreach my $line (@matching_samples){
		my @linea=split "\t", $line;
		push @SDnumbers, $linea[0];
	}
	foreach my $sample(@SDnumbers){
		my $matching_line= (grep { $_ =~ /\b$sample\b/ } @RD)[0];	
		if (defined $matching_line and length $matching_line > 0){
			(undef, my $rDNA, my $cup1, my $mito, my $ty1, undef, undef, undef, undef, undef, my $telo, undef, undef)= split "\t", $matching_line;
			push @OUT, $telo;		
		}
	}
	my $string=join "\t", @OUT;
	printf "$lab\t$string\n";
}

