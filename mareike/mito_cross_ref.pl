#!/usr/bin/env perl

# Author:       mh23
# Maintainer:   mh23
# Created: 		Sep 2016
# Name: mito_cross_ref.pl


use warnings;
use strict;

#my $SGD = shift;
my $genelist = shift;

#Open genelist 
# Open the Conversion file and make hashes for the systematic gene names:
open(F, $genelist ) or die ("Unable to open file $genelist: $!\n" );
        #go through file line by line
        while ( my $line = <F>) {
           next if ($line !~ /Del/);
           #print $line;
           chomp $line;
           my $gene = 'NA';
           my @s = split '_', $line;
           #print "$s[1]\n";
           if( $s[1] =~ /([A-Z]+[0-9]+(,\d)?(\w)?)/ig)
			{
  				$gene = $1;
			}
           
           #check its existence in the mito_cand_genes_SGD.txt file  
           my $command = "cat mito_cand_genes_SGD.txt | grep $gene";
           my $result = `$command`;       
           if ($result eq ''){
           	print "$gene\tNot_on_SGD\n";
           }
           else{
           	print "$result";
           } 
           #exit if $. == 10;     
                
         
        }# close while
close F or die "Cannot close $genelist: $!\n";

