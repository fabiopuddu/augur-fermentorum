
# Author:       mh23
# Maintainer:   mh23
# Created: 		Dec 2016
# Name: rDNA_cross_ref.pl


use warnings;
use strict;

#my $SGD = shift;
#input should be the lower/higher sig genes
my $genelist = shift;

#Open genelist 
# Open the Conversion file and make hashes for the systematic gene names:
open(F, $genelist ) or die ("Unable to open file $genelist: $!\n" );
        #go through file line by line
        while ( my $line = <F>) {
           #print "$.\n";
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
           #my $command = "cat /mnt/scratch/jackson/fp305/results/rDNA_cand_genes_kobayashi.txt | grep $gene -w";
           my $command = "cat /mnt/scratch/jackson/fp305/results/rDNA_cand_genes_kobayashi_all.txt   | grep $gene -w";
           #my $command = "cat /mnt/scratch/jackson/fp305/results/rep.txt | grep $gene -w";
           my $result = `$command`;       
           if ($result eq ''){
           	print "$gene\tNot_on_Kobayashi\n";
           }
           else{
           	my @p = split "\n", $result;
           	#print "$result";
           	print "$p[0]\n";
           } 
           #exit if $. == 10;     
                
         
        }# close while
close F or die "Cannot close $genelist: $!\n";

