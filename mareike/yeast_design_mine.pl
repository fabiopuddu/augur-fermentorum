


#!/usr/bin/perl
# Author:       mh23
# Maintainer:   mh23

# Name: yeast_KO_oligo_designer.pl

use strict;
use warnings;

use lib '/Users/mh23/Documents/Data/2017_01_10-Package_writing';
use Design_KO_Oligos_Mine('Design_KO_Oligos_Mine');

my $input = shift;

my %results = Design_KO_Oligos_Mine($input);


print "$results{'Standard_Name'}-F1\t$results{'F1'}\n";
print "$results{'Standard_Name'}-R1\t$results{'R1'}\n";
print "$results{'Standard_Name'}.3\t$results{'.3'}\n";  
print "$results{'Standard_Name'}.4\t$results{'.4'}\n";  

#print "PCR length\t$results{'PCR_length'}\n";  
#print "Systematic\t$results{'Systematic_Name'}\n"; 




