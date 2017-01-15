

#!/usr/bin/perl
# Author:       mh23
# Maintainer:   mh23

# Name: yeast_KO_oligo_designer.pl

use strict;
use warnings;

use lib '/Users/mh23/Documents/Data/2017_01_10-Package_writing';
use Design_KO_Oligos('Design_KO_Oligos');

my $input = shift;

my %results = Design_KO_Oligos($input);


#print "Gene: $results{'Standard'}: $results{'Name'}, $results{'Location'}, $results{'Strand'}\n";
print "$results{'Standard'}-F1\t$results{'F1'}\n";
print "$results{'Standard'}-R1\t$results{'R1'}\n";
print "$results{'Standard'}.3\t$results{'.3'}\n";  #$results{'Location3'}bp from Stop\n";
print "$results{'Standard'}.4\t$results{'.4'}\n";  #$results{'Location4'}bp from Stop\n";
#print "Expected PCR length: $results{'PCR_length'}\n";



