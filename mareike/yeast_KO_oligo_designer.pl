

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


foreach my $g (sort keys %results){
	print "$g\t$results{$g}\n";
}






