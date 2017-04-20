#!/usr/bin/env perl

use strict;
use warnings;
my @FILES=`ls */ploidy_data/*_plstats.txt`;
my @CHROMOSOMES=('chr01', 'chr02', 'chr03', 'chr04', 'chr05', 'chr06', 'chr07', 'chr08', 'chr09', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16');
my $delta;
printf "Sample name\tchr01\tchr02\tchr03\tchr04\tchr05\tchr06\tchr07\tchr08\tchr09\tchr10\tchr11\tchr12\tchr13\tchr14\tchr15\tchr16\n";
for my $filename (@FILES){ 
	$delta=0;
	chomp ($filename);
	my $fn= (split "/", $filename)[-1]; 
	$fn =~ s/_plstats.txt//;
	open (my $fh, '<', $filename);
	chomp(my @CCN = <$fh>);
	close ($fh);
	printf "$fn\t";
        for my $chr (@CHROMOSOMES){
                my $matching_line = (grep { /$chr/ } @CCN)[0];
		my $ccn=(split "\t", $matching_line)[1];
	 	$delta += abs($ccn-2);
                printf "$ccn\t";
        }
        printf "$delta\n";
}

