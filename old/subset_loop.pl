#!/usr/bin/env perl
# 
# Author:       mh23	
# Maintainer:   mh23
# Created: 

#Test: perl /nfs/users/nfs_m/mh23/Scripts/subset_loop_no_conversion.pl Sae2.merged.vcf
     

use Carp;
use strict;
use warnings;
use Getopt::Long;

my ($input);

GetOptions
(
'i|input=s'         => \$input,
);

( $input && -f $input ) or die qq[Usage: $0 -i <input vcf>\n];



#open the vcf file
my $ifh;
if( $input =~ /\.gz/ ){open($ifh, qq[gunzip -c $input|]);}else{open($ifh, $input ) or die $!;}

#Get an array of all sample names
my @samples;
while( my $l = <$ifh> )
{
    chomp( $l );
    
    next if( $l =~ /^#/ && $l !~ /^#CHROM/);
    my @s = split( /\t/, $l );
    
    if( $l =~ /#CHROM/ ){for(my $i=0;$i<@s;$i++){$samples[$i]=$s[$i];}next;}
}
close( $ifh );

for(my $count=0;$count<9;$count++) {
    my $removed = shift @samples;
}

# for each sample execute the systems command
foreach my $samp (@samples) {

	my $new_name = $samp .'.vcf';
	 my $command = "vcf-subset -c $samp $input -e > $new_name";
    #my $cmd = 'bsub -q long -R "select[type==X86_64 && mem > 10000] rusage[mem=10000]" -M10000 -o sub.o -e sub.e -J mouse_sub " '.$command.'"';
	 print "$command\n";
	 system($command);
}

