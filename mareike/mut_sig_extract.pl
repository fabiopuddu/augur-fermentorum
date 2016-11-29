#!/usr/bin/env perl
# 
# Author:       tk2
# Maintainer:   tk2
# Created: 

#Test: perl /nfs/users/nfs_m/mh23/Scripts/mut_sig_extract.pl -i tof1_merge.vcf
# SampleID	Chr	Locus	Nucleotide_Change
#	SampleA	1	1001	A>T
     

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



my $ifh;
if( $input =~ /\.gz/ ){open($ifh, qq[gunzip -c $input|]);}else{open($ifh, $input ) or die $!;}

my @samples;
my $sampleID;
my $chr;
my $locus;
my $nuc_change;

# SampleID	Chr	Locus	Nucleotide_Change
#	SampleA	1	1001	A>T

print qq[SampleID\tChr\tLocus\tNucleotide_Change\n];
while( my $l = <$ifh> )
{
    chomp( $l );
    next if( $l =~ /^#/ && $l !~ /^#CHROM/);
    my @s = split( /\t/, $l );
    
    if( $l =~ /#CHROM/ ){for(my $i=0;$i<@s;$i++){$samples[$i]=$s[$i];}next;}
    next if( $l =~ /INDEL/);
	$chr = $s[0];
	$locus = $s[1];
	$nuc_change = $s[3].'>'.$s[4];
	for(my $i=9;$i<@s;$i++){
		if($s[$i]=~/1\/1:/ ||  $s[$i]=~/0\/1:/){
			$sampleID = $samples[$i];
			print("$sampleID\t$chr\t$locus\t$nuc_change\n");
		}
	}

}
close( $ifh );