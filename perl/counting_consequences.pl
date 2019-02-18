#!/usr/bin/env perl
# 
# Author:       mh23
# Maintainer:   mh23
# Created: 08.01.2014
# Run command: /nfs/users/nfs_m/mh23/Scripts/counting_consequences.pl -i <vcf>
#I am sooo happy     

use Carp;
use strict;
use warnings;
use Getopt::Long;

#the categories of consequences we are looking at
my @csqs = ('stop_gained','stop_lost','frameshift_variant','missense_variant','transcript_ablation','splice_donor_variant','splice_acceptor_variant','initiator_codon_variant','inframe_insertion','inframe_deletion','transcript_amplification','splice_region_variant','incomplete_terminal_codon_variant','5_prime_UTR_variant','3_prime_UTR_variant', 'synonymous_variant' );



my $input;
GetOptions
(
'i|input=s'         => \$input,
);

( $input && -f $input ) or die qq[Usage: $0 -i <input vcf>\n];


my $counter=0;
my $counter2 = 0;
my %consequences_snp;
my %consequences_indel;
$consequences_snp{up_downstream}=0;
$consequences_indel{up_downstream}=0;
foreach my $c (@csqs){
	$consequences_snp{$c}=0;
    $consequences_indel{$c}=0;
}


my $ifh;
if( $input =~ /\.gz/ ){open($ifh, qq[gunzip -c $input|]);}else{open($ifh, $input ) or die $!;}

while( my $l = <$ifh> ){
    $counter=0;
    if ($l =~ /^#/ || $l =~ /^#CHROM/ || $l =~ /\.\/\./) {}
    else { $counter2++;}
    chomp( $l );
    next if( $l =~ /^#/ || $l =~ /^#CHROM/);
    next if( $l =~ /\.\/\./);
    foreach my $cons (@csqs){
    	next if $counter == 1;
    	if ($l =~ /$cons/i){
            if ($l =~ /INDEL/i){
                $consequences_indel{$cons}=$consequences_indel{$cons}+1;
                $counter=1; }
            else{
                $consequences_snp{$cons}=$consequences_snp{$cons}+1;
                $counter=1;}
    	}
    }
    	if ($counter == 0){
    		if ($l =~ /INDEL/i){
                $consequences_indel{up_downstream}=$consequences_indel{up_downstream}+1;}
            else{
            $consequences_snp{up_downstream}=$consequences_snp{up_downstream}+1;}
            
    	}
}
close( $ifh );

foreach my $t  ( keys %consequences_snp) {
	print ("$consequences_snp{$t}\tSNP\t$t\n");
    print ("$consequences_indel{$t}\tINDEL\t$t\n");
}
print("\n$counter2\n\n");