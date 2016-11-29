#!/usr/bin/env perl
# 
# Author:       tk2
# Maintainer:   tk2
# Created: 

#Test: perl /nfs/users/nfs_m/mh23/Scripts/vcf_to_gene_list.pl -i Sae2.merged.vcf > Sae2.gene_list.txt
     

use Carp;
use strict;
use warnings;
use Getopt::Long;

my @csqs = ('stop_gained', 'missense_variant','STOP_LOST','frameshift_variant','initiator_codon_variant','synonymous_variant');

my ($input);
my $aminoacid;

#the conversion file
#S. cerevisiae:
my $path = __FILE__;
my @path_components = split( /\//, $path );
my $rem = pop @path_components;
$path = (join "/",  @path_components);
my $conversion_file="$path".'/yeast_genelist_nameconv.tsv';
# Declare variables:
my %name_conversion; # a hash with the systematic name as the key and the common name as the entry
my %description; # Build a hash with the systematic name as key and the description as entry

# Open the Conversion file and make hashes for the systematic gene names:
open(F, $conversion_file ) or die ("Unable to open file $conversion_file: $!\n" );
        #go through file line by line
        while ( my $line = <F>) {
                next if $. == 1; #first line skipped
                #split the line into its columns
                my($f1, $f2, $f3, $f4, $f5, $f6, $f7, $f8) = split '\t', $line;
        		# $f2 is the systematic gene name; $f3 is the common name and $f4 is the description
                # Build a hash with the systematic name as the key and the common name as the entry
                $name_conversion{$f2}=$f3;
                # Build a hash with the systematic name as key and the description as entry
                $description{$f2}=$f4;
                #Test print:
        }# close while
close F or die "Cannot close $conversion_file: $!\n";


GetOptions
(
    'i|input=s'         => \$input,
);

( $input && -f $input ) or die qq[Usage: $0 -i <input vcf>\n];

my $ifh;
if( $input =~ /\.gz/ ){open($ifh, qq[gunzip -c $input|]);}else{open($ifh, $input ) or die $!;}

my @samples;

print qq[TYPE\tCHROM\tPOS\tGENE\tGENE\tDESCRIPTION\tCSQ\tHETS\tHET Samples\n];
while( my $l = <$ifh> )
{
    my $m = $l;
    chomp( $l );
    $aminoacid = 'NA';
    next if( $l =~ /^#/ && $l !~ /^#CHROM/);
    my @s = split( /\t/, $l );
    
    if( $l =~ /#CHROM/ ){for(my $i=0;$i<@s;$i++){$samples[$i]=$s[$i];}next;}
    
    foreach my $csq(@csqs)
    {
        if( $l =~ /$csq/i ) 
        {
            my @s1 = split( /;/, $s[ 7 ] );
            
            foreach my $tag(@s1)
            {
                next unless $tag =~ /^(CSQ=)(.+)/;
                
                my @s2 = split( /\+/, $2 );
                foreach my $csq(@s2)
                {
                    foreach my $csq_del(@csqs)
                    {
                        if( $csq =~ /$csq_del/i )
                        {
                            #get number of samples affected
                            my $homs = 0;my $hets = 0;
                            my $samplesHet='';
                            my $samplesHom='';
                            my $c = 0;for(my $i=9;$i<@s;$i++){if($s[$i]=~/1\/1:/){$homs++;$samplesHom.=qq[$samples[$i];];}elsif($s[$i]=~/1\/0:/||$s[$i]=~/0\/1:/){$hets++;$samplesHet.=qq[$samples[$i];];}}
                            #the consequence line will be split on the :
                            #YLR442C:YLR442C:missense_variant:2732:911:C>Y
                            #that makes [2] the consequence
                            my @s3 = split( /\:/, $csq );
                            if ( $csq =~ /missense/i )
                        	{
                        		#wanted: C911Y;
                        		#YLR442C[0]:YLR442C[1]:missense_variant[2]:2732[3]:911[4]:C>Y[5]
                                my @s5 = split( /\>/, $s3[5]);
                                $aminoacid = '£x'.$s5[0].$s3[4].$s5[1].'£';
                            }
                            if ( $csq =~ /initiator_codon_variant/i )
                        	{
                        		#wanted: M1I;
                        		#YIL129C[0]:YIL129C[1]:initiator_codon_variant[2]:3[3]:1[4]:M>I[5]
                        		my @s5 = split( /\>/, $s3[5]);
                        		$aminoacid = '£'.$s5[0].$s3[4].$s5[1].'£';	
                            }
                            if ( $csq =~ /frameshift_variant/i )
                        	{
                        		#wanted: range
                        		#YLL066W-B[0]:YLL066W-B[1]:frameshift_variant,feature_truncation[2]:59[3]:20[4]
                                $aminoacid = '£FS@'.$s3[4].'£';
                            }
                            if ( $csq =~ /stop_gained/i )
                            {
                                if ($m =~ /INDEL/) {
                                    #delta(number-1)
                                    #YLR442C[0]:YLR442C[1]:stop_gained[2]:2733-2775[3]:911-20[4]
                                    my @new_num = split ( /\-/, $s3[4] );
                                    my $number=$new_num[0]-1;
                                    $aminoacid = '£'.'Δ(FS)'.$number.'£';}
                                
                                else{
                                    #delta(number-1)
                                    #YLR442C[0]:YLR442C[1]:stop_gained[2]:2733[3]:911[4]
                                    my $number=$s3[4]-1;
                                    $aminoacid = '£'.'Δ'.$number.'£';}
                            }
                            if ( $csq =~ /synonymous/i )
                        	{
                        		#wanted £P91:C>T:2£
                        		#YDR420W[0]:YDR420W[1]:synonymous_variant[2]:1533[3]:511[4]:P>P[5]
                        		my @s5 = split( /\>/, $s3[5]);
                        		my $remainder = $s3[3] % 3;
                        		my $codon;
                        		if ($remainder == 0){$codon=3;}
                        		if ($remainder == 1){$codon=1;}
                        		if ($remainder == 2){$codon=2;}
                        		$aminoacid = '£'.$s5[0].$s3[4].':'.$s[3].'>'.$s[4].':'.$codon.'£';
                            }
                            if ( $csq =~ /missense/i )
                            {
                                #wanted: C911Y;
                                #YLR442C[0]:YLR442C[1]:missense_variant[2]:2732[3]:911[4]:C>Y[5]
                                my @s5 = split( /\>/, $s3[5]);
                                $aminoacid = '£x'.$s5[0].$s3[4].$s5[1].'£';
                            }

                            if( $homs > 0 || $hets > 0)
                            {
                            if( $l=~/INDEL/){print qq[INDEL\t];}else{print qq[SNP\t];}
                            if(defined($name_conversion{$s3[1]})){
        					print qq[$s[0]\t$s[1]\t$aminoacid\t$s3[1]\t$name_conversion{$s3[1]}\t$description{$s3[1]}\t$s3[2]\t$homs\t$samplesHom\t%\t$hets\t$samplesHet\n];
                            }
                            else{
                            print qq[$s[0]\t$s[1]\t$samplesHet\t$s3[1]\t$s3[2]\t$s3[2]\t$s3[2]\t$hets\n];
                            }
                            }
                        }
                    }
                }
            }
        }
    }
}
close( $ifh );