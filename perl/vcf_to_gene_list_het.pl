#!/usr/bin/env perl
#
# Author:       tk2
# Maintainer:   Fabio Puddu & Mareike Herzog
# Created:
use Carp;
use strict;
use warnings;
use Getopt::Long;

my @csqs = ('stop_gained', 'missense_variant','stop_lost','frameshift_variant','initiator_codon_variant','synonymous_variant', 'inframe_insertion', 'inframe_deletion');
my ($input);
my $aminoacid;

my $path = __FILE__;
my @path_components = split( /\//, $path );
my $rem = pop @path_components;
$path = (join "/",  @path_components);
my $conversion_file="$path".'/yeast_genelist_nameconv.tsv';     #the conversion file
# Declare variables:
my %name_conversion;                                            # a hash with the systematic name as the key and the common name as the entry
my %description;                                                # Build a hash with the systematic name as key and the description as entry

# Open the Conversion file and make hashes for the systematic gene names:
open(F, $conversion_file ) or die ("Unable to open file $conversion_file: $!\n" );
while ( my $line = <F>) {                               #go through file line by line
        next if $. == 1;                                #first line skipped
        my($f1, $f2, $f3, $f4, $f5, $f6, $f7, $f8) = split '\t', $line; #split the line into its columns
	# $f2 is the systematic gene name; $f3 is the common name and $f4 is the description
        $name_conversion{$f2}=$f3;      # Build a hash with the systematic name as the key and the common name as the entry
        $description{$f2}=$f4;          # Build a hash with the systematic name as key and the description as entry
}
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
while( my $l = <$ifh> ){                                #for every line in the file
        my $m = $l;;
        chomp( $l );
        $aminoacid = 'NA';
        next if( $l =~ /^#/ && $l !~ /^#CHROM/);            # skip the header
        my @s = split( /\t/, $l );                          # split the string
        if( $l =~ /#CHROM/ ){
                for(my $i=0;$i<@s;$i++){
                        $samples[$i]=$s[$i];                       #store the line header in a list
                }
        next;                                              #and go to the next line
        }
        next unless $l =~ /PASS/;                          #skip mutations not passing filters
        foreach my $csq_del(@csqs){                        #for every consequence that we are checking
                if ($l =~ /$csq_del/i ){                   #does the line contain that consequence ?
                        my @s1 = split( /;/, $s[ 7 ] );    #if so we need to isolate the INFO field $s1="VDB=5.5115e-07;SGB=-0.688148;RPB=2.32974e-05;MQB=0.0148624;MQSB=0.734517;BQB=0.677425;MQ0F=0.0105263;ICB=1;HOB=0.5;MQ=54;CSQ=Y....."
                        foreach my $tag(@s1){              #go through all the tags
                                next unless $tag =~ /^(CSQ=)(.+)/;      #find the one starting with CSQ=
                                my @s2 = split( /\+/, $2 );#split on plus
                                foreach my $csq(@s2){      #run through all the consequences
                                        my @s3 = split( /\:/, $csq );           #the consequence line will be split on the : YLR442C:YLR442C:missense_variant:2732:911:C>Y
                                        if ( $csq =~ /missense/i ){                         #YLR442C[0]:YLR442C[1]:missense_variant[2]:2732[3]:911[4]:C>Y[5];  #wanted: C911Y;
                                                my @s5 = split( /\>/, $s3[5]);
                                                $aminoacid = '£x'.$s5[0].$s3[4].$s5[1].'£';
                                        }
                                        elsif ( $csq_del  =~ /initiator_codon_variant/i ){
                        		        if ($m =~ /INDEL/) {
                                                        $aminoacid = '£FS@'.$s3[4].'£';
        				        }
        				        else {
        					        my @s5 = split( /\>/, $s3[5]);                     #YIL129C[0]:YIL129C[1]:initiator_codon_variant[2]:3[3]:1[4]:M>I[5]
                        			        $aminoacid = '£'.$s5[0].$s3[4].$s5[1].'£';         #wanted: M1I;
                            		        }
        			        }
                                        elsif ( $csq =~ /frameshift_variant/i ){
                                                $aminoacid = '£FS@'.$s3[4].'£';                 #YLL066W-B[0]:YLL066W-B[1]:frameshift_variant,feature_truncation[2]:59[3]:20[4];  #wanted: range
                                        }
                                        elsif ( $csq =~ /inframe_insertion/i ){
                                                $aminoacid = '£II@'.$s3[4].":".$s3[5].'£';                 #YLL066W-B[0]:YLL066W-B[1]:frameshift_variant,feature_truncation[2]:59[3]:20[4] #wanted: range
                                        }
                                        elsif ( $csq =~ /inframe_deletion/i ){
                                                $aminoacid = '£ID@'.$s3[4].":".$s3[5].'£';             #YLL066W-B[0]:YLL066W-B[1]:frameshift_variant,feature_truncation[2]:59[3]:20[4] #wanted: range
                                        }
                                        elsif ( $csq =~ /stop_gained/i ){
                                                if ($m =~ /INDEL/) {
                                                        my @new_num = split ( /\-/, $s3[4] );
                                                        my $number=$new_num[0]-1;               #YLR442C[0]:YLR442C[1]:stop_gained[2]:2733-2775[3]:911-20[4]        #delta(number-1)
                                                        $aminoacid = '£'.'Δ(FS)'.$number.'£';}
                                                else{
                                                        my $number=$s3[4]-1;                    #YLR442C[0]:YLR442C[1]:stop_gained[2]:2733[3]:911[4]                #delta(number-1)
                                                        $aminoacid = '£'.'Δ'.$number.'£';}
                                                }
                                        elsif ( $csq =~ /synonymous/i ){
                        		        my @s5 = split( /\>/, $s3[5]);
                                        	my $remainder = $s3[3] % 3;
                        		        my $codon;
                                		if ($remainder == 0){$codon=3;}
                                		if ($remainder == 1){$codon=1;}
                                		if ($remainder == 2){$codon=2;}
                                		$aminoacid = '£'.$s5[0].$s3[4].':'.$s[3].'>'.$s[4].':'.$codon.'£';    #YDR420W[0]:YDR420W[1]:synonymous_variant[2]:1533[3]:511[4]:P>P[5]    #wanted £P91:C>T:2£
                                        }
                                        elsif ( $csq =~ /stop_lost/i ){
                                                $aminoacid = '£st_lost£';
                                        }
                                        #NOW WE NEED TO COUNT HOW MANY SAMPLES ARE AFFECTED BY EACH MUTATION AND STORE THE SAMPLE NAMES
                                        if( $csq =~ /$csq_del/i ){      #and isolate the first one matching the one that we are checking
                                                my $homs = 0;my $hets = 0;  #get number of samples affected
                                                my $samplesHet='';
                                                my $samplesHom='';
                                                my $c = 0;
                                                for(my $i=9;$i<@s;$i++){    #run through all the sample fields(starting from index 9)
                                                    if($s[$i]=~/\//){
                                                            my @gty = split( /\:/, $s[$i]);     #extract the genotype
                                                            my @g = split( /\//, $gty[0]);      #extract the genotype components
                                                            my $ploidy = scalar @g;
                                                            if ($ploidy eq '2'){
                                                                    if ($g[0] ne '.'){
                                                                            if($g[0]==$g[1]){            #1/1, 2/2, etc
                                                                                    $homs++;$samplesHom.=qq[$samples[$i];];
                                                                            }
                                                                            else{                       #0/1, 1/2, etc
                                                                                    $hets++;$samplesHet.=qq[$samples[$i];];
                                                                            }
                                                                    }

                                                            }

                                                            if ($ploidy eq '4'){
                                                                    if ($g[0] ne '.'){
                                                                            if(@g == grep { $_ == $g[0] } @g ){         #if all the values in the g array are the same e.g. 1/1/1/1 or 2/2/2/2
                                                                                    $homs++;$samplesHom.=qq[$samples[$i];];
                                                                            }
                                                                            else{
                                                                                    $hets++;$samplesHet.=qq[$samples[$i];];
                                                                            }
                                                                    }
                                                            }
                                                    }
                                                }
                                                if( $homs > 0 || $hets > 0){
                                                        if( $l=~/INDEL/){
                                                            print qq[INDEL\t];
                                                        }
                                                        else{
                                                            print qq[SNP\t];
                                                        }
                                                        if(defined($name_conversion{$s3[1]})){
                                                            print qq[$s[0]\t$s[1]\t$aminoacid\t$s3[1]\t$name_conversion{$s3[1]}\t$description{$s3[1]}\t$csq_del\t$homs\t$samplesHom\t%\t$hets\t$samplesHet\n];
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
close( $ifh );
