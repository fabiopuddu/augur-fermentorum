#!/usr/bin/env perl
use strict;
use warnings;
# Depending on which thing you are analyzing you will uncomment one of the below
# Each has two sample names (input and IP) and corresponding normalisation factor
###########################################################
#notag
# my $ip='SRR1101856';	my $ip_norm_factor=60;
# my $input='???????'; my $input_norm_factor=153;
###########################################################
#H4K16
#my $ip='SRR1105630';	my $ip_norm_factor=222;
#my $input='SRR1105631'; my $input_norm_factor=98;
###########################################################
#H3
# my $ip='SRR1106079';	my $ip_norm_factor=157;
# my $input='SRR1106080'; my $input_norm_factor=152;
###########################################################
#GFP 
# my $ip='SRR1009212';	my $ip_norm_factor=49;
# my $input='SRR1009213'; my $input_norm_factor=139;
###########################################################
#H3-mn
# my $ip='SRR1105777';	my $ip_norm_factor=188;
# my $input='SRR1105781'; my $input_norm_factor=117;
###########################################################
#sir2
# my $ip='SRR1068439';	my $ip_norm_factor=39;
# my $input='SRR1068442'; my $input_norm_factor=75;
############################################################
#sir3
# my $ip='SRR1068450';	my $ip_norm_factor=31;
# my $input='SRR1068454'; my $input_norm_factor=42;
############################################################
#sir3 replicate 2
my $ip='SRR1104091';	my $ip_norm_factor=33;
my $input='SRR1068454'; my $input_norm_factor=42;
############################################################
#sir3 in Sir2 mutant
# my $ip='SRR1105632';	my $ip_norm_factor=73;
# my $input='SRR1105633'; my $input_norm_factor=163;
############################################################
#sir1
# my $ip='SRR2889257';	my $ip_norm_factor=132;
# my $input='SRR2889258'; my $input_norm_factor=91;
############################################################
#sir4
# my $ip='SRR3355027';	my $ip_norm_factor=25; 
# my $input='SRR1068587'; my $input_norm_factor=57; 
############################################################
#sir2_from Li et al
# my $ip='SRR583850';	my $ip_norm_factor=38; 
# my $input='SRR583856'; my $input_norm_factor=176; 
############################################################


my $enrichment; 
#Use the cov subroutine to get arrays for Chromosome, Position and Values  
my ($wg_chr_i, $wg_pos_i, $wg_INPUT) = cov("$input");
my ($wg_chr_ip, $wg_pos_ip, $wg_IP) = cov("$ip");

my @wg_INPUT = @$wg_INPUT;
my @wg_IP = @$wg_IP;
my @wg_chr = @$wg_chr_i;
my @wg_pos = @$wg_pos_i;

#Open results file to write 
open( my $fh, '>', 'wg-result.txt');
#Loop through the arrays
my $len = @wg_INPUT;
for (my $i=0; $i < $len; $i++){
			#Normalise the read value for INPUT and IP
			$wg_INPUT[$i] = $wg_INPUT[$i] / $input_norm_factor;
			$wg_IP[$i] = $wg_IP[$i] / $input_norm_factor;
			#Calculate the enrichment
			if ( $wg_INPUT[$i] != 0){
						$enrichment = $wg_IP[$i] / $wg_INPUT[$i];						
						}
			else {
			$enrichment = 0;
			}	
	#Print results to file 		
	print $fh "$wg_chr[$i]\t$wg_pos[$i]\t$wg_INPUT[$i]\t$wg_IP[$i]\t$enrichment\n";				
}
close $fh;

#reformat output to allow sorting by Chromosome number
system('cat wg-result.txt |grep -v Mito |  sed "s|^I	|1	I	|g" |  sed "s|^II	|2	II	|g" | sed "s|^III	|3	III	|g" |  sed "s|^IV	|4	IV	|g" |  sed "s|^V	|5	V	|g" |  sed "s|^VI	|6	VI	|g" |  sed "s|^VII	|7	VII	|g" |  sed "s|^VIII	|8	VIII	|g" |  sed "s|^IX	|9	IX	|g" |  sed "s|^X	|10	X	|g" |  sed "s|^XI	|11	XI	|g" |  sed "s|^XII	|12	XII	|g" |  sed "s|^XIII	|13	XIII	|g" |  sed "s|^XIV	|14	XIV	|g" |  sed "s|^XV	|15	XV	|g" |  sed "s|^XVI	|16	XVI	|g" >wg-resu.txt; sort -n -k1 -n -k3 wg-resu.txt > wg-results.txt');


# my @ty1_INPUT = cov('YDRWTy1-5:1-11000',"$input");
# my @ty1_IP = cov('YDRWTy1-5:1-11000',"$ip");
# 
# my @ty2_INPUT = cov('YLRWTy2-1:1-11000',"$input");
# my @ty2_IP = cov('YLRWTy2-1:1-11000',"$ip");
# 
# my @ty3_INPUT = cov('YILWTy3-1:1-11000',"$input");
# my @ty3_IP = cov('YILWTy3-1:1-11000',"$ip");
# 
# my @ty4_INPUT = cov('YHLWTy4-1:1-11000',"$input");
# my @ty4_IP = cov('YHLWTy4-1:1-11000',"$ip");
# 
# my @ty5_INPUT = cov('YCLWTy5-1:1-11000',"$input");
# my @ty5_IP = cov('YCLWTy5-1:1-11000',"$ip");
# 
# open(my $fh, '>', 'ty1-results.txt');
# my $len = @ty1_INPUT;
# for (my $i=0; $i < $len; $i++){
# 			print $fh "$ty1_INPUT[$i]\t$ty1_IP[$i]\n"
# }
# close $fh;
# 
# open( $fh, '>', 'ty2-results.txt');
# $len = @ty2_INPUT;
# for (my $i=0; $i < $len; $i++){
# 			print $fh "$ty2_INPUT[$i]\t$ty2_IP[$i]\n"
# }
# close $fh;
# 
# open( $fh, '>', 'ty3-results.txt');
# $len = @ty3_INPUT;
# for (my $i=0; $i < $len; $i++){
# 			print $fh "$ty3_INPUT[$i]\t$ty3_IP[$i]\n"
# }
# close $fh;
# 
# open( $fh, '>', 'ty4-results.txt');
# $len = @ty4_INPUT;
# for (my $i=0; $i < $len; $i++){
# 			print $fh "$ty4_INPUT[$i]\t$ty4_IP[$i]\n"
# }
# close $fh;
# 
# open( $fh, '>', 'ty5-results.txt');
# $len = @ty5_INPUT;
# for (my $i=0; $i < $len; $i++){
# 			print $fh "$ty5_INPUT[$i]\t$ty5_IP[$i]\n"
# }
# close $fh;

sub cov{
	my $file= $_[0];
	my @values;
	my @chromosome;
	my @position;
	my $ifh;
	open($ifh,"samtools view -q10 -b $file"."_marked.bam | genomeCoverageBed -d -ibam stdin -g |") || die "Failed: $!\n";
	while (my $l = <$ifh> )
	{
  	chomp( $l );    
  	my @s = split( /\t/, $l );
  	push(@values, $s[2]);
  	push(@chromosome, $s[0]);
  	push(@position, $s[1]);
	}
    return (\@chromosome, \@position, \@values);
}
