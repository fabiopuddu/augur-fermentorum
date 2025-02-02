#!/usr/bin/env perl
#
# Author:       Mareike Herzog
# Maintainer:   Fabio Puddu
# Created: 		23.02.2015

use Carp;
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use List::Util qw/first/;

#Declare variables
my @s; my $sample; my $num;
my $indel_count;
my $snp_count;
$indel_count='0';$snp_count='0';
my $check=0;
my $transitions=0; my $transversions=0;
my @number; my %result;my $number_of_samples;
my $A_C=0;my $A_G=0;my $A_T=0;
my $C_A=0;my $C_G=0;my $C_T=0;
my $G_A=0;my $G_C=0;my $G_T=0;
my $T_A=0;my $T_C=0;my $T_G=0;
#INDELs
my $one=0;my $two=0;my $three=0;my $more=0;
my $minone=0;my $mintwo=0;my $minthree=0;my $less=0;


###################################  CHANGES TO ACCEPT POLYPLOIDS  ###################################
##      1. Identify ploidy by amount of '/' fileds in GT of vcf file.
##      2. According to the ploidy, grep for homozygous and heterozygous accordingly
##      3. Replace vcf-stats with bcftools stats
##      4. Change the output parser, as the output of bcftools stats is different.
##         (Should be easier for global statistics)
######################################################################################################

#go through file line by line
#open my $fh, "<", $input or die $!;
my $input;
GetOptions
(
    'i|input=s'         => \$input,
);

( $input && -f $input ) or die qq[Usage: $0 -i <input vcf>\n];

my $line;my $ploidy = 0;
# Get DelN from the input string (which should contain the DelFolder)
my @samplein=split( /\//, $input);
my $deln=first { $_ =~ 'Del[0-9]+_.+' } @samplein;

# Determine the ploidy of the sample by looking at the GT field of the vcf file
my $ifh;
if( $input =~ /\.gz/ ){open($ifh, qq[gunzip -c $input|]);}else{open($ifh, $input ) or die $!;}
my @FILE=<$ifh>;
chomp @FILE;
my $number_of_columns = scalar(split("\t", $FILE[-1])); #get the number of columns by counting how many fields in the last line of the vcf file
for (my $column=9; $column<$number_of_columns; $column++){
	for my $line(@FILE) {
    		next if ($line =~ /#/);
    		next unless ($line =~ /\//);
    		my @fields = split("\t", $line);
    		# Check in the columns corresponding to the first 5 samples, to make sure that at least one of them has some variant (controls may not have any)
    		my $select = $fields[$column] unless $fields[$column] eq '.' || $fields[$column] =~ /\.:/;
	        if (defined $select){
			my @gt = split(':', $select);
    			my @gtsep = split('/', $gt[0]);
    			$ploidy = scalar @gtsep;
		}
	}
last if ($ploidy > 0);
}
close($ifh);

my $command='';
if ($ploidy eq 2){
    #$command =  "cat $input"." | grep ".'"'.'1/1" -v'." | vcf-stats";
    $command = "pigz -dc $input". '| grep "[0-9]/[0-9]\|#" | grep "PASS\|#" | bcftools stats'; #edited by FABIO to resolve bug that does eliminates rows where one of the mutation is masked
}
if ($ploidy eq 4){
    $command = "pigz -dc $input".' | grep "[0-9]/[0-9]/[0-9]/[0-9]\|#" | grep "PASS\|#" | bcftools stats'; #edited by FABIO to resolve bug that does eliminates rows where one of the mutation is masked
}
my @stats_out =  readpipe("$command");
print("INDEL\tSNP\tTs\tTv\tSampleName\tC>T\tA>G\tA>T\tC>G\tG>T\tA>C\t<-3\t-3\t-2\t-1\t1\t2\t3\t>3\tNumber of samples\n");
foreach my $l (@stats_out){
    my $k = $l;
    chomp( $l );
    next if ($l =~ /#/);
        if ($l =~ /SN/){
            if ($l =~ /number of indels/){
                @s = split("\t", $l); $indel_count=$s[3];
            }
            if ($l =~ /number of SNP/){
                @s = split("\t", $l); $snp_count=$s[3];
            }
            if ($l =~ /number of samples/){
                @s = split("\t", $l); $number_of_samples=$s[3];
            }
    }
    if ($l =~ /TSTV/){
        @s = split("\t", $l); $transitions=$s[2]; $transversions=$s[3];
    }
    if ($l =~ /ST/){
        if ($l =~ /A>C/){
        @s = split("\t", $l); $A_C=$s[3];
        }
        if ($l =~ /A>G/){
            @s = split("\t", $l); $A_G=$s[3];
        }
        if ($l =~ /A>T/){
            @s = split("\t", $l); $A_T=$s[3];
        }
        if ($l =~ /C>A/){
            @s = split("\t", $l); $C_A=$s[3];
        }
        if ($l =~ /C>G/){
            @s = split("\t", $l); $C_G=$s[3];
        }
        if ($l =~ /C>T/){
            @s = split("\t", $l); $C_T=$s[3];
        }
        if ($l =~ /G>A/){
            @s = split("\t", $l); $G_A=$s[3];
        }
        if ($l =~ /G>C/){
            @s = split("\t", $l); $G_C=$s[3];
        }
        if ($l =~ /G>T/){
            @s = split("\t", $l); $G_T=$s[3];
        }
        if ($l =~ /T>A/){
            @s = split("\t", $l); $T_A=$s[3];
        }
        if ($l =~ /T>C/){
            @s = split("\t", $l); $T_C=$s[3];
        }
        if ($l =~ /T>G/){
            @s = split("\t", $l); $T_G=$s[3];
        }
    }
    if ($l =~ /IDD/){
        @s = split("\t", $l);
        $num = $s[2];
        if ($num == 1){$one = $s[3];}
        if ($num == 2){$two = $s[3];}
        if ($num == 3){$three = $s[3];}
        if ($num == -1){$minone = $s[3];}
        if ($num == -2){$mintwo = $s[3];}
        if ($num == -3 ){$minthree = $s[3];}
        if ($num > 3) {$more = $more + $s[3];}
        if ($num < -3) {$less = $less + $s[3];}
    }
}# close

#transitions
$C_T =    $C_T + $G_A;
$A_G =    $A_G + $T_C;
#transversions
$A_T =    $A_T + $T_A;
$C_G =    $C_G + $G_C;
$G_T =  $G_T + $C_A;
$A_C =    $A_C + $T_G;

#$transitions=$A_G+$C_T;
#$transversions=$A_T+$A_C+$C_G+$G_T;
print("$indel_count\t$snp_count\t$transitions\t$transversions\t$deln\t$C_T\t$A_G\t$A_T\t$C_G\t$G_T\t$A_C\t");
print("$less\t$minthree\t$mintwo\t$minone\t$one\t$two\t$three\t$more\t$number_of_samples\n");
