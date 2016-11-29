#!/usr/bin/env perl

# Author:       mh23
# Maintainer:   mh23
# Created: 		Oktober 2015
# Name:			Ty1-5Estimator.pl
# Test: 		perl /lustre/scratch109/sanger/mh23/Transposons/Ty1-5Estimator.pl -i SC_MFY5616297.merged.bam.Ty.bam



use Carp;
use strict;
use warnings;
use Getopt::Long;

###############
## Get input ##
###############


#Get the bam file as input
my ($input);

GetOptions
(
'i|input=s'         => \$input,
);

( $input && -f $input ) or die qq[Usage: $0 -i <input vcf>\n];


##########################
## Get name of Bam file ##
##########################
my $bam_file='';
my $path = $input;

my @n = split( /\t/, $path );

foreach my $frag (@n){
	if ($frag =~ /.bam/){
		$bam_file = $frag;
	}
}


######################
## Define variables ##
######################
my $ty1; my $ty2;
my $ty3; my $ty4;
my $ty5;



#########
## TY1 ##
#########

#>YDRWTy1-5  Chr 4   from 1203704 to 1215621  
#Control:   YDRWTy1-5:100-2699
#Ty1:       YDRWTy1-5:4000-6999

my $ty1_ctrl = 'YDRWTy1-5:100-2699';
my $ty1_reg = 'YDRWTy1-5:4000-6999';

my @control1 = mpile($ty1_ctrl);
my @sample1 = mpile($ty1_reg);

my $av_c1 = average(\@control1);
my $av_s1 = average(\@sample1);



#########
## TY2 ##
#########

#>YLRWTy2-1  Chr 12   from 938192 to 950150  
#Control:   YLRWTy2-1:100-2699
#Ty2:       YLRWTy2-1:4000-6999

my $ty2_ctrl = 'YLRWTy2-1:100-2699';
my $ty2_reg = 'YLRWTy2-1:4000-6999';

my @control2 = mpile($ty2_ctrl);
my @sample2 = mpile($ty2_reg);

my $av_c2 = average(\@control2);
my $av_s2 = average(\@sample2);



#########
## TY3 ##
#########

#>YILWTy3-1  Chr 9   from 202220 to 213647  
#Control:   YILWTy3-1:100-2699
#Ty3:       YILWTy3-1:4000-6999

my $ty3_ctrl = 'YILWTy3-1:100-2699';
my $ty3_reg = 'YILWTy3-1:4000-6999';

my @control3 = mpile($ty3_ctrl);
my @sample3 = mpile($ty3_reg);

my $av_c3 = average(\@control3);
my $av_s3 = average(\@sample3);



#########
## TY4 ##
#########

#>YHLWTy4-1  Chr 8   from 82539 to 94761  
#Control: YHLWTy4-1:100-2499
#Ty4:      YHLWTy4-1:4000-6999

my $ty4_ctrl = 'YHLWTy4-1:100-2699';
my $ty4_reg = 'YHLWTy4-1:4000-6999';

my @control4 = mpile($ty4_ctrl);
my @sample4 = mpile($ty4_reg);

my $av_c4 = average(\@control4);
my $av_s4 = average(\@sample4);




#########
## TY5 ##
#########

#>YCLWTy5-1  Chr 3   from 679 to 7322 
#Control:  YCLWTy5-1:2000-2999
#Ty5:      YCLWTy5-1:5000-5999

my $ty5_reg = 'YCLWTy5-1:2000-2999';

my $ty5_ctrl = 'YCLWTy5-1:5000-5999';

my @control5 = mpile($ty5_ctrl);
my @sample5 = mpile($ty5_reg);

my $av_c5 = average(\@control5);
my $av_s5 = average(\@sample5);



########################
## Calculate averages ##
########################

my @controls;

if ($av_c1 != 0) {push @controls, $av_c1; }
if ($av_c2 != 0) {push @controls, $av_c2; }
if ($av_c3 != 0) {push @controls, $av_c3; }
if ($av_c4 != 0) {push @controls, $av_c4; }


my $control_mean = average(\@controls);

if ($control_mean != 0){
	$ty1 = sprintf("%.0f",$av_s1/$control_mean);
	$ty2 = sprintf("%.0f",$av_s2/$control_mean);
	$ty3 = sprintf("%.0f",$av_s3/$control_mean);
	$ty4 = sprintf("%.0f",$av_s4/$control_mean);
	$ty5 = sprintf("%.0f",$av_s5/$control_mean); }
else {$ty1 = 0;$ty2 = 0;$ty3 = 0;$ty4 = 0;$ty5 = 0;}


######################
## print the result ##
######################
#print ("@sample5\n");

#print ("Controls\t $av_c1  $av_c2  $av_c3  $av_c4 $av_c5\n");
#print ("Averages\t $av_s1  $av_s2  $av_s3  $av_s4 $av_s5\n");

print ("\t$ty1\t$ty2\t$ty3\t$ty4\t$ty5\n");






###########################################
## subroutine to calculate mean coverage ##
###########################################

sub average{
		my @numbers = @{$_[0]};
		my $data = scalar @numbers;
        if (not $data) {return 0;}
        my $total = 0;
        foreach my $n (@numbers) {
        	#print "$n\n";
        	$total = $total + $n;
        }
       my $average = $total / $data;
       return $average;
}


##############################################
## use samtools mpileup to get the coverage ##
##############################################

sub mpile {
	my $region = $_[0];
	my @values;
	my $ifh;
	open($ifh,"samtools mpileup  -A -Q 0 -r $region -f ~/sw/bin/PF/mpileup_defaults/Ty_ref/Ty1-5.fa $input 2>/dev/null |") || die "Failed: $!\n";
	while (my $l = <$ifh> )
	{
  	chomp( $l );    
  	my @s = split( /\t/, $l );
  	push(@values, $s[3]);
	}
 return @values;
}



























