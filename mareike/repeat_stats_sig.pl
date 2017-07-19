#!/usr/bin/env perl

# Author:       mh23
# Maintainer:   mh23
# Created: 		Nov 2016
# Name:			repeat_stats_sig.pl

use autodie;
use utf8;
use Carp;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Pod::Usage;
use Cwd 'abs_path';
use POSIX qw(ceil);
use POSIX qw(floor);
use IPC::Open3;
use IO::File;
use Math::Complex qw(sqrt);



#######################################
##    								 ##
##				USAGE				 ##
##									 ##
#######################################

#Get the input and print usage if help of wrong input

#DEFAULTS
my $input = '0';
my $help = 0;
#Get the location of the reference genome
my $script_location = abs_path($0);


## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions('help|?' => \$help, 
		   'i|input=s' => \$input,
		   ) or pod2usage(2);
pod2usage(1) if $help;
## If no input argument were given, then allow STDIN to be used only
## if it's not connected to a terminal (otherwise print usage)
pod2usage("$0: No input given.")  if (($input eq 0) && (-t STDIN));

#Check that input exists and has a size bigger than 0		
pod2usage("$0: File $input does not exist.")  unless ( -e $input);
pod2usage("$0: File $input is empty.")  if ( -z $input);



#######################################
##    								 ##
##				WT STATS			 ##
##									 ##
#######################################

###################### Hashes ##########################

my %averages;
my %Stderr;
my %StDev;
my %CI_lower;
my %CI_higher;


### Step1

# Get wild type stats and numbers


# Sample	rDNA	CUP1	Ty1		Ty2		Ty3		Ty4		Ty5		Genome_wide_median	 ERS	 Del

##Get the average
# open(PS,"cat $input | grep 'WT' |") || die "Failed: $!\n";
# while ( <PS> )
# {
#  my $l = $_;
#  chomp $l;
#   
# }

###################### AVERAGE ##########################

my $command = "cat $input | grep 'WT-' | awk '".'{ sum += $2; n++ } END { if (n > 0) print'." sum / n; }'";
my $rDNA_av = `$command`;

$command = "cat $input | grep 'WT-' | awk '".'{ sum += $3; n++ } END { if (n > 0) print'." sum / n; }'";
my $Cup1_av = `$command`;

$command = "cat $input | grep 'WT-' | awk '".'{ sum += $4; n++ } END { if (n > 0) print'." sum / n; }'";
my $Mito_av = `$command`; 
 
$command = "cat $input | grep 'WT-' | awk '".'{ sum += $5; n++ } END { if (n > 0) print'." sum / n; }'";
my $Ty1_av = `$command`;

$command = "cat $input | grep 'WT-' | awk '".'{ sum += $6; n++ } END { if (n > 0) print'." sum / n; }'";
my $Ty2_av = `$command`;

$command = "cat $input | grep 'WT-' | awk '".'{ sum += $7; n++ } END { if (n > 0) print'." sum / n; }'";
my $Ty3_av = `$command`;

$command = "cat $input | grep 'WT-' | awk '".'{ sum += $8; n++ } END { if (n > 0) print'." sum / n; }'";
my $Ty4_av = `$command`;

$command = "cat $input | grep 'WT-' | awk '".'{ sum += $9; n++ } END { if (n > 0) print'." sum / n; }'";
my $Ty5_av = `$command`;

$command = "cat $input | grep 'WT-' | awk '".'{ sum += $10; n++ } END { if (n > 0) print'." sum / n; }'";
my $gwm_av = `$command`;

$command = "cat $input | grep 'WT-' | awk '".'{ sum += $11; n++ } END { if (n > 0) print'." sum / n; }'";
my $tel_av = `$command`;

#chomp all of those
chomp $rDNA_av; chomp $Cup1_av; chomp $Mito_av; chomp $Ty1_av; chomp $Ty2_av; chomp $Ty3_av; chomp $Ty4_av; chomp $Ty5_av; chomp $gwm_av; chomp $tel_av;

$averages{'rDNA'}=$rDNA_av;
$averages{'CUP1'}=$Cup1_av;
$averages{'Mito'}=$Mito_av;
$averages{'Ty1'}=$Ty1_av;
$averages{'Ty2'}=$Ty2_av;
$averages{'Ty3'}=$Ty3_av;
$averages{'Ty4'}=$Ty4_av;
$averages{'Ty5'}=$Ty5_av;
$averages{'gwm'}=$gwm_av;
$averages{'Tel'}=$tel_av;


###################### STD-Err ##########################

#rDNA
$command = "cat $input | grep 'WT-' | awk '". '{sum+=$2; array[NR]=$2} END '."{for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'"; 
my $stdev = `$command`; chomp $stdev;
$command = "cat $input | grep 'WT-' | wc -l";
my $no_samp = `$command`; chomp $no_samp;
my $root = sqrt($no_samp);

$StDev{'rDNA'}=$stdev;
my $rDNA_stderr = $stdev / $root;

#CUP1
$command = "cat $input | grep 'WT-' | awk '". '{sum+=$3; array[NR]=$3} END '."{for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'"; 
$stdev = `$command`; chomp $stdev;
$command = "cat $input | grep 'WT-' | wc -l";
$no_samp = `$command`; chomp $no_samp;
$root = sqrt($no_samp);

$StDev{'CUP1'}=$stdev;
my $Cup1_stderr = $stdev / $root;


#Mito
$command = "cat $input | grep 'WT-' | awk '". '{sum+=$4; array[NR]=$4} END '."{for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'"; 
$stdev = `$command`; chomp $stdev;
$command = "cat $input | grep 'WT-' | wc -l";
$no_samp = `$command`; chomp $no_samp;
$root = sqrt($no_samp);

$StDev{'Mito'}=$stdev;
my $Mito_stderr = $stdev / $root;



#Ty1
$command = "cat $input | grep 'WT-' | awk '". '{sum+=$5; array[NR]=$5} END '."{for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'"; 
$stdev = `$command`; chomp $stdev;
$command = "cat $input | grep 'WT-' | wc -l";
$no_samp = `$command`; chomp $no_samp;
$root = sqrt($no_samp);

$StDev{'Ty1'}=$stdev;
my $Ty1_stderr = $stdev / $root;

#Ty2
$command = "cat $input | grep 'WT-' | awk '". '{sum+=$6; array[NR]=$6} END '."{for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'"; 
$stdev = `$command`; chomp $stdev;
$command = "cat $input | grep 'WT-' | wc -l";
$no_samp = `$command`; chomp $no_samp;
$root = sqrt($no_samp);

$StDev{'Ty2'}=$stdev;
my $Ty2_stderr = $stdev / $root;

#Ty3
$command = "cat $input | grep 'WT-' | awk '". '{sum+=$7; array[NR]=$7} END '."{for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'"; 
$stdev = `$command`; chomp $stdev;
$command = "cat $input | grep 'WT-' | wc -l";
$no_samp = `$command`; chomp $no_samp;
$root = sqrt($no_samp);

$StDev{'Ty3'}=$stdev;
my $Ty3_stderr = $stdev / $root;

#Ty4
$command = "cat $input | grep 'WT-' | awk '". '{sum+=$8; array[NR]=$8} END '."{for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'"; 
$stdev = `$command`; chomp $stdev;
$command = "cat $input | grep 'WT-' | wc -l";
$no_samp = `$command`; chomp $no_samp;
$root = sqrt($no_samp);

$StDev{'Ty4'}=$stdev;
my $Ty4_stderr = $stdev / $root;

#Ty5
$command = "cat $input | grep 'WT-' | awk '". '{sum+=$9; array[NR]=$9} END '."{for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'"; 
$stdev = `$command`; chomp $stdev;
$command = "cat $input | grep 'WT-' | wc -l";
$no_samp = `$command`; chomp $no_samp;
$root = sqrt($no_samp);

$StDev{'Ty5'}=$stdev;
my $Ty5_stderr = $stdev / $root;

#Genome-wide median
$command = "cat $input | grep 'WT-' | awk '". '{sum+=$10; array[NR]=$10} END '."{for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'"; 
$stdev = `$command`; chomp $stdev;
$command = "cat $input | grep 'WT-' | wc -l";
$no_samp = `$command`; chomp $no_samp;
$root = sqrt($no_samp);

$StDev{'gwm'}=$stdev;
my $gwm_stderr = $stdev / $root;

#Telomeres
$command = "cat $input | grep 'WT-' | awk '". '{sum+=$11; array[NR]=$11} END '."{for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'"; 
#$command = "cat $input | grep 'WT-' | awk '{sum+=$11; sumsq+=$11*$11} END {print sqrt(sumsq/NR - (sum/NR)**2)}'";
$stdev = `$command`; chomp $stdev;
$command = "cat $input | grep 'WT-' | wc -l";
$no_samp = `$command`; chomp $no_samp;
$root = sqrt($no_samp);

$StDev{'Tel'}=$stdev;
my $Tel_stderr = $stdev / $root;

#print "$no_samp, $stdev and $Tel_stderr\n";  

#$Stderr{'rDNA'} =  $rDNA_stderr;
#$Stderr{'CUP1'} =  $Cup1_stderr;
#$Stderr{'Mito'} =  $Mito_stderr;
#$Stderr{'Ty1'} =   $Ty1_stderr;
#$Stderr{'Ty2'} =   $Ty2_stderr;
#$Stderr{'Ty3'} =   $Ty3_stderr;
#$Stderr{'Ty4'} =   $Ty4_stderr;
#$Stderr{'Ty5'} =   $Ty5_stderr;
#$Stderr{'gwm'} =   $gwm_stderr;
#$Stderr{'Tel'} =   $Tel_stderr;



my @replist=('rDNA','CUP1', 'Mito', 'Ty1', 'Ty2', 'Ty3', 'Ty4', 'Ty5', 'gwm', 'Tel');
print "Repeat\tWTave\tWTStDev\tConfidence intervals\n";
foreach my $rep (@replist){
		printf ("%s\t%.2f\t%.2f\t(%.1f - %.1f)\n", $rep,$averages{$rep},$StDev{$rep},$averages{$rep}-3*$StDev{$rep},$averages{$rep}+3*$StDev{$rep});
}








###################### CI ##########################

my $multiplier = 5;

### rDNA

my $CI_rDNA = $multiplier * $rDNA_stderr;
my $CI_rDNA_lower = floor($rDNA_av-$CI_rDNA);
my $CI_rDNA_higher = ceil($rDNA_av+$CI_rDNA);
#print "CI: $CI_rDNA; lower: $CI_rDNA_lower; higher: $CI_rDNA_higher\n ";


### CUP1

my $CI_Cup1 = $multiplier * $Cup1_stderr;
my $CI_Cup1_lower = floor($Cup1_av-$CI_Cup1);
my $CI_Cup1_higher = ceil($Cup1_av+$CI_Cup1);
#print "CI: $CI_Cup1; lower: $CI_Cup1_lower; higher: $CI_Cup1_higher\n ";

### Mito

my $CI_Mito = $multiplier * $Mito_stderr;
my $CI_Mito_lower = floor($Mito_av-$CI_Cup1);
my $CI_Mito_higher = ceil($Mito_av+$CI_Cup1);
#print "CI: $CI_Mito; lower: $CI_Mito_lower; higher: $CI_Mito_higher\n ";


### Ty1

my $CI_Ty1 = $multiplier * $Ty1_stderr;
my $CI_Ty1_lower = floor($Ty1_av-$CI_Ty1);
my $CI_Ty1_higher = ceil($Ty1_av+$CI_Ty1);
#print "CI: $CI_Ty1; lower: $CI_Ty1_lower; higher: $CI_Ty1_higher\n ";


### Ty2

my $CI_Ty2 = $multiplier * $Ty2_stderr;
my $CI_Ty2_lower = floor($Ty2_av-$CI_Ty2);
my $CI_Ty2_higher = ceil($Ty2_av+$CI_Ty2);
#print "CI: $CI_Ty2; lower: $CI_Ty2_lower; higher: $CI_Ty2_higher\n ";

### Ty3

my $CI_Ty3 = $multiplier * $Ty3_stderr;
my $CI_Ty3_lower = floor($Ty3_av-$CI_Ty3);
my $CI_Ty3_higher = ceil($Ty3_av+$CI_Ty3);
#print "CI: $CI_Ty3; lower: $CI_Ty3_lower; higher: $CI_Ty3_higher\n ";

### Ty4

my $CI_Ty4 = $multiplier * $Ty4_stderr;
my $CI_Ty4_lower = floor($Ty4_av-$CI_Ty4);
my $CI_Ty4_higher = ceil($Ty4_av+$CI_Ty4);
#print "CI: $CI_Ty4; lower: $CI_Ty4_lower; higher: $CI_Ty4_higher\n ";

### Ty5

my $CI_Ty5 = $multiplier * $Ty5_stderr;
my $CI_Ty5_lower = floor($Ty5_av-$CI_Ty5);
my $CI_Ty5_higher = ceil($Ty5_av+$CI_Ty5);
#print "CI: $CI_Ty5; lower: $CI_Ty5_lower; higher: $CI_Ty5_higher\n ";

### Genome wide median

my $CI_gwm = $multiplier * $gwm_stderr;
my $CI_gwm_lower = floor($gwm_av-$CI_gwm);
my $CI_gwm_higher = ceil($gwm_av+$CI_gwm);
#print "CI: $CI_gwm; lower: $CI_gwm_lower; higher: $CI_gwm_higher\n ";

### Telomeres

my $CI_Tel = $multiplier * $Tel_stderr;
my $CI_Tel_lower = floor($tel_av-$CI_Tel);
my $CI_Tel_higher = ceil($tel_av+$CI_Tel);
#print "CI: $CI_Tel; lower: $CI_Tel_lower; higher: $CI_Tel_higher\n ";




#### If you instead want the CI to be 3* the standard deviation
$CI_rDNA_lower = floor($rDNA_av-($StDev{'rDNA'}*3));
$CI_rDNA_higher = ceil($rDNA_av+($StDev{'rDNA'}*3));
$CI_Cup1_lower = floor($Cup1_av-($StDev{'CUP1'}*3));
$CI_Cup1_higher = ceil($Cup1_av+($StDev{'CUP1'}*3));
$CI_Mito_lower = floor($Mito_av-($StDev{'Mito'}*3));
$CI_Mito_higher = ceil($Mito_av+($StDev{'Mito'}*3));
$CI_Ty1_lower = floor($Ty1_av-($StDev{'Ty1'}*3));
$CI_Ty1_higher = ceil($Ty1_av+($StDev{'Ty1'}*3));
$CI_Ty2_lower = floor($Ty2_av-($StDev{'Ty2'}*3));
$CI_Ty2_higher = ceil($Ty2_av+($StDev{'Ty2'}*3));
$CI_Ty3_lower = floor($Ty3_av-($StDev{'Ty3'}*3));
$CI_Ty3_higher = ceil($Ty3_av+($StDev{'Ty3'}*3));
$CI_Ty4_lower = floor($Ty4_av-($StDev{'Ty4'}*3));
$CI_Ty4_higher = ceil($Ty4_av+($StDev{'Ty4'}*3));
$CI_Ty5_lower = floor($Ty5_av-($StDev{'Ty5'}*3));
$CI_Ty5_higher = ceil($Ty5_av+($StDev{'Ty5'}*3));
$CI_gwm_lower = floor($gwm_av-($StDev{'gwm'}*3));
$CI_gwm_higher = ceil($gwm_av+($StDev{'gwm'}*3));
$CI_Tel_lower = floor($tel_av-($StDev{'Tel'}*3));
$CI_Tel_higher = ceil($tel_av+($StDev{'Tel'}*3));


#### Save in CI hashes
$CI_lower{'rDNA'}  = $CI_rDNA_lower;
$CI_higher{'rDNA'} = $CI_rDNA_higher;
$CI_lower{'CUP1'}  = $CI_Cup1_lower;
$CI_higher{'CUP1'} = $CI_Cup1_higher;
$CI_lower{'Mito'}  = $CI_Mito_lower;
$CI_higher{'Mito'} = $CI_Mito_higher;
$CI_lower{'Ty1'}  = $CI_Ty1_lower;
$CI_higher{'Ty1'} = $CI_Ty1_higher;
$CI_lower{'Ty2'}  = $CI_Ty2_lower;
$CI_higher{'Ty2'} = $CI_Ty2_higher;
$CI_lower{'Ty3'}  = $CI_Ty3_lower;
$CI_higher{'Ty3'} = $CI_Ty3_higher;
$CI_lower{'Ty4'}  = $CI_Ty4_lower;
$CI_higher{'Ty4'} = $CI_Ty4_higher;
$CI_lower{'Ty5'}  = $CI_Ty5_lower;
$CI_higher{'Ty5'} = $CI_Ty5_higher;
$CI_lower{'gwm'}  = $CI_gwm_lower;
$CI_higher{'gwm'} = $CI_gwm_higher;
$CI_lower{'Tel'}  = $CI_Tel_lower;
$CI_higher{'Tel'} = $CI_Tel_higher;



#######################################
##    								 ##
##			LOOP THRU FILE			 ##
##									 ##
#######################################

###################### Hashes ##########################

my %samples;
my %rDNA_Sig_low; my %rDNA_Sig_high;
my %CUP1_Sig_low; my %CUP1_Sig_high;
my %Mito_Sig_low; my %Mito_Sig_high;
my %Ty1_Sig_low; my %Ty1_Sig_high;
my %Ty2_Sig_low; my %Ty2_Sig_high;
my %Ty3_Sig_low; my %Ty3_Sig_high;
my %Ty4_Sig_low; my %Ty4_Sig_high;
my %Ty5_Sig_low; my %Ty5_Sig_high;
my %Tel_Sig_low; my %Tel_Sig_high;

my $ifh;
if( $input =~ /\.gz/ ){open($ifh, qq[gunzip -c $input|]);}else{open($ifh, $input ) or die $!;}

#print "Lower threshold: $CI_lower{'rDNA'}\n";

while( my $l = <$ifh> ){
	#read the cover into a hash
	chomp( $l );
    my @s = split( /\t/, $l );
	$samples{$s[0]}=$s[12];
	###################### rDNA ##########################
	if ($s[1] < $CI_lower{'rDNA'})  { 
		$rDNA_Sig_low{$s[0]} =  $s[12];
		#print "$s[0]\t$s[11]\n";
	}
	if ($s[1] > $CI_higher{'rDNA'}) { 
		$rDNA_Sig_high{$s[0]} =  $s[12];
	}
	###################### CUP1 ##########################
	if ($s[2] < $CI_lower{'CUP1'})  { 
		$CUP1_Sig_low{$s[0]} =  $s[12];
		#print "$s[0]\t$s[11]\n";
	}
	if ($s[2] > $CI_higher{'CUP1'}) { 
		$CUP1_Sig_high{$s[0]} =  $s[12];
	} 
	###################### Mito ##########################
	if ($s[3] < $CI_lower{'Mito'})  { 
		$Mito_Sig_low{$s[0]} =  $s[12];
		#print "$s[0]\t$s[11]\n";
	}
	if ($s[3] > $CI_higher{'Mito'}) { 
		$Mito_Sig_high{$s[0]} =  $s[12];
	} 
	###################### Ty1 ##########################
	#print ("$s[3] thresh: $CI_lower{'Ty1'}\n");
	if ($s[4] < $CI_lower{'Ty1'})  { 
		$Ty1_Sig_low{$s[0]} =  $s[12];
		#print "$s[0]\t$s[11]\n";
	}
	if ($s[4] > $CI_higher{'Ty1'}) { 
		$Ty1_Sig_high{$s[0]} =  $s[12];
	}
	###################### Ty2 ##########################
	if ($s[5] < $CI_lower{'Ty2'})  { 
		$Ty2_Sig_low{$s[0]} =  $s[12];
		#print "$s[0]\t$s[11]\n";
	}
	if ($s[5] > $CI_higher{'Ty2'}) { 
		$Ty2_Sig_high{$s[0]} =  $s[12];
	}
	###################### Ty3 ##########################
	if ($s[6] < $CI_lower{'Ty3'})  { 
		$Ty3_Sig_low{$s[0]} =  $s[12];
		#print "$s[0]\t$s[11]\n";
	}
	if ($s[6] > $CI_higher{'Ty3'}) { 
		$Ty3_Sig_high{$s[0]} =  $s[12];
	}
	###################### Ty4 ##########################
	if ($s[7] < $CI_lower{'Ty4'})  { 
		$Ty4_Sig_low{$s[0]} =  $s[12];
		#print "$s[0]\t$s[11]\n";
	}
	if ($s[7] > $CI_higher{'Ty4'}) { 
		$Ty4_Sig_high{$s[0]} =  $s[12];
	}
	###################### Ty5 ##########################
	if ($s[8] < $CI_lower{'Ty5'})  { 
		$Ty5_Sig_low{$s[0]} =  $s[12];
		#print "$s[0]\t$s[11]\n";
	}
	if ($s[8] > $CI_higher{'Ty5'}) { 
		$Ty5_Sig_high{$s[0]} =  $s[12];
	}
	###################### Tel ##########################
	if ($s[10] < $CI_lower{'Tel'})  { 
		$Tel_Sig_low{$s[0]} =  $s[12];
		#print "$s[0]\t$s[11]\n";
	}
	if ($s[10] > $CI_higher{'Tel'}) { 
		$Tel_Sig_high{$s[0]} =  $s[12];
	}
			
}
close( $ifh );

#######################################
##    								 ##
##		  Check sample number		 ##
##									 ##
#######################################


###################### rDNA ##########################

########## LOW

my %rDNA_results_low;
foreach my $SD (keys %rDNA_Sig_low){
	my $gene = $rDNA_Sig_low{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$rDNA_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $rDNA_Sig_low{$SD_2}){$rDNA_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$rDNA_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $rDNA_Sig_low{$SD_2}){$rDNA_results_low{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$rDNA_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $rDNA_Sig_low{$SD_2}){$rDNA_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$rDNA_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $rDNA_Sig_low{$SD_2}){$rDNA_results_low{$gene}='sig_two_samples';}}
	}	
}

########### HIGH

my %rDNA_results_high;
foreach my $SD (keys %rDNA_Sig_high){
	my $gene = $rDNA_Sig_high{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$rDNA_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $rDNA_Sig_high{$SD_2}){$rDNA_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$rDNA_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $rDNA_Sig_high{$SD_2}){$rDNA_results_high{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$rDNA_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $rDNA_Sig_high{$SD_2}){$rDNA_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$rDNA_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $rDNA_Sig_high{$SD_2}){$rDNA_results_high{$gene}='sig_two_samples';}}
	}	
}


###################### CUP1 ##########################

########## LOW

my %CUP1_results_low;
foreach my $SD (keys %CUP1_Sig_low){
	my $gene = $CUP1_Sig_low{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$CUP1_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $CUP1_Sig_low{$SD_2}){$CUP1_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$CUP1_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $CUP1_Sig_low{$SD_2}){$CUP1_results_low{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$CUP1_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $CUP1_Sig_low{$SD_2}){$CUP1_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$CUP1_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $CUP1_Sig_low{$SD_2}){$CUP1_results_low{$gene}='sig_two_samples';}}
	}	
}

########### HIGH

my %CUP1_results_high;
foreach my $SD (keys %CUP1_Sig_high){
	my $gene = $CUP1_Sig_high{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$CUP1_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $CUP1_Sig_high{$SD_2}){$CUP1_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$CUP1_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $CUP1_Sig_high{$SD_2}){$CUP1_results_high{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$CUP1_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $CUP1_Sig_high{$SD_2}){$CUP1_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$CUP1_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $CUP1_Sig_high{$SD_2}){$CUP1_results_high{$gene}='sig_two_samples';}}
	}	
}



###################### Mito ##########################

########## LOW

my %Mito_results_low;
foreach my $SD (keys %Mito_Sig_low){
	my $gene = $Mito_Sig_low{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Mito_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Mito_Sig_low{$SD_2}){$Mito_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Mito_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Mito_Sig_low{$SD_2}){$Mito_results_low{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Mito_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Mito_Sig_low{$SD_2}){$Mito_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Mito_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Mito_Sig_low{$SD_2}){$Mito_results_low{$gene}='sig_two_samples';}}
	}	
}

########### HIGH

my %Mito_results_high;
foreach my $SD (keys %Mito_Sig_high){
	my $gene = $Mito_Sig_high{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Mito_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Mito_Sig_high{$SD_2}){$Mito_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Mito_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Mito_Sig_high{$SD_2}){$Mito_results_high{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Mito_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Mito_Sig_high{$SD_2}){$Mito_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Mito_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Mito_Sig_high{$SD_2}){$Mito_results_high{$gene}='sig_two_samples';}}
	}	
}








###################### Ty1 ##########################

########## LOW

my %Ty1_results_low;
foreach my $SD (keys %Ty1_Sig_low){
	my $gene = $Ty1_Sig_low{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Ty1_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty1_Sig_low{$SD_2}){$Ty1_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Ty1_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty1_Sig_low{$SD_2}){$Ty1_results_low{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Ty1_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty1_Sig_low{$SD_2}){$Ty1_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Ty1_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty1_Sig_low{$SD_2}){$Ty1_results_low{$gene}='sig_two_samples';}}
	}	
}

########### HIGH

my %Ty1_results_high;
foreach my $SD (keys %Ty1_Sig_high){
	my $gene = $Ty1_Sig_high{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Ty1_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty1_Sig_high{$SD_2}){$Ty1_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Ty1_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty1_Sig_high{$SD_2}){$Ty1_results_high{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Ty1_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty1_Sig_high{$SD_2}){$Ty1_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Ty1_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty1_Sig_high{$SD_2}){$Ty1_results_high{$gene}='sig_two_samples';}}
	}	
}



###################### Ty2 ##########################

########## LOW

my %Ty2_results_low;
foreach my $SD (keys %Ty2_Sig_low){
	my $gene = $Ty2_Sig_low{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Ty2_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty2_Sig_low{$SD_2}){$Ty2_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Ty2_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty2_Sig_low{$SD_2}){$Ty2_results_low{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Ty2_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty2_Sig_low{$SD_2}){$Ty2_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Ty2_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty2_Sig_low{$SD_2}){$Ty2_results_low{$gene}='sig_two_samples';}}
	}	
}

########### HIGH

my %Ty2_results_high;
foreach my $SD (keys %Ty2_Sig_high){
	my $gene = $Ty2_Sig_high{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Ty2_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty2_Sig_high{$SD_2}){$Ty2_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Ty2_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty2_Sig_high{$SD_2}){$Ty2_results_high{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Ty2_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty2_Sig_high{$SD_2}){$Ty2_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Ty2_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty2_Sig_high{$SD_2}){$Ty2_results_high{$gene}='sig_two_samples';}}
	}	
}


###################### Ty3 ##########################

########## LOW

my %Ty3_results_low;
foreach my $SD (keys %Ty3_Sig_low){
	my $gene = $Ty3_Sig_low{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Ty3_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty3_Sig_low{$SD_2}){$Ty3_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Ty3_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty3_Sig_low{$SD_2}){$Ty3_results_low{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Ty3_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty3_Sig_low{$SD_2}){$Ty3_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Ty3_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty3_Sig_low{$SD_2}){$Ty3_results_low{$gene}='sig_two_samples';}}
	}	
}

########### HIGH

my %Ty3_results_high;
foreach my $SD (keys %Ty3_Sig_high){
	my $gene = $Ty3_Sig_high{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Ty3_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty3_Sig_high{$SD_2}){$Ty3_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Ty3_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty3_Sig_high{$SD_2}){$Ty3_results_high{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Ty3_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty3_Sig_high{$SD_2}){$Ty3_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Ty3_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty3_Sig_high{$SD_2}){$Ty3_results_high{$gene}='sig_two_samples';}}
	}	
}



###################### Ty4 ##########################

########## LOW

my %Ty4_results_low;
foreach my $SD (keys %Ty4_Sig_low){
	my $gene = $Ty4_Sig_low{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Ty4_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty4_Sig_low{$SD_2}){$Ty4_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Ty4_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty4_Sig_low{$SD_2}){$Ty4_results_low{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Ty4_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty4_Sig_low{$SD_2}){$Ty4_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Ty4_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty4_Sig_low{$SD_2}){$Ty4_results_low{$gene}='sig_two_samples';}}
	}	
}

########### HIGH

my %Ty4_results_high;
foreach my $SD (keys %Ty4_Sig_high){
	my $gene = $Ty4_Sig_high{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Ty4_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty4_Sig_high{$SD_2}){$Ty4_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Ty4_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty4_Sig_high{$SD_2}){$Ty4_results_high{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Ty4_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty4_Sig_high{$SD_2}){$Ty4_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Ty4_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty4_Sig_high{$SD_2}){$Ty4_results_high{$gene}='sig_two_samples';}}
	}	
}



###################### Ty5 ##########################

########## LOW

my %Ty5_results_low;
foreach my $SD (keys %Ty5_Sig_low){
	my $gene = $Ty5_Sig_low{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Ty5_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty5_Sig_low{$SD_2}){$Ty5_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Ty5_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty5_Sig_low{$SD_2}){$Ty5_results_low{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Ty5_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty5_Sig_low{$SD_2}){$Ty5_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Ty5_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty5_Sig_low{$SD_2}){$Ty5_results_low{$gene}='sig_two_samples';}}
	}	
}

########### HIGH

my %Ty5_results_high;
foreach my $SD (keys %Ty5_Sig_high){
	my $gene = $Ty5_Sig_high{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Ty5_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty5_Sig_high{$SD_2}){$Ty5_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Ty5_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty5_Sig_high{$SD_2}){$Ty5_results_high{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Ty5_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty5_Sig_high{$SD_2}){$Ty5_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Ty5_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Ty5_Sig_high{$SD_2}){$Ty5_results_high{$gene}='sig_two_samples';}}
	}	
}



###################### Tel ##########################

########## LOW

my %Tel_results_low;
foreach my $SD (keys %Tel_Sig_low){
	my $gene = $Tel_Sig_low{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Tel_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Tel_Sig_low{$SD_2}){$Tel_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Tel_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Tel_Sig_low{$SD_2}){$Tel_results_low{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Tel_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Tel_Sig_low{$SD_2}){$Tel_results_low{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Tel_results_low{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Tel_Sig_low{$SD_2}){$Tel_results_low{$gene}='sig_two_samples';}}
	}	
}

########### HIGH

my %Tel_results_high;
foreach my $SD (keys %Tel_Sig_high){
	my $gene = $Tel_Sig_high{$SD};
	my $SD_2 = '';
	if ($SD =~ /b2/) { 
		$SD_2 = substr ($SD, 0, -1); 
		if (!defined $samples{$SD_2}){$Tel_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Tel_Sig_high{$SD_2}){$Tel_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b$/) {
		$SD_2 = $SD.'2'; 
		if (!defined $samples{$SD_2}){$Tel_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Tel_Sig_high{$SD_2}){$Tel_results_high{$gene}='sig_two_samples';}}
	}	
	if ($SD =~ /b3/) {
		$SD_2 = substr ($SD, 0, -1); $SD_2 = $SD.'4'; 
		if (!defined $samples{$SD_2}){$Tel_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Tel_Sig_high{$SD_2}){$Tel_results_high{$gene}='sig_two_samples';}}
	}
	elsif ($SD =~ /b4/) {
		$SD_2 = substr ($SD, 0, -1);$SD_2 = $SD.'3'; 
		if (!defined $samples{$SD_2}){$Tel_results_high{$gene}='sig_one_sample_only';}
		elsif (defined $samples{$SD_2}){ if (defined $Tel_Sig_high{$SD_2}){$Tel_results_high{$gene}='sig_two_samples';}}
	}	
}





#######################################
##    								 ##
##			PRINT OUTPUT			 ##
##									 ##
#######################################


###################### rDNA ##########################
#LOWER rDNA - print to file
open(my $fh, '>', 'rDNA_-.txt');
for my $out (sort keys %rDNA_results_low){
	print $fh "$out\t$rDNA_results_low{$out}\n";
}
close $fh;

#LOWER rDNA - print to file
open(my $fh2, '>', 'rDNA_+.txt');
for my $out (sort keys %rDNA_results_high){
	print $fh2 "$out\t$rDNA_results_high{$out}\n";
}
close $fh2;

###################### CUP1 ##########################
#LOWER Cup1 - print to file
open(my $fh3, '>', 'CUP1_-.txt');
for my $out (sort keys %CUP1_results_low){
	print $fh3 "$out\t$CUP1_results_low{$out}\n";
}
close $fh3;

#LOWER Cup1 - print to file
open(my $fh4, '>', 'CUP1_+.txt');
for my $out (sort keys %CUP1_results_high){
	print $fh4 "$out\t$CUP1_results_high{$out}\n";
}
close $fh4;

###################### Mito ##########################
#LOWER Mito - print to file
open(my $fh31, '>', 'M_-.txt');
for my $out (sort keys %Mito_results_low){
	print $fh31 "$out\t$Mito_results_low{$out}\n";
}
close $fh31;

#LOWER Mito - print to file
open(my $fh41, '>', 'M_+.txt');
for my $out (sort keys %Mito_results_high){
	print $fh41 "$out\t$Mito_results_high{$out}\n";
}
close $fh41;



###################### Ty1 ##########################
#LOWER Ty1 - print to file
open(my $fh5, '>', 'T1_-.txt');
for my $out (sort keys %Ty1_results_low){
	print $fh5 "$out\t$Ty1_results_low{$out}\n";
}
close $fh5;

#LOWER Ty1 - print to file
open(my $fh6, '>', 'T1_+.txt');
for my $out (sort keys %Ty1_results_high){
	print $fh6 "$out\t$Ty1_results_high{$out}\n";
}
close $fh6;

###################### Ty2 ##########################
#LOWER Ty2 - print to file
open(my $fh7, '>', 'T2_-.txt');
for my $out (sort keys %Ty2_results_low){
	print $fh7 "$out\t$Ty2_results_low{$out}\n";
}
close $fh7;

#LOWER Ty2 - print to file
open(my $fh8, '>', 'T2_+.txt');
for my $out (sort keys %Ty2_results_high){
	print $fh8 "$out\t$Ty2_results_high{$out}\n";
}
close $fh8;

###################### Ty3 ##########################
#LOWER Ty3 - print to file
open(my $fh9, '>', 'T3_-.txt');
for my $out (sort keys %Ty3_results_low){
	print $fh9 "$out\t$Ty3_results_low{$out}\n";
}
close $fh9;

#LOWER Ty3 - print to file
open(my $fh10, '>', 'T3_+.txt');
for my $out (sort keys %Ty3_results_high){
	print $fh10 "$out\t$Ty3_results_high{$out}\n";
}
close $fh10;


###################### Ty4 ##########################
#LOWER Ty4 - print to file
open(my $fh11, '>', 'T4_-.txt');
for my $out (sort keys %Ty4_results_low){
	print $fh11 "$out\t$Ty4_results_low{$out}\n";
}
close $fh11;

#LOWER Ty4 - print to file
open(my $fh12, '>', 'T4_+.txt');
for my $out (sort keys %Ty4_results_high){
	print $fh12 "$out\t$Ty4_results_high{$out}\n";
}
close $fh12;


###################### Ty5 ##########################
#LOWER Ty5 - print to file
open(my $fh13, '>', 'T5_-.txt');
for my $out (sort keys %Ty5_results_low){
	print $fh13 "$out\t$Ty5_results_low{$out}\n";
}
close $fh13;

#LOWER Ty5 - print to file
open(my $fh14, '>', 'T5_+.txt');
for my $out (sort keys %Ty5_results_high){
	print $fh14 "$out\t$Ty5_results_high{$out}\n";
}
close $fh14;

###################### Tel ##########################
#LOWER Tel - print to file
open(my $fh15, '>', 'TEL_-.txt');
for my $out (sort keys %Tel_results_low){
	print $fh15 "$out\t$Tel_results_low{$out}\n";
}
close $fh15;

#LOWER Tel - print to file
open(my $fh16, '>', 'TEL_+.txt');
for my $out (sort keys %Tel_results_high){
	print $fh16 "$out\t$Tel_results_high{$out}\n";
}
close $fh16;






#print "$SD and $SD_2\n";

##################
##				##
## 	USAGE POD	##
##				##
##################


__END__


=head1 SYNOPSIS

repeat_stats_sig.pl [options] -i <filename>
 
 Options:
   -help	brief help message
   -i		input (repeat file)
   
   

=head1 OPTIONS

=over 4

=item B<-help>

Prints a brief help message and exits.

=item B<-i>

Accepts a path to file with results for DNA repeat estimates.

=back


=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something useful with the contents thereof.

=cut

