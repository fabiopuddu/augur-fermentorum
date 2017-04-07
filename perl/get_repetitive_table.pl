#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Term::ANSIColor;
use Getopt::Long;

#################################################
#                                               #
#  CALCULATING/DISPLAYING REP DNA TABLE	        #
#                                               #
#################################################

my ($input);
GetOptions
(
'i|input=s'	  => \$input,
);
( $input && -f $input ) or die qq[Usage: $0 \n
					 	-i <input file>\n
						];

#open input file
open (my $repDNA_file, '<', $input) or die "cannot open input file";
chomp(my @REPDNA = <$repDNA_file>);
close ($repDNA_file);
my %rDNA;
my %CUP1;
my %mito;
my %ty1;
my %ty2;
my %ty3;
my %ty4;
my %ty5; 
my %gwm;
my %telo;
my @SAMPLES;

foreach my $line (@REPDNA) {
	next if $line =~ "Genome_wide_median";
	my @LINE=split("\t", $line);
	my $ERSNO=$LINE[11];
	$rDNA{$ERSNO}= (defined $LINE[1] and $LINE[1] ne "-1" and $LINE[1] ne "") ? $LINE[1] : "N/A" ;
	$CUP1{$ERSNO}=(defined $LINE[2] and $LINE[2] ne "-1" and $LINE[2] ne "") ? $LINE[2] : "N/A" ;
	$mito{$ERSNO}=(defined $LINE[3] and $LINE[3] ne "-1" and $LINE[3] ne "") ? $LINE[3] : "N/A" ;
	$ty1{$ERSNO}=(defined $LINE[4] and $LINE[4] ne "-1" and $LINE[4] ne "") ? $LINE[4] : "N/A" ;
	$ty2{$ERSNO}=(defined $LINE[5] and $LINE[5] ne "-1" and $LINE[5] ne "") ? $LINE[5] : "N/A" ;;
	$ty3{$ERSNO}=(defined $LINE[6] and $LINE[6] ne "-1" and $LINE[6] ne "") ? $LINE[6] : "N/A" ;;
	$ty4{$ERSNO}=(defined $LINE[7] and $LINE[7] ne "-1" and $LINE[7] ne "") ? $LINE[7] : "N/A" ;;
	$ty5{$ERSNO}=(defined $LINE[8] and $LINE[8] ne "-1" and $LINE[8] ne "") ? $LINE[8] : "N/A" ;;
	$gwm{$ERSNO}=(defined $LINE[9] and $LINE[9] ne "-1" and $LINE[9] ne "") ? $LINE[9]  : "N/A" ;;
	$telo{$ERSNO}=(defined $LINE[10] and $LINE[10] ne "-1" and $LINE[10] ne "") ? $LINE[10] : "N/A" ;;
	push @SAMPLES, $ERSNO;
}


open (my $nameconversion, '<', '../../name conversion.tsv') or die "cannot open name conversion file";
chomp(my @NC = <$nameconversion>);
close ($nameconversion);
my %nc;
for my $line (@NC){
		(my $barcodeID, my $delname, my $plate, my $aka, my $filename, undef, my $ERSNO)=split("\t", $line);
		$nc{$ERSNO}=$aka;
}

print "REPETITIVE DNA QUANTIFICATION";
print "\n=========================================================================================================================================================================\n";
print "ERS.NO\tSample name\t\trDNA\tCUP1\tTy1\tTy2\tTy3\tTy4\tTy5\tTelomeres\tMitochondrial DNA\n";
print "=========================================================================================================================================================================\n";
#run through the has sorted by keys name
foreach my $ERSNO (@SAMPLES){
	print "$ERSNO\t$nc{$ERSNO}\t\t$rDNA{$ERSNO}\t$CUP1{$ERSNO}\t$ty1{$ERSNO}\t$ty2{$ERSNO}\t$ty3{$ERSNO}\t$ty4{$ERSNO}\t$ty5{$ERSNO}\t$telo{$ERSNO}\t$mito{$ERSNO}\t";#print the key (current sample ERS number)
	
	print "\n"
}