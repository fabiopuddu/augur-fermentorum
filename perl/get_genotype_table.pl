#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Term::ANSIColor;
use Getopt::Long;
use List::Util 1.33 'any';

#################################################
#                                               #
#  CALCULATING/DISPLAYING GENOTYPE TABLE        #
#                                               #
#################################################

##convert output into genotypes
#create genotype tables with both common names and systematic names; add 'NONAME' whenever a non standard gene name is found

my ($s);
my ($input);
my ($pubmed);
GetOptions
(
's|show_synonymous'   => \$s,
'i|input=s'	  => \$input,
'p|pubmed'	  => \$pubmed,
);
( $input && -f $input ) or die qq[Usage: $0 \n
					 	-i <input file>\n
					 	-s <show synonymous mutations (default: no)> \n
						];

my %mutation = (
    "rad5-535"  => "YLR032W:missense_variant:1603:535:G>R",
    "rad50S" => "YNL250W:missense_variant:242:81:K>I",
    "sae2-F276S"  => "YGL175C:missense_variant:827:276:F>S",
);

my %mutations;

#open input file

open (my $mutation_file, '<', $input) or die "cannot open input file";
chomp(my @mutation_list = <$mutation_file>);
close ($mutation_file);
#get a list of all the ERS files; each one is a sample
my @files = <ERS*.vcf.gz>;
#run through all the samples
foreach my $file (@files) {
	my $sample_name = $file =~ s/.vcf.gz//r; #remove extensions from ERS filename
	my @genotype;
	#run through every line of the input file 
	foreach my $line (@mutation_list){
		next if $line !~ m/$sample_name/; #skip if the line does not contain the current sample name
		(undef, undef, undef, my $mut, my $system_name, my $common_name,undef,undef,undef,undef)=split("\t", $line);
		# push all the genotypes in the genotype list; use either the common name if present or the systematic name
		if ((defined $common_name) and (not $common_name eq '""') ){
			push (@genotype, $common_name."-".$mut);
		}
		else{
			push (@genotype, $system_name."-".$mut);
		}
	}
	#put the reference to the current genotype list into a hash that has current sample name as key
	$mutations{$sample_name}=\@genotype;
}
#print Dumper(\%mutations);

open (my $nameconversion, '<', '../../name conversion.tsv') or die "cannot open name conversion file";
chomp(my @NC = <$nameconversion>);
close ($nameconversion);
my %nc;
for my $line (@NC){
		(my $barcodeID, my $delname, my $plate, my $aka, my $filename, undef, my $ERSNO)=split("\t", $line);
		$nc{$ERSNO}=$aka;
}

print "Differential Genotype";
print "\n=========================================================================================================================================================\n";
#run through the has sorted by keys name
foreach my $gt (sort keys %mutations){
	if (any { /TOP1/ } @{$mutations{$gt}}) {
	print color("red"), "$gt\t$nc{$gt}\t", color("reset");#print the key (current sample ERS number)
	}
	else{
	print "$gt\t$nc{$gt}\t";#print the key (current sample ERS number)
	}
	for my $mut(@{$mutations{$gt}}){
		if (($mut =~ m/FS@/) or ($mut =~ m/£Δ/)){
		$mut=~ s/£//g;
		print color("red"), "$mut\t", color("reset");
		}
		elsif ($mut =~ m/£x/) {
		$mut=~ s/£//g;
		$mut=~ s/x//;
		my $PMCID_NUM='';
		if ($pubmed){
			(my $a, my $b) = split("-", $mut);
			chop $b;
			$PMCID_NUM=`scraper.py $a $b | wc -l`; #${gen::-1} returns the genotipe minus the last character
			chop $PMCID_NUM;
			$PMCID_NUM='('.$PMCID_NUM.')';
		}
		print color("yellow"), "$mut$PMCID_NUM\t", color("reset"), ;
		}
		elsif ($mut =~ m/>/) {
		$mut=~ s/£//g;
		print color("blue"), "$mut\t", color("reset") if ($s);
		}
		elsif (($mut =~ m/£II/) or ($mut =~ m/£ID/)){
		$mut=~ s/£//g;
		print color("green"), "$mut\t", color("reset");
		}
		else{
		print "error!\t"
		}
	}
	print "\n"
}
