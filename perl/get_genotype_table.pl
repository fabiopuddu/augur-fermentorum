#!/usr/bin/env perl

# Author:		Fabio Puddu
# Maintainer:	Fabio Puddu
# Created:	Apr 2017
# Description:

use strict;
use warnings;
use Data::Dumper;
use Term::ANSIColor;
use Getopt::Long;
use List::Util 1.33 'any';
use Sort::Naturally;
use Cwd;
#################################################
#                                               #
#  CALCULATING/DISPLAYING GENOTYPE TABLE        #
#                                               #
#################################################
##convert output into genotypes
my ($s);my ($input);my ($pubmed);my ($nc);
GetOptions
(
's|show_synonymous'   => \$s,
'i|input=s'	  => \$input,
'n|nameconversion=s'	  => \$nc,
'p|pubmed'	  => \$pubmed,
);
( $input && -f $input && $nc && -f $nc) or die qq[Usage: $0 \n
					 	-i <input file>\n
						-n <name conversion file>\n
					 	-s <show synonymous mutations (default: no)> \n
						];

my %mutation = (
    "rad5-535"  => "YLR032W:missense_variant:1603:535:G>R",
    "rad50S" => "YNL250W:missense_variant:242:81:K>I",
    "sae2-F276S"  => "YGL175C:missense_variant:827:276:F>S",
);
my %NC;
my $cwd=(split "/", $input)[0];
open (my $nc_file, '<', $nc) or die "cannot open input file";
while (my $line=<$nc_file>){
				chomp $line;
				my ($bcID, $Del, $plate, $aka, $fname, $ers, $pl)=split "\t", $line;
				next unless $Del eq $cwd;
				$NC{'bcID'}{$bcID}=$bcID;
				$NC{'Del'}{$bcID}=$Del;
				$NC{'plate'}{$bcID}=$plate;
				$NC{'aka'}{$bcID}=$aka;
				$NC{'fname'}{$bcID}=$fname;
				$NC{'ers'}{$bcID}=$ers;
				$NC{'pl'}{$bcID}=$pl;
}
close($nc_file);

my %mutations;
#open input file
open (my $mutation_file, '<', $input) or die "cannot open input file";
chomp(my @mutation_list = <$mutation_file>);
close ($mutation_file);
#get a list of all the ERS files; each one is a sample
my @files= values %{$NC{'bcID'}};
my %gene_counts;
#first we find which genes are most commonly mutated
foreach (@mutation_list){
	next if $_ =~ 'CHROM';
  my (undef, undef, undef, undef, $system_name, $common_name,undef,undef,$hom_count,undef, undef, $het_count, undef)=split("\t", $_);
	$gene_counts{$system_name}+=$hom_count;
	$gene_counts{$system_name}+=$het_count;
}
#then we sort the list of mutations based on that
my @sorted_mutation_list;
foreach my $gene (sort {$gene_counts{$b} <=> $gene_counts{$a}} keys %gene_counts ){
	push @sorted_mutation_list, grep {$_ =~ /\t$gene\t/} @mutation_list;
}
#print Dumper \@sorted_mutation_list;
#finally we go through the list and create the genotypes
foreach my $sample_name (@files) {				#run through all the samples
	my @genotype;
	my @genotype_het;
	foreach my $line (@sorted_mutation_list){			#run through every line of the input file
		next if $line !~ m/$sample_name/; 		#skip if the line does not contain the current sample name
		my (undef, undef, undef, $mut, $system_name, $common_name,undef,undef,undef,$ERS_hom, undef, undef, $ERS_het)=split("\t", $line);
		my $gtype;

		if ((defined $common_name) and (not $common_name eq '""') ){ # use either the common name if present or the systematic name
			$gtype=$common_name."-".$mut;
		}
		else{
			$gtype=$system_name."-".$mut
		}

		if ( $ERS_hom =~ $sample_name){
			push (@genotype, $gtype);			#push the genotype if mutation is hom
		}
		elsif($ERS_het  =~ $sample_name){
			push (@genotype_het, '('.$gtype.')') #push the genotype in brackets if mutation is het
		}
		else {
			die ("This should never happen\n")
		};
	}
	#put the reference to the current genotype list into a hash that has current sample name as key
	push @genotype, @genotype_het;
	$mutations{$sample_name}=\@genotype;
}
print "Differential Genotype (mutations in parentheses are heterozygous or subclonal)";
print "\n=========================================================================================================================================================\n";
#run through the has sorted by keys name
open(my $ofh,'>',"$cwd/analysis/predicted_genotypes.txt");
foreach my $gt (sort {ncmp($NC{'aka'}{$a},$NC{'aka'}{$b})} keys %mutations){
	print $ofh "$gt\t$NC{'aka'}{$gt}\t";
	if (any { /TOP1/ } @{$mutations{$gt}}) {
	print color("red"), "$gt\t$NC{'aka'}{$gt}\t", color("reset");#print the key (current sample ERS number)
	}
	else{
	print "$gt\t$NC{'aka'}{$gt}\t";#print the key (current sample ERS number)
	}
	for my $mut(@{$mutations{$gt}}){
		if (($mut =~ m/FS@/) or ($mut =~ m/£Δ/) or ($mut =~ m/£st_lo/)){
			$mut=~ s/£//g;
			print color("red"), "$mut\t", color("reset");
			print $ofh "$mut\t";
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
		print $ofh "$mut\t";
		}
		elsif ($mut =~ m/>/) {
		$mut=~ s/£//g;
		print color("blue"), "$mut\t", color("reset") if ($s);
		print $ofh "$mut\t";
		}
		elsif (($mut =~ m/£II/) or ($mut =~ m/£ID/)){
		$mut=~ s/£//g;
		print color("green"), "$mut\t", color("reset");
		print $ofh "$mut\t";
		}
		else{
		print "error!\t";
		print $ofh "error!\t";
		}
	}
	print "\n";
	print $ofh "\n";
}

close $ofh;
