#!/usr/bin/perl
use strict;
use warnings;
#This script will go through a list of genes and attempt to plot on their 
#respective chromosomes them with gnuplot

my $filename_prefix="TEL";
open (my $up, '<', $filename_prefix.'_+.txt');
open (my $down, '<', $filename_prefix.'_-.txt');
open (my $gene_database, '<', '/mnt/home1/jackson/fp305/sw/bin/PF/defaults/all_yeast_genes.txt');

chomp(my @UP = <$up>);
chomp(my @DOWN = <$down>);

my @UPDOWN = (@UP, @DOWN);

chomp(my @GENES = <$gene_database>);

close ($up);
close ($down);
close ($gene_database);
open(my $out, '>', 'tempdata.tsv');

foreach my $row (@UPDOWN){
	my $gene=(split "\t", $row)[0];
	$gene =~ s/Del.*_//;
	$gene =~ s/\"//g;
	#print "$gene\n";
	my @match= split "\t", (grep { $_ =~ $gene } @GENES)[0];	
	my $match_chr=$match[1];
	my $match_pos=($match[2]+$match[3])/2;
	my $cn;
	if ($match_chr eq 'I'){$cn='01'}
	elsif ($match_chr eq 'II'){ $cn='02'}
	elsif ($match_chr eq 'III'){ $cn='03'}
	elsif ($match_chr eq 'IV'){ $cn='04'}
	elsif ($match_chr eq 'V'){ $cn='05'}
	elsif ($match_chr eq 'VI'){ $cn='06'}
	elsif ($match_chr eq 'VII'){ $cn='07'}
	elsif ($match_chr eq 'VIII'){ $cn='08'}
	elsif ($match_chr eq 'IX'){ $cn='09'}
	elsif ($match_chr eq 'X'){ $cn='10'}
	elsif ($match_chr eq 'XI'){ $cn='11'}
	elsif ($match_chr eq 'XII'){ $cn='12'}
	elsif ($match_chr eq 'XIII'){ $cn='13'}
	elsif ($match_chr eq 'XIV'){ $cn='14'}
	elsif ($match_chr eq 'XV'){ $cn='15'}
	elsif ($match_chr eq 'XVI'){ $cn='16'}
	printf  $out "$match_pos\t$cn\t-1\n" if /$gene/i ~~ @UP;
	printf  $out "$match_pos\t-1\t$cn\n" if /$gene/i ~~ @DOWN;
}
close ($out);
system("gnuplot /mnt/home1/jackson/fp305/sw/bin/PF/gnuplot/results_genomeplot.gpl");
system("rm tempdata.tsv");