#!/usr/bin/perl
use strict;
use warnings;
use IO::Zlib;
#This script will go through a vcf file and feed to gnuplot a list of positions to plot
my $unfiltered=shift;
my $u = new IO::Zlib;

if ($u->open($unfiltered, "rb")){
	while(my $line = <$u>){
		chomp $line;
		next if $line =~ /^#/;
		my ($chr, $pos, @rest) = split "\t", $line;
		my $c = get_as_number($chr);
	  if ($line =~ /INDEL/){
			 printf  "0\t0\t$c\t$pos\n";
		}
		else{
			printf  "$c\t$pos\t0\t0\n";
		}
	}
}
close($u);

sub get_as_number{
	my $cn;
	my $match_chr = shift;
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
	else {die "Error on chromosome name"}
	return $cn;
}
