#!/usr/bin/env perl
 
# Author:       Fabio Puddu
# Maintainer:   Fabio Puddu
# Created:      Jun 2017
# Description: This script runs through a file containing absolute paths to all the bam files to be processed and writes a script that can be run in IGV to obtain snapshots of the  area around a gene knockout
#	       to do so it reads the coordinates of each gene from a database file

use strict;
use warnings;
use Cwd 'abs_path';
#just get the directory this script is run from as a reference
my $script_location = abs_path($0);
my @path = split ('/', $script_location);
pop @path; pop @path;
my $dir = join ('/',@path);

my $outdir="snapshots"; #define the name of the output directory for the IGV snapshots
my $window_span=5000; # define the window size in bp before/after start/stop

my $cwd=`pwd`;
chomp($cwd);
$outdir=$cwd."/".$outdir;
#define where to find the files to be analysed
my $bamlist="all_bams.txt";
my $genepositions=$dir."/defaults/all_yeast_genes.txt";



open (my $position_file, '<', $genepositions) or die "cannot open gene position file";
my %positions;
while (my $line = <$position_file> ){
	next if $.==1;
	chomp($line);
	(my $sysname, my $chr, my $start, my $end, my $commname) = split("\t",$line);	
	my $window_start=$start-$window_span;
	my $window_end=$end+$window_span;
	#print "$commname\t$sysname\n";
	if (defined $commname && $commname ne ""){
	$positions{$commname}="chr".$chr.":".$window_start."-"."$window_end";
	}
	else{
	$positions{$sysname}="chr".$chr.":".$window_start."-"."$window_end";
	}
}
close ($position_file);
$positions{'WT-1'}='chrI:1-1000';
$positions{'WT-2'}='chrI:1-1000';
$positions{'WT-3'}='chrI:1-1000';
$positions{'WT-4'}='chrI:1-1000';

open (my $bamlist_file, '<', $bamlist) or die "cannot open input file";

while (my $file =<$bamlist_file>){
	my @pathlist=split("/",$file);
	my $SD=$pathlist[-1];
	chomp($SD);
	$SD =~ s/.bam//g;
	my $gene_name=$pathlist[-3];
	$gene_name =~ s/Del.*_//g;
	print "new\n";
	print "genome sacCer3\n";
	print "snapshotDirectory ",$outdir,"\n";
	print "load  ".$cwd."/".$file;
	print "goto ".$positions{$gene_name}."\n";
	print "viewaspairs\n";
	print "snapshot ".$SD.".png\n";
}
close ($bamlist_file);
