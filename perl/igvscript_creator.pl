#!/usr/bin/env perl

# Author:       Fabio Puddu
# Maintainer:   Fabio Puddu
# Created:      Jun 2017
# Description: This script runs through a file containing absolute paths to all the bam files to be processed and writes a script that can be run in IGV to obtain snapshots of the  area around a gene knockout
#	       to do so it reads the coordinates of each gene from a database file

use strict;
use warnings;
use Cwd 'abs_path';
my $bamlist=shift;
my @bam_list = split ",", $bamlist;
#just get the directory this script is run from as a reference
my $script_location = abs_path($0);
my @path = split ('/', $script_location);
pop @path; pop @path;
my $dir = join ('/',@path);

my $outdir="snapshots"; #define the name of the output directory for the IGV snapshots
my $window_span=5000; # define the window size in bp before/after start/stop

my $cwd=`pwd`;
chomp($cwd);
my @pathlist=split "/", $cwd;
$outdir=$cwd."/".$outdir;
#define where to find the files to be analysed

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
	$positions{uc($commname)}="chr".$chr.":".$window_start."-"."$window_end";
	}
	else{
	$positions{uc($sysname)}="chr".$chr.":".$window_start."-"."$window_end";
	}
}
close ($position_file);
$positions{'WT-1'}='chrI:1-1000';
$positions{'WT-2'}='chrI:1-1000';
$positions{'WT-3'}='chrI:1-1000';
$positions{'WT-4'}='chrI:1-1000';

open(my $ofh, ">", "snapshots/script.igv");
foreach my $SD (@bam_list){
	chomp($SD);
	$SD =~ s/.cram//g;
	my $gene_name=$pathlist[-1];
	$gene_name =~ s/Del.*_//g;
  $gene_name =~ s/BAM\///g;
	print $ofh "new\n";
	print $ofh "genome sacCer3\n";
	print $ofh "snapshotDirectory ",$outdir,"\n";
	print $ofh "load  ".$cwd."/BAM/".$SD.".cram\n";
	print $ofh "goto ".$positions{uc($gene_name)}."\n";
	print $ofh "viewaspairs\n";
	print $ofh "snapshot ".$SD.".png\n";
}
close ($ofh);
