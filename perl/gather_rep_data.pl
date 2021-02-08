#!/usr/bin/env perl

# Author:		Fabio Puddu
# Maintainer:	Fabio Puddu
# Created:	Sep 2016
# Description:	Once the analysis has been run, run this program from the main analysis folder to gather all the repetitive DNA data in one fole called and redirect the output to a file

use strict;
use warnings;
use Data::Dumper;
my @FILES=`ls */ploidy/*_plstats.txt`;
my @SAMPLES=`cat "name_conversion.tsv" | cut -f1,2,6,7 | sort | uniq`;
printf "#bcID\trDNA\tCUP1\tmito\t2micron\tTy1\tTy2\tTy3\tTy4\tTy5\tGWM\tMatscore\tTelo\tERSno\tDeletion\tchr1\tchr2\tchr3\tchr4\tchr5\tchr6\tchr7\tchr8\tchr9\tchr10\tchr11\tchr12\tchr13\tchr14\tchr15\tchr16\tAneup\tAneupChrom\tCR\n";
my $rep_group;
my $telo;
my @CHROMOSOMES=('chr01', 'chr02', 'chr03', 'chr04', 'chr05', 'chr06', 'chr07', 'chr08', 'chr09', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16');
for my $line (@SAMPLES){
	chomp($line);
	(my $sample_code, my $delname, my $ERS,  my $pl_exp)=split("\t", $line);
	if (-e "$delname/repDNA/$sample_code.txt"){
		$rep_group=`cat "$delname/repDNA/$sample_code.txt" | grep -v Sample`;
		chomp($rep_group);
        	$telo=`cat "$delname/repDNA/$sample_code.tel"`;
		chomp($telo);
	}
	else{
		$rep_group=undef;
		$telo=undef;
	}
  my $ploidy_ref=&get_ploidy($sample_code,$delname, $pl_exp);
	my $gcr=get_gcr($sample_code,$delname,$ploidy_ref);
	if (defined $ploidy_ref and defined $rep_group and defined $telo and defined $gcr){
		my $ploidy=join("\t", @$ploidy_ref);
		print "$rep_group\t$telo\t$ERS\t$delname\t$ploidy\t$gcr\n";
	}
}


sub get_gcr{
	my $fn=shift;
	my $dname=shift;
	my $pref=shift;
	my %ploidy_by_chr;
	#read the chromosome ploidy of that sample in a hash
	for (my $i=0; $i<16; $i++){
		$ploidy_by_chr{$CHROMOSOMES[$i]}=${$pref}[$i];
	}
	if ( -e "$dname/ploidy/".$fn."_breakpoints.txt" && !-z "$dname/ploidy/".$fn."_breakpoints.txt"){
		open (my $fh, '<', "$dname/ploidy/$fn"."_breakpoints.txt");
		my @BKP=<$fh>;
		chomp @BKP;
		close $fh;
		my $bkp_count=0;
		#go through the breakpoint file
		foreach my $line (@BKP){
			next if $line =~ /^$/;
			(my $chr, my $start, my $end, my $pl) = split( "\t",$line );
			my $bkp_size=$end-$start;
			my $ploidy_deviation=abs($pl-$ploidy_by_chr{$chr});
			#only count breakpoints if the size is larger than 2kb and if the deviation from the ploidy of that particular chromosome is larger than 0.5
			$bkp_count++ if $bkp_size>20000 and $ploidy_deviation>0.5
		}
		return $bkp_count;
	}
	else{
	return 0;
	}
}

sub get_ploidy(){
	my $fn=$_[0];
	my $dname=$_[1];
	my $pl_exp=$_[2];
	my $delta;
	$delta=0;
	my $cc_aff=0;
	my $filename="$dname/ploidy/$fn"."_plstats.txt";
	if (-e $filename){
		open (my $fh, '<', $filename);
		chomp(my @CCN = <$fh>);
		close ($fh);
		my @output;
		#Go through each chromosome
        	for my $chr (@CHROMOSOMES){
                	my $matching_line = (grep { /$chr/ } @CCN)[0];
			my $ccn=(split "\t", $matching_line)[1];
			push @output, $ccn;
			#increase total aneuploidy and number of chromosomes affected;
			#if ($ccn <=1 or $ccn >=3){
			 $delta += abs($ccn-$pl_exp);
			$cc_aff++ if $ccn != $pl_exp;
			#}
        	}
        	push @output, $delta;
        	push @output, $cc_aff;
		return \@output;
	}
	else{
	return undef;
	}
}
