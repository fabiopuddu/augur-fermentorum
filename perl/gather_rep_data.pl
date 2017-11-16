#!/usr/bin/env perl

use strict;
use warnings;
my @FILES=`ls */ploidy_data/*_plstats.txt`;
my @SAMPLES=`cat "name conversion.tsv"`;
printf "#SDname\trDNA\tCUP1\tmitochondria\tTy1\tTy2\tTy3\tTy4\tTy5\tGWM\tTelomeres\tERSno\tDeletion\tchr01\tchr02\tchr03\tchr04\tchr05\tchr06\tchr07\tchr08\tchr09\tchr10\tchr11\tchr12\tchr13\tchr14\tchr15\tchr16\tAneupNumber\n";
my $rep_group;
my $telo;
for my $line (@SAMPLES){
	chomp($line);
	(my $sample_code, my $delname, my $plate, my $aka, my $fname, undef, my $ERS)=split("\t", $line);
	if (-e "$delname/repDNA/$fname.txt"){ 
		$rep_group=`cat $delname/repDNA/$fname.txt | grep -v Sample`;
		chomp($rep_group);
        	$telo=`cat $delname/repDNA/$fname.tel`;
		chomp($telo);
	}
	else{
		$rep_group=undef;
		$telo=undef;
	}
	my $ploidy_ref=&get_ploidy($fname,$delname);
	if (defined $ploidy_ref and defined $rep_group and defined $telo){
	my $ploidy=join("\t", @$ploidy_ref);
	print "$rep_group\t$telo\t$ERS\t$delname\t$ploidy\n";
	}
}

sub get_ploidy(){
	my $fn=$_[0];
	my $dname=$_[1];
	my @CHROMOSOMES=('chr01', 'chr02', 'chr03', 'chr04', 'chr05', 'chr06', 'chr07', 'chr08', 'chr09', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16');
	my $delta;
	$delta=0;
	my $filename="$dname/ploidy_data/$fn"."_plstats.txt";
	if (-e $filename){
		open (my $fh, '<', $filename);
		chomp(my @CCN = <$fh>);
		close ($fh);
		my @output;
        	for my $chr (@CHROMOSOMES){
                	my $matching_line = (grep { /$chr/ } @CCN)[0];
			my $ccn=(split "\t", $matching_line)[1];
		 	$delta += abs($ccn-2);
                	push @output, $ccn;
        	}
        	push @output, $delta;
		return \@output;
	}
	else{
	return undef;
	}
}
