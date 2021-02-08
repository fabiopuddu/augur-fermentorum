#!/usr/bin/env perl

# Author:       	Mareike Herzog
# Maintainer:   	Fabio Puddu
# Created: 	Sep 2016
# Description:	This script calculates copy number of different repetitive elements and mating type using coverage data
use autodie;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Pod::Usage;
use Cwd 'abs_path';
use POSIX qw(ceil);
use IPC::Open3;
use IO::File;
use Parallel::Subs;
use Data::Dumper;
#######################################
#######################################
##                                   ##
##	USAGE                        ##
##                                   ##
#######################################
#######################################
#Get the input and print usage if help of wrong input
#DEFAULTS
my $input = '0';
my $help = 0;
my $ploidy = 0; 					#default ploidy
#Get the location of the reference genome
my $script_location = abs_path($0);
my @path = split ('/', $script_location);
pop @path; pop @path; 					#this is the equivalent of .. from where the script is
my $dir = join ('/',@path);
my $ref_genome;
my $mt_only=0;
my $qPCR=0;
## Parse options and print usage if there is a syntax error, or if usage was explicitly requested.
GetOptions('help|?' => \$help,
		   'i|input=s' => \$input,
		   'p|ploidy=s' => \$ploidy,
		   'r|reference=s' => \$ref_genome,
		   'M|mitochondrial_data_only=s' => \$mt_only,
		   'Q|qPCR_data_only=s' => \$qPCR) or pod2usage(2);
pod2usage(1) if $help;
pod2usage("$0: No input given.")  if (($input eq 0) && (-t STDIN));		## If no input argument were given, then allow STDIN to be used only if it's not connected to a terminal (otherwise print usage)
pod2usage("$0: No ploidy given.")  if ($ploidy eq 0);
pod2usage("$0: No .bam file given.")  if ($input !~ /\.cram$|\.bam$/i); 	#Check that input is a bam file (at least in name: checks that file ends in .bam)

pod2usage("$0: File $ref_genome does not exist.")  unless ( -e $ref_genome);	#Check that reference files exist and have a size bigger than 0
pod2usage("$0: File $ref_genome is empty.")  if ( -z $ref_genome);

pod2usage("$0: File $input does not exist.")  unless ( -e $input);		#Check that input exists and has a size bigger than 0
pod2usage("$0: File $input is empty.")  if ( -z $input);

my $index_file1 = $input.'.bai';
my $index_file2 = $input.'.crai';							#Check that there is a .bai index file
pod2usage("$0: File $input is not indexed.")  unless ( -e $index_file1 or -e $index_file2);
pod2usage("$0: File $input is not indexed properly.")  if ( -z $index_file1 or -z $index_file2);

my $check = 0;									#Check that there is a chromosome XII in the bam file check that running 'samtools view' does not create an error
my $in = '';
local *CATCHERR = IO::File->new_tmpfile;
my $pid = open3($in, \*CATCHOUT, ">&CATCHERR", "samtools view -H $input");
while( <CATCHOUT> ) {if ($_ =~ /SN:XII/){$check=1; last;}}
waitpid($pid, 0);
seek CATCHERR, 0, 0;
while( <CATCHERR> ) {pod2usage("$0: Samtools is not installed properly.");}
waitpid($pid, 0);
pod2usage("$0: File $input does not seem to be S. cerevisiae. It does not have a ChrXII.")  if ( $check == 0);

#check that either of the genome coverage commands works
my $bed_command = '';
my @bed_out =  `genomeCoverageBed 2>&1 1>/dev/null`;		#genomeCoverageBed
if ($bed_out[0] !~ /.+/) {$bed_command = 'genomeCoverageBed';}
@bed_out =  `bedtools genomecov 2>&1 1>/dev/null`;		#bedtools genomecov
if ($bed_out[0] !~ /.+/) {$bed_command = 'bedtools genomecov';}
pod2usage("$0: BedTools does not seem to work at start.") if ($bed_command !~ /.+/);

#######################################
#######################################
##    			             ##
##		 Main program	     ##
##			             ##
#######################################
#######################################
my @s = split ('/', $input);
my $sa = $s[-1];
my @samp = split (/\./, $sa);
my $gen_cov_ref = genome_cov($input);			#read in the data from the bam/cram file
my @g_cov;
my @chromosomes=qw[I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI];
foreach (@chromosomes){					#|
	push @g_cov, values %{%{$gen_cov_ref}{$_}}	#|calculate the genomewide median
}							#|
my $GWM = median(@g_cov);

if ( ! $qPCR and ! $mt_only){ 				#In the standard case where we want all data
	print "Sample\trDNA\tCUP1\tMito\t2-micron\tTy1\tTy2\tTy3\tTy4\tTy5\tGenome_wide_median\tMatType\n";	#header
	print "$samp[0]\t";				#print sample name
	my $p = Parallel::Subs->new();			#run all the copy number estimations in parallel
	$p->add(
    		sub {
			#####    Ribosomal DNA    #######
			#   Average number of copies per haploid genome  ####
			my $location = 'XII:452000-459000';
			my $rDNA_estimate = (local_coverage($gen_cov_ref, $location)/$GWM)*2;
			#$rDNA_estimate = sprintf("%.3f", $rDNA_estimate);
		});
	$p->add(
                sub {
			#####    CUP1     #######
			#   Average number of copies per haploid genome  ####
			my $location = 'VIII:212986-213525';
			my $cup1_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM)*2;
			#$cup1_estimate = sprintf("%.3f", $cup1_estimate);
	});
	$p->add(
                sub {
			#####    mitochondrial DNA     #######
			#   Average number of copies per cell  ####
			my $location = 'Mito:14000-20000'; #COX1
			my $mito_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM)*$ploidy;
			#$mito_estimate = sprintf("%.3f", $mito_estimate);
	});
	$p->add(
                sub {	#####    2 micron plasmid     #######
			#   Average number of copies per cell  ####
			my $location = '2-micron:2000-4500';
			my $twom_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM)*$ploidy;
			#$mito_estimate = sprintf("%.3f", $mito_estimate);
	});
	$p->add(
                sub {	#####    Ty1 transposon     #######
			#   Average number of copies per haploid genome  ####
			my $location = 'YDRWTy1-5:1000-3999';
			my $ty1 =  local_coverage($gen_cov_ref, $location)/$GWM;
	});
	$p->add(
                sub {	#####    Ty2 transposon     #######
			#   Average number of copies per haploid genome  ####
			my $location = 'YLRWTy2-1:1000-3999';
			my $ty2 =  local_coverage($gen_cov_ref, $location)/$GWM;
	});
	$p->add(
                sub {	#####    Ty3 transposon     #######
			#   Average number of copies per haploid genome  ####
			my $location = 'YILWTy3-1:1000-3999';
			my $ty3 =  local_coverage($gen_cov_ref, $location)/$GWM;
	});
	$p->add(
                sub {	#####    Ty4 transposon     #######
			#   Average number of copies per haploid genome  ####
			my $location = 'YHLWTy4-1:1000-3999';
			my $ty4 =  local_coverage($gen_cov_ref, $location)/$GWM;
	});
	$p->add(
                sub {	#####    Ty5 transposon     #######
			#   Average number of copies per haploid genome  ####
			my $location = 'YCLWTy5-1:0-999';
			my $ty5 =  local_coverage($gen_cov_ref, $location)/$GWM;
	});
	$p->add(
                sub {	#####    MATa/MATalpha   #######
			#  MAT_type estimation   ####
			my $MATa = 'MATa_HMR:1400-2000';
			my $MATalpha=  'MATalpha_HML:1700-2700';
			my $MATa_estimate = local_coverage($gen_cov_ref,$MATa)/$GWM;
			my $MATalpha_estimate = local_coverage($gen_cov_ref,$MATalpha)/$GWM;
			my $sex_estimate=log2($MATa_estimate/$MATalpha_estimate);
	});
	$p->wait_for_all();		#wait for all parallel estimates to finish
	my @res=$p->results();		#and gather the results
	my $rDNA_estimate=$res[0][0];	#|
	my $cup1_estimate=$res[0][1];	#|
	my $mito_estimate=$res[0][2];	#|
	my $twom_estimate=$res[0][3];	#|
	my $ty1=$res[0][4];		#|sort the results in appropriate variables
	my $ty2=$res[0][5];		#|
	my $ty3=$res[0][6];		#|
	my $ty4=$res[0][7];		#|
	my $ty5=$res[0][8];		#|
	my $sex_estimate=$res[0][9];	#|
	#and print out
	print "$rDNA_estimate\t$cup1_estimate\t$mito_estimate\t$twom_estimate\t$ty1\t$ty2\t$ty3\t$ty4\t$ty5\t$GWM\t$sex_estimate";
}
if ($mt_only){ 				#In the case were I want only estimates about mtDNA, with both COX1 and COX3
	print "Sample\tCOX1\tCOX3\tC1C3 ratio\n";
	print "$samp[0]\t";
	#####    mitochondrial DNA COX1     #######
	#   Average number of copies per cell estimated via COX1 ####
	my $location = 'Mito:14000-20000'; 		#COX1
	my $cox1_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM)*$ploidy;
	$cox1_estimate = sprintf("%.3f", $cox1_estimate);
	#####    mitochondrial DNA COX3     #######
	#   Average number of copies per cell estimated via COX3 ####
	$location = 'Mito:79213-80022'; #COX3
	my $cox3_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM)*$ploidy;
	$cox3_estimate = sprintf("%.3f", $cox3_estimate);

	my $ratio = $cox1_estimate/$cox3_estimate;
	$ratio = sprintf("%.3f", $ratio);
	print "$cox1_estimate\t$cox3_estimate\t$ratio\n";
}
if ($qPCR){   				#If I want only estimates of CUP1 and mtDNA using the same region of qPCR amplicons
	print "Sample\tCOX1_qPCR\tGAL1_qPCR\tCUP1_qPCR\tmtDNA_qPCR_ratio\tCUP1_qPCR_ratio\n";
	print "$samp[0]\t";

	####qPCR amplicon on COX1 (Oligos Cox1-qF2 + Cox1-qR2)
	my $location = 'Mito:25574-25686';
	my $mito_qpcr_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM)*$ploidy;
	$mito_qpcr_estimate = sprintf("%.3f", $mito_qpcr_estimate);
	print "$mito_qpcr_estimate\t";

	#####qPCR amplicon on GAL1 (GAL1(+1442)up + GAL1(+1442)low)
	$location = 'II:280382-280459';
	my $gal1_qpcr_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM);
	$gal1_qpcr_estimate = sprintf("%.3f", $gal1_qpcr_estimate);
	print "$gal1_qpcr_estimate\t";

	####qPCR amplicon on CUP1( Oligos CUP1-qPCR-F1 + CUP1-qPCR-R1/R2 )
	$location = 'VIII:212551-212669';
	my $cup1_qpcr_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM)*2;
	$cup1_qpcr_estimate = sprintf("%.3f", $cup1_qpcr_estimate);
	print "$cup1_qpcr_estimate\t";

	my $ratio = sprintf("%.3f", $mito_qpcr_estimate/$gal1_qpcr_estimate);
	print "$ratio\t";
	$ratio = sprintf("%.3f", $cup1_qpcr_estimate/$gal1_qpcr_estimate);
	print "$ratio\t";
	print "\n";
}

#################
## Subroutines ##
#################
sub log2 {
        my $n = shift;
        return log($n)/log(2);
    }
#################
sub median {
  return sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
}
#################
sub genome_cov{
        my $file= $_[0];
        my $CH;
	my %genome_cov;
	$in = '';
	local *CATCHERR = IO::File->new_tmpfile;
		my $pid = open3($in, \*CATCHOUT, ">&CATCHERR", "samtools view -@ 8 -F 0x0400 -b $file | $bed_command -d -ibam stdin -g");
		while( <CATCHOUT> ) {
		  my $line = $_;
		  chomp $line;
		  my @s = split(/\t/,$line);
		  $genome_cov{$s[0]}{$s[1]}=$s[2];			  #if the output is as expected then $s[0] is the chr, [1] is the position, [2] is the coverage
		}
		waitpid($pid, 0);
		seek CATCHERR, 0, 0;
		while( <CATCHERR> ) {pod2usage("$0: BedTools does not seem to work. File: $file")}
    	return \%genome_cov;
}
#################
sub local_coverage{
	my $coverage_hash=shift;
	my %coverage=%{$coverage_hash};
	my $position=shift;
	my ($CHR,$RANGE)=split ":", $position;
	my ($START, $END)=split "-", $RANGE;
	my @loc_cov;
	foreach my $k (sort {$a <=> $b} keys %{$coverage{$CHR}}){
		last if $k>$END;
		if ($k>=$START){
			push @loc_cov, $coverage{$CHR}{$k};
		}
	}
	my $out=median(@loc_cov);
	return $out;
}
##################
##		##
## USAGE POD	##
##		##
##################
__END__
=head1 SYNOPSIS

rDNA_cnv_estimate.pl [options] -i F<filename.bam>

 Options:
   -help	brief help message
   -i		input (bam file)
   -p		ploidy (default: 1)
   -r		a S. cerevisiae reference genome, FASTA
   -t		a Ty element custom reference genome, FASTA
   -M=1		only report mtDNA comparisons
   -Q=1		only report qPCR comparisons
This program requires that samtools and bedtools are installed and in the path.

=head1 OPTIONS

=over 4

=item B<-help>

Prints a brief help message and exits.

=item B<-i>

Accepts a path to file for a S. cerevisiae whole genome sequencing bam file.

=item B<-p>

Gives the ploidy of the sequenced strain. (Default: haploid, p=1)

=back

=item B<-r>

a S. cerevisiae reference genome, FASTA. (default: ../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa)

=item B<-t>

a S. cerevisiae ty element custom reference genome, FASTA. (default: ../mpileup_defaults/mpileup_defaults/Ty_ref/Ty1-5.fa)


=back


=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something useful with the contents thereof.

=cut
