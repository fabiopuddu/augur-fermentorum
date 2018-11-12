#!/usr/bin/env perl

# Author:       	Mareike Herzog
# Maintainer:   	Fabio Puddu
# Created: 	Sep 2016
# Description:	This script calculates copy number of different repetitive elements and mating type using coverage data


use autodie;
#use utf8;
#use Carp;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Pod::Usage;
use Cwd 'abs_path';
use POSIX qw(ceil);
use IPC::Open3;
use IO::File;
#use Data::Dumper;


#######################################
#######################################
##    			       ##
##	USAGE		       ##
##			       ##
#######################################
#######################################

#Get the input and print usage if help of wrong input

#DEFAULTS
my $input = '0';
my $help = 0;
my $ploidy = 0; #default ploidy
#Get the location of the reference genome
my $script_location = abs_path($0);
my @path = split ('/', $script_location);
pop @path; pop @path; #this is the equivalent of .. from where the script is
my $dir = join ('/',@path);
my $ref_genome = $dir.'/mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa'; #default reference genome
my $ty_ref = $dir.'/mpileup_defaults/repDNA_ref/repDNA.fa'; # default Ty reference genome
my $mat_ref = $dir.'/mpileup_defaults/repDNA_ref/repDNA.fa';
my $twom_ref = $dir.'/mpileup_defaults/repDNA_ref/repDNA.fa';
my $mt_only=0;
my $qPCR=0;
## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions('help|?' => \$help, 
		   'i|input=s' => \$input,
		   'p|ploidy=s' => \$ploidy,
		   'r|reference=s' => \$ref_genome,
		   'M|mitochondrial_data_only=s' => \$mt_only,	
		   'Q|qPCR_data_only=s' => \$qPCR,
		   't|ty=s' => \$ty_ref,) or pod2usage(2);

pod2usage(1) if $help;
## If no input argument were given, then allow STDIN to be used only
## if it's not connected to a terminal (otherwise print usage)
pod2usage("$0: No input given.")  if (($input eq 0) && (-t STDIN));
pod2usage("$0: No ploidy given.")  if ($ploidy eq 0);
#Check that input is a bam file (at least in name: checks that file ends in .bam)
pod2usage("$0: No .bam file given.")  if ($input !~ /\.bam$/i); 

#Check that reference files exist and have a size bigger than 0		
pod2usage("$0: File $ref_genome does not exist.")  unless ( -e $ref_genome);
pod2usage("$0: File $ref_genome is empty.")  if ( -z $ref_genome);
pod2usage("$0: File $ty_ref does not exist.")  unless ( -e $ty_ref);
pod2usage("$0: File $ty_ref is empty.")  if ( -z $ty_ref);
pod2usage("$0: File $twom_ref does not exist.")  unless ( -e $twom_ref);
pod2usage("$0: File $twom_ref is empty.")  if ( -z $twom_ref);

#Check that input exists and has a size bigger than 0		
pod2usage("$0: File $input does not exist.")  unless ( -e $input);
pod2usage("$0: File $input is empty.")  if ( -z $input);

#Check that there is a .bai index file
my $index_file = $input.'.bai';
pod2usage("$0: File $input is not indexed.")  unless ( -e $index_file);
pod2usage("$0: File $input is not indexed properly.")  if ( -z $index_file);

#Check that there is a chromosome XII in the bam file
#check that running 'samtools view' does not create an error
my $check = 0;
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
#genomeCoverageBed
my @bed_out =  `genomeCoverageBed 2>&1 1>/dev/null`;
if ($bed_out[0] !~ /.+/) {$bed_command = 'genomeCoverageBed';}
#bedtools genomecov
@bed_out =  `bedtools genomecov 2>&1 1>/dev/null`;
if ($bed_out[0] !~ /.+/) {$bed_command = 'bedtools genomecov';}
pod2usage("$0: BedTools does not seem to work.") if ($bed_command !~ /.+/);


#######################################
#######################################
##    			       ##
##	 Main programme   	       ##
##			       ##
#######################################
#######################################

 #For regions of the genome with two copies in the genome multiply by two (eg rDNA, CUP1)
 
# check that the TY bam is there 

#split the path 
my @b = split ('/', $input);
#get rid of the last three elements
my $bam = pop @b; #pop @b;
my $sample_name = substr($bam, 0, -4);
#join and add the others
my $bam_dir = join ('/',@b);
my $rep_bam = $bam_dir.'/../TR_BAMS/'.$sample_name.'.Ty.bam';
my @s = split ('/', $input);
my $sa = $s[-1];
my @samp = split (/\./, $sa);

my $ty_bam = $bam_dir.'/../TR_BAMS/'.$sample_name.'.Ty.bam';
my $mat_bam = $bam_dir.'/../TR_BAMS/'.$sample_name.'.Ty.bam';
my $twom_bam = $bam_dir.'/../TR_BAMS/'.$sample_name.'.Ty.bam';

my $gen_cov_ref = genome_cov($input);
my @g_cov;
foreach (keys %{$gen_cov_ref}){
	push @g_cov, values %{%{$gen_cov_ref}{$_}}
}
my $GWM = median(@g_cov);

my $rep_gen_cov_ref = genome_cov($rep_bam);


if ( ! $qPCR and ! $mt_only){  
	print "Sample\trDNA\tCUP1\tMito\t2-micron\tTy1\tTy2\tTy3\tTy4\tTy5\tGenome_wide_median\tMatType\n";
	print "$samp[0]\t";
	#####    Ribosomal DNA    #######
	#   Average number of copies per haploid genome  ####
	my $location = 'XII:452000-459000';
	my $rDNA_estimate = (local_coverage($gen_cov_ref, $location)/$GWM)*2;
	#$rDNA_estimate = sprintf("%.3f", $rDNA_estimate);
	print "$rDNA_estimate\t";
	
	#####    CUP1     #######
	#   Average number of copies per haploid genome  ####
	$location = 'VIII:212986-213525';
	my $cup1_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM)*2;
	#$cup1_estimate = sprintf("%.3f", $cup1_estimate);
	print "$cup1_estimate\t";
	
	#####    mitochondrial DNA     #######
	#   Average number of copies per cell  ####
	$location = 'Mito:14000-20000'; #COX1
	my $mito_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM)*$ploidy;
	#$mito_estimate = sprintf("%.3f", $mito_estimate);
	print "$mito_estimate\t";
	
	#####    2 micron plasmid     #######
	#   Average number of copies per cell  ####
	$location = '2-micron:2000-4500';
	my $twom_estimate =  (local_coverage($rep_gen_cov_ref, $location)/$GWM)*$ploidy;
	#$mito_estimate = sprintf("%.3f", $mito_estimate);
	print "$twom_estimate\t";

}

if ($mt_only){
	print "Sample\tCOX1\tCOX3\tC1C3 ratio\n";
	print "$samp[0]\t";
	#####    mitochondrial DNA     #######
	#   Average number of copies per cell estimated via COX1 ####
	my $location = 'Mito:14000-20000'; #COX1
	my $mito_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM)*$ploidy;
	$mito_estimate = sprintf("%.3f", $mito_estimate);
	print "$mito_estimate\t";
	
	#####    mitochondrial DNA     #######
	#   Average number of copies per cell estimated via COX3 ####
	$location = 'Mito:79213-80022'; #COX3
	my $mito2_estimate =  (local_coverage($gen_cov_ref, $location)/$GWM)*$ploidy;
	$mito2_estimate = sprintf("%.3f", $mito2_estimate);
	print "$mito2_estimate\t";
	my $ratio = $mito_estimate/$mito2_estimate;
	$ratio = sprintf("%.3f", $ratio);
	print "$ratio\t";	
	print "\n";
}

if ($qPCR){
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


#######################################
#######################################
##      			       ##
##	Ty elements	       ##
##			       ##
#######################################
#######################################
if ( ! $qPCR and ! $mt_only){  


	my $ty1 = 0;
	my $ty2 = 0;
	my $ty3 = 0;
	my $ty4 = 0;
	my $ty5 = 0;

	my $genome_cov_file = $bam.'.genomewide_median';

	#Save Genome Wide Median in a File 
	open(my $FH, '>', $genome_cov_file);
	print $FH "$bam\t$GWM\n";
	close $FH;

	#print "TY bam: $ty_bam\n";
	#If the Ty bam file does not exist set them all to -1
	if (! -e $ty_bam || -z $ty_bam){
		$ty1 = -1;
		$ty2 = -1;
		$ty3 = -1;
		$ty4 = -1;
		$ty5 = -1;	
	}
	else {
		#Check that input exists and has a size bigger than 0		
		pod2usage("$0: File $ty_bam is empty.")  if ( -z $ty_bam);
		#Check that there is a .bai index file
		$index_file = $ty_bam.'.bai';
		pod2usage("$0: File $ty_bam is not indexed.")  unless ( -e $index_file);
		pod2usage("$0: File $ty_bam is not indexed properly.")  if ( -z $index_file);
	
		#get locations of the ty elements
		my $ty1_ctrl = 'YDRWTy1-5:100-2699';
		my $ty1_reg = 'YDRWTy1-5:4000-6999';
		my $ty2_ctrl = 'YLRWTy2-1:100-2699';
		my $ty2_reg = 'YLRWTy2-1:4000-6999';
		my $ty3_ctrl = 'YILWTy3-1:100-2699';
		my $ty3_reg = 'YILWTy3-1:4000-6999';
		my $ty4_ctrl = 'YHLWTy4-1:100-2699';
		my $ty4_reg = 'YHLWTy4-1:4000-6999';
		my $ty5_reg = 'YCLWTy5-1:2000-2999';
		my $ty5_ctrl = 'YCLWTy5-1:5000-5999';

		$ty1 =  local_coverage($rep_gen_cov_ref, $ty1_reg)/$GWM;
		$ty2 =  local_coverage($rep_gen_cov_ref, $ty2_reg)/$GWM;
		$ty3 =  local_coverage($rep_gen_cov_ref, $ty3_reg)/$GWM;
		$ty4 =  local_coverage($rep_gen_cov_ref, $ty4_reg)/$GWM;
		$ty5 =  local_coverage($rep_gen_cov_ref, $ty5_reg)/$GWM;
	
		#$ty1 = sprintf("%.3f", $ty1);
		#$ty2 = sprintf("%.3f", $ty2);
		#$ty3 = sprintf("%.3f", $ty3);
		#$ty4 = sprintf("%.3f", $ty4);
		#$ty5 = sprintf("%.3f", $ty5);

	}

		
		########################################
		###	MATING TYPE ESTIMATION
		#########################################

		my $MATa = 'MATa_HMR:1400-2000';
		my $MATalpha=  'MATalpha_HML:1700-2700';
		#print "$mat_bam,$MATa,$mat_ref\n";
		my $MATa_estimate = local_coverage($rep_gen_cov_ref,$MATa)/$GWM;
		#exit;
		my $MATalpha_estimate = local_coverage($rep_gen_cov_ref,$MATalpha)/$GWM;
		my $sex_estimate=log2($MATa_estimate/$MATalpha_estimate);
		my $sex;
		if    ($sex_estimate <= -0.35){$sex="alpha"}
		elsif ($sex_estimate >=  0.35){$sex="a"}
		elsif ($sex_estimate > -0.35 and $sex_estimate < 0.35){$sex="a/alpha"}
		else { die "Cannot safely determine mating type and ploidy" }
		
		print ("$ty1\t");
		print ("$ty2\t");
		print ("$ty3\t");
		print ("$ty4\t");
		print ("$ty5\t");
		print ("$GWM\t");
		printf ("%s (%.1f)\n",$sex,$sex_estimate);

}

#################
## Subroutines ##
#################

sub log2 {
        my $n = shift;
        return log($n)/log(2);
    }

sub median {
  return sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
}

sub mean {
 return sum(@_)/@_;
}

sub genome_cov{
        my $file= $_[0];
        my $CH;
	my %genome_cov;
		#choose either of the open commands below to exclude low quality reads (2nd option) or not
		#open($CH, "$bed_command -d -ibam $input -g 2>/dev/null |") || die "Failed: $!\n";
	$in = '';
	local *CATCHERR = IO::File->new_tmpfile;
		my $pid = open3($in, \*CATCHOUT, ">&CATCHERR", "samtools view -@ 8 -F 0x0400 -b $file | $bed_command -d -ibam stdin -g");
		while( <CATCHOUT> ) {
		  my $line = $_;
		  chomp $line;
		  my @s = split(/\t/,$line);
		  #if the output is as expected then $s[0] is the chr, [1] is the position, [2] is the coverage
		  $genome_cov{$s[0]}{$s[1]}=$s[2];
		}
		waitpid($pid, 0);
		seek CATCHERR, 0, 0;
		while( <CATCHERR> ) {pod2usage("$0: BedTools does not seem to work.")}
    	return \%genome_cov;
}

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
##				##
## 	USAGE POD	##
##				##
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

