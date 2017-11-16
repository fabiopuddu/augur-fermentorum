#!/usr/bin/env perl

# Author:       mh23
# Maintainer:   mh23
# Created: 		Sep 2016
# Name:			rDNA_cnv_estimate.pl


# test bam: /mnt/home3/jackson/fp305/data/Retta_di_taratura/Del1_Fob1/BAM/SC_MFY5971838/SC_MFY5971838.bam

#Another test bam: /mnt/home3/jackson/fp305/data/Fabio_Mrc1_nico/Del1_Mrc1/BAM/SC_MFY6097876/SC_MFY6097876.bam
# the ty bam is here: /mnt/home3/jackson/fp305/data/Fabio_Mrc1_nico/Del1_Mrc1/TR_BAMS/SC_MFY6097876.bam.Ty.bam 


use autodie;
use utf8;
use Carp;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Pod::Usage;
use Cwd 'abs_path';
use POSIX qw(ceil);
use IPC::Open3;
use IO::File;


#######################################
#######################################
##    								 ##
##				USAGE				 ##
##									 ##
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

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions('help|?' => \$help, 
		   'i|input=s' => \$input,
		   'p|ploidy=s' => \$ploidy,
		   'r|reference=s' => \$ref_genome,
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
##    								 ##
##	 		Main programme			 ##
##									 ##
#######################################
#######################################

 #For regions of the genome with two copies in the genome multiply by two (eg rDNA, CUP1)
 
my $gen_cov_ref = genome_cov($input);
my $genome_wide_median_coverage = median(@$gen_cov_ref);

########################################
###	MATING TYPE AND PLOIDY ESTIMATION
#########################################

# check that the TY bam is there 

#split the path 
my @b = split ('/', $input);
#get rid of the last three elements
my $bam = pop @b; #pop @b;
my $sample_name = substr($bam, 0, -4);
#join and add the others
my $bam_dir = join ('/',@b);
my $ty_bam = $bam_dir.'/../TR_BAMS/'.$sample_name.'.Ty.bam';
my $mat_bam = $bam_dir.'/../TR_BAMS/'.$sample_name.'.Ty.bam';

my $twom_bam = $bam_dir.'/../TWOMICRON_BAMS/'.$sample_name.'.2m.bam';


my $MATa = 'MATa_HMR:1400-2000';
my $MATalpha=  'MATalpha_HML:1700-2700';
#print "$mat_bam,$MATa,$mat_ref\n";
my $MATa_estimate = repeat_estimate($mat_bam,$MATa,$mat_ref);
#exit;
my $MATalpha_estimate = repeat_estimate($mat_bam,$MATalpha,$mat_ref);
my $sex_estimate=log2($MATa_estimate/$MATalpha_estimate);
my $sex;
if    ($sex_estimate <= -0.6){$sex="alpha"}
elsif ($sex_estimate >=  0.6){$sex="a"}
elsif ($sex_estimate > -0.6 and $sex_estimate < 0.6){$sex="a/alpha"}
else { die "Cannot safely determine mating type and ploidy" }


my @s = split ('/', $input);
my $sa = $s[-1];
my @samp = split (/\./, $sa);

print "Sample\trDNA\tCUP1\tMito\t2-micron\tTy1\tTy2\tTy3\tTy4\tTy5\tGenome_wide_median\n";
print "$samp[0]\t";
my $rDNA_loc = 'XII:452000-459000';
my $rDNA_estimate = repeat_estimate($input,$rDNA_loc,$ref_genome)*2;
#$rDNA_estimate = sprintf("%.3f", $rDNA_estimate);
print "$rDNA_estimate\t";

my $cup1_loc = 'VIII:212986-213525';
my $cup1_estimate = repeat_estimate($input,$cup1_loc,$ref_genome)*2;
#$cup1_estimate = sprintf("%.3f", $cup1_estimate);
print "$cup1_estimate\t";

my $mito_loc = 'Mito:14000-20000';
my $mito_estimate = repeat_estimate($input,$mito_loc,$ref_genome)*$ploidy;
#$mito_estimate = sprintf("%.3f", $mito_estimate);
print "$mito_estimate\t";

my $twom_loc = '2-micron:2000-4500';
my $twom_estimate = repeat_estimate($twom_bam,$twom_loc,$twom_ref)*$ploidy;
#$mito_estimate = sprintf("%.3f", $mito_estimate);
print "$twom_estimate\t";


#######################################
#######################################
##    								 ##
##	 		 Ty elements			 ##
##									 ##
#######################################
#######################################

my $ty1 = 0;
my $ty2 = 0;
my $ty3 = 0;
my $ty4 = 0;
my $ty5 = 0;

my $genome_cov_file = $bam.'.genomewide_median';

#Save Genome Wide Median in a File 
open(my $FH, '>', $genome_cov_file);
print $FH "$bam\t$genome_wide_median_coverage\n";
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

	$ty1 =  repeat_estimate($ty_bam, $ty1_reg, $ty_ref);
	$ty2 =  repeat_estimate($ty_bam, $ty2_reg, $ty_ref);
	$ty3 =  repeat_estimate($ty_bam, $ty3_reg, $ty_ref);
	$ty4 =  repeat_estimate($ty_bam, $ty4_reg, $ty_ref);
	$ty5 =  repeat_estimate($ty_bam, $ty5_reg, $ty_ref);
	
	#$ty1 = sprintf("%.3f", $ty1);
	#$ty2 = sprintf("%.3f", $ty2);
	#$ty3 = sprintf("%.3f", $ty3);
	#$ty4 = sprintf("%.3f", $ty4);
	#$ty5 = sprintf("%.3f", $ty5);

}

	print ("$ty1\t");
	print ("$ty2\t");
	print ("$ty3\t");
	print ("$ty4\t");
	print ("$ty5\t");
	print ("$genome_wide_median_coverage\t");
	printf ("%s (%.1f)\n",$sex,$sex_estimate);


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
		my @genome_cov;
		#choose either of the open commands below to exclude low quality reads (2nd option) or not
		#open($CH, "$bed_command -d -ibam $input -g 2>/dev/null |") || die "Failed: $!\n";
		$in = '';
		local *CATCHERR = IO::File->new_tmpfile;
		my $pid = open3($in, \*CATCHOUT, ">&CATCHERR", "samtools view -q10 -b $file | $bed_command -d -ibam stdin -g");
		while( <CATCHOUT> ) {
		  my $line = $_;
		  chomp $line;
		  my @s = split(/\t/,$line);
		  #if the output is as expected then $s[0] is the chr, [1] is the position, [2] is the coverage
		  push (@genome_cov, $s[2]);
		}
		waitpid($pid, 0);
		seek CATCHERR, 0, 0;
		while( <CATCHERR> ) {pod2usage("$0: BedTools does not seem to work.")}
    	return (\@genome_cov);
}

sub repeat_estimate {
	#This subroutine will be given a location string like: "XII:123000-789000"
	#and the variable containing the appropriate reference
	my $in_file = $_[0];
	my $loc = $_[1];
	my $ref_g = $_[2];
	#split this to get the numbers
	my @t = split (':', $loc);
	my @n = split ('-', $t[1]);
	my %Cov_hash;
	#Set the coverage to 0 at each position
	foreach my $i ($n[0]..$n[1]) { $Cov_hash{$i}=0;}

	#check that the samtools mpileup command works
	#Ignore the standard STDERR
	$in = '';
	$check = 0;
	my $err = '';
	local *CATCHERR = IO::File->new_tmpfile;
	$pid = open3($in, \*CATCHOUT, ">&CATCHERR", "samtools mpileup  -A -Q 0 -r $loc -f $ref_g $in_file");
	while( <CATCHOUT> ) {
	my $l = $_;
	chomp $l;
	my @s = split("\t",$l);
	#if the output is as expected then $s[0] is the chr, [1] is the position, [3] is the coverage
	$Cov_hash{$s[1]}=$s[3];
	#print "$s[1]\t$s[3]\n";
	}
	waitpid($pid, 0);
	seek CATCHERR, 0, 0;
	while( <CATCHERR> ) {if ($_ =~ /[0-9]+\s+samples\s+in\s+[0-9]+\s+input/i) {$check = 1;} else {$err = $_;}}
	pod2usage("$0: samtools mpileup did not run correctly.\nError: $_\n")  if ( $check == 0);
	
	#Get the median of all rDNA coverage values extracted
	my $med_cov = median(values(%Cov_hash));
	#print "$med_cov\n";	

	#Divide query coverage by genome wide median coverage and round
	my $estimate = $med_cov/$genome_wide_median_coverage;
	
	return $estimate;
	
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

