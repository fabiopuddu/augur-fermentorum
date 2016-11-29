#!/usr/bin/env perl

# Author:       mh23
# Maintainer:   mh23
# Created: 		Sep 2016
# Name:			Ty1-5_estimate.pl


#A test bam file: 
# /mnt/home3/jackson/fp305/data/Retta_di_taratura/Del1_Fob1/BAM/SC_MFY5971838/SC_MFY5971838.bam

# /mnt/home3/jackson/fp305/data/Fabio_Mec1_HU/Del1_Mec1/bams_for_mpileup

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
my $ploidy = 1; #default ploidy
#Get the location of the reference genome
my $script_location = abs_path($0);
my @path = split ('/', $script_location);
pop @path; pop @path; #this is the equivalent of .. from where the script is
my $dir = join ('/',@path);
my $ref_genome = $dir.'/mpileup_defaults/Ty_ref/Ty1-5.fa'; #default reference genome

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions('help|?' => \$help, 
		   'i|input=s' => \$input,
		   'p|ploidy=s' => \$ploidy,
		   'r|reference=s' => \$ref_genome,) or pod2usage(2);
pod2usage(1) if $help;
## If no input argument were given, then allow STDIN to be used only
## if it's not connected to a terminal (otherwise print usage)
pod2usage("$0: No input given.")  if (($input eq 0) && (-t STDIN));

#Check that input exists and has a size bigger than 0		
pod2usage("$0: File $input does not exist.")  unless ( -e $input);
pod2usage("$0: File $input is empty.")  if ( -z $input);

#The programs that are required are samtools and bwa
#Check that bwa doesn't throw an error
my @bwa_check = `bwa mem 2>&1 1>/dev/null`;
chomp $bwa_check[0];
pod2usage("$0: bwa mem is not installed properly.") if ($bwa_check[0] ne '');
#Check that samtools doesn't throw an error
my @sam_check = `samtools 2>&1 1>/dev/null`;
chomp $sam_check[0];
pod2usage("$0: samtools is not installed properly.") if ($sam_check[0] ne '');

#Check that there are bam files in the file and that they exist
my @files_err = `cat $input 2>&1 1>/dev/null`;
pod2usage("$0: the cat command does not work.") if (defined $files_err[0]);
#Get the bam files
my @b_files = `cat $input 2>/dev/null`; 








##################################

#########

#########  UNFINISHED!!!!!

#########

###################################
















##################
##				##
## 	USAGE POD	##
##				##
##################


__END__


=head1 SYNOPSIS

Ty1-5_estimate.pl [options] -i F<filename.bam>
 
 Options:
   -help	brief help message
   -i		input (text file with a list of bams)
   -p		ploidy (default: 1)
   -r		a Ty reference genome, FASTA 
   -f		path to fq sequencing files 
   
This program requires that samtools and bwa are installed and in the path.

=head1 OPTIONS

=over 4

=item B<-help>

Prints a brief help message and exits.

=item B<-i>

Accepts a path to file containing paths to S. cerevisiae whole genome sequencing bam files.

=item B<-p>

Gives the ploidy of the sequenced strains. (Default: haploid, p=1)

=back

=item B<-r>

a S. cerevisiae Ty Element reference genome, FASTA. (default: ..//mpileup_defaults/Ty_ref/Ty1-5.fa)

=item B<-f>

path to the fq files containing the sequencing reads. (default: same directory as bam file)

=back


=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something useful with the contents thereof.

=cut