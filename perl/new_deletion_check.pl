#!/usr/bin/env perl
# Author:		Mareike Herzog & Fabio Puddu
# Maintainer:	Fabio Puddu
# Created:	Jul 2019
# Description: Short script to detect presence of specific point mutations in bam files

use Carp;
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use List::MoreUtils qw(uniq);
use List::Util qw(sum);
use Data::Dumper;

use DBI;
# MySQL database configuration
my $host = 'spj-fp1.gurdon.private.cam.ac.uk';
my $db = 'LaboratoryDatabase';
my $dsn="DBI:mysql:$db:$host";
my $username = "cb-head";
my $password = 'charles';

my %attr = (PrintError=>0,RaiseError=>1 );
my $dbh = DBI->connect($dsn,$username,$password,\%attr);
# Test file (P301R): /mnt/scratch/jackson/fp305/from_CGP/Pol_BAM/SD0824b.bam
# Test command: perl /mnt/home1/jackson/fp305/PF/perl/new_deletion_check.pl -i /mnt/scratch/jackson/fp305/from_CGP/Pol_BAM/SD0824b.bam -d Del012_pol2-P301R -o yeast

#####################################
# Get command line options and path #
#####################################

my ($input,$del_string, $organism, $result);
GetOptions
(
'i|input=s'  => \$input,
'd|del_string=s' => \$del_string,
'o|organism=s' => \$organism,
);
( $del_string && $input && -f $input && $organism) or die qq[Usage: $0 \n
					 	-i <input .bam file> \n
						-d <Deletion string e.g. Del012_pol2-P301R,...> \n
						-o (organism: yeast, mouse, human")\n];
my $sample_name = (split "/",$input)[-1];
$sample_name =~ s/.cram|.bam//g;
my @path = split( '/' , abs_path($0));
pop(@path);
my $local_folder = join('/',@path);
my $dirscript = dirname(__FILE__); # Directory where the script is stored]

my $ref = $dirscript."/../mpileup_defaults/new_reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm_MASKED+REPDNA.toplevel.fa";

########################################
# Coordinates -> may pass a file later #
########################################
my %sysname_locations;

my %sysname;
my %cname;
open (my $fh, '<', $dirscript."/../defaults/all_yeast_genes.txt");
while (my $line=<$fh>){
	chomp $line;
	my @riga=split "\t", $line;
	$sysname_locations{uc($riga[0])}{'chr'}=$riga[1];
	$sysname_locations{uc($riga[0])}{'start'}=$riga[2];
	$sysname_locations{uc($riga[0])}{'end'}=$riga[3];

	$sysname{uc($riga[4])}=uc($riga[0]) if defined $riga[4];
	$cname{uc($riga[0])}= defined $riga[4] ? uc($riga[4]) : uc($riga[0]);
}
close ($fh);

my $sql = "SELECT * FROM alleles ";
my $stmt = $dbh->prepare($sql);
$stmt->execute();
my %complex_mutations;
while (my @row = $stmt->fetchrow_array) {
	$complex_mutations{$row[1]}=$row[2];
}
####################
# Get the mutation #
####################
my %result;
#split Del change
my $genotype = (split "_",$del_string)[1];
my @mutations;
my $logic;
if ($genotype =~ /\+/ ){
    @mutations = split /\+/, $genotype;
    $logic='and';
}
elsif ($genotype =~ /\^/ ){
    @mutations = split /\^/, $genotype;
    $logic='or';
}
else{
    $mutations[0] = $genotype;
}
foreach (@mutations){
	if ($_ !~ /WildtypeControl/){
		my @allele=split( '-', $_);
		if (scalar @allele == 2 ){
			if ($allele[1] eq "A" or $allele[1] eq "B" or $allele[1] eq "C" or ($allele[1] eq "6" and $allele[0] eq "YRF1")){
				#WE HAVE A YLR432C-A type of gene name	and we need to check for a deletion
				my $gene="$allele[0]-$allele[1]"; #stitch it back together
				$result{$gene}=check_deletion($gene, $input);
			}
			else{
				#WE NEED TO CHECK FOR A POINT MUTATION
				my $mutations=$allele[1];
				my $gene=$allele[0];
				$result{$gene."-".$mutations}=check_point_mut($input, $gene, $mutations);
			}
		}
		elsif (scalar @allele == 1){
			my $gene=$allele[0];
			$result{$gene}=check_deletion($gene, $input)#we need to check for a deletion
			}
		else{
			die ("What do you think you're doing!?!")
		}
	}
	else{
		$result{"WildtypeControl"}{"answer"}="Y";
	}
}
#print Dumper \%result;
my @r;
foreach (keys %result){
	push @r, $result{$_}{'answer'};
}
if ("FAILED" ~~ @r){
	print "$input\tFAILED\n";
	exit 1;
}
else{
	my $answer;
	if ($logic eq 'and' or not defined $logic){
		$answer = (scalar (uniq(map {uc($_)} @r)) == 1 and (uniq(map {uc($_)} @r))[0] eq "Y") ? "Y" : "N";
	}
	elsif ($logic eq 'or'){
		my @R=map {uc($_)} @r;
		$answer = ("Y" ~~ @R ) ? "Y" : "N";
	}
	print "$input\tGlobal>\tCheck>\tResult>\tPass?\tAnswer:$answer\n";
}

foreach (keys %result){
	if ($_ !~ /WildtypeControl/) {
		my $gene_name=""; my $gene_oi="";
		my $gene=uc((split "-", $_)[0]);
		if ($gene =~ /Y[ABCDEFGHIJKLMNOP][LR][0-9]+[WC].*/){
			$gene_name=$gene; $gene_oi=$cname{$gene};
		}
		else{
			$gene_name=$sysname{$gene}; $gene_oi=$gene;
		}
		my @allele=split( '-', $_);
		if (scalar @allele == 2 ){
			print "$input\t$gene_oi\t$_\tPointMut\tAR:$result{$_}{'perc'};DP:$result{$_}{'DP'};CN:$result{$_}{'CN'}\tAnswer:$result{$_}{'answer'}\n";

		}
		elsif (scalar @allele == 1){
			print "$input\t$gene_oi\t$gene_name\t$result{$_}{'cov'}:$result{$_}{'perc'}\tDeletion\tAnswer:$result{$_}{'answer'}\n";
		}
		else{
			die ("What do you think you're doing!?!");
		}
	}
	else{
		print "$input\t$_\tno check enforced\tOK\tOK\tAnswer:$result{$_}{'answer'}\n";
	}

}

sub check_point_mut{
	my $input=$_[0];
	my $gene=uc($_[1]);
	my $sysgene = $gene =~ /Y[ABCDEFGHIJKLMNOP][LR][0-9]+[WC].*/ ? $gene : $sysname{uc($gene)};
	my $mutation=$_[2];
	my @changes;
	my @all_f;
	my @DP;
	#extract info on what exactly we are looking for
	if ($mutation =~ /^[A-Z][0-9]+[A-Z]$/ ){
		push @changes, $mutation; #an aminoacid change was specified
	}
	else{
		#a complex mutation was specified : look it up in dictionary;
		@changes=split /\+/, $complex_mutations{"$_"} or die ("Cannot find the mutation $_ in the database")
	}
	my @out;
	foreach my $change(@changes){
		my $location;
		my ($position_in_gene) = $change =~ /[0-9]+/g;
		my ($strand) = $sysgene =~ m/(^Y[ABCDEFGHIJKLMNOP][LR][0-9]+[WC]$)/g;
		$strand =~ s/Y[ABCDEFGHIJKLMNOP][LR][0-9]+//g;
		if ($strand eq 'W') {
			my $start=$sysname_locations{$sysgene}{'start'}+3*($position_in_gene-1);
			my $end=$start+2;
			$location="$sysname_locations{$sysgene}{'chr'}:$start-$end";
		}
		elsif ($strand eq 'C') {
			my $end=$sysname_locations{$sysgene}{'end'}-3*($position_in_gene-1);
			my $start=$end-2;
			$location="$sysname_locations{$sysgene}{'chr'}:$start-$end";
		}
		else {die("Which strand is this gene located on again?");}
		#my $vep_command = "samtools mpileup -f $ref -r $location -g -t DP,DV -C0 -pm3 -F0.2 -d10000 $input 2>/dev/null| bcftools call -vm -f GQ | vep --species saccharomyces_cerevisiae --format vcf -o /dev/stdout --no_stats --no_progress --force_overwrite --offline";
		#print ("samtools mpileup -f $ref -r $location $input \n");
		open (my $f, "samtools mpileup -f $ref -r $location $input |");
		my @ref_codon;
		my @obs_codon;
		my @allele_freq;
		while (my $line =<$f>){
			chomp $line;
			my @riga=split "\t", $line;
			push @ref_codon, $riga[2];
			my $depth=$riga[3];
			my @observed_line=split "",$riga[4];
			my %observed_counts;
			++$observed_counts{uc($_)} for @observed_line;		#count unique
			$observed_counts{$_}=$depth>0 ? $observed_counts{$_}/$depth : 0  for keys %observed_counts;	#make percentage
			#print Dumper \%observed_counts;
			my $obs_nt;
			my $obs_ct;
			foreach (sort {$observed_counts{$a} <=> $observed_counts{$b}} keys %observed_counts){
				next if $observed_counts{$_} <0.2; #at least 20% of the reads suppporting a mutation in that position
				if (!( $_ eq "." or  $_ eq ",")){
					$obs_nt=$_ ;
					$obs_ct=$observed_counts{$_} ;
				}

			}
			$obs_nt=$riga[2] if ! defined $obs_nt;
			push @obs_codon, $obs_nt;
			push @allele_freq, $obs_ct if defined $obs_ct;
			push @DP, $depth;
		}
		close $f;
		my %aminoacid=gen_code();
		my $obs_mutation;
		#print Dumper \@allele_freq;
		if (@allele_freq > 0) {
			push @all_f, sum(@allele_freq)/@allele_freq;
			if ($strand eq "W"){
				$obs_mutation=$aminoacid{(join "", @ref_codon)}.$position_in_gene.$aminoacid{(join "", @obs_codon)};
			}
			elsif ($strand eq "C"){
				$obs_mutation=$aminoacid{revcomp(join "", @ref_codon)}.$position_in_gene.$aminoacid{revcomp(join "", @obs_codon)};
			}
			if ($obs_mutation eq $change){
				push @out, "Y";
			}
			else{
				push @out, "N";
			}
		}
		else {
			push @all_f, 0;
			push @out, "N";
		}
	}
	my $answer;
	if (scalar(uniq(@out)) == 1 and (uniq(@out))[0] eq "Y"){
		$answer="Y";
	}
	else{
		$answer="N";
	}
	my %output;
	$output{'answer'}=$answer;
	$output{'perc'}= sum(@all_f)/@all_f;
	my $average_depth;
	if (@DP > 0){
		$average_depth=sum(@DP)/@DP;
	}
	else{
		$average_depth=0;
	}
	my %gene_data=%{get_gene_info($gene)};
	$output{'DP'}=$average_depth;
	$output{'CN'}=$output{'DP'}/get_cwm($input, $gene_data{'chrom'});
	return (\%output);
}


sub check_deletion{
	######Set here the threshold for median gene coverage
	my $T=15;   # (15%)
	#Set here the threshold for small gaps coverage
	my $t=0.05; #(5%)
	#########
	use Carp;use strict;use warnings;
	use Getopt::Long;use List::Util qw(sum);
	my $gene_name = shift; my $input = shift; #Get the bam file as input
	if ( ! -e $input) {print "Could not find $input\n";  die}; #Exit if the input file does not exist
	my %gene_data=%{get_gene_info($gene_name)};
	my %res;
	if (defined $gene_data{'chrom'} and defined $gene_data{'st'} and $gene_data{'en'}){
		my $CWM=get_cwm($input, $gene_data{'chrom'});
		my $threshold=$CWM * $t; #Set the $threshold
		#Get the coverage of the gene of interest:
		my $range_goi = "$gene_data{'chrom'}".':'."$gene_data{'st'}".'-'."$gene_data{'en'}";
		my @mp_out = `samtools view -@ 8 -b $input -F 0x0400 \'$range_goi\'| genomeCoverageBed -dz -ibam stdin -g  `;
		my %cov; #Get the coverage of the gene of interest:
		for my $line(@mp_out){
		  	chomp $line;
		  	my($s0, $s1, $s2) = split '\t', $line;
		  	$cov{$s1}=$s2;
		}
		my @coverage;my $cons_bas=0;my $longest_cons_bas=0;my $prev_cov=1000;
		for (my $i=$gene_data{'st'}; $i<=$gene_data{'en'}; $i++){
		  	my $s2=$cov{$i} || 0; # || 0 is required because genomeCoverageBed -dz does not output positions with zero coverage
		  	push @coverage, $s2;
		  	#calculate if there are gaps in the coverage that are smaller than the whole gene
		  	if ($s2 < $threshold and $prev_cov < $threshold){
		  		$cons_bas++;
		  		$longest_cons_bas=$cons_bas if $cons_bas>$longest_cons_bas;
		  	}
		  	else{
		  		$cons_bas=0;
		  	}
		  		$prev_cov=$s2;
		}
		#Calculate statistics
		my $cov_of_interest = median(@coverage);
		my $perc_cov = sprintf("%.1f", $cov_of_interest/$CWM*100);
		$cov_of_interest = sprintf("%.0f", $cov_of_interest);

		#Output the result
		my $answer;
		if ($perc_cov < $T) {$answer="Y"} elsif ($longest_cons_bas > 9) {$answer="y"} else {$answer="N"}
			#return ("$input\t$gene_oi\t$gene_name\t$cov_of_interest\t$perc_cov %\tDeleted:$answer\n");
			$perc_cov = "$perc_cov%";
			%res=(
			  	"answer" => $answer,
				"cov" => $cov_of_interest,
				"perc" => $perc_cov
			  );
	}
	else{
		%res=("answer" => "FAILED");
	}
	return (\%res);
}

sub get_cwm{
	my $in=shift;
	my $chr=shift;
	my @mp_out = `samtools view -b $in -F 0x0400 \'$chr\' | genomeCoverageBed -dz -ibam stdin -g  `;
	my @cwm_coverage;
	for my $line(@mp_out){
		chomp $line;
		my($s0, $s1, $s2) = split '\t', $line;
		push @cwm_coverage, $s2;
	}
	return median(@cwm_coverage);
}

sub get_gene_info{
	my $gene_name=shift;
	#Get the path of the current script to localise gene_name list file
	my $path = __FILE__; my @path_components = split( /\//, $path );
	my $rem = pop @path_components; $rem = pop @path_components;
	$path = (join "/",  @path_components);
	my $gene_file="$path".'/defaults/all_yeast_genes.txt';
	my $ifi;
	# Open the Conversion file and make an array for the location of the gene:
	#       Ensembl Gene ID	  Chromosome Name	Gene Start (bp)	Gene End (bp)	Associated Gene Name
	#       YHR055C		         VIII		        214533		      214718		    CUP1-2
	if( $gene_file =~ /\.gz/ ){
	    open($ifi, qq[gunzip -c $gene_file|]);
	  }
	else{
	    open($ifi, $gene_file ) or die $!;
	  }
	my %results; my $gene_oi; my $chrom; my $st; my $en;
	while( my $line = <$ifi> ) {
	#        	chomp $line;
		  next if $line =~ /^#/ ; #header lines skipped
		  if ($line =~ /$gene_name\s/i){
		    my($f1, $f2, $f3, $f4, $f5) = split '\t', $line; #split the line into its columns
		    my $f6 = $f4 - $f3 +1;
		    $results{'gene_oi'} = $f1; $results{'chrom'} = $f2; $results{'st'} = $f3; $results{'en'} = $f4;
		  }
	}
	#print Dumper \%results;
	close( $ifi );
	return \%results;
}
sub median{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
	return $vals[int($len/2)];
    }
    else #even
    {
	return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub revcomp{
	my $revcomp = reverse $_[0];
	$revcomp =~ tr/ATGCatgc/TACGtacg/;
	return $revcomp;
}

sub gen_code{
my (%genetic_code) = (
'TCA' => 'S', # Serine
'TCC' => 'S', # Serine
'TCG' => 'S', # Serine
'TCT' => 'S', # Serine
'TTC' => 'F', # Phenylalanine
'TTT' => 'F', # Phenylalanine
'TTA' => 'L', # LeTcine
'TTG' => 'L', # LeTcine
'TAC' => 'Y', # Tyrosine
'TAT' => 'Y', # Tyrosine
'TAA' => '_', # Stop
'TAG' => '_', # Stop
'TGC' => 'C', # Cysteine
'TGT' => 'C', # Cysteine
'TGA' => '_', # Stop
'TGG' => 'W', # Tryptophan
'CTA' => 'L', # LeTcine
'CTC' => 'L', # LeTcine
'CTG' => 'L', # LeTcine
'CTT' => 'L', # LeTcine
'CCA' => 'P', # Proline
'CAT' => 'H', # Histidine
'CAA' => 'Q', # GlTtamine
'CAG' => 'Q', # GlTtamine
'CGA' => 'R', # Arginine
'CGC' => 'R', # Arginine
'CGG' => 'R', # Arginine
'CGT' => 'R', # Arginine
'ATA' => 'I', # IsoleTcine
'ATC' => 'I', # IsoleTcine
'ATT' => 'I', # IsoleTcine
'ATG' => 'M', # Methionine
'ACA' => 'T', # Threonine
'ACC' => 'T', # Threonine
'ACG' => 'T', # Threonine
'ACT' => 'T', # Threonine
'AAC' => 'N', # Asparagine
'AAT' => 'N', # Asparagine
'AAA' => 'K', # Lysine
'AAG' => 'K', # Lysine
'AGC' => 'S', # Serine
'AGT' => 'S', # Serine
'AGA' => 'R', # Arginine
'AGG' => 'R', # Arginine
'CCC' => 'P', # Proline
'CCG' => 'P', # Proline
'CCT' => 'P', # Proline
'CAC' => 'H', # Histidine
'GTA' => 'V', # Valine
'GTC' => 'V', # Valine
'GTG' => 'V', # Valine
'GTT' => 'V', # Valine
'GCA' => 'A', # Alanine
'GCC' => 'A', # Alanine
'GCG' => 'A', # Alanine
'GCT' => 'A', # Alanine
'GAC' => 'D', # Aspartic Acid
'GAT' => 'D', # Aspartic Acid
'GAA' => 'E', # GlTtamic Acid
'GAG' => 'E', # GlTtamic Acid
'GGA' => 'G', # Glycine
'GGC' => 'G', # Glycine
'GGG' => 'G', # Glycine
'GGT' => 'G'  # Glycine
);
return %genetic_code;
}
