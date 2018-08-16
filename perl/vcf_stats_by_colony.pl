#!/usr/bin/env perl

use strict; 
use warnings;
use Cwd;
use Data::Dumper;

#Declare variables
 
my $pwd = (split /\//, cwd() )[-2];

my @indel_count; my @snp_count;
my @transitions; my @transversions;
my @A_C;my @A_G;my @A_T;my @C_A;my @C_G;my @C_T;my @G_A;my @G_C;my @G_T;my @T_A;my @T_C;my @T_G;
my @ERS;
#INDELs
my @one;my @two;my @three;my @more;
my @minone;my @mintwo;my @minthree;my @less;
my %mutations;
my %null; my $totlines=0;
open (my $mergehandle, '<', 'experiment_merge.vcf');
  my @samples;
while (my$line=<$mergehandle>){
	chomp $line;
	next if $line =~ /^##/;
	if ($line =~ /^#/){
			#read header
			my @row=split "\t", $line;
			@samples=splice @row, 9;
			for my $sample (@samples){
				$null{$sample}=0;
			}
			next
	}
	(my $chrom, my $pos, my $ID, my $REF, my $ALT, my $QUAL, my $FILT, my $INFO, my $FORMAT, my @mutations) = split "\t", $line;
	next if $FILT ne 'PASS';
	$totlines++;
	for (my $c=0; $c<scalar @samples; $c++){
		if ($mutations[$c] =~ /[0-9]\/[0-9]/){
#			print "$line\n";
			if ($INFO =~ /INDEL/){
				my $size = length($ALT) - length($REF);
			     	$mutations{$samples[$c]}{'IND_COUNT'}++;
				if (abs($size)<4){
					$mutations{$samples[$c]}{$size}++
				}
				elsif ($size > 3){
					$mutations{$samples[$c]}{'>3'}++
				}
				 elsif ($size < -3){
                                	$mutations{$samples[$c]}{'<-3'}++
	                        }
				else {die "Something odd happened, investigate!!!\n"}
			}
			elsif (length $REF  == 1 and length $ALT == 1){
				$mutations{$samples[$c]}{$REF.">".$ALT}++;
			        $mutations{$samples[$c]}{'SNP_COUNT'}++;
			}
			else {
				die "Bad line in vcf file\n";
			}
		}	
		elsif ($mutations[$c] =~ /^\.$/){
			$null{$samples[$c]}++;
		}
	}	
}

close $mergehandle;

#Identify control sample
foreach my $k (@samples){
	$mutations{$k}{'SAMPLES'}=$k;
	$mutations{$k} = 'control' if $null{$k} == $totlines;
}
my $s=0;
foreach my $k (@samples){
	next if $mutations{$k} eq 'control';
	$mutations{$k}{'G>T'} = ( $mutations{$k}{'G>T'} || 0 )+ ( $mutations{$k}{'C>A'} || 0);
	$mutations{$k}{'C>T'} = ( $mutations{$k}{'C>T'} || 0 )+ ( $mutations{$k}{'G>A'} || 0);
	$mutations{$k}{'C>G'} = ( $mutations{$k}{'C>G'} || 0 )+ ( $mutations{$k}{'G>C'} || 0);
	$mutations{$k}{'A>T'} = ( $mutations{$k}{'A>T'} || 0 )+ ( $mutations{$k}{'T>A'} || 0);
	$mutations{$k}{'A>G'} = ( $mutations{$k}{'A>G'} || 0 )+ ( $mutations{$k}{'T>C'} || 0);
	$mutations{$k}{'A>C'} = ( $mutations{$k}{'A>C'} || 0 )+ ( $mutations{$k}{'T>G'} || 0);
	$mutations{$k}{'TS'} =  $mutations{$k}{'C>T'} +  $mutations{$k}{'A>G'};
	$mutations{$k}{'TV'} =  $mutations{$k}{'G>T'} +  $mutations{$k}{'C>G'} +  $mutations{$k}{'A>T'} +  $mutations{$k}{'A>C'};
	$s++;	
}


my @fields=qw[IND_COUNT SNP_COUNT TS TV C>T A>G A>T C>G G>T A>C <-3 -3 -2 -1 1 2 3 >3 SAMPLES];
my %out;
my @output;
for my $field (@fields){
	for my $sample(@samples){
		push @{$out{$field}}, ($mutations{$sample}{$field} || 0) unless  $mutations{$sample} eq 'control';

	}
	push @output , (join ':', @{$out{$field}});
}

splice @output, 4, 0, $pwd;
splice @output, -1, 0, $s;
#print Dumper \@output;

print("INDEL\tSNP\tTs\tTv\tSample    \tC>T\tA>G\tA>T\tC>G\tG>T\tA>C\t<-3\t-3\t-2\t-1\t1\t2\t3\t>3\tNumber of samples\n");
print join "\t", @output;
print "\n";

