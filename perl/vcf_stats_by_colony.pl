#!/usr/bin/env perl

use strict; 
use warnings;
use Cwd;
use Data::Dumper;
no warnings 'experimental::smartmatch';

#Declare variables
 
my $pwd = (split /\//, cwd() )[-2];
my $p=0;
$p=1 if shift eq '-p';
my @indel_count; my @snp_count;
my @transitions; my @transversions;
my @A_C;my @A_G;my @A_T;my @C_A;my @C_G;my @C_T;my @G_A;my @G_C;my @G_T;my @T_A;my @T_C;my @T_G;
my @ERS;
#INDELs
my @one;my @two;my @three;my @more;
my @minone;my @mintwo;my @minthree;my @less;
my %mutations;
my %null; my $totlines=0;
my @VCF;
if (open (my $mergehandle, '<', 'experiment_merge.vcf')) {
	@VCF=<$mergehandle>;
	close $mergehandle;
} ;
@VCF=<STDIN> unless @VCF;
my @samples;
foreach my $line (@VCF){
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
	next if $FILT ne 'PASS'; next if $QUAL<180;
	$totlines++;
	for (my $c=0; $c<scalar @samples; $c++){
		if ($mutations[$c] =~ /[0-9]\/[0-9]/){
#			print "$line\n";
			if ($INFO =~ /INDEL/){
				push @{$mutations{$samples[$c]}{'ID_INDEL'}},  "$chrom-$pos";
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
				push @{$mutations{$samples[$c]}{'ID_SNV'}},  "$chrom-$pos";
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



#Identify control sample
foreach my $k (@samples){
	$mutations{$k}{'SAMPLES'}=$k;
	$mutations{$k} = 'control' if $null{$k} == $totlines;
}

@samples= grep {$mutations{$_} ne 'control'} @samples;

#Find how many unique mutations in the file in each sample
#For every sample that is not a control
foreach my $sample1 (@samples){
	#and again for every sample
	foreach my $sample2 (@samples){
		# skip when the two sample name match, we are not interested 
		next if $sample1 eq $sample2;
		#for both the ID of SNV and INDELS
		foreach my $elem ("ID_SNV", "ID_INDEL"){
			#for every mutation in the list of mutation ID of sample 1
			foreach my $mutation (@{$mutations{$sample1}{$elem}}){
#				print "$mutation\n";
				#increase the counter only if the mutation currently examined is not found in sample2
				push @{$mutations{$sample1}{"uniq".$elem}}, $mutation unless $mutation ~~ @{$mutations{$sample2}{$elem}} 
			}
		
#		print "Sample1: $sample1 ; Sample2: $sample2; ";
#		print Dumper \@{$mutations{$sample1}{'uniqID_SNV'}};
#		print  "\n";
		}
	}
}
#

my $s=0;
foreach my $k (@samples){
	$mutations{$k}{'G>T'} = ( $mutations{$k}{'G>T'} || 0 )+ ( $mutations{$k}{'C>A'} || 0);
	$mutations{$k}{'C>T'} = ( $mutations{$k}{'C>T'} || 0 )+ ( $mutations{$k}{'G>A'} || 0);
	$mutations{$k}{'C>G'} = ( $mutations{$k}{'C>G'} || 0 )+ ( $mutations{$k}{'G>C'} || 0);
	$mutations{$k}{'A>T'} = ( $mutations{$k}{'A>T'} || 0 )+ ( $mutations{$k}{'T>A'} || 0);
	$mutations{$k}{'A>G'} = ( $mutations{$k}{'A>G'} || 0 )+ ( $mutations{$k}{'T>C'} || 0);
	$mutations{$k}{'A>C'} = ( $mutations{$k}{'A>C'} || 0 )+ ( $mutations{$k}{'T>G'} || 0);
	$mutations{$k}{'TS'} =  $mutations{$k}{'C>T'} +  $mutations{$k}{'A>G'};
	$mutations{$k}{'TV'} =  $mutations{$k}{'G>T'} +  $mutations{$k}{'C>G'} +  $mutations{$k}{'A>T'} +  $mutations{$k}{'A>C'};
#	$mutations{$k}{'uniqID_SNV_NO'} = scalar(@{$mutations{$k}{"uniqID_SNV"}});
#       $mutations{$k}{'uniqID_INDEL_NO'} = scalar(@{$mutations{$k}{"uniqID_INDEL"}});
	$mutations{$k}{'uniqID_SNV'} = uniq(@{$mutations{$k}{"uniqID_SNV"}});
	$mutations{$k}{'uniqID_INDEL'} = uniq(@{$mutations{$k}{"uniqID_INDEL"}});
	$mutations{$k}{'uniqID_SNV_no'} = scalar @{$mutations{$k}{'uniqID_SNV'}}; 
        $mutations{$k}{'uniqID_INDEL_no'} = scalar @{ $mutations{$k}{'uniqID_INDEL'}};
	$s++;	
}

#print Dumper \%mutations;



my @fields=qw[IND_COUNT SNP_COUNT TS TV C>T A>G A>T C>G G>T A>C <-3 -3 -2 -1 1 2 3 >3 uniqID_SNV_no uniqID_INDEL_no SAMPLES];
my %out;
my @output;
for my $field (@fields){
	for my $sample(@samples){
		if ($p){
			if ($field =~ /[ACTG]>[ACTG]/){
				my $percentage=sprintf ("%.2f", ( ( $mutations{$sample}{$field} || 0)/$mutations{$sample}{'SNP_COUNT'}) * 100);
				push @{$out{$field}}, $percentage."%" unless  $mutations{$sample} eq 'control';
			}
			elsif ($field =~ /[0-9]/ || $field =~ /[<>][0-9]/ ){
                                my $percentage=sprintf ("%.2f", ( ( $mutations{$sample}{$field} || 0)/$mutations{$sample}{'IND_COUNT'}) * 100);
                                push @{$out{$field}}, $percentage."%" unless  $mutations{$sample} eq 'control';
                        }
			else {
                        	push @{$out{$field}}, ($mutations{$sample}{$field} || 0) unless  $mutations{$sample} eq 'control';
                	}

		}
		else {
			push @{$out{$field}}, ($mutations{$sample}{$field} || 0) unless  $mutations{$sample} eq 'control';
		}
	}
	if (exists $out{$field}){
		push @output , (join ':', @{$out{$field}});
	}
	else{
		push @output ,"0";
	}
}	

splice @output, 4, 0, $pwd;
splice @output, -1, 0, $s;
#print Dumper \@output;

print("INDEL\tSNP\tTs\tTv\tSample \tC>T\tA>G\tA>T\tC>G\tG>T\tA>C\t<-3\t-3\t-2\t-1\t1\t2\t3\t>3\tuqSNV\tuqIND\tNumber of samples\n");
print join "\t", @output;
print "\n";

sub uniq {
  my %seen;
  my @uniq=grep { !$seen{$_}++ } @_;
  return \@uniq;
}

