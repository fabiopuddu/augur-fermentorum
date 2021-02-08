#!/usr/bin/env perl
# Author:       Inigo Ayestaran
# Maintainer:   Inigo Ayestaran
# Created:      Sep 2017
# Description:

use strict;
use Cwd;
use warnings;
use Parallel::Loops;
use List::Util qw(sum uniq reduce );
use Data::Dumper;



#Get the path to reference files
my $path = __FILE__;
my @path_components = split( /\//, $path );
my $rem = pop @path_components;
my $rem2 = pop @path_components;
$path = (join "/",  @path_components);
my $tag_seq_ref="$path".'/defaults/strain_heterozygous_diploid.txt';
my $gene_name_ref="$path".'/defaults/all_yeast_genes.txt';
my $common_seq_ref="$path".'/defaults/common_tag_seqs.txt';

# Import sequences that are common to all samples and next to the gene-specific 20mer
open (my $ifh, $common_seq_ref) or die $!;
chomp(my @seqs = <$ifh>);
close ($ifh);

(undef, my $upbef, my $upaft) = split("\t", $seqs[1]);
(undef, my $dnbef, my $dnaft) = split("\t", $seqs[2]);

# Also check reverse complementary sequences
my $revupbef= reverse($upbef);
$revupbef =~ tr/ACGT/TGCA/;
my $revupaft= reverse($upaft);
$revupaft =~ tr/ACGT/TGCA/;
my $revdnbef= reverse($dnbef);
$revdnbef =~ tr/ACGT/TGCA/;
my $revdnaft= reverse($dnaft);
$revdnaft =~ tr/ACGT/TGCA/;


# Define all Fastq files
# my $fqdir = "$dir"."/BAM";
# opendir(my $dh, $fqdir);
# my @fqfiles = grep(/\.fq[12]\.gz/, readdir($dh));
# closedir($dh);
# @fqfiles = sort @fqfiles;
# my @sdnames; my $fq;
# for (my $fi=0; $fi<$#fqfiles; $fi++){
#     $fq = $fqfiles[$fi];
#     $fq =~ s/\.fq[12]\.gz//;
#     $sdnames[$fi] = $fq;
# }
#
# @sdnames = uniq @sdnames;
my $SDname=shift;
my $DELdir=shift;
$DELdir = $DELdir =~ /Del[0-9]+_.*/ ? $DELdir."/": "";
# Find occurrences of the tag sequences in all fastq files
# Run in parallel
# my $maxProcs = 20;
# my $pl = Parallel::Loops->new($maxProcs);
my %returnValues;
# $pl->share( \%returnValues );
# $pl->foreach( [ 0..($#sdnames) ], sub {
        # my $l = $_;
        my $prefix = $DELdir."BAM/"."$SDname";
        #my @opfiles = ( "$prefix".".fq1.gz" , "$prefix".".fq2.gz");
        # Parameters derived from which sequence is found
        my $up_down; my $reverse; my $bef_aft; my $seq;
        # Dictionary for genewise results
        my %genewise;
        #foreach my $opf (@opfiles){
            open($ifh, qq[samtools view -F 2048 "$prefix.cram" | cut -f10 |]) or die $!;
            while ( my $line = <$ifh> ){
                chomp $line;
                if ($line =~ $upbef){
                    $seq = $upbef;
                    $up_down = 'up';
                    $reverse = 0;
                    $bef_aft = 1;
                }
                if ($line =~ $upaft){
                    $seq = $upaft;
                    $up_down = 'up';
                    $reverse = 0;
                    $bef_aft = 0;
                }
                if ($line =~ $revupbef){
                    $seq = $revupbef;
                    $up_down = 'up';
                    $reverse = 1;
                    $bef_aft = 1;
                }
                if ($line =~ $revupaft){
                    $seq = $revupaft;
                    $up_down = 'up';
                    $reverse = 1;
                    $bef_aft = 0;
                }
                if ($line =~ $dnbef){
                    $seq = $dnbef;
                    $up_down = 'dn';
                    $reverse = 0;
                    $bef_aft = 1;
                }
                if ($line =~ $dnaft){
                    $seq = $dnaft;
                    $up_down = 'dn';
                    $reverse = 0;
                    $bef_aft = 0;
                }
                if ($line =~ $revdnbef){
                    $seq = $revdnbef;
                    $up_down = 'dn';
                    $reverse = 1;
                    $bef_aft = 1;
                }
                if ($line =~ $revdnaft){
                    $seq = $revdnaft;
                    $up_down = 'dn';
                    $reverse = 1;
                    $bef_aft = 0;
                }
                next unless ($up_down); # Skip line if up_down is not defined (aka no match)

                # Based on the parameters, select the sequence that belongs to the unique 20mer
                my @splitres = split($seq, $line);
                next unless ((scalar @splitres) == 2); # To avoid cases where $seq is in one of the ends, leaving one single subsequence
                my $index = abs($bef_aft - $reverse);
                my $s20mer = $splitres[$index];

                # Lookup the 20mer in the reference table and get what gene it belongs to
                open (my $iftag, $tag_seq_ref) or die $!;
                my $tagseq; my $matchcounter = 0; my $matchgene;
                while ( my $tagref = <$iftag> ){
                    chomp $tagref;
                    if ($up_down eq 'up'){
                        $tagseq = uc((split "\t", $tagref)[22]);
                    }
                    if ($up_down eq 'dn'){
                        $tagseq = uc((split "\t", $tagref)[23]);
                    }
                    $tagseq =~ s/^\s+|\s+$//g;
                    next unless ($s20mer =~ $tagseq);
                    $matchcounter++;
                    if ($matchcounter > 1){  # Avoid multiple matches
                        last;
                    }
                    $matchgene = uc((split "\t", $tagref)[1]);
                    $matchgene =~ s/^\s+|\s+$//g;

                }
                close ($iftag);
                next if ($matchcounter != 1);
                $genewise{$matchgene}++;
            }
            close($ifh);
        #}

        $returnValues{$SDname} = \%genewise;
  #  });


# Get common name of expected gene
my $gene = uc((split '_', $DELdir)[-1]);
$gene =~ s/^\s+|\s+$//g;
$gene = substr $gene, 0, -1;
# Print results to output file

# for (my $i=0; $i<=$#sdnames; $i++) {
    open (my $outfile, '>', $DELdir.'del+bcd/'.$SDname.'.bcd');
    my $reshash = $returnValues{$SDname};
    my %reshash2 = %$reshash;
    # printf "$sdnames[$i]\n";
    # print Dumper(\%reshash2);
    print $outfile "$gene\t$SDname\t";
    # Look up common name from systematic name
    my %commonname;
    open (my $ifgene, $gene_name_ref) or die $!;
    while ( my $geneline = <$ifgene> ){
            chomp $geneline;
            my ($systname, $aka) = (split "\t", $geneline)[0, 4];
	          if (! $aka){
                    $aka = $systname;
	          }
            $commonname{$systname} = $aka;
    }
    close ($ifgene);
    my $tot=sum values %reshash2;
    foreach my $key (sort {$reshash2{$b}<=>$reshash2{$a}} keys %reshash2) {
            my $printcomname = $commonname{$key};
	          my $perc = sprintf ("%.1f", $reshash2{$key}/$tot*100);
            print $outfile "$printcomname $perc% ($reshash2{$key});";
    }
    print $outfile "\n";
    close ($outfile);
# }
