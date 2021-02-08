#!/usr/bin/env perl
use strict;
use warnings;

my $rep_file=shift;

print STDERR "Reading repDNA file\n";
my $fh;
open $fh, "<", $rep_file or die $!;  #open provided repDNA file
chomp(my @REP_FILE=<$fh>);
close $fh;
my @normalised;
my @header;
foreach my $line (@REP_FILE){
  if ($line =~ /^#/){
    @header=split "\t", $line;

    #print "$header[1]\t$header[22]\t$header[14]\n";
  }
  else{
    my @riga=split "\t", $line;
    #print "$riga[1]\t$riga[26]\t$riga[14]\n" if $riga[14] =~ "RTT109";
    #print "$riga[1]\t$riga[22]\t$riga[14]\n";
    $riga[1]=($riga[1]/$riga[26])*2; #normalise rDNA on chr12
    $riga[2]=($riga[2]/$riga[22])*2; #normalise CUP1 on chr8
    push @normalised, (join "\t", @riga);
  }
}
print join "\t", @header;
print "\n";
print "$_\n" foreach @normalised;
