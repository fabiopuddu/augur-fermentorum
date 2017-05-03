#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(ceil floor);
use Data::Dumper;
#This script will go through a list of genes and attempt to plot on their 
#respective chromosomes them with gnuplot

my $bin_size=15;
my $protein_length=768;
open (my $mutationlist, '<', 'mutations.txt');

chomp(my @MUT = <$mutationlist>);
close ($mutationlist);
my @nonsense;
my @missense;
my @frameshift;
foreach my $row (@MUT){ 
	(my $gene, my $mutation)=(split "-", $row);
	my $AA;
	if ($mutation =~ /^Δ/){
		#print 'c';	
		$AA=substr $mutation, 2;
		push @nonsense, $AA;
	}
	elsif ($mutation =~ /^FS/){
		$AA=(substr $mutation, 3);
		push @frameshift, $AA;
	}
	elsif ($mutation =~ /[a-zA-Z][0-9]*[a-zA-Z]/){
		$AA=substr $mutation, 1;
		chop($AA);
		push @missense, $AA;
	}
	else{
		next;
	}
	my $previous_gene=$gene;
}
my $file=histogram($bin_size , @missense);
open(my $out, '>', 'prot_tempdata.tsv');
print $out join "\n",@$file;
close ($out);
my $command="gnuplot -e \"plen='$protein_length"."'\" /mnt/home1/jackson/fp305/sw/bin/PF/gnuplot/plot_protein_map.gpl";
printf "$command\n";
system($command);
# system("rm tempdata.tsv");





sub round{
    my($number) = shift;
    return int($number + .5 * ($number <=> 0));
}

sub histogram{
  my ($bin_width, @list, $protein_length) = @_;
  my %histogram;
  my $line;
  my @out;
  $histogram{ceil(($_ + 1) / $bin_width) -1}++ for @list;
  for my $k (sort { $a <=> $b } keys %histogram){
	my $x=$k*$bin_size;
	my $y=$histogram{$k};
	if ($y>1){
		for (my $n=1; $n<=$y; $n++){
			$line= "$x\t$n";
			push @out, $line;
		}
	}
	else{
		$line= "$x\t$y";
		push @out, $line;
	} 	
   }  
   return \@out
}



