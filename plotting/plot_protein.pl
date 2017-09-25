#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(ceil floor);
use Data::Dumper;
use File::Basename;
#get the directory the script is in
my $dirname = dirname(__FILE__);

# This script will go through a list of mutations in the following format
# TOP1-X123Y => missense
# TOP1-∆123  => nonsense
# TOP1-FS@22-24 => frameshift
# and will create a graphical represenation of their position using a gnuplot script 
# This script also allows to plot two different sets of mutations (eg. spontaneous and EMS induced)
# The mutations belonging to the second set should be prepended by an asterisk '*' in the input file
#
# Usage: plot_protein.pl input.txt [missense|nonsense|frameshift] protein_length

my $bin_size=10;
my $filename=shift;
my $mutation_analysed = shift;
my $protein_length=shift;

open (my $mutationlist, '<', $filename);
chomp(my @MUT = <$mutationlist>);
close ($mutationlist);
my @nonsense;
my @missense;
my @frameshift;
my $gene='';
foreach my $row (@MUT){
	($gene, my $mutation)=(split "-", $row);
	my $AA;
	my $EMS='';
	if ($gene =~ /\*/){
		$gene = substr $gene,1;
		$EMS='*' 
	};
	if ($mutation =~ /^Δ/){
		#print 'c';	
		$AA=substr $mutation, 2;
		push @nonsense, $AA.$EMS;
	}
	elsif ($mutation =~ /^FS/){
		$AA=(substr $mutation, 3);
		push @frameshift, $AA.$EMS;
	}
	elsif ($mutation =~ /[a-zA-Z][0-9]*[a-zA-Z]/){
		$AA=substr $mutation, 1;
		chop($AA);
		push @missense, $AA.$EMS;
	}
	else{
		next;
	}
	my $previous_gene=$gene;
}
my $file;
if ($mutation_analysed eq 'missense'){
	$file=histogram($bin_size , @missense);
}
elsif ($mutation_analysed eq 'nonsense'){
	$file=histogram($bin_size , @nonsense);
}
elsif ($mutation_analysed eq 'frameshift'){
	$file=histogram($bin_size , @frameshift);
}
else{
	exit 1
}

open(my $out, '>', 'protein_data.tsv');
print $out join "\n", @$file;
close ($out);

my $command="gnuplot -e \"plen='$protein_length"."'; outfile='$gene"."_"."$mutation_analysed.png"."'\" $dirname/plot_protein_map.gpl";
printf "$command\n";
system($command);
system("rm protein_data.tsv");

sub histogram{
  my ($bin_width, @list, $protein_length) = @_;
  my %spontaneous;
  my %EMS;
  my %all_mut;
  my $line;
  my @out;
  for my $mut(@list){
  	if ($mut =~ /\*/){
  		chop($mut); 
  	 	$EMS{ceil(($mut + 1) / $bin_width) -1}++;
  		$all_mut{ceil(($mut + 1) / $bin_width) -1}++
  	}
  	else{
  		$spontaneous{ceil(($mut + 1) / $bin_width) -1}++ ;
  		$all_mut{ceil(($mut + 1) / $bin_width) -1}++
  	}
  }
  for my $k (sort { $a <=> $b } keys %all_mut){
	my $x=$k*$bin_size;
	my $y=$all_mut{$k};
	if ($y>1){
		my $ys = 0;
		$ys = $spontaneous{$k} if defined $spontaneous{$k}  ;
		my $yE = 0;
		$yE = $EMS{$k} if defined $EMS{$k}  ;
		for (my $n=1; $n<=$ys; $n++){
			$line= "$x\t$n\tNA";
			push @out, $line;
		}
		for (my $m=$ys+1; $m<=$y; $m++){
			$line= "$x\tNA\t$m";
			push @out, $line;
		}
	}
	else{
		if (defined $spontaneous{$k}){
			$line= "$x\t$y\tNA";
		}
		elsif (defined $EMS{$k}){		
			$line= "$x\tNA\t$y";
		}
		push @out, $line;
	} 	
   }  
   return \@out
}



