
#!/usr/bin/perl -w

use strict;
use warnings;

my $i=1;
#my $tresh1= shift @ARGV; #get first threshold: minimum number of GT/AC repeats
#my $tresh2= shift@ARGV; #get second threshold: minimum GT% in the read
my $infastq= shift @ARGV; #get path to fastq file; can be omitted, it will read from STDIN
my %base_numbers; #counter to get total number of base in each read
my %hist;

my $fh;
my $is_stdin = 0; #boolean value to determine wheter STDIN should be used
if (defined $infastq){
	open my $fh, "<", $infastq or die $!;  #open provided fastq file
} 
else {
  $fh = *STDIN; #open STDIN 
  $is_stdin++; #STDIN is being used
}

my $total_reads=0;

while(my $row=<$fh>){ #read file or STDIN line by line
        next unless $row =~ /^[TAGCN]+$/; #skip line if it does not contain DNA
	next unless $row =~ /GT{1,4}|A{1,4}C/; #skip the line if the DNA does not contain any telomeric sequence
	chomp $row;
	$total_reads++;
	my %base_numbers= (
  			  	"A" => "0",	#
    				"C" => "0",	#
    				"G" => "0",	# initialise counters
				"T" => "0",	#
				"N" => "0",	#
	);
	my @sequence = split '', $row;
	foreach my $base (@sequence){		
		   $base_numbers{"$base"}++; #count each base		
	}
	my $gt_perc = ($base_numbers{'G'}+$base_numbers{'T'})/@sequence*100; #calculate GT percent
	my $count=0; #initialise counter of repeat instances
	while ($row =~ /GT{1,4}|A{1,4}C/g){
	$count++	#count repeat instances
	}
	$hist{$count}++;        

#printf "$row\n" if ($count>=$tresh1 && ($gt_perc > $tresh2 || $gt_perc < (100-$tresh2))); #output the reads if the number of matches is higher than thresh1 and the GT or AC percent is higher than tresh2
}

#Print historgram
foreach my $num (sort {$a <=> $b} keys %hist){
	$hist{$num}=$hist{$num}/$total_reads;
	print "$num\t$hist{$num}\n";

}
#print "Total:\t$total_reads\n";
