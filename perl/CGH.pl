	use Carp;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Data::Dumper;
#Get the bam file as input
my ($input);
my ($ploidy);
GetOptions
(
'i|input=s'         => \$input,
'p|ploidy=s'         => \$ploidy,
);
( $ploidy && $input && -f $input ) or die qq[Usage: $0 -i <input .bam file> -p <ploidy 1,2,...>\n];

my $bin_size=200; #define the size of the bin to make averages

my %genome;#define a hash ;genome'  chromosomes names as keys and the coverage hash as value.

my @val;
foreach my $chrom ('I','II','III','IV','V', 'VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV', 'XVI'){
	print "processing chromosome $chrom\n";
	my %c_c = chrom_cov($chrom);
	$genome{$chrom}=\%c_c;
	foreach my $elem (values %c_c){
		push (@val,$elem); 
	}
}	
my $gw_median = median(@val);
print "Genome wide median: $gw_median\n";

foreach my $chrom ('I','II','III','IV','V', 'VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV', 'XVI'){
	my %c_c = %{$genome{$chrom}};
	foreach my $k (keys %c_c){
			$c_c{$k}=$c_c{$k}/$gw_median*$ploidy; 
		
	}
	open( my $fh, '>>', "ploidy_data.txt");
	my $cn;
	if ($chrom eq 'I'){$cn='01'}
	elsif ($chrom eq 'II'){ $cn='02'}
	elsif ($chrom eq 'III'){ $cn='03'}
	elsif ($chrom eq 'IV'){ $cn='04'}
	elsif ($chrom eq 'V'){ $cn='05'}
	elsif ($chrom eq 'VI'){ $cn='06'}
	elsif ($chrom eq 'VII'){ $cn='07'}
	elsif ($chrom eq 'VIII'){ $cn='08'}
	elsif ($chrom eq 'IX'){ $cn='09'}
	elsif ($chrom eq 'X'){ $cn='10'}
	elsif ($chrom eq 'XI'){ $cn='11'}
	elsif ($chrom eq 'XII'){ $cn='12'}
	elsif ($chrom eq 'XIII'){ $cn='13'}
	elsif ($chrom eq 'XIV'){ $cn='14'}
	elsif ($chrom eq 'XV'){ $cn='15'}
	elsif ($chrom eq 'XVI'){ $cn='16'}
	foreach my $pos (sort {$a <=> $b} keys %c_c) {
		my $end = $pos+($bin_size-1);
    	printf $fh "chr$cn\t$pos\t$end\t$c_c{$pos}\n";	
	}
	close ($fh);
}

open (my $fh, '<', "ploidy_data.txt");
open (my $out, '>>', "highlights.txt");
while (my $line=<$fh>){
	chomp $line;
	my @linea=split("\t",$line);
	print join("\t",@linea),"\n";
	if ($linea[3]<0.15){
		printf $out "$linea[0]\t$linea[1]\t$linea[2]\tfill_color=blue\n";
	}
} 
close ($fh);
close ($out);

sub median {
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

sub mean {
    return sum(@_)/@_;
}
sub chrom_cov{
	# execute samtools mpileup as a systems command and store the output in an array
	# define a subroutine to calculate the mean
    my $c = $_[0];
	my @mpileup_out =  `samtools view -b $input \'$c\'| genomeCoverageBed -d -ibam stdin -g | grep -w \'$c\'`; 
	#print @mpileup_out;
	# declare variables
	my %mpileup_hash;
	my @positions;
	my @values;
	#go through the mpileup array and extract the position and coverage 
	foreach my $l (@mpileup_out){
        chomp( $l );    
        my @s = split( /\t/, $l );
        #the position is $s[1] & coverage is $s[3]
        #make it into a hash and arrays
        $mpileup_hash{$s[1]} = $s[2];
        push(@positions, $s[1]); 
        push(@values, $s[2]); 
 	}
	#go through the hash and take the mean of each 25	
	#get length of array to determine how long the for loop should be
	my $length_ar = scalar @positions;
	my $length = $length_ar/$bin_size;
		#declare variable 
	my %output;
	#calculate the average every 25 
	for (my $i=0; $i < $length; $i++) {
       		 #get a subset of 25 values at a time
       		 my @mean_values=();
       		 @mean_values = splice (@values, 0, $bin_size);
       		 my @mean_positions=();
        	 @mean_positions = splice (@positions, 0, $bin_size); 
        	 #get the first position from the list of 25
        	 my $pos = shift @mean_positions;
       		 #get the mean of the 25 coverage values
       		 my $mean = eval(join("+", @mean_values)) / @mean_values; 
       	 	 #put the output into a hash
        	 $output{$pos}=$mean;
	}
	return %output;	  
}