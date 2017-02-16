#!/usr/bin/env perl
use Carp;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Data::Dumper;
use Cwd 'abs_path';

my @path = split( '/' , abs_path($0));
pop(@path);
my $local_folder = join('/',@path);

system("if [ -e ploidy_data.txt ]; then rm ploidy_data.txt; fi");

#Get the bam file as input
my ($input);
my ($ploidy);
my $fil;
my $label;
my @labels;
GetOptions
(
'f|filter'   => \$fil,
'i|input=s'  => \$input,
'p|ploidy=s' => \$ploidy,
'l|label=s' => \$label,
);
( $ploidy && $input && -f $input ) or die qq[Usage: $0 \n
					 	-i <input .bam file> \n
						-p <ploidy 1,2,...> \n
						-f (optional:filter LRT, transposons and telomeres) \n
						-l (optional: label circos plots with strain name in the form of yfg1∆:Del1234:SD1234b")\n];

my $sample_name = (split "/",$input)[-1];
$sample_name =~ s/.bam//;
if (defined $label and length $label>0) {
	 @labels=split(":",$label);
	 (scalar @labels == 3) or die qq"Not enough arguments in label";
}
my $bin_size=200; #define the size of the bin to make averages
my %genome;#define a hash ;genome'  chromosomes names as keys and the coverage hash as value.
my %filter;
#Load the regions to be filtered from file into a hash
if ($fil){
	open (my $ff, '<', "$local_folder/../defaults/region-filter.txt") or die;
	while (my $line=<$ff>){
		chomp $line;
		my @linea=split("\t",$line);
		my $chr=$linea[0];
		my @insert_val=($linea[1],$linea[2]);
		push(@{$filter{$chr}}, \@insert_val);
	} 
	close ($ff);
#	print Dumper(\%filter);
}
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
	open( my $fh, '>>', $sample_name."_ploidy_data.txt");
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
    	#print "chr$cn\t$pos\t$end\t$c_c{$pos}\n";
    	printf $fh "chr$cn\t$pos\t$end\t$c_c{$pos}\n";	
	}
	close ($fh);
}


open (my $fplo, '<', $sample_name."_ploidy_data.txt");
open (my $out, '>', $sample_name."_highlights.txt");
chomp(my @PLO = <$fplo>);
my @highlight_block;
foreach my $line (@PLO){
	my @linea=split("\t",$line); 
	#if the line is a "hit" push the line into an highlight array
	if (($linea[3] < $ploidy-0.5) or ($linea[3] > $ploidy+0.5)){
		push @highlight_block, "$linea[0]\t$linea[1]\t$linea[2]\tfill_color=red\n";
	}
	#if the line is not a hit we need to print out the previous block of highlights
	else { 
	#but only if the size of the highlight array is greater than 3
		if (scalar (@highlight_block) >=3 ){
			#print the block of highlights
			for my $output_line (@highlight_block){
				printf $out $output_line;
				#print "$output_line\n";
			}
		}
		#always re-initialise the @highlight_block array
		@highlight_block=(); 
	}
	
} 
close ($fplo);
close ($out);

print "Executing circos";
#Execute circos silently 
system ("circos -silent -conf $local_folder/../defaults/circos_aneuploidy.conf -param highlights/highlight/file=SC_MFY5784096_highlights.txt -param plots/plot/file=SC_MFY5784096_ploidy_data.txt -outputfile ".$sample_name);
#If labels have been defined write annotation on the png file
if (scalar @labels>0){
	my $command="convert ".$sample_name.".png -font Helvetica -weight 70  -gravity center -pointsize 60 -annotate 0 \"$labels[0]\n\n \"  -pointsize 30 -annotate 0 \"$labels[1]   $labels[2]\" out.".$sample_name.".png";
	system($command);
}

system ("mv out.".$sample_name.".png  ".$sample_name.".png");

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
	#Get current chromosome
	my $c = $_[0];
	# execute samtools mpileup as a systems command and store the output in an array
	# define a subroutine to calculate the mean
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
	#first of all, check wether we have filtering on and if so import the coordinates of the regions to be excluded
       	my @curr_filt;
       	if ($fil){
       		 @curr_filt=@{$filter{$c}};
       	}
	my %output;
	#calculate the average every 25 
	for (my $i=0; $i < $length; $i++) {
		#get a subset of 25 values at a time
       		my @mean_values=();
       		@mean_values = splice (@values, 0, $bin_size);
       		my @mean_positions=();
        	 	@mean_positions = splice (@positions, 0, $bin_size); 
        		#check wether filtering is on
	 	my $skip=0;
	 	if ($fil){
 			foreach my $range(@curr_filt){
				#check wether the current slice overlaps any of the ranges
		 	 	if (($mean_positions[0] > $range->[0] and $mean_positions[0] < $range->[1] ) or 
		 	 	    ($mean_positions[-1] > $range->[0] and $mean_positions[-1] < $range->[1]) or
		 	 	    ($mean_positions[0] < $range->[0] and $mean_positions[-1] > $range->[0] ) or
		 	 	    ($mean_positions[0] < $range->[1] and $mean_positions[-1] > $range->[1])
		 	 	    ){
					$skip=1;	 	 	    
		 	 		#print "Filter triggered:: Start: $mean_positions[0] End: $mean_positions[-1] Filter start: $range->[0] Filter ends: $range->[1] \n";  
		 	 		last;
		 	 	}
		 	 }
		 } 
		 #print "$skip\n";
		 next if $skip;
        		 #get the first position from the list of 25
        		 my $pos = shift @mean_positions;
       		 #get the mean of the 25 coverage values
       		 my $mean = eval(join("+", @mean_values)) / @mean_values; 
       	 	 #put the output into a hash
       	 	 #print "Position = $pos\n";
        		 $output{$pos}=$mean;
	}
	#print Dumper(\%output);
	return %output;	  
}
