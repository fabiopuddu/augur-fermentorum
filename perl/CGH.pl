#!/usr/bin/env perl
use Carp;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Data::Dumper;
use Cwd 'abs_path';
use Math::Round;
#####
my $bin_size=200; #define the size of the bin to make averages in base pairs
my $min_span_highlight=2000; #define the minimum lenght of a jump in ploidy to be reported in highlights
my $threshold=int($min_span_highlight / $bin_size);
#####
#define a hash of centromere coordinates(+5000 bp on both sides)
my %centromere=(
	chr01    => [ "146465", "156582" ],
    	chr02    => [ "233207", "243323" ],
    	chr03    => [ "109385", "119501" ],
   	chr04    => [ "444711", "454821" ],
        chr05    => [ "146987", "157104" ],
        chr06    => [ "143510", "153627" ],
   	chr07    => [ "491920", "502038" ],
        chr08    => [ "100586", "110703" ],
        chr09    => [ "350629", "360745" ],
   	chr10    => [ "431307", "441425" ],
        chr11    => [ "435129", "445246" ],
        chr12    => [ "145828", "155947" ],
   	chr13    => [ "263031", "273149" ],
        chr14    => [ "623758", "633875" ],
        chr15    => [ "321584", "331702" ],
   	chr16    => [ "550957", "561073" ],
);

#####
my @path = split( '/' , abs_path($0));
pop(@path);
my $local_folder = join('/',@path);
#system("if [ -e ploidy_data.txt ]; then rm ploidy_data.txt; fi");
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
open (my $statout, '>', $sample_name."_plstats.txt");#
open (my $firstdev, '>', $sample_name."_deriv.txt");

chomp(my @PLO = <$fplo>);
close ($fplo);
my @highlight_block;
my @chromosome_block;
my @centromere_block;
my $chr_len=10;
my $prev_chr='chr01';
#get starting ploidy data from the first line of the file
my @firstline=split("\t",$PLO[0]);
my $prev_plo=$firstline[3];
my $count=0;
#initialise the value of the previous chromosome
my $chrom_end=0;
my $aneup_fract;
my $prev_pos=0;
my @average_array;
my $i=0;
my $average;
my $flag_up=0;
my $flag_down=0;
my $flag_end=0;
my $brkpt1='';
my $brkpt2='';
my @breakpoints;
print ($statout "Chromosome\tPred. ploidy\tFr.deviating\n");
foreach my $line (@PLO){
	my @linea=split("\t",$line);		#split the line on tabs 0:chr 1:start 2:end 3:ploidy 
	if ( $prev_chr ne $linea[0]){		# if the current chromosome doesnt match the old we need to calculate per chromosome information
		
		my $mean_ploidy;
		$aneup_fract=$count/$chr_len;
		if ($aneup_fract <0.25 or $aneup_fract > 0.75){
			$mean_ploidy=ploidy_mode(@chromosome_block);
		}
		else{
			$mean_ploidy="-1"
		}
		$chr_len=0;
		$count=0;
		my $centromere_ploidy=ploidy_mode(@centromere_block);
		print ($statout "$prev_chr\t$centromere_ploidy\t$aneup_fract\n");
		@chromosome_block=();
		@centromere_block=();
	}
	$chr_len++;	
	# if the current slice matches one of the centromeres push the information into a list.
	if ((($linea[1] >$centromere{$linea[0]}[0]) and  #if the start of the slice is after the start of the centromere and
	     ($linea[1] <$centromere{$linea[0]}[1])) or  # the start of the slice is before the end of the centromere or
	    (($linea[2] >$centromere{$linea[0]}[0]) and #the end of the slice is after the start of the centromere and 
	     ($linea[2] >$centromere{$linea[0]}[0]))){   # the end of the slice is before the end of the centromere
		push @centromere_block, "$linea[3]"
	}
#HERE WE CREATE A HIGHLIGHT FILE THAT CIRCOS CAN USE TO HIGHLIGHT REGIONS OF DIFFERENT PLOIDY
	
	push @chromosome_block, "$linea[3]";	

	#if the line is a "hit" and we have not reached the end of the chromosome
		#then push the line into an highlight array	
	if ((($linea[3] < $ploidy-0.5) or ($linea[3] > $ploidy+0.5)) and ($prev_chr eq $linea[0])){
		push @highlight_block, "$linea[0]\t$linea[1]\t$linea[2]\tfill_color=red\n";
		$count++;
	}
	#if the line is not a hit or we reached the end of the chromosome 
	else { 
	#but only if the size of the highlight array is greater than 3
		if (scalar (@highlight_block) >= $threshold ){
			#print the block of highlights
			for my $output_line (@highlight_block){
				printf $out $output_line;		
			}
			#Detect  breakpoint coordinates
			my @highlight_start=split"\t", $highlight_block[0];
			my @highlight_end=split"\t", $highlight_block[-1];
			my $breakpoint_start=$highlight_start[1];
			my $breakpoint_end=$highlight_end[2];
			if ($prev_chr eq $linea[0]){ #if we finished because the ploidy has come back to normal > print the current chromosome
				push @breakpoints, "$linea[0]\t$breakpoint_start\t$breakpoint_end\n";
			}
			else {		#if we finished because the chromosome has finished it means linea[0] now contains the next line; therefore print the previous choromosome
				push @breakpoints, "$prev_chr\t$breakpoint_start\t$breakpoint_end\n"
			}
		}
		#always re-initialise the @highlight_block array
		@highlight_block=(); 
	}
	$prev_chr=$linea[0];
} 

#HERE WE USE THE HIGHLIGHT FILE TO DETECT BREAKPOINTS

print Dumper \@breakpoints;
push @breakpoints, '-\t-\t-'; #put in an extra record for comparison purposes
my $endcycle=scalar(@breakpoints)-1;
#print $endcycle;
for (my $i=0; $i < $endcycle; $i++){	
	my @line=split "\t", $breakpoints[$i];
#	printf "Next line =  $breakpoints[$i+1]";
#	printf "Ciclo $i; Endcycle: $endcycle";
	my @next_line=split "\t", $breakpoints[$i+1];
	chomp @line; chomp @next_line;
	#print ("$line[0] is $next_line[0] & $next_line[1]-$line[2]<2000\n");	
	#If the chr on the next line equals the current one and the difference between start and end is less than 2kb
	if (($line[0] eq $next_line[0]) and ($next_line[1]-$line[2]<15000)){
		my $new_end=$next_line[2]; 				#the new end will be the end of the next block
		splice @breakpoints, $i+1, 1; 			#remove the next line
		$breakpoints[$i]="$line[0]\t$line[1]\t$new_end";		#assign to the current line the new end
		$endcycle=$endcycle-1;				#reduce the number of iterations by 1 since we removed an element
		$i=$i-1;						#decrease the counter by 1 since we want to work on the same element for the next cycle	
	}
#	print "$i\n";
	#print "$line[0]\t$next_line[0]\t$i\n";
}
pop @breakpoints; #remove the extra record put in before
print Dumper \@breakpoints;

my $mean_ploidy;
$aneup_fract=$count/$chr_len;
if ($aneup_fract <0.25 or $aneup_fract > 0.75){
			$mean_ploidy=ploidy_mode(@chromosome_block);
		}
		else{
			$mean_ploidy="."
		}
my $centromere_ploidy=ploidy_mode(@centromere_block);
print ($statout "$prev_chr\t$centromere_ploidy\t$aneup_fract\n");
close ($statout);

close ($out);
close ($firstdev);
sleep 10;
system("mkdir -p png; mkdir -p svg");
print "Executing circos";
#Execute circos silently 
system ("circos -silent -conf $local_folder/../defaults/circos_aneuploidy.conf -param highlights/highlight/file=".$sample_name."_highlights.txt -param plots/plot/file=".$sample_name."_ploidy_data.txt -outputfile ".$sample_name);
#If labels have been defined write annotation on the png file
if (scalar @labels>0){
	my $command="convert ".$sample_name.".png -font Helvetica -weight 70  -gravity center -pointsize 60 -annotate 0 \"$labels[0]\n\n \"  -pointsize 30 -annotate 0 \"$labels[1]   $labels[2]\" out.".$sample_name.".png";
	system($command);
}

system ("mv out.".$sample_name.".png  ".$sample_name.".png");
system("mv ".$sample_name.".png png/; mv ".$sample_name.".svg svg/");


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

sub geomean{
my $product=1;
foreach my $number (@_){
	$product=$number*$product;
}
my $log_e = log($product);
return exp($log_e/@_);
}

sub ploidy_mode{
        my %counts;
        foreach my $element (@_){
                $element=round($element);
                $counts{$element}++;
        }
	my @sorted_array = sort { $counts{$a} <=> $counts{$b} } keys %counts;
        return $sorted_array[-1];
}

sub chrom_cov{
	#Get current chromosome
	my $c = $_[0];
	# execute samtools mpileup as a systems command and store the output in an array
	# define a subroutine to calculate the mean
	my @mpileup_out =  `samtools view -b $input \'$c\'| genomeCoverageBed -dz -ibam stdin -g  `; 
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
