#!/usr/bin/env perl
use Carp;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Data::Dumper;
use Cwd 'abs_path';
#####
my $bin_size=200; #define the size of the bin to make averages in base pairs
my $min_span_highlight=2000; #define the minimum lenght of a jump in ploidy to be reported in highlights
my $threshold=int($min_span_highlight / $bin_size);
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
open (my $firstdev, '>', $sample_name."_deriv.txt");

chomp(my @PLO = <$fplo>);
my @highlight_block;
#get starting ploidy data from the first line of the file
my @firstline=split("\t",$PLO[0]);
my $prev_plo=$firstline[3];
#initialise the value of the previous chromosome
my $chrom_end=0;

my $prev_pos=0;
my @average_array;
my $i=0;
my $average;
my $flag_up=0;
my $flag_down=0;
my $flag_end=0;
my $brkpt1='';
my $brkpt2='';
foreach my $line (@PLO){

	#THE FIRST BLOCK OF CODE IDENTIFIES REGIONS WITH "SIGNIFICANT" PLOIDY CHANGES 
	#AND HIGHLIGHTS THEM IN A FILE

	my @linea=split("\t",$line); 
	#if the line is a "hit" push the line into an highlight array
	if (($linea[3] < $ploidy-0.5) or ($linea[3] > $ploidy+0.5)){
		push @highlight_block, "$linea[0]\t$linea[1]\t$linea[2]\tfill_color=red\n";
	}
	#if the line is not a hit we need to print out the previous block of highlights
	else { 
	#but only if the size of the highlight array is greater than 3
		if (scalar (@highlight_block) >= $threshold ){
			#print the block of highlights
			for my $output_line (@highlight_block){
				printf $out $output_line;
				#print "$output_line\n";
			}
		}
		#always re-initialise the @highlight_block array
		@highlight_block=(); 
	}

	#THE SECOND BLOCK OF CODE TRIES TO DETECT THE COORDINATES OF THE BREAKPOINTS
	
	my @next_line = split("\t",$PLO[$i+1]);
	my $next_line_chr = $next_line[0];
	print "$PLO[$i+1]\n";
	printf "cur line $linea[1] cur chr $linea[0] next chr $next_line_chr\n";
	if ($linea[0] ne $next_line_chr){
		$chrom_end=1
	}
	#push into average array
	push @average_array, $linea[3];
	#if n=20 average the array and report
	if (($i==40) or $chrom_end){
		$i=0;
		$average=sum(@average_array)/scalar(@average_array);
		#printf "$average\n";
		@average_array=();
		#calculate dx and report
		#my $dx=-$prev_pos;
		my $dy=($average-$prev_plo);
		printf $firstdev "$linea[0]\t$dy\n";
		#Determine if dy is outside a threshold
		if ($dy > 0.4){$flag_up = 1; $brkpt1=$linea[1];print"up $linea[0] $linea[1]\n";}
		elsif ($dy < -0.4) {$flag_down = 1; $brkpt2=$linea[1];print"down $linea[0] $linea[1]\n"}
		elsif ($chrom_end){$flag_end=1}
		if (($flag_up && $flag_down ) || (($flag_down || $flag_up) && $flag_end) ) {
			print "Curr_line $linea[1] Flags: up $flag_up down $flag_down end $flag_end\n";
			$brkpt1=$linea[1] unless $flag_up;
			$brkpt2=$linea[1] unless $flag_down;
			if ($brkpt1>$brkpt2){
				print "Breakpoint: $linea[0]\t$brkpt2\t$brkpt1\n";
			}
			else {
				print "Breakpoint: $linea[0]\t$brkpt1\t$brkpt2\n";
			}
			$flag_up = 0;
			$flag_down = 0;
			$brkpt1 = 0;
			$brkpt2 = 0;
		}

		$prev_pos=$linea[1];
		$prev_plo=$average;
	}
	$i++;

} 
close ($fplo);
close ($out);
close ($firstdev);
sleep 10;
system("mkdir png; mkdir svg");
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
