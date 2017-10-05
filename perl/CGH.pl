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
	I     => [ "146465", "156582" ],
    	II    => [ "233207", "243323" ],
    	III   => [ "109385", "119501" ],
   	IV    => [ "444711", "454821" ],
          V     => [ "146987", "157104" ],
          VI    => [ "143510", "153627" ],
   	VII   => [ "491920", "502038" ],
          VIII  => [ "100586", "110703" ],
          IX    => [ "350629", "360745" ],
   	X     => [ "431307", "441425" ],
          XI    => [ "435129", "445246" ],
          XII   => [ "145828", "155947" ],
   	XIII  => [ "263031", "273149" ],
          XIV   => [ "623758", "633875" ],
          XV    => [ "321584", "331702" ],
   	XVI   => [ "550957", "561073" ],
);


my %chr_ends=(
	I     => [ "15000", "215218" ],
    	II    => [ "15000", "798184" ],
    	III   => [ "15000", "301620" ],
   	IV    => [ "15000", "1516933" ],
          V     => [ "15000", "561874" ],
          VI    => [ "15000", "255161" ],
   	VII   => [ "15000", "1075940" ],
          VIII  => [ "15000", "547643" ],
          IX    => [ "15000", "424888" ],
   	X     => [ "15000", "730751" ],
          XI    => [ "15000", "651816" ],
          XII   => [ "15000", "1063177" ],
   	XIII  => [ "15000", "909431" ],
          XIV   => [ "15000", "769333" ],
          XV    => [ "15000", "1076291" ],
   	XVI   => [ "15000", "933066" ],
);




my @chromosomes=('I','II','III','IV','V', 'VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV', 'XVI');

@chromosomes=('V',, 'VII');



###################################################
##### INPUT DATA
###################################################
my @path = split( '/' , abs_path($0));
pop(@path);
my $local_folder = join('/',@path);
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



###################################################
###READ THE FILTER IF FILTERING HAS BEEN CHOSEN
###################################################
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
#       print Dumper(\%filter);
}




###################################################
#RUN THROUGH ALL THE CHROMOSOMES AND EXTRACT COVERAGE INFORMATION FROM BAM FILES
###################################################

my %genome;#define a hash ;genome'  chromosomes names as keys and the coverage hash as value.
my @val;
foreach my $chrom (@chromosomes){
	print "processing chromosome $chrom\n";
	my %c_c = chrom_cov($chrom);
	$genome{$chrom}=\%c_c;
	foreach my $elem (values %c_c){
		push (@val,$elem); 
	}
}	
my $gw_median = median(@val);
print "Genome wide median: $gw_median\n";




###################################################
#RUN THROUGH ALL THE CHROMOSOMES AGAIN, NORMALISE DATA,
#EXTRACT CHROMOSOME PLOIDY INFO FROM CENTROMERES AND OUTPUT COVERAGE IN FORMAT READABLE BY CIRCOS
###################################################
my %ploidy_by_chr; #<- this will contain the ploidy by chromosome for future reference
foreach my $chrom (@chromosomes){
	my @centromere_block=();;
	my %c_c = %{$genome{$chrom}};
	#normalise all the coverage data using the genome wide median and the expected ploidy
	foreach my $k (keys %c_c){
			$c_c{$k}=$c_c{$k}/$gw_median*$ploidy; 		
			if (($k>$centromere{$chrom}[0]) and ($k<$centromere{$chrom}[1])){
   		             push @centromere_block, $c_c{$k}
        			}
          }
	$ploidy_by_chr{$chrom}=ploidy_mode(@centromere_block);
	open( my $fh, '>>', $sample_name."_ploidy_data.txt");
	my $cn=get_as_chr($chrom);
	foreach my $pos (sort {$a <=> $b} keys %c_c) {
		my $end = $pos+($bin_size-1);
    	#print "chr$cn\t$pos\t$end\t$c_c{$pos}\n";
    	printf $fh "chr$cn\t$pos\t$end\t$c_c{$pos}\n";	
	}
	close ($fh);
}

###write the chromosome ploidy info in a file
open (my $statout, '>', $sample_name."_plstats.txt");

print ($statout "Chromosome\tPred.ploidy\n");
foreach my $chrom (@chromosomes){
	my $cn=get_as_chr($chrom);		
	print ($statout "chr$cn\t$ploidy_by_chr{$chrom}\n");
}

close ($statout);

#NOW WE OPEN THE FILE JUST CREATED AND START WORKING ON IT

open (my $fplo, '<', $sample_name."_ploidy_data.txt");
chomp(my @PLO = <$fplo>);
close ($fplo);


my $chr_len=10;

#get starting ploidy data from the first line of the file
my $count=0;
#initialise the value of the previous chromosome
my $i=0;


my @highlight_block;
my @breakpoint_block;
my @breakpoints;
my $prev_chr='chr01';


open (my $out, '>', $sample_name."_highlights.txt");
foreach my $line (@PLO){
	my @linea=split("\t",$line);		#split the line on tabs 0:chr 1:start 2:end 3:ploidy 
	
	#>>>>HERE WE CREATE A HIGHLIGHT FILE THAT CIRCOS CAN USE TO HIGHLIGHT REGIONS OF DIFFERENT PLOIDY
	#if the line is a "hit" then push the line into an highlight array	
	if ((($linea[3] < $ploidy-0.5) or ($linea[3] > $ploidy+0.5))){
		push @highlight_block, "$linea[0]\t$linea[1]\t$linea[2]\tfill_color=red\n";
		$count++;
	}
	#if the line is not a hit we need to process the highlight block
	else { 
		#but only if the size of the highlight array is greater than 3
		if (scalar (@highlight_block) >= $threshold ){
			#print the block of highlights
			for my $output_line (@highlight_block){
				printf $out $output_line;		
			}
		}
		#always re-initialise the @highlight_block array
		@highlight_block=(); 
	}
	#>>>>AND HERE WE CREATE A BREAKPOINT FILE
	#check wether the ploidy of the current line falls outside of what is expected for that chromosome
	my $chr_name=get_as_rom($linea[0]);
	
	if (($linea[3] < $ploidy_by_chr{$chr_name}-0.5) or ($linea[3] > $ploidy_by_chr{$chr_name}+0.5)){
		push @breakpoint_block, "$linea[0]\t$linea[1]\t$linea[2]\t$linea[3]\n";
	}
	else { 
		#but only if the size of the highlight array is greater than 3
		if (scalar (@breakpoint_block) >= $threshold ){
			#print the block of highlights
			my $block=collapse_region(@breakpoint_block);
			push @breakpoints, @$block;
			
		}
		#always re-initialise the @highlight_block array
		@breakpoint_block=(); 
	}
	
	
	$prev_chr=$linea[0];
} 


my $final_breakpoints=collapse_region(@breakpoints);
print Dumper $final_breakpoints;




close ($out);

sleep 10;
system("mkdir -p png; mkdir -p svg");
print "Executing circos\n";
#Execute circos silently 
system ("circos -silent -conf $local_folder/../defaults/circos_aneuploidy.conf -param highlights/highlight/file=".$sample_name."_highlights.txt -param plots/plot/file=".$sample_name."_ploidy_data.txt -outputfile ".$sample_name);
#If labels have been defined write annotation on the png file
if (scalar @labels>0){
	my $command="convert ".$sample_name.".png -font Helvetica -weight 70  -gravity center -pointsize 60 -annotate 0 \"$labels[0]\n\n \"  -pointsize 30 -annotate 0 \"$labels[1]   $labels[2]\" out.".$sample_name.".png";
	system($command);
}

#system ("mv out.".$sample_name.".png  ".$sample_name.".png");
system("mv ".$sample_name.".png png/; mv ".$sample_name.".svg svg/");


#############################################################
#SUBROUTINES
#############################################################
sub collapse_region{
	my @input=@_;
	#print Dumper \@input;
	push @input, '-\t-\t-'; #put in an extra record for comparison purposes
	my $endcycle=scalar(@input)-1;
	my @br_ploidy;
	#print $endcycle;
	for (my $i=0; $i < $endcycle; $i++){	
		my @line=split "\t", $input[$i];
		my @next_line=split "\t", $input[$i+1];
		chomp @line; chomp @next_line;
		#If the chr on the next line equals the current one and the difference between start and end is less than 2kb
		if (($line[0] eq $next_line[0]) and ($next_line[1]-$line[2]<20000)){
			push @br_ploidy, $next_line[3];
			my $new_end=$next_line[2]; 				#the new end will be the end of the next block
			splice @input, $i+1, 1; 			#remove the next line
			my $mean_ploidy=mean(@br_ploidy);
			$input[$i]="$line[0]\t$line[1]\t$new_end\t$mean_ploidy";		#assign to the current line the new end
			$endcycle=$endcycle-1;				#reduce the number of iterations by 1 since we removed an element
			$i=$i-1;						#decrease the counter by 1 since we want to work on the same element for the next cycle	
		}
		else{
		@br_ploidy=();
		}
	#	print "$i\n";
		#print "$line[0]\t$next_line[0]\t$i\n";
	}
	pop @input; #remove the extra record put in before
	#print Dumper \@input;
	return \@input;
}



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
#############################################################
sub mean {
    return sum(@_)/@_;
}
#############################################################
sub geomean{
my $product=1;
foreach my $number (@_){
	$product=$number*$product;
}
my $log_e = log($product);
return exp($log_e/@_);
}
#############################################################
sub get_as_chr{
	my $chrom=$_[0];
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
	return $cn;
	
}
################################################################
sub get_as_rom{
	my $chrom=$_[0];
	my $cn;
	if ($chrom eq 'chr01'){$cn='I'}
	elsif ($chrom eq 'chr02'){ $cn='II'}
	elsif ($chrom eq 'chr03'){ $cn='III'}
	elsif ($chrom eq 'chr04'){ $cn='IV'}
	elsif ($chrom eq 'chr05'){ $cn='V'}
	elsif ($chrom eq 'chr06'){ $cn='VI'}
	elsif ($chrom eq 'chr07'){ $cn='VII'}
	elsif ($chrom eq 'chr08'){ $cn='VIII'}
	elsif ($chrom eq 'chr09'){ $cn='IX'}
	elsif ($chrom eq 'chr10'){ $cn='X'}
	elsif ($chrom eq 'chr11'){ $cn='XI'}
	elsif ($chrom eq 'chr12'){ $cn='XII'}
	elsif ($chrom eq 'chr13'){ $cn='XIII'}
	elsif ($chrom eq 'chr14'){ $cn='XIV'}
	elsif ($chrom eq 'chr15'){ $cn='XV'}
	elsif ($chrom eq 'chr16'){ $cn='XVI'}
	return $cn;
	
}



#############################################################
sub ploidy_mode{
        my %counts;
        foreach my $element (@_){
                $element=round($element);
                $counts{$element}++;
        }
	my @sorted_array = sort { $counts{$a} <=> $counts{$b} } keys %counts;
        return $sorted_array[-1];
}

#############################################################

#THIS SUBROUTINE TAKES IN INPUT A CHROMOSOME NAME EG 'IV' AND RETURNS A HASH WITH POSITIONS AS KEYS AND COVERAGES AS VALUES
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
