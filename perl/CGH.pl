#!/usr/bin/env perl

# Author:		Fabio Puddu  
# Maintainer:	Fabio Puddu
# Created:	Sep 2016
# Description:

use Carp;
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);
use Data::Dumper;
use Cwd 'abs_path';
use AutoLoader qw/AUTOLOAD/;
use Math::Round;
use Parallel::Loops;

#####
my $bin_size=400; #define the size of the bin to make averages in base pairs
my $min_span_highlight=2500; #define the minimum lenght of a jump in ploidy to be reported in highlights
my $threshold=$min_span_highlight / $bin_size;
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

my @chromosomes=('I','II','III','IV','V', 'VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV', 'XVI');

#@chromosomes=('V','VII');

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
#read the data from bedtools genomeCoverageBed
my %raw_data=read_data($input);
#define a hash 'genome' with chromosomes names and average coverage across

my @val;
my %threads;
my $maxProcs = 16;
my $pl = Parallel::Loops->new($maxProcs);
my %genome;
$pl->share( \%genome );

$pl->foreach( \@chromosomes, sub{
	print "processing chromosome $_\n";
	$genome{$_} = chrom_cov( $_, \%raw_data);
});
	#my %c_c = chrom_cov($chrom, \%raw_data);
foreach my $chrom (@chromosomes){
	#$genome{$chrom}=\%c_c;
	foreach my $elem (values %{$genome{$chrom}}){
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
my %centromere_ploidy;
foreach my $chrom (@chromosomes){
	my @centromere_block=();
	my %c_c = %{$genome{$chrom}};
	#normalise all the coverage data using the genome wide median and the expected ploidy
	foreach my $k (keys %c_c){
			$c_c{$k}=$c_c{$k}/$gw_median*$ploidy; 		
			if (($k>$centromere{$chrom}[0]) and ($k<$centromere{$chrom}[1])){
   		             push @centromere_block, $c_c{$k}
        			}
          }
	$ploidy_by_chr{$chrom}=ploidy_mode(@centromere_block);
	$centromere_ploidy{$chrom}=median(@centromere_block);
	open( my $fh, '>>', $sample_name."_ploidy_data.txt");
	my $cn=get_as_chr($chrom);
	foreach my $pos (sort {$a <=> $b} keys %c_c) {
		my $start=$pos-($bin_size/2); 
		my $end = $pos+(($bin_size/2)-1);
    		printf $fh "chr$cn\t$start\t$end\t$c_c{$pos}\n";	
	}
	close ($fh);
}

###write the chromosome ploidy info in a file
open (my $statout, '>', $sample_name."_plstats.txt");
print ($statout "Chromosome\tPred.ploidy\n");
my $tot_aneup=0;
foreach my $chrom (@chromosomes){
	my $cn=get_as_chr($chrom);		
        #assign an integer number to ploidy if the mean ploidy of the chromosome does not diverge too much from the predicted
        my $chr_ploidy;
        if (abs($centromere_ploidy{$chrom}-$ploidy_by_chr{$chrom})<0.20){
                $chr_ploidy=$ploidy_by_chr{$chrom};
        }
        else{
                $chr_ploidy=$centromere_ploidy{$chrom};
        }
	$tot_aneup += abs($chr_ploidy-$ploidy);
	print ($statout "chr$cn\t$chr_ploidy\n");
}
print ($statout "Total.aneuploidies\t$tot_aneup\n");
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
my $prev_plo=$ploidy;

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
	#print "$linea[0]\t$linea[1]\tabs($linea[3]-$prev_plo)\n";
	if (abs ($linea[3]-$centromere_ploidy{$chr_name}) > 0.4 and abs($linea[3]-$prev_plo)<0.8){
		push @breakpoint_block, "$linea[0]\t$linea[1]\t$linea[2]\t$linea[3]\t$centromere_ploidy{$chr_name}";
		#print "positive\n";
	}
	else { 
		#but only if the size of the highlight array is greater than 3
		if (scalar (@breakpoint_block) >= ($threshold*4) ){
			#print the block of highlights
			#print "$threshold*4\n";
			#print Dumper \@breakpoint_block;
			my @block=collapse_region(@breakpoint_block);
			push @breakpoints, @block;
			
		}
		#always re-initialise the @highlight_block array
		@breakpoint_block=(); 
	}
	
	
	$prev_chr=$linea[0];
	$prev_plo=$linea[3];
} 
#
#print Dumper \@breakpoints;
my @final_breakpoints=collapse_region(@breakpoints);

open (my $out_brkp, '>', $sample_name."_breakpoints.txt");
print $out_brkp join("\n", @final_breakpoints);
close ($out_brkp);


sleep 2;
system("mkdir -p png; mkdir -p svg");
print "Executing circos\n";
#Execute circos silently 
system ("circos -silent -conf $local_folder/../defaults/circos_aneuploidy.conf -param highlights/highlight/file=".$sample_name."_highlights.txt -param plots/plot/file=".$sample_name."_ploidy_data.txt -outputfile ".$sample_name);
#If labels have been defined write annotation on the png file
if (scalar @labels>0){
	my $command="convert ".$sample_name.".png -font Helvetica -weight 70  -gravity center -pointsize 60 -annotate 0 \"$labels[0]\n\n \"  -pointsize 30 -annotate 0 \"$labels[1]   $labels[2]\" out.".$sample_name.".png";
	system($command);
	system ("mv out.".$sample_name.".png  ".$sample_name.".png");
}

system("mv ".$sample_name.".png png/; mv ".$sample_name.".svg svg/");
system("convert png/".$sample_name.".png -quality 96 -resize 500x500  png/".$sample_name."_web.jpg");

#############################################################
#SUBROUTINES
#############################################################

#################

sub collapse_region{
	my @input=@_;
	my @output;
	my @block;
	my $c=0;
	my $prev_chr;
	my $prev_start;
	my $prev_end;
 	my $prev_plo;
	
	#print Dumper \@input;
	#push @input, '-\t-\t-'; #put in an extra record for comparison purposes
	foreach my $line (@input){	
		my ($chr, $start, $end, $plo, $cen_plo)=split "\t", $line;	
		if ($c>0){
			if (   ($chr ne $prev_chr) or ($start-$prev_end>20000) or (abs($prev_plo-$plo)>0.8) ){
				my $block_chr= (split "\t", $block[0])[0];#process block
				my $block_start= (split "\t", $block[0])[1];#process block
				my $block_end= (split "\t", $block[-1])[2];#process block
				my @plo;
				foreach (@block){
					push @plo, (split "\t", $_)[3] 
				}
				my $block_ploidy = sprintf('%.1f', mean(@plo));
				my $chr_plo=(split "\t", $block[-1])[4];#process block
				#and push the information
				push @output, "$block_chr\t$block_start\t$block_end\t$block_ploidy\t$chr_plo";
				#reinitialise the block
				@block=();	
			}	
		}
		push @block, $line;
		$prev_chr = $chr;
		$prev_start = $start;
		$prev_end = $end;
		$prev_plo = $plo;
		$c++;
	}
	if (scalar @block > 0){
		my $block_chr= (split "\t", $block[0])[0];#process block
		my $block_start= (split "\t", $block[0])[1];#process block
		my $block_end= (split "\t", $block[-1])[2];#process block
		my @plo;
		foreach (@block){
			push @plo, (split "\t", $_)[3] 
		}
		my $block_ploidy = sprintf('%.1f', mean(@plo));
		my $chr_plo=(split "\t", $block[-1])[4];#process block
		#and push the information
		push @output, "$block_chr\t$block_start\t$block_end\t$block_ploidy\t$chr_plo";
	}
	return @output;
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
                my $rounded_element=round($element);
                $counts{$rounded_element}++;
        }
	my @sorted_array = sort { $counts{$a} <=> $counts{$b} } keys %counts;
        return $sorted_array[-1];
}


#############################################################

#THIS SUBROUTINE TAKES IN INPUT A CHROMOSOME NAME EG 'IV' AND RETURNS A HASH WITH POSITIONS AS KEYS AND COVERAGES AS VALUES

sub read_data{
	my $in=shift;
	# execute samtools mpileup as a systems command and store the output in an array
	my $command=  "samtools view -@ 8 -b $in -F 0x0400 | genomeCoverageBed -d -ibam stdin -g |" ; 
	open (my $mpileup_out, $command);
	#print @mpileup_out;
	# declare variables
	my %mpileup_hash;
	#go through the mpileup array and extract the position and coverage 
	while (my $l=<$mpileup_out>){
	        chomp( $l );    
	        my @s = split( /\t/, $l );
	        #the position is $s[1] & coverage is $s[3]
	        #make it into a hash and arrays
	        $mpileup_hash{$s[0]}{$s[1]} = $s[2];
 	}
 	close $mpileup_out;
	return %mpileup_hash;
}




sub chrom_cov{
	my $c=shift; # get the chromosome we are working on
	my $data_ref=shift; #get the data structure containing the coverage data
	my %data=%{$data_ref};
	#let's load the filter for the correct chromosome
	my @curr_filt;
       	if ($fil){
       		 @curr_filt=@{$filter{$c}};
       	}
	my %output;
	for (my $i=1; $i< scalar keys %{$data{$c}}; $i+=$bin_size){
		#go through the array in steps of bin size
		my @values=();
		my $skip=0;
		#within each bin go in step of one nucleotide
		for (my $j=1; $j<$bin_size; $j++){
			#if the filtering was chosen
			if ($fil){
	 			#we need to check if the current position falls within a region to be filtered
	 			foreach my $range(@curr_filt){
					if ((($i+$j) > $range->[0] )and (($i+$j) < $range->[1])){
						#if it does, let's raise a flag and stop searching
						$skip=1;	 	 	    
		 		 		last;
		 	 		}
		 	 	}
		 	}
		 	last if $skip; 			
			push @values, $data{$c}{$i+$j} if exists $data{$c}{$i+$j};				
		}
		next if $skip;
		$output{$i+($bin_size/2)}=mean(@values) if scalar @values>0;
	}
	return \%output;
}	

