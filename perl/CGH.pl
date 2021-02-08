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
use Parallel::ForkManager;

####
#Inputs
####
my ($input,$ploidy, $fil, $organism);
GetOptions
(
'f|filter'   => \$fil,		#are we filtering telomeres, tys, LTRs rDNA, etc or not?
'i|input=s'  => \$input,	#what file are we working on?
'p|ploidy=s' => \$ploidy,	#what is the expected ploidy?
'o|organism=s' => \$organism,	#what organism are we working on?
);
( $ploidy && $input && -f $input && $organism) or die qq[Usage: $0 \n
					 	-i <input .bam file> \n
						-p <ploidy 1,2,...> \n
						-f (optional:filter LRT, transposons and telomeres) \n
						-o (organism: yeast, mouse, man")\n];


########################
#Get sample info from the input file
my $sample_name = (split "/",$input)[-1];
$sample_name =~ s/.cram|.bam//g;
my $delname = (split "/",$input)[0];
my @path = split( '/' , abs_path($0));
pop(@path);
my $local_folder = join('/',@path);
########################
#Declare some variables
my ($bin_size, $min_span_highlight,$threshold, %centromere, @chromosomes, %filter);
my $forks;
########################
#auto-configure the script depending on the organism
if ($organism =~ /yeast/){
	$bin_size=400; 					#define the size of the bin to make averages in base pairs
	$min_span_highlight=2500; 			#define the minimum lenght of a jump in ploidy to be reported in highlights
	$threshold=$min_span_highlight / $bin_size;
	if ($fil){					#Load the regions to be filtered from file into a hash if filtering has been chosen
	open (my $ff, '<', "$local_folder/../defaults/region-filter.txt") or die;
	        while (my $line=<$ff>){
	                chomp $line;
	                my @linea=split("\t",$line);
	                my $chr=$linea[0];
	                my @insert_val=($linea[1],$linea[2]);
	                push(@{$filter{$chr}}, \@insert_val);
        	}
	        close ($ff);
	}
	#define the position of the centromeres to calculate chromosome ploidy
	%centromere=(
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
	#define the chromosome names
	@chromosomes=('I','II','III','IV','V', 'VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV', 'XVI');
	$forks=8;					#define the number of parallel processes
}
elsif ($organism =~ /mouse/){
	$bin_size=10000; 				#define the size of the bin to make averages in base pairs
	$min_span_highlight=2500; 			#define the minimum lenght of a jump in ploidy to be reported in highlights
	$threshold=$min_span_highlight / $bin_size;
	#define the chromosome names
	@chromosomes=('1','2','3','4','5', '6','7','8','9','10','11','12','13','14','15', '16', '17', '18', '19', 'X', 'Y');
	$forks=3;					#define the number of parallel processes
}
elsif ($organism =~ /man/){
	$bin_size=10000; 				#define the size of the bin to make averages in base pairs
	$min_span_highlight=2500; 			#define the minimum lenght of a jump in ploidy to be reported in highlights
	$threshold=$min_span_highlight / $bin_size;
	@chromosomes=('1','2','3','4','5', '6','7','8','9','10','11','12','13','14','15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y');
	$forks=3;					#define the number of parallel processes
}

######################################################################################################
########RUN THROUGH ALL THE CHROMOSOMES AND EXTRACT COVERAGE INFORMATION FROM BAM FILES
######################################################################################################
#start the fork manager
my $pm = Parallel::ForkManager->new($forks);
my %genome;							#this will contain all the coverage data
$pm->run_on_finish( sub {
    			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
    			my $r=$data_structure_reference->{result};	#extract results info from thread datastore
    			my $c=$data_structure_reference->{chr};		#extract which chromosomome the thread was working on
			$genome{$c}=$r;					#store data in permanent variable
  		});
#Run through all the chromosomes
foreach my $c (@chromosomes){
	my $pid = $pm->start and next;				#start the fork for the code in this loop
	my %data_for_thread;					#initialise a datastore for this thread
	#TESTING COMMAND
	#my $command=  "samtools view -@ 2 -b -F 0x0400 $input $c:0-12000000 | samtools depth  -a -m0 - | head -n 12000000 |" ;
	#REAL COMMAND
	my $command=  "samtools view -b -F 0x0400 $input $c | samtools depth  -a -m0 - |" ;
	open (my $mpileup_out, $command);			#run the command
	print "processing chromosome $c\n";
	my $array_ref;
	my $n=1;
	my $count=1;
	while (my $l=<$mpileup_out>){				#loop through the output
		chomp( $l );
		my ($chr, $pos, $cov, @res)=split "\t", $l;	#extract information on chromosome, position and coverage
		$data_for_thread{$chr}{$pos} = $cov;		#store data in thread datastore
		if ($count >= 1000000*$n){			#print some stats to reassure the user
			$n++;
			print "Done ".($count/1000000)." megabases on chr $c \n";
		}
		$count++;
	}
  	close $mpileup_out;
	$array_ref = chrom_cov( $c, %data_for_thread);		#process the data from thread datastore to filter and make bin averages and store the reference of the resulting hash
	$pm->finish(0, { result => $array_ref, chr => $c });	#return to the fork manager the reference to the results and to the chromosome we were working on.
}
$pm->wait_all_children;					#wait for all threads before proceeding
my @val;
push @val, values %{$_} foreach values %genome;		#calculate genomewide median (from binned data)
my $gw_median = median(@val);
@val=(); undef @val;
print "Genome wide median: $gw_median\n";

###################################################
#RUN THROUGH ALL THE CHROMOSOMES AGAIN, NORMALISE DATA,
#EXTRACT CHROMOSOME PLOIDY INFO FROM CENTROMERES AND OUTPUT COVERAGE IN FORMAT READABLE BY CIRCOS
###################################################
my %ploidy_by_chr; 							#<- this will contain the ploidy by chromosome for future reference
my %centromere_ploidy;
my $cn;
#for biological reason we only use centromere coverage to determine a chromosome's ploidy
foreach my $chrom (@chromosomes){					#foreach chromosome
			#do this block only if we are working with yeast
	if ($organism eq 'yeast'){
		my @centromere_block=();
		foreach my $k (keys %{$genome{$chrom}}){		#for every position in the binned genome
			$genome{$chrom}{$k}=$genome{$chrom}{$k}/$gw_median*$ploidy;		#normalise coverage by genomewide median
			if (($k>$centromere{$chrom}[0]) and ($k<$centromere{$chrom}[1])){	#if we are looking at a centromeric region
   	   			push @centromere_block, $genome{$chrom}{$k}			#save the data in a temporary array
    			}
  		}
	$ploidy_by_chr{$chrom}=ploidy_mode(@centromere_block);		#store the "normal" ploidy of the centromere (not considering for subclonal variation; see subroutine)
	$centromere_ploidy{$chrom}=median(@centromere_block);		#store the "median" ploidy of the centromere, which includes subclonal variation
	$cn="chr".(get_as_chr($chrom));					#convert roman chromosome to chromosome name number for circos
	}
			#do this block only if we are working with mouse
	elsif ($organism =~ /mouse/){
		foreach my $k (keys %{$genome{$chrom}}){
			$genome{$chrom}{$k}=$genome{$chrom}{$k}/$gw_median*$ploidy;		#normalise coverage by genomewide median
		}
		$cn="mm".$chrom;					#convert chromosome number to chromosome name for circos
	}
			#do this block only if we are working with man
	elsif ($organism =~ /man/){
                foreach my $k (keys %{$genome{$chrom}}){
                        $genome{$chrom}{$k}=$genome{$chrom}{$k}/$gw_median*$ploidy;		#normalise coverage by genomewide median
                }
                $cn="hs".$chrom;					#convert chromosome number to chromosome name for circos
        }
	#OUTPUT FIRST ROUND OF RESULTS TO FILE
	open( my $fh, '>>', $delname."/ploidy/".$sample_name."_ploidy_data.txt");	#open file for appending
	foreach my $pos (sort {$a <=> $b} keys %{$genome{$chrom}}) {			#sort the hash by chromosome position
		my $start=$pos-($bin_size/2);						#calculate start (the key in the hash is the center of the bin)
		my $end = $pos+(($bin_size/2)-1);					#calculate end (the key in the hash is the center of the bin)
    		printf $fh "$cn\t$start\t$end\t$genome{$chrom}{$pos}\n";		#write to file in circos format
	}
	close ($fh);
}
#LOOK FOR BREAKPOINTS IN COVERAGE DEPTH TO REPORT AND HIGHLIGHT IN CIRCOS, BUT ONLY IF WE ARE WORKING WITH YEAST
if ($organism eq 'yeast'){
	#FIRST WE WRITE CHROMOSOME PLOIDY STASTS TO A FILE
	open (my $statout, '>', $delname."/ploidy/".$sample_name."_plstats.txt");	#write the chromosome ploidy info in a file
	print ($statout "Chromosome\tPred.ploidy\n");

	open (my $statout1, '>', $delname."/ploidy/".$sample_name."_plstats_raw.txt");	#the raw file will contain raw data, i.e. with ploidy < unit.2 o >unit.8 not rounded to the closest unit
	print ($statout1 "\tChr1\tChr2\tChr3\tChr4\tChr5\tChr6\tChr7\tChr8\tChr9\tChr10\tChr11\tChr12\tChr13\3tChr14\tChr15\tChr16\n");
	print ($statout1 "$sample_name\t");
	my @STTOUT1;
	my $tot_aneup=0;
	foreach my $chrom (@chromosomes){						#for each chromosome
		my $cn=get_as_chr($chrom);
	        my $chr_ploidy;
	        if (abs($centromere_ploidy{$chrom}-$ploidy_by_chr{$chrom})<0.20){	#if the mean ploidy of the chromosome does not diverge too much from its closest unit
			$chr_ploidy=$ploidy_by_chr{$chrom};				#assign an integer number to ploidy
		}
	        else{
			$chr_ploidy=$centromere_ploidy{$chrom};				#otherwise keep the rational figure
		}
		$tot_aneup += abs($chr_ploidy-$ploidy);					#keep a counter for total aneuploidies
		print ($statout "chr$cn\t$chr_ploidy\n");			#write the chromosome ploidy info in a file
		push @STTOUT1,$centromere_ploidy{$chrom};			#push the raw data (non rounded) in a list
	}
	print ($statout "Total.aneuploidies\t$tot_aneup\n");			#write the total number of aneuploidies
	print $statout join "\t", @STTOUT1 ;					#write the joined list to the raw data
	close ($statout);							#and close the files
	close ($statout1);

	#NOW WE OPEN THE FILE JUST CREATED AND START WORKING ON IT TO SEARCH FOR BREAKPOINTS FOR HIGHLIGHTING AND REPORTING
	open (my $fplo, '<', $delname."/ploidy/".$sample_name."_ploidy_data.txt");
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

	#PROCESS THE PLOIDY DATA FOR HIGHLIGHTING AND BREAKPOINT REARRANGEMENT
	open (my $out, '>', $delname."/ploidy/".$sample_name."_highlights.txt");	#open a highlight file that circos will use to highlight regions with deviating ploidy
	foreach my $line (@PLO){
	# PROCESS THE HIGHLIGHTING
		my @linea=split("\t",$line);								#split the line on tabs 0:chr 1:start 2:end 3:ploidy
		if ((($linea[3] < $ploidy-0.5) or ($linea[3] > $ploidy+0.5))){ 				#if the line is a "hit", i.e. is > 0.5 units different from the expected ploidy (like 2 for a diploid)
			push @highlight_block, "$linea[0]\t$linea[1]\t$linea[2]\tfill_color=red\n"; 	#then push the line into an highlight array
			$count++;									#and count the occurrences
		}
		else { 											#if the line is not a hit we need to process the highlight array
			if (scalar (@highlight_block) >= $threshold ){ 					#but only if the size of the highlight array is greater than the threshold i.e. 2.5 kb default
				for my $output_line (@highlight_block){
					printf $out $output_line;					#print the block of highlights to the circos highlight file
				}
			}
			@highlight_block=(); 								# and we clear the temp array
		}
	# PROCESS THE BREAKPOINTS
		my $chr_name=get_as_rom($linea[0]);							# get the chromosome name
		if (abs ($linea[3]-$centromere_ploidy{$chr_name}) > 0.4 				# if the line is a "hit", i.e. is > 0.4 units different from the raw ploidy of the centromere of that chromsome
			and abs($linea[3]-$prev_plo)<0.8){						# and if the ploidy is similar to the ploidy of the previous block (i.e. <0.8 units of difference)
			push @breakpoint_block, "$linea[0]\t$linea[1]\t$linea[2]\t$linea[3]\t$centromere_ploidy{$chr_name}";	#we push that slice in a breakpoint block
		}
		else {											#if the line is not a hit we need to process the breakpoint block
			if (scalar (@breakpoint_block) >= ($threshold*4) ){				# but only if the length of the block is > 4 times the threshold (i.e. 10kb by default)
				my @block=collapse_region(@breakpoint_block);				# we collapes the block (see subroutine)
				push @breakpoints, @block;						#and then we push it to a list of breakpoints
			}
			@breakpoint_block=();								#and we clear the temp array
		}
		$prev_chr=$linea[0];									#we update the previous chromosome
		$prev_plo=$linea[3];									#and the previous ploidy for the next loop
	}
	close $out;									#close the circos file
	#print Dumper \@breakpoints;
	my @final_breakpoints=collapse_region(@breakpoints);				#finally we collapse once again the breakpoint list, this is to avoid that masked regions such as transposons could split one single rearrangement in two
	#FIANLLY WRITE THE BREAKPOINTS TO A FILE
	open (my $out_brkp, '>', $delname."/ploidy/".$sample_name."_breakpoints.txt");
	print $out_brkp join("\n", @final_breakpoints);
	close ($out_brkp);
}
#############################################################
#SUBROUTINES
#############################################################
#This subroutine takes in input a list of strings with the following format
#chr	start1	end1	ploidy1	centromere_ploidy
#chr	start2	end2	ploidy2	centromere_ploidy
#and return a shorter list, merging all the data when the blocks are contiguous or nearly contiguous
##chr	start1	end2	(ploidy1+ploidy2)/2	centromere_ploidy
sub collapse_region{
	my @input=@_;
	my @output;
	my @block;
	my $c=0;
	my $prev_chr;
	my $prev_start;
	my $prev_end;
 	my $prev_plo;
	foreach my $line (@input){							#for every line
		my ($chr, $start, $end, $plo, $cen_plo)=split "\t", $line;		#split to extract variables
		if ($c>0){								#skip the first line
			if (   ($chr ne $prev_chr) 					#if the chromosome is different from the previous chromsome
				or ($start-$prev_end>20000) 				#or if this block is not contiguous* to the previous (*20 kb avoid  masked regions such as two transposons  spliting one single rearrangement in two)
				or (abs($prev_plo-$plo)>0.8) ){				#or if the ploidy of this block is substantially different from the previous
											#we need to process the block (see below how the block is defined)
					my $block_chr= (split "\t", $block[0])[0];	#the chromosome of the block will be equal to the chromosome of the first record of the block
					my $block_start= (split "\t", $block[0])[1];	#the start of the block will be equal to the start of the first record of the block
					my $block_end= (split "\t", $block[-1])[2];	#the end of the block will be equal to the end of the last record of the block
					my $chr_plo=(split "\t", $block[-1])[4];	#the chromosome ploidy of that block will be extracted from the last record as well
					my @plo;
					foreach (@block){				#|
						push @plo, (split "\t", $_)[3]		#|
					}						#|the ploidy of the block will be equal to the average ploidy of all the records
					my $block_ploidy = sprintf('%.1f', mean(@plo));	#|
					push @output, "$block_chr\t$block_start\t$block_end\t$block_ploidy\t$chr_plo";	#push all the information in our return variable
					@block=();					#flush data from block
				}
			}
		push @block, $line; 		#always push the current line to the block
		$prev_chr = $chr;		#|
		$prev_start = $start;		#|
		$prev_end = $end;		#|update the previous* variables for the next loop
		$prev_plo = $plo;		#|
		$c++;				#and increment the counter
	}
	#after running to all the data, we still need to process the last datablock as we have done before
	if (scalar @block > 0){					#but only if it contains something
		my $block_chr= (split "\t", $block[0])[0];	#the chromosome of the block will be equal to the chromosome of the first record of the block
		my $block_start= (split "\t", $block[0])[1];	#the start of the block will be equal to the start of the first record of the block
		my $block_end= (split "\t", $block[-1])[2];	#the end of the block will be equal to the end of the last record of the block
		my $chr_plo=(split "\t", $block[-1])[4];	#the chromosome ploidy of that block will be extracted from the last record as well
		my @plo;
		foreach (@block){				#|
			push @plo, (split "\t", $_)[3]		#|
		}						#|the ploidy of the block will be equal to the average ploidy of all the records
		my $block_ploidy = sprintf('%.1f', mean(@plo));	#|
		push @output, "$block_chr\t$block_start\t$block_end\t$block_ploidy\t$chr_plo";	#push the remaining information in our return variable
	}
	return @output;			#and return
}
#############################################################
sub qselect_median{
    my ($list, $k) = @_;
    my $pivot = @$list[int rand @{ $list } - 1];
    my @left  = grep { $_ < $pivot } @$list;
    my @right = grep { $_ > $pivot } @$list;
    if ($k <= @left)
    {
        return qselect(\@left, $k);
    }
    elsif ($k > @left + 1)
    {
        return qselect(\@right, $k - @left - 1);
    }
    else { $pivot }
}
#############################################################
sub quick_sort {
  my @a = @_;
  return @a if @a < 2;
  my $p = pop @a;
  quick_sort(grep $_ < $p, @a), $p, quick_sort(grep $_ >= $p, @a);
}
#############################################################
sub median {
    my @vals = quick_sort @_;
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
#This subroutine takes as input a chromsome name and a dictionary of locations to coverage
#and returns a dictionary of binned and filtered locations to coverage
sub chrom_cov{
	my $c=shift;							# get the chromosome we are working on
	my %data=@_;
	my @curr_filt;
       	if ($fil){							#let's load the filter for the correct chromosome if filtering was chosen
       		 @curr_filt=@{$filter{$c}};
       	}
	my %output;
	for (my $i=1; $i< scalar keys %{$data{$c}}; $i+=$bin_size){	#go through the array in steps of bin size
		my @values=();
		my $skip=0;
		for (my $j=1; $j<$bin_size; $j++){			#within each bin go in step of one nucleotide
			if ($fil){					#if the filtering was chosen
	 			foreach my $range(@curr_filt){		#we need to check if the current position falls within a region to be filtered
					if ((($i+$j) > $range->[0] )and (($i+$j) < $range->[1])){
						$skip=1;		#if it does, let's raise a flag and stop searching
		 		 		last;
		 	 		}
		 	 	}
		 	}
		 	last if $skip;					#if the skip flag was set, we do not want to store coverage info for this position, go to the next position
			push @values, $data{$c}{$i+$j} if exists $data{$c}{$i+$j};	#otherwise push the coverage of the current nucleotide in a temp array
		}
		next if $skip;						#if the skip flag was set for one nucleotide of this bin, no need to to anything else for this bin, go straight to the next bin
		$output{$i+($bin_size/2)}=mean(@values) if scalar @values>0;	#otherwise we calculate the mean value of the bin and store it in the returning hash
	}
	return \%output;
}
#############################################################
