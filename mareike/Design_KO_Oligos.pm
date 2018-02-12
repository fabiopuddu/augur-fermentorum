#!/usr/bin/perl
# Author:       mh23
# Maintainer:   mh23

# Name: yeast_KO_oligo_designer.pl

package Design_KO_Oligos;
use Exporter;
@ISA = ('Exporter');
@EXPORT = ('Design_KO_Oligos');


############

# RULES

# The F1 oligo includes 40bp immediately upstream of the ATG of the gene 
# followed by 20bp that are complementary to the Kan Cassette containing plasmid

# The R1 oligo includes 40bp immediately downstream of the stop codon of the gene 
# followed by 20bp that are complementary to the Kan Cassette containing plasmid

# The oligos for checking the deletion should be made like this 
# The .3 is upstream of the stop codon
# The .4 is downstream of stop codon
# They should be ~20bp long and end in two G/C (not more than two)
# They should be no more than 400bp away from stop codon



sub Design_KO_Oligos {
	
	my $stable_id = shift;
	
	#Tails used
	my $F1_tail = ' CGGATCCCCGGGTTAATTAA';
	my $R1_tail = ' GAATTCGAGCTCGTTTAAAC';
	my $F2_tail = $F1_tail;
#	$F2_tail = ' CGTACGCTGCAGGTCGAC';	

	my $F1_oligo='';
	my $R1_oligo='';
	my $F2_oligo='';
	my $three_oligo='';
	my $four_oligo='';


	##########################################################################
	# Step 1: Design a function that gets the coordinates for the input gene #
	##########################################################################

	#########################
	##  Connect to API     ##   
	#########################

	#First the coordinates and strand of a query gene is needed
	#This is done by using the ensembl API

	#First connect to the API


	use Bio::EnsEMBL::Registry;

	my $registry = 'Bio::EnsEMBL::Registry';

	$registry->load_registry_from_db(
		-host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
		-user => 'anonymous'
	);

	#Retrieve gene coordinated by a yeast gene name

	#First get adaptors
	my $gene_adaptor  = $registry->get_adaptor( 'Saccharomyces_cerevisiae', 'Core', 'Gene' );
	my $slice_adaptor = $registry->get_adaptor( 'Saccharomyces_cerevisiae', 'Core', 'Slice' );

	#Get gene
	my $gene = $gene_adaptor->fetch_by_stable_id($stable_id);

	#Get the coordinates and the strand of the gene in question
	my $chr = $gene->slice->seq_region_name();
	my $start = $gene->start();
	my $end = $gene->end();
	my $strand = $gene->strand();
	my $gene_name = $gene->external_name();


	####################################
	# Step 2: Get the F1 and F2 and R1 oligos #
	####################################

	#Get the gene sequence from coordinates
	my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start, $end );
	my $gene_sequence = $slice->seq();

	#For the For and the Rev oligos for a KO cassette we want 40nts up-/downstream 

	my $offset = $start - 40;
	$slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $offset, $start-1 );
	if ( $strand == 1 ) {
	$F1_oligo =  $slice->seq();}
	if ( $strand == -1 ) {
	$R1_oligo =  $slice->seq();}

	$offset = $end + 40;
	$slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $end+1, $offset );
	my $temp_seq = $slice->seq();
	if ( $strand == 1 ) {
	$R1_oligo = RevComp($temp_seq)}
	if ( $strand == -1 ) {
	$F1_oligo = RevComp($temp_seq)}

	#F2 oligo 
	if ( $strand == 1 ) {
	$offset = $end - 42;
        $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $offset, $end-3 );
	$F2_oligo =  $slice->seq();
	}
	if ( $strand == -1 ) {
	$offset = $start + 42;
	$slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start+3, $offset );
        $temp_seq = $slice->seq();
	$F2_oligo = RevComp($temp_seq);
	}

	#######################################
	# Step 3: Design the .3 and .4 oligos #
	#######################################

	#Get a 400bp either side of stop codon sequence in the coding strand

	my $sequence_stop = '';

	if ( $strand == 1 ) { 
		$slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $end-400, $end+400 );
		$sequence_stop = $slice->seq();
	}

	if ( $strand == -1 ) { 
		$slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start-400, $start+400 );
		$sequence_stop = $slice->seq();
		$sequence_stop = RevComp($sequence_stop);
	}

	#print ("\n\n$sequence_stop\n\n");

	my $three_location = '';
	my $four_location = '';

	# Start with the .3 oligo
	($three_oligo, $three_location) = Search_Check_Oligo($sequence_stop);


	# Then with the .4 oligo
	#This bit searches from outside towards the stop codon for .4
	my $right_flank = substr ($sequence_stop,400);
	$right_flank = RevComp($right_flank);
	($four_oligo, $four_location) = Search_Check_Oligo($right_flank);
	$four_location = 800-$four_location;


	#This bit searches from stop codon outwards for .4
	# for (my $i=500; $i <= 750; $i++)  {
	# 	my $length = 20;
	# 	my $oligo_t = substr ($sequence_stop,$i,$length);
	# 	
	# 	#Check conditions
	# 	my $beginning_two = substr ($oligo_t,0,2);
	# 	my $third_from_beg = substr ($oligo_t,2,1);
	# 	#print("$oligo_t\n$beginning_two\n$third_from_beg\n\n");
	# 		if ($beginning_two =~ /^[CG]+$/){
	# 		if ($third_from_beg ne "C" and $third_from_beg ne "G") {
	# 			$four_oligo = $oligo_t;
	# 			$four_oligo = RevComp($four_oligo);
	# 			$four_location = $i;
	# 			#print "$i $four_oligo\n";
	# 			last;
	# 		}
	# 	}
	# 
	# }


	##########################
	## PCR product estimate ##
	##########################

	my $product_length = $four_location - $three_location;

	my $three_from_stop = 400 - $three_location;
	my $four_from_stop = $four_location - 400; 



	#####################
	##	    Output	   ##
	#####################

	my %output;

	#Add the tails to F1 and R1
	$F1_oligo = $F1_oligo . $F1_tail ;
	$R1_oligo = $R1_oligo . $R1_tail ;
	$F2_oligo = $F2_oligo . $F2_tail ;	#The tail is the same as of the F1	

	#Fill Output
	$output{'Name'}=$stable_id;
	$output{'Standard'}=$gene_name;
	$output{'Location'}=$chr.':'.$start.'-'.$end;
	$output{'Strand'}=$strand;
	$output{'F1'}=$F1_oligo;
	$output{'R1'}=$R1_oligo;
	$output{'F2'}=$F2_oligo;
	$output{'.3'}=$three_oligo;
	$output{'.4'}=$four_oligo;
	$output{'Location3'}='-'.$three_from_stop;
	$output{'Location4'}='+'.$four_from_stop;
	$output{'PCR_length'}=$product_length;
	# print ("Name: $stable_id Location: $chr $start - $end Strand: $strand\n");
	# print ("F1: $F1_oligo\nR1: $R1_oligo\n");
	# print (".3: $three_oligo Location: -$three_from_stop from Stop Codon\n.4: $four_oligo Location: +$four_from_stop from Stop Codon\n");
	 print ("Estimated PCR length: $product_length\n"); 

	return %output;
}


#######################
## Other Subroutines ##
#######################


sub RevComp {
	my $seq = shift;
	my $revseq = scalar reverse($seq);
	$revseq =~ tr/actgACTG/tgacTGAC/;
	return $revseq;
}

sub Search_Check_Oligo {
	my $oligo; my $location;
	my $seq = shift;
	for (my $i=100; $i <= 300; $i++)  {
	my $length = 20;
	my $oligo_t = substr ($seq,$i,$length);  
	#print ("$oligo_t\n");
	
	#Check conditions
	my $end_two = substr ($oligo_t,18);
	my $third_from_last = substr ($oligo_t,17,1);
	
	if ($end_two =~ /^[CG]+$/){
		if ($third_from_last ne "C" and $third_from_last ne "G") {
			$oligo = $oligo_t;
			$location = $i;
			#print "$i $three_oligo\n";
			last;
		}
	}
}
return($oligo, $location);
}

1;
