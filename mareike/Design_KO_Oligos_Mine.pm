
#!/usr/bin/perl


#################


package Design_KO_Oligos_Mine;
use Exporter;
@ISA = ('Exporter');
@EXPORT = ('Design_KO_Oligos_Mine');



#################

# Basic Package Info

# This module returns a hash that contains information that you might find useful

# See below for the types of Keys you can use for the hash and what their output is

# Standard_Name => Returns the standard yeast gene name

# Systematic_Name => Returns the systematic yeast gene name

# Whole_seq => Returns the genomic sequence of the gene +1kb up- and downstream

# F1 => Returns the F1 oligo

# R1 => Returns the R1 oligo

# .3 => Returns the .3 oligo

# .4 => Returns the .4 oligo

# PCR_length => Returns the length of the expected PCR product of the .3 and .4 reaction 

#################



#################

# OLIGO DESIGN RULES

# The F1 oligo includes 40bp immediately upstream of the ATG of the gene 
# followed by 20bp that are complementary to the Kan Cassette containing plasmid

# The R1 oligo includes 40bp immediately downstream of the stop codon of the gene 
# followed by 20bp that are complementary to the Kan Cassette containing plasmid

# The oligos for checking the deletion should be made like this 
# The .3 is upstream of the stop codon
# The .4 is downstream of stop codon
# They should be ~20bp long and end in two G/C (not more than two)
# They should be no more than 400bp away from stop codon





######################################################################
# YEASTMINE INFORMATION
# Part of this is an automatically generated script to run your query.
# To use it you will require the InterMine Perl client libraries.
# These can be installed from CPAN, using your preferred client, eg:
#
#    sudo cpan Webservice::InterMine
#
# For help using these modules, please see these resources:
#
#  * https://metacpan.org/pod/Webservice::InterMine
#       - API reference
#  * https://metacpan.org/pod/Webservice::InterMine::Cookbook
#       - A How-To manual
#  * http://intermine.readthedocs.org/en/latest/web-services
#       - General Usage
#  * http://iodoc.labs.intermine.org
#       - Reference documentation for the underlying REST API
#
# Obtained from: http://yeastmine.yeastgenome.org/yeastmine/template.do?name=Gene_Flanking_Sequence&scope=all
#
######################################################################

sub Design_KO_Oligos_Mine {


my $gene_name = shift;
my $stand_name = '';

#Tails used
my $F1_tail = ' CGGATCCCCGGGTTAATTAA';
my $R1_tail = ' GAATTCGAGCTCGTTTAAAC';

# Step 1: Get the sequence
########### YEASTMINE CODE #######################

# Set the output field separator as tab
$, = "\t";
# Print unicode to standard out
binmode(STDOUT, 'utf8');
# Silence warnings when printing null fields
no warnings ('uninitialized');

# This code makes use of the Webservice::InterMine library.
# The following import statement sets YeastMine as your default
use Webservice::InterMine 0.9904 'http://yeastmine.yeastgenome.org/yeastmine';

# Description: For a given gene(s), get a chosen length of upstream and/or
# downstream sequence along with the gene sequence. Choose to retrieve 0.5 kb,
# 1.0 kb, or 2.0 kb of upstream or downstream flanking regions, or both. To
# retrieve only flanking regions without the coding sequence, select “false” for
# the “Gene Flanking Region > Include Gene” option (for this option, you must
# select only “upstream” or “downstream” flanking regions, not “both”).

my $template = Webservice::InterMine->template('Gene_Flanking_Sequence')
    or die 'Could not find a template called Gene_Flanking_Sequence';

# Use an iterator to avoid having all rows in memory at once.
my $it = $template->results_iterator_with(
    # B:  Gene
    opB    => 'LOOKUP',
    valueB => $gene_name,
    extra_valueB => 'S. cerevisiae',
    # C:  Gene.flankingRegions.direction
    opC    => '=',
    valueC => 'both',
    # A:  Gene.flankingRegions.distance
    opA    => '=',
    valueA => '1.0kb',
    # D:  Gene.flankingRegions.includeGene
    opD    => '=',
    valueD => 'true',
);

my $whole_sequence = '';


while (my $row = <$it>) {
 #   print $row->{'secondaryIdentifier'}, $row->{'symbol'}, $row->{'length'},
 #     $row->{'flankingRegions.direction'}, $row->{'flankingRegions.sequence.length'},
 #       $row->{'flankingRegions.sequence.residues'}, "\n";
	$whole_sequence = $row->{'flankingRegions.sequence.residues'};
	$stand_name = $row->{'symbol'}; 

}

#print "$stand_name\n\n";

# print "Whole sequence\n$whole_sequence\n";

# With the yeastmine seq returned we can get all the pertinent information

# Get the flanking sequences

my $upstream_seq = substr $whole_sequence, 0, 1000;

my $downstream_seq = substr $whole_sequence, -1000;

#print "$downstream_seq\n";
##############################################################

####################################
# Step 2: Get the F1 and R1 oligos #
####################################


# The F1 oligo includes 40bp immediately upstream of the ATG of the gene 
# followed by 20bp that are complementary to the Kan Cassette containing plasmid

# That means the 40bp can be subsetted from $upstream_seq

my $F1_oligo = substr $upstream_seq, -40;

$F1_oligo .= $F1_tail;

#print "$F1_oligo\n"; 

# The R1 oligo includes 40bp immediately downstream of the stop codon of the gene 
# followed by 20bp that are complementary to the Kan Cassette containing plasmid

my $R1_oligo = substr $downstream_seq, 0, 40;

#The R1 oligo requires the rev-comp before adding the tail

$R1_oligo = RevComp($R1_oligo);

$R1_oligo .= $R1_tail;

#print "$R1_oligo\n";


# The F2 oligo includes 40bp immediately upstream of the stop codon of the gene 
# followed by 20bp that are complementary to the Kan Cassette containing plasmid

# That means the 40bp can be subsetted from $upstream_seq

my $F2_oligo = substr $whole_sequence, -1043, 40;


$F2_oligo .= $F1_tail;

#print "$F2_oligo\n";



############################################################


#######################################
# Step 3: Design the .3 and .4 oligos #
#######################################



# The oligos for checking the deletion should be made like this 
# The .3 is upstream of the stop codon
# The .4 is downstream of stop codon
# They should be ~20bp long and end in two G/C (not more than two)
# They should be no more than 400bp away from stop codon

## Get +/-400bp from End site 
# The Stop Codon is 1000 from the end
# That means the 800bp start at -1400

my $sequence_stop  = substr $whole_sequence, -1400, 800;

#print "$check_sequence\n";

my $four_oligo = '';
my $three_oligo = '';
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


#print ".3 = $three_oligo\n.4 = $four_oligo\n";


##########################
## PCR product estimate ##
##########################

my $product_length = $four_location - $three_location;

my $three_from_stop = 400 - $three_location;
my $four_from_stop = $four_location - 400; 




#####################
##	Output	   ##
#####################


my %output;

#Fill Output
$output{'Systematic_Name'}=$gene_name;
$output{'Standard_Name'}=$stand_name;
#$output{'Location'}=$chr.':'.$start.'-'.$end;
#$output{'Strand'}=$strand;
$output{'F1'}=$F1_oligo;
$output{'F2'}=$F2_oligo;
$output{'R1'}=$R1_oligo;
$output{'.3'}=$three_oligo;
$output{'.4'}=$four_oligo;
#$output{'Location3'}='-'.$three_from_stop;
#$output{'Location4'}='+'.$four_from_stop;
$output{'PCR_length'}=$product_length;
$output{'Whole_seq'}=$whole_sequence;


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






