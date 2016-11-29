#! /usr/bin/perl
import os;
my @aminoacids = ('', 'A', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'F', 'Y', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'G');
my $gene='MRE11'; 
my $residue_n=110;
my $residue=P;

my $sl; my @output;
print "Protein\tPos.\tREF\tALT\tTitle\tURL\n";
for my $aa (@aminoacids){
    my $search_string="DRS2 K876";
	$search_string = $search_string.$aa;
	print ("$search_string\n");
	@output=`python scholar.py --csv -c10 --cookie-file=cookies.txt --no-patents -A \"$search_string\" -s \"yeast cerevisiae\"`; 
	for my $line(@output){
		chomp ($line);
		my @result = split('\|',$line);
		print ("$gene\t$residue_n\t$residue\t$aa\t$result[0]\t$result[1]\n");
	}
	$sl=rand()*20+5;
	print "sleeping $sl\n";
	sleep($sl);
}