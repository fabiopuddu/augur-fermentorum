#!/usr/bin/env perl

# Author:       	Fabio Puddu
# Maintainer:   	Fabio Puddu
# Created: 	Apr 2017
# Description:	

use strict;
use warnings;
#This scripts runs through the deletion check log file and the name conversion file to determine
#How did undeleted and failed samples distributes across the different studies
#Two known bugs: WT samples and deletion of MF(ALPHA)1 and MF(ALPHA)2 are reported as not deleted but are good 
#at least in the samples before the propagation
open (my $rep, '<', 'rep.txt');
open (my $name_conversion, '<', 'name conversion.tsv');

chomp(my @RD = <$rep>);
chomp(my @NC = <$name_conversion>);

my @OUT;

my @plate_list=('C001','C002','C003','C004','C005','C006','C007','C008','C009','C010','C011','C012','C013','C014','C015','C016','C017','C018','C019','C020','C021','C022','C023','C024','C025','C026','C027','C028','C029','C030','C031','C032','C033','C034','C035','C036','C037','C038','C039','C040','C041','C042','C043','C044','C045','C046','C047','C048','C049','C050','C051','C052','C053','C054','C055','C056','C057','C058','C059','C060','C061','C062','C063','C064','C065','C066','C067','C068','C069','C070','C071','C072','C073','C074','C075','C076','C077','C078','C079','C080','C081','C082','C083','C084','C085','C086','C087','C088','C089','C090','C091','C092','C093','C094','C095','C096','C097'); 
my @SDnumbers;
foreach my $plate (@plate_list){
	@SDnumbers=();
	@OUT=();
	my @matching_samples= grep { $_ =~ $plate } @NC;
	foreach my $line (@matching_samples){
		my @linea=split "\t", $line;
		push @SDnumbers, $linea[4];
	}
 	foreach my $sample(@SDnumbers){
                my $matching_line= (grep { $_ =~ /\b$sample\b/ } @RD)[0];
                if (defined $matching_line and length $matching_line > 0){
                        (undef, my $rDNA, my $cup1, my $mito, my $ty1, undef, undef, undef, undef, undef, my $telo, undef, undef)= split "\t", $matching_line;
                        push @OUT, $rDNA;
                }
        }
        my $string=join "\t", @OUT;
        printf "$plate\t$string\n";
}
