#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Term::ANSIColor;
use Getopt::Long;

#################################################
#                                               #
#  CALCULATING/DISPLAYING GENOTYPE TABLE        #
#                                               #
#################################################

##convert output into genotypes
#create genotype tables with both common names and systematic names; add 'NONAME' whenever a non standard gene name is found

my ($s);
my ($input);
GetOptions
(
's|show_synonymous'   => \$s,
'i|input=s'	  => \$input,
);
( $input && -f $input ) or die qq[Usage: $0 \n
					 	-i <input file>\n
					 	-s <show synonymous mutations (default: no)> \n
						];

my $rad5='YLR032W:missense_variant:1603:535:G>R';  #check for the presence of rad5-535
my $rad50S='YNL250W:missense_variant:242:81:K>I';  #check for the presence of rad50S mutation
my $sae2_F276S= 'YGL175C:missense_variant:827:276:F>S';  #check for the presence of sae2-F276S mutation

my %mutations;

open (my $mutation_file, '<', 'hom.table.file') or die "cannot open file";
chomp(my @mutation_list = <$mutation_file>);
close ($mutation_file);
my @files = <ERS*.vcf.gz>;
foreach my $file (@files) {
	my $sample_name = $file =~ s/.vcf.gz//r;
	my @genotype;
	foreach my $line (@mutation_list){
		next if $line !~ m/$sample_name/;
		(undef, undef, undef, my $mut, my $system_name, my $common_name,undef,undef,undef,undef)=split("\t", $line);
		if ((defined $common_name) and (not $common_name eq '""') ){
			push (@genotype, $common_name."-".$mut);
		}
		else{
			push (@genotype, $system_name."-".$mut);
		}
	}
	$mutations{$sample_name}=\@genotype;
}
#print Dumper(\%mutations);
print "Differential Genotype";
print "\n=========================================================================================================================================================\n";
foreach my $gt (sort keys %mutations){
	print "$gt\t";
	for my $mut(@{$mutations{$gt}}){
		if (($mut =~ m/FS@/) or ($mut =~ m/£Δ/)){
		$mut=~ s/£//g;
		print color("red"), "$mut\t", color("reset");
		}
		elsif ($mut =~ m/£x/) {
		$mut=~ s/£//g;
		$mut=~ s/x//;
		print color("yellow"), "$mut\t", color("reset");
		}
		elsif ($mut =~ m/>/) {
		$mut=~ s/£//g;
		print color("blue"), "$mut\t", color("reset") if ($s);
		}
		elsif (($mut =~ m/£II/) or ($mut =~ m/£ID/)){
		$mut=~ s/£//g;
		print color("green"), "$mut\t", color("reset");
		}
		else{
		print "error!\t"
		}
	}
	print "\n"
}





#shell code follows


#                             then printf "\e[0m$tab\t" | tr "\n" "\t"
#                     fi
#                     if [[ $c -gt 3 ]]
#                         then    gen="`sed -n "$l"p mutation.file | cut -f$column | tr "\n" "\t" | tr -d "\t"`"
#                                 if [[ "$tab" == 'NONAME' ]] #print in colors: red for frameshift and stop gained, yellow point mutations, gray synonymous
#                                     then    sist=`sed -n "$l"p systematic_name.file | cut -f$column | tr "\n" "\t" | tr -d "\t"`
#                                             if [[ $gen =~ 'Δ' ]] || [[ $gen =~ 'FS@' ]]
#                                                 then printf "\e[32m$sist-$gen\t"
#                                             elif [[ $gen =~ ':' ]]
#                                                 then if [[ "$show_syno" == '1' ]]    
#                                                         then printf "\e[34m$sist-$gen\t"
#                                                      fi    
#                                             elif [[ $gen =~ 'x' ]]
#                                                 then PMCID_NUM='0' #PMCID_NUM=`scraper.py ${sist} ${gen:1:-1} | wc -l` #${gen::-1} returns the genotipe minus the last character
#                                                 	 printf "\e[33m$sist-$gen(${PMCID_NUM})\t" | tr -d 'x' ######################
#                                                 else printf "\e[97m$sist$gen\t"
#                                             fi
#                                     else    if [[ $gen =~ 'Δ' ]] || [[ $gen =~ 'FS@' ]]
#                                                 then printf "\e[32m$tab-$gen\t" | tr "\n" "\t"
#                                             elif [[ $gen =~ ':' ]]
#                                                 then if [[ "$show_syno" == '1' ]]    
#                                                          then printf "\e[34m$tab-$gen\t"
#                                                      fi
#                                             elif [[ $gen =~ 'x' ]]
#                                                 then PMCID_NUM='0' #PMCID_NUM=`scraper.py ${tab} ${gen:1:-1} | wc -l` #${gen::-1} returns the genotipe minus the last character
#                                                 	 printf "\e[33m$tab-$gen(${PMCID_NUM})\t" | tr -d 'x'
#                                                 else printf "\e[97m$tab-$gen\t" | tr "\n" "\t"
#                                             fi
#                                 fi
#                     fi
#                 #sleep 1
#                 done
#     done