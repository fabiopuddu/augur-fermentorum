#!/usr/bin/env perl
use strict;
use warnings;

#################################################
#                                               #
#  CALCULATING/DISPLAYING GENOTYPE TABLE        #
#                                               #
#################################################

##convert output into genotypes
#create genotype tables with both common names and systematic names; add 'NONAME' whenever a non standard gene name is found

rad5='YLR032W:missense_variant:1603:535:G>R'  #check for the presence of rad5-535
rad50S='YNL250W:missense_variant:242:81:K>I'  #check for the presence of rad50S mutation
sae2-F276S= 'YGL175C:missense_variant:827:276:F>S'  #check for the presence of sae2-F276S mutation




@files = <ERS*.vcf.gz>;
foreach $file (@files) {
	open (my $input, '<', "gunzip -c $file | ") or die "cannot open $file";
	while( my $line = <$input>){
		next if $line =~ /^#/
		
	}
	close ($input)
  
  
  }

#shell code follows


for x in ERS*.vcf.gz
    do  n=$(echo $x | sed 's/.vcf.gz//g')
        samp_name=`cat "../../name conversion.tsv" | grep -w $n | cut -f2`;    
        num=$(echo $n | sed 's/ERS//g')
        geno_common=`grep $n hom.table.file |grep Y.[LR]........... -o |  tr "\t" ";" | sed 's/Y[ABCDEFGHIJKLMNOP][LR][0123456789][0123456789][0123456789][WC]-*[ABCD]*;//g' | sed 's/;.*$//' | sed 's/""/NONAME/g' | sed 's/-[ABCDEF]""/NONAME/g' | tr '\n' '\t'`
        geno_systematic=`grep $n hom.table.file | grep Y.[LR].......... -o | tr "\t" ";" |  sed -e 's/(//g' |sed -e 's/[[:space:]][A-Z][A-Z][A-Z][0123456789]*//g' | sed 's/""//g' | sed 's/;.*$//g' | tr '\n' '\t'`
        geno_details=`cat hom.table.file | grep $n | grep -o £.*£ | sed 's/£//g' | tr '\n' '\t'`
        if [[ "$control" =~ "$num" ]]
            then    if [[ $rad5 != "0" ]] #print rad5-535 if the mutation is present
                        then add_gen="\e[31mrad5-535\e[0m"
                        else add_gen=''
                    fi
                    if [[ $mut != "0" ]] #print the mutation to be checked, if the mutation is present
                        then add_gen="$add_gen\t\e[96m$mut_check\e[0m"
                        else add_gen="$add_gen"
                    fi
                    printf "$n\t$samp_name\t++++\t$add_gen\n" >> common_name.file
            else    printf "$n\t$samp_name\t.\t$geno_common\n" >> common_name.file
        fi
        printf "$n\t$samp_name\t$geno_systematic\n" >>systematic_name.file
        printf "$n\t$samp_name\t$geno_details\n" >>mutation.file
    done
# identify the samples with known resistant genes
top1=`cat hom.table.file | grep TOP1 | grep -o ERS...... | tr '\n' '\t'`
pdr=`cat hom.table.file | grep PDR | grep -o ERS...... | tr '\n' '\t'` #pleiotropic drug resistance
#zds=`cat hom.table.file | grep ZDS | grep -o ERS...... | tr '\n' '\t'` #zillion different screens
#print headers
printf "\nDIFFERENTIAL GENOTYPE\n"
if [[ $ploidy == "2" ]]
        then printf "ERS NUMBER\tSAMPLE NAME\tREF\tHOMOZYGOUS GENOTYPE"
        else printf "ERS NUMBER\tSAMPLE NAME\tREF\tGENOTYPE"
fi
printf "\n=========================================================================================================================================================\n"
l=0
r=0
#parse common name genotype table and substitute the corresponding systematic name when 'NONAME' is found
cat common_name.file | while read line
    do
        l=$(($l+1))
        c=0
        if [[ "$l" != '1' ]]
            then printf "\n"
        fi
        for tab in $line
                do
                    c=$(($c+1))
                    column=$(($c-1))

                    if [[ $c -eq 1 ]]
                            then    if [[ $top1 =~ $tab || $pdr =~ $tab || $zds =~ $tab ]] #if the sample contains mutations in known suppressors (top1, pdr1, pdr3,zds1,zds2) write the sample name in grey
                                    then    printf "\e[31m$l $tab\t\e[0m" | tr "\n" "\t"
                                    else    printf "\e[0m$l $tab\t" | tr "\n" "\t"
                                fi
                    fi
                    if [[ $c -eq 2 ]]
                            then printf "\e[0m$tab\t" | tr "\n" "\t"
                    fi
                    if [[ $c -eq 3 ]]
                            then printf "\e[0m$tab\t" | tr "\n" "\t"
                    fi
                    if [[ $c -gt 3 ]]
                        then    gen="`sed -n "$l"p mutation.file | cut -f$column | tr "\n" "\t" | tr -d "\t"`"
                                if [[ "$tab" == 'NONAME' ]] #print in colors: red for frameshift and stop gained, yellow point mutations, gray synonymous
                                    then    sist=`sed -n "$l"p systematic_name.file | cut -f$column | tr "\n" "\t" | tr -d "\t"`
                                            if [[ $gen =~ 'Δ' ]] || [[ $gen =~ 'FS@' ]]
                                                then printf "\e[32m$sist-$gen\t"
                                            elif [[ $gen =~ ':' ]]
                                                then if [[ "$show_syno" == '1' ]]    
                                                        then printf "\e[34m$sist-$gen\t"
                                                     fi    
                                            elif [[ $gen =~ 'x' ]]
                                                then PMCID_NUM='0' #PMCID_NUM=`scraper.py ${sist} ${gen:1:-1} | wc -l` #${gen::-1} returns the genotipe minus the last character
                                                	 printf "\e[33m$sist-$gen(${PMCID_NUM})\t" | tr -d 'x' ######################
                                                else printf "\e[97m$sist$gen\t"
                                            fi
                                    else    if [[ $gen =~ 'Δ' ]] || [[ $gen =~ 'FS@' ]]
                                                then printf "\e[32m$tab-$gen\t" | tr "\n" "\t"
                                            elif [[ $gen =~ ':' ]]
                                                then if [[ "$show_syno" == '1' ]]    
                                                         then printf "\e[34m$tab-$gen\t"
                                                     fi
                                            elif [[ $gen =~ 'x' ]]
                                                then PMCID_NUM='0' #PMCID_NUM=`scraper.py ${tab} ${gen:1:-1} | wc -l` #${gen::-1} returns the genotipe minus the last character
                                                	 printf "\e[33m$tab-$gen(${PMCID_NUM})\t" | tr -d 'x'
                                                else printf "\e[97m$tab-$gen\t" | tr "\n" "\t"
                                            fi
                                fi
                    fi
                #sleep 1
                done
    done