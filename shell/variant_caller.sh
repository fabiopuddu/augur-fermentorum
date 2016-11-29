#!/bin/bash
OPTIND=1         # Reset in case getopts has been used previously in the shell.
force_rewrite=0; v=0; check_deletions=0;DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ); f=0
#################################
#                               #
#   GET COMMAND LINE OPTIONS    #
#                               #
#################################
while getopts "fvkh" opt
    do  case "$opt" in
                        h)  printf "############   HELP   ###############\nOPTIONS\n"
                            printf "\t-h\tThis Help\n"	
                            printf "\t-f\tForce Rewrite ***NOT IMPLEMENTED YET***\n"
                            printf "\t-v\tverbose: print out detailed analysis progression ***NOT IMPLEMENTED YET\n"
                            printf "\t-k\tcheck deletions in BAM files before calling mutations\n"
                            exit 0
                        ;;
                        f)
                        force_rewrite=1
                        ;;
                        v)
                        v=1
                        ;;
                        k)
                        check_deletions=1
                        ;;
        esac
    done
#######################
#######PROGRAM STARTS##
#############################
waitforcompletion(){
printf "Waiting for process to complete"
finito=`squeue -u fp305 | wc -l | tr -d "\t"`
while [[ $finito != '1' ]]    
    do  finito=`squeue -u fp305 | wc -l | tr -d "\t"`
        printf '.'
        sleep 5
    done  
printf "\n"    
}

###########################################################################################
###### GET LIST OF BAM FILES TO PROCESS AND INDEX THEM IF THE BAI FILES ARE MISSING########
###########################################################################################
#cleanup
if [[ -a bams_for_mpileup ]]
    then rm bams_for_mpileup
fi
#create a list of BAM files to analyse
ls BAM/*/*/*/*.bam >> bams_for_mpileup
ls BAM/*/*/*.bam >> bams_for_mpileup
ls BAM/*/*.bam >> bams_for_mpileup
ls BAM/*.bam >> bams_for_mpileup
sort bams_for_mpileup > sorted_bams_for_mpileup
cat sorted_bams_for_mpileup | grep -v 'sorted' > sieved_bams_for_mpileup #remove'sorted' files, which come from re-alignments for ty and rDNA analysis
#index files if the index does not exist
cat sieved_bams_for_mpileup | parallel -j20 --no-notice ' x={};
													 n=`echo {} | tr '/' "\n" | tail -n1 | sed 's/.bam//g'`; 
	    											 printf "\r\033[K  Indexing $n"
													 if [[ ! -a $x.bai || $force_rewrite == 1  ]] #if the file does not exist or if force_rewrite is called
                											then    samtools index $x
        											 fi'
printf "\n"
cat sieved_bams_for_mpileup | grep 'merged' > merged_bams_for_mpileup #generate a list of already merged BAM files
#remove from sieved bams for mpileup any reference to files that have already been merged
cat sieved_bams_for_mpileup | grep 'merged' | grep -o "SC_MFY.......\|SD....." | sort -u | while read line
	do cat sieved_bams_for_mpileup | sed "/$line/d" > new_bams_for_mpileup
	   mv new_bams_for_mpileup sieved_bams_for_mpileup
	done	
#combine the lists of unmerged samples with the list of already merged samples	
cat sieved_bams_for_mpileup merged_bams_for_mpileup > bams_for_mpileup	
rm sieved_bams_for_mpileup merged_bams_for_mpileup sorted_bams_for_mpileup	
#I don't understand what the following lines are for: should have commented
# cat bams_for_mpileup | while read line
#     do if echo $line | grep -q 'merged.bam'
#             echo $line >> new_bams_for_mpileup
#         fi
#     done
#####################################################################
##### MERGE MULTIPLE BAM FILES IF SAMPLES WERE SEQUENCED TWICE#######
#####################################################################
#temporarily commented these out as they are not relevant anymore
#sort bams_for_mpileup > sorted_bams_for_mpileup
#merge_bams_samtools.pl sorted_bams_for_mpileup 
#waitforcompletion
#rm bams_for_mpileup
#mv sorted_bams_for_mpileup bams_for_mpileup
#cat bams_for_mpileup | grep merged | parallel --no-notice 'samtools index {}' #re-index the merged files

#generate igv.txt
sed 's|^BAM|load "/Volumes/LuCia/BAM files/Yeast|' bams_for_mpileup > igv.txt 
sed 's|$|"|' igv.txt > igv_loader.txt
rm igv.txt
###################################################################################################################
#IF DELETION CHECK IS ENFORCED, CHECK THAT THE STRAIN IS DELETED IN THE GENE, AND EXCLUDE THE BAM FILE IF NOT
###################################################################################################################
if [[ $check_deletions == 1 ]] 
	then	touch bams_for_mpileup_filtered
			printf "Checking deletions\n"
			gene_name=`echo $PWD | tr '/' "\n" | grep Del | sed 's/Del[01234567.89]*_//g' | tr '[:lower:]' '[:upper:]'` 
				cat bams_for_mpileup | while read line
					do
						name=`echo $line | grep -o "SC_MFY.......\|SD......" | tr -d '.'`
						printf "Gene: $gene_name\tProcessing...$name\t" >>../deletion_check_log
						check=`detect_deletion_chr_region.pl $gene_name $line | tr "\t" "\n" | grep Deleted: | sed 's/Deleted://g'`
						#check=`detect_deletion_chr_region.pl $gene_name $line`
						printf "$check\n" >>../deletion_check_log
						if [[ $check == 'Y' || $check == 'Yp (p=partial)' ]] 
							then echo $line >> bams_for_mpileup_filtered
						fi	
					done		
			mv bams_for_mpileup_filtered bams_for_mpileup
fi
#################################
#####MUTATION CALLING STARTS#####
#################################
mkdir -p calling
cd calling
cat ../bams_for_mpileup |  while read line 
	do	n=$(echo $line |  tr '/' "\n" | grep "SD.*\|SC.*" | tail -n1 | sed "s/\.bam//g" | sed "s/\.merged//g")
		ers=`grep -w $n ../../name\ conversion.tsv | tr "\t" "\n" | tail -n 1`
		printf "\r\033[K  Processing $ers"
#		if [[ ! -a $m.vcf.gz ]] 
		 	 sbatch -n1 -N1 -p1604 --mem 100 --exclude="cb-node17" -o slurm.%N.%j.out.txt -e slurm.%N.%j.err.txt --wrap="samtools mpileup -f $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa -g -t DP,DV -C0 -pm3 -F0.2 -d10000 ../$line | bcftools call -vm -f GQ | bgzip -f > $ers.vcf.gz "  > /dev/null
		 	#I set the -C flag from 50 to 0 because when many artificially introduced mutantions fall within the same read, it downgraded mapping quality too much resulting in those 			mutations not being called.					
 #		fi
 	done
 printf "\nDone\n"	
#  Created by Fabio on 25/01/2015.
#
