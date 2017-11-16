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
finito=`squeue -u ia327 | wc -l | tr -d "\t"`
while [[ $finito != '1' ]]    
    do  finito=`squeue -u ia327 | wc -l | tr -d "\t"`
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
#ls BAM/*/*/*/*.bam >> bams_for_mpileup
#ls BAM/*/*/*.bam >> bams_for_mpileup
#	ls BAM/*/*.bam >> bams_for_mpileup
ls BAM/*.bam >> bams_for_mpileup
sort bams_for_mpileup > sorted_bams_for_mpileup
#cat sorted_bams_for_mpileup | grep -v 'sorted' > sieved_bams_for_mpileup #remove'sorted' files, which come from re-alignments for ty and rDNA analysis
#index files if the index does not exist
cat bams_for_mpileup | parallel -j20 --no-notice ' x={};
				n=`echo {} | tr '/' "\n" | tail -n1 | sed 's/.bam//g'`; 
	    			printf "\r\033[K  Indexing $n"
		 		if [[ ! -a $x.bai || $force_rewrite == 1  ]] #if the file does not exist or if force_rewrite is called
                				then    samtools index $x
        				fi'
printf "\n"
#cat sieved_bams_for_mpileup | grep 'merged' > merged_bams_for_mpileup #generate a list of already merged BAM files
#remove from sieved bams for mpileup any reference to files that have already been merged
#cat sieved_bams_for_mpileup | grep 'merged' | grep -o "SC_MFY.......\|SD....." | sort -u | while read line
#	do cat sieved_bams_for_mpileup | sed "/$line/d" > new_bams_for_mpileup
#	   mv new_bams_for_mpileup sieved_bams_for_mpileup
#	done	
#combine the lists of unmerged samples with the list of already merged samples	
#cat sieved_bams_for_mpileup merged_bams_for_mpileup > bams_for_mpileup	
#rm sieved_bams_for_mpileup merged_bams_for_mpileup sorted_bams_for_mpileup	
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
						name=`echo $line | grep -oE '\bSD[^ ]*\.|\bSC_MFY[^ ]*\.' | tr -d '.'`
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
cat ../bams_for_mpileup | while read line
	do	n=$(echo $line |  tr '/' "\n" | tail -n1 | sed "s/\.bam//g" | sed "s/\.merged//g")
		ers=`grep -w $n ../../name\ conversion.tsv | tr "\t" "\n" | tail -n 1 | sed -e 's/[[:space:]]*$//'`
        check_pool=`echo $n | grep 'pool'`
        if [ -z $check_pool ]
            then ploidy=2
            else ploidy=4
        fi
		printf "\r\033[K  Processing $ers\t"
#		if [[ ! -a $m.vcf.gz ]]

        # Set some parameters according to ploidy
        minfrac=`bc -l <<< "scale=1; 2/$ploidy/5"` # Bash doesn't do floating point arithmetic
        minfrac=`echo "0$minfrac"`

        ### Check if bam contains 2-micron reads
        hits2m=`samtools view ../$line | awk '$3 == "2-micron"' | wc -l` # Optimise? (stop checking if one 2-micron is found)
            if [ "$hits2m" -gt 0 ]
                then
                    samtools view -h ../$line | awk '!($3 == "2-micron")' > $n.no2m.sam
                    command2="samtools view -Sb $n.no2m.sam | samtools sort -o $n.no2m.bam"
                    command3="freebayes -f $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa -p $ploidy -m0 -C3 -F $minfrac --max-coverage 10000 -J -= $n.no2m.bam | bgzip -f > $ers.vcf.gz"
                    command4="rm $n.no2m.*"
                    PROC2=$(sbatch --partition=LONG --wrap="${command2}" | sed 's/Submitted batch job //g')
                    PROC3=$(sbatch --partition=LONG -o slurm.%N.%j.out.txt -e slurm.%N.%j.err.txt --dependency=afterok:${PROC2} --wrap="${command3}" | sed 's/Submitted batch job //g')
                    sbatch --partition=LONG --dependency=afterok:${PROC3} --wrap="${command4}"

                else
                    sbatch --partition=LONG -o slurm.%N.%j.out.txt -e slurm.%N.%j.err.txt --wrap="freebayes -f $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa -p $ploidy -m0 -C3 -F $minfrac --max-coverage 10000 -J -= ../$line | bgzip -f > $ers.vcf.gz"
            fi

# freebayes outputs to vcf instead of bcf (no equivalent to -g mpileup tag)
# freebayes can't choose just two tags for output (-t DP,AD/DV)

			sleep 1
		 	#I set the -C flag from 50 to 0 because when many artificially introduced mutantions fall within the same read, it downgraded mapping quality too much resulting in those 			mutations not being called.					
 #		fi
 	done
 printf "\nDone\n"	
#  Created by Fabio on 25/01/2015.
#
