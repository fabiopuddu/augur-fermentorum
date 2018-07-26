#!/bin/bash
# Author:       Fabio Puddu  
# Maintainer:   Fabio Puddu
# Created:      Jan 2015
# Description:


OPTIND=1         # Reset in case getopts has been used previously in the shell.
force_rewrite=0; v=0; check_deletions=0;DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ); f=0;
used_caller=''
#################################
#                               #
#   GET COMMAND LINE OPTIONS    #
#                               #
#################################
while getopts "fvkhSF" opt
    do  case "$opt" in
                        h)  printf "############   HELP   ###############\nOPTIONS\n"
                            printf "\t-h\tThis Help\n"	
                            printf "\t-f\tForce Rewrite ***NOT IMPLEMENTED YET***\n"
                            printf "\t-v\tverbose: print out detailed analysis progression ***NOT IMPLEMENTED YET\n"
                            printf "\t-k\tcheck deletions in BAM files before calling mutations\n"
                            printf "\t-c\tVariant caller used:\n\t\t\tS: samtools mpileup (default)\n\t\t\tF: freebayes\n"
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
                        S)
                        used_caller="S"
                        ;;
                        F)
                        used_caller="F"
                        ;;
        esac
    done

## Default variant caller and check
if [ "$used_caller" == '' ]
    then used_caller="S"
    printf "samtools mpileup chosen as variant caller by default\n"
elif [ "$used_caller" == "S" ]
    then echo "samtools mpileup chosen as variant caller"
elif [ "$used_caller" == "F" ]
    then echo "freebayes chosen as variant caller"
fi

#######################
#######PROGRAM STARTS##
#############################

###########################################################################################
###### GET LIST OF BAM FILES TO PROCESS AND INDEX THEM IF THE BAI FILES ARE MISSING########
###########################################################################################
#cleanup
rm -f bams_for_mpileup
#create a list of BAM files to analyse
ls BAM/*.bam >> bams_for_mpileup
#index files if the index does not exist
#cat bams_for_mpileup | while read line;
#	do n=`echo $line | tr '/' "\n" | tail -n1 | sed 's/.bam//g'`; 
#	printf "\r\033[K  Indexing $n"
#	if [[ ! -a $line.bai || $force_rewrite == 1  ]] #if the file does not exist or if force_rewrite is called
#        	then sbatch  --wrap="samtools index $line"
#        fi
#printf "\n"
#generate igv.txt
#sed 's|^BAM|load "/Volumes/LuCia/BAM files/Yeast|' bams_for_mpileup > igv.txt 
#sed 's|$|"|' igv.txt > igv_loader.txt
#rm igv.txt
###################################################################################################################
#IF DELETION CHECK IS ENFORCED, CHECK THAT THE STRAIN IS DELETED IN THE GENE, AND EXCLUDE THE BAM FILE IF NOT
###################################################################################################################
if [[ $check_deletions == 1 ]] 
	then	touch bams_for_mpileup_filtered
			printf "Checking deletions\n"
			gene_name=`echo $PWD | tr '/' "\n" | grep Del | sed 's/Del[01234567.89]*_//g' | sed 's/_.*$//g' | tr '[:lower:]' '[:upper:]'` 
				cat bams_for_mpileup | while read line
					do
						name=`echo $line | grep -oE '\bSD[^ ]*\.|\bSC_MFY[^ ]*\.' | tr -d '.'`
						printf "Gene: $gene_name\tProcessing...$name\t" >>deletion_check_log
						check=`detect_deletion_chr_region.pl $gene_name $line | tr "\t" "\n" | grep Deleted: | sed 's/Deleted://g'`
						printf "$check\n" >>deletion_check_log
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
		ers=`grep -w $n ../../name\ conversion.tsv | awk '{print $6}'`
        	check_pool=`echo $n | grep 'pool'`
        if [ -z $check_pool ]
            	then ploidy=2
            	else ploidy=4
        fi
	printf "\r\033[K  Processing $ers\t"
	#What caller are we using?
	if [ "$used_caller" == "S" ]
		then export SBATCH_CMD="samtools mpileup -f $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa -g -t DP,DV -C0 -pm3 -F0.2 -d10000 ../$line | bcftools call -vm -f GQ | bgzip -f > $ers.vcf.gz"
	 elif [ "$used_caller" == "F" ]
                then
			# Set some parameters according to ploidy
			minfrac=`bc -l <<< "scale=1; 2/$ploidy/5"` # Bash doesn't do floating point arithmetic
        		minfrac=`echo "0$minfrac"`
			export SBATCH_CMD="freebayes -f $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa -p $ploidy -m0 -C3 -F $minfrac --max-coverage 10000 -E-1 -J -= $n.no2m.bam | bgzip -f > $ers.vcf.gz"
        fi
        ### Check if bam contains 2-micron reads
        hits2m=`samtools idxstats ../$line | cut -f1 | grep "2-micron" | wc -l` # Optimise? (stop checking if one 2-micron is found)
        if [ "$hits2m" -gt 0 ]
                then
                    #samtools view -h ../$line | awk '!($3 == "2-micron")' > $n.no2m.sam
                    command2="samtools view -Sb $n.no2m.sam | samtools sort -o $n.no2m.bam"
                    command4="rm $n.no2m.*"
                    PROC2=$(sbatch --partition=LONG --wrap="${command2}" | sed 's/Submitted batch job //g')
                    PROC3=$(sbatch --dependency=afterok:${PROC2} submit_sbatch_caller.sh | sed 's/Submitted batch job //g')
                    sbatch --partition=LONG --dependency=afterok:${PROC3} --wrap="${command4}"
		    export SBATCH_CMD=''

                else
		   sbatch  submit_sbatch_caller.sh
		   export SBATCH_CMD=''
            fi
# freebayes outputs to vcf instead of bcf (no equivalent to -g mpileup tag)
# freebayes can't choose just two tags for output (-t DP,AD/DV)
	sleep 0.3
#I set the -C flag from 50 to 0 because when many artificially introduced mutantions fall within the same read, it downgraded mapping quality too much resulting in those mutations not being called.					
 	done
printf "\nDone\n"	
