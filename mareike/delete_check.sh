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
printf "$1\n"
printf "Waiting for jobs $1  to complete"
list="ciao${1}"
finito=`squeue -u fp305 | grep "${list}" |wc -l | tr -d "\t"`
while [[ $finito != '0' ]]    
    do  finito=`squeue -u fp305 | grep "${list}" |wc -l | tr -d "\t"`
        printf ".${finito}"
        sleep 5
    done  
printf "\n"    
}

###########################################################################################
###### GET LIST OF BAM FILES TO PROCESS AND INDEX THEM IF THE BAI FILES ARE MISSING########
###########################################################################################
###################################################################################################################
#IF DELETION CHECK IS ENFORCED, CHECK THAT THE STRAIN IS DELETED IN THE GENE, AND EXCLUDE THE BAM FILE IF NOT
###################################################################################################################
PROCLIST=''
ls -d Del* |  while read fol
do  echo $fol
    cd "$fol"
    printf "Checking deletions\n"
    gene_name=`echo $PWD | tr '/' "\n" | grep Del | sed 's/Del[01234567.89]*_//g' | tr '[:lower:]' '[:upper:]'` 
    proclist=''
    ls BAM/*.bam |  while read line
					do	#get the name of the sample
						name=`echo $line | grep -o "SC_MFY.......\|SD......" | tr -d '.'`
						#print the start of the line into a sample specific temp file
						printf "\n Gene: $gene_name\tProcessing...$name\t" >>../$name.temp
						#run the command to get the deletion checked and the information printed into the file 
						#command1="echo 'test' >> ../$name.temp"
						command1="detect_deletion_chr_region.pl $gene_name $line >> ../$name.temp"
						sbatch --partition=LONG --error slurm-%j.err --wrap="${command1}" 
					done
    cd ..
done







