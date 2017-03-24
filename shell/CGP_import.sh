#!/bin/bash
#This scripts imports the CGP bam files from the raw_from_CGP folder
#and creates an appropriate structure. It has to be run from the main "from_CGP" directory.
#It assumes the presence of a directory "raw_from_CGP" containing raw data downloaded from globus
#It also assumes the presence of a directory "Gpgkey" containing decryption keys for each raw file

## FUNCTIONS ##
#This function submits commands to SLURM workload manager using sbatch
#Argument 1 = the command you want to run
#Argument 2 = any additional argument you want to pass to sbatch
submit_sbatch(){
sbatch --partition=LONG --output /dev/null --error /dev/null ${1} --wrap="${2}" 
}

#This checks the slurm queue and sleeps until the queue is empty
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

##MAIN##
#cd raw_from_CGP
mkdir -p /mnt/scratch/jackson/fp305/from_CGP/CGP_temp  
printf "Decrypting...\n"

for line in *SD*.gpg
    do passfile=`echo $line | sed 's/.tar.gz.gpg/.gpgkey/g' `
       strain=`echo $line | sed 's/.tar.gz.gpg//g' |sed 's/^.*_//g' `
       passphrase=`cat ../Gpgkey/$passfile | cut -f2`
       command="gpg --passphrase $passphrase --decrypt $line > /mnt/scratch/jackson/fp305/data/CGP_temp/$strain.tar.gz"
       if [[ -d ../CGP_BAM/$strain ]]
               then printf "Skipping $strain..."
               else submit_sbatch " " "${command}" > /dev/null
       fi        
    done
    
waitforcompletion
cd /mnt/scratch/jackson/fp305/from_CGP/CGP_temp
printf "Expanding...\n"
for line in SD*.tar.gz
    do submit_sbatch " " "tar xvf $line"  > /dev/null
    done
waitforcompletion
mkdir  -p ../CGP_BAM
rm -rf SD*.tar.gz
printf "Moving in place...\n"
for line in  *SD*
    do strain=`echo $line | sed 's/.tar.gz.gpg//g' |sed 's/^.*_//g' `
    mkdir -p ../CGP_BAM/$strain
    mv $line/mapped_sample/CEREVISIAE_S288c_R64-1-1_genomic_$strain.dupmarked.bam  ../CGP_BAM/$strain/$strain.bam
    mv $line/mapped_sample/CEREVISIAE_S288c_R64-1-1_genomic_$strain.dupmarked.bam.bai  ../CGP_BAM/$strain/$strain.bam.bai
    done
cd ../CGP_BAM
ls -d SD* | while read fld
    do     if [[ -a ${fld}/${fld}.fq1.gz ]]
           then printf "Skipping ${fld}"
           else command1="samtools sort -n -T ${fld}/temp.bam -O BAM -o ${fld}/${fld}_sorted.bam ${fld}/${fld}.bam"
                command2="bamToFastq -i ${fld}/${fld}_sorted.bam -fq ${fld}/${fld}.fq1 -fq2 ${fld}/${fld}.fq2"
                command3="gzip -9 -f ${fld}/${fld}.fq1 ${fld}/${fld}.fq2"
                command4="rm ${fld}/${fld}_sorted.bam"
                PROC1=$(submit_sbatch " " "${command1}" | sed 's/Submitted batch job //g') 
                PROC2=$(submit_sbatch "--dependency=afterok:${PROC1}" "${command2}"  | sed 's/Submitted batch job //g')
                PROC3=$(submit_sbatch "--dependency=afterok:${PROC2}" "${command3}"  | sed 's/Submitted batch job //g')
                submit_sbatch "--dependency=afterok:${PROC3}" "${command4}"    
        fi        
    done
    
    
        
