#!/bin/bash
#Declare function that checks if the current jobs are done
waitforcompletion(){
if [[ ${v} -eq '1' ]]; then printf "Waiting for jobs $1  to complete"; else printf "Waiting for jobs to complete"; fi
list="ciao${1}"
finito=`squeue | grep "${list}" |wc -l | tr -d "\t"`
while [[ $finito != '0' ]]    
    do  finito=`squeue  | grep "${list}" |wc -l | tr -d "\t"`
        if [[ ${v} -eq '1' ]]; then printf ".${finito}"; else printf "."; fi
        sleep 5
    done  
printf "\n"    
}

#Declare function to submit jobs
submit_sbatch(){
sbatch --partition=LONG --output /dev/null --error %j.err ${1} --wrap="${2}" 
}



DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) #get the directory this program is stored in
proclist=''
while read -r line
	do 
 	ERSNO=`cat ../experiments/name\ conversion.tsv | grep $line | head -n1 | cut -f7`
	SNAME=`cat ../experiments/name\ conversion.tsv | grep $line | head -n1 | cut -f4`
	command1="bwa mem -R \"@RG\tID:${ERSNO}\tSM:${SNAME}\" $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa $line.fq1.gz $line.fq2.gz | samtools view -bS - | samtools sort -o $line.bam  -O bam -T $line.temp"
	command2="samtools index $line.bam"
	echo $command
	PROC1=$(submit_sbatch "" "$command1" |  sed 's/Submitted batch job //g')
	PROC2=$(submit_sbatch "--dependency=afterok:${PROC1}" "$command2" |  sed 's/Submitted batch job //g')
	proclist="${proclist}\|${PROC2}"
	done  < <(cat ../experiments/name\ conversion.tsv | cut -f5 | sort | uniq)

waitforcompletion "$proclist"

#proclist=''
#while read -r line
#	do command="samtools index $f"
#	PROC1=$(submit_sbatch "" "$command" |  sed 's/Submitted batch job //g')
#	done < <(ls *.bam)
#waitforcompletion "$proclist"



cat ../experiments/name\ conversion.tsv | cut -f5 | sort | uniq | while read line; 
	do mkdir $line;
	mv $line.* $line; 
	done
