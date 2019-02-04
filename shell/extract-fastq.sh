#!/bin/bash

#submit_sbatch(){
#export SBATCH_CMD="$2"
#echo $SBATCH_CMD
#sbatch {$1} submit_sbatch_cramtobam.sh | sed 's/Submitted batch job //g'
#export SBATCH_CMD=""
#}

waitforcompletion(){
if [[ ${v} -eq '1' ]]; then printf "Waiting for jobs $1  to complete"; else printf "Waiting for jobs to complete"; fi
list="2146328${1}"
finito=`squeue | grep "${list}" |wc -l | tr -d "\t"`
while [[ $finito != '0' ]]    
    do  finito=`squeue  | grep "${list}" |wc -l | tr -d "\t"`
        if [[ ${v} -eq '1' ]]; then printf ".${finito}"; else printf "."; fi
        sleep 5
    done  
printf "\n"    
}

if ls *.cram 1> /dev/null 2>&1
then 
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

#cramTObam preprocessing

proclist=''
while read -r f
	do 	fname=`echo $f | sed "s/.cram//g"`
		export SBATCH_CMD="samtools view  -@ 8 -b -o $fname.bam -T $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa  $f"
		PROC1=$(sbatch submit_sbatch_cramtobam.sh | sed 's/Submitted batch job //g')
		echo $SBATCH_CMD
		export SBATCH_CMD=""
		#$command
		proclist="${proclist}\|${PROC1}"
		sleep 0.5
	 done  < <(ls *.cram)

waitforcompletion "${proclist}"
fi 


#extraction
proclist=''
while read -r line
 do              fl=`echo $line | sed "s/.bam//g"`
                echo $fl
                 export SBATCH_CMD="samtools sort -m 500M -@ 8 -n -T ${fl}_temp.bam -O BAM -o ${fl}_sorted.bam ${fl}.bam"
                PROC1=$(sbatch submit_sbatch_cramtobam.sh | sed 's/Submitted batch job //g')
		 export SBATCH_CMD=""
		export SBATCH_CMD="bamToFastq -i ${fl}_sorted.bam -fq ${fl}.fq1 -fq2 ${fl}.fq2"
                PROC2=$(sbatch --dependency=afterok:${PROC1} submit_sbatch_cramtobam.sh | sed 's/Submitted batch job //g')
		export SBATCH_CMD="" 
		export SBATCH_CMD="pigz -9 -f ${fl}.fq1 ${fl}.fq2"
                PROC3=$(sbatch --dependency=afterok:${PROC2} submit_sbatch_cramtobam.sh| sed 's/Submitted batch job //g')
		export SBATCH_CMD=""
		export SBATCH_CMD="rm ${fl}_sorted.bam"
                PROC4=$(sbatch --dependency=afterok:${PROC3} submit_sbatch_cramtobam.sh | sed 's/Submitted batch job //g')
		export SBATCH_CMD="" 
         	proclist="${proclist}\|${PROC4}"
		 sleep 1
 done  < <(ls *.bam)

waitforcompletion "${proclist}"


for f in *.bam; 
do 	name=`echo "${f}" | sed 's/.bam//g'` ;
	mkdir "${name}"; 
	mv ${name}.* ${name};
done

for f in */*.bam; 
	do sbatch --wrap "samtools index $f" ; 
	done
