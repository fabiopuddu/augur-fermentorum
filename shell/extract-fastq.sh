#!/bin/bash

submit_sbatch(){
sbatch --partition=LONG --output /dev/null  ${1} --wrap="${2}" 
}

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


DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

#cramTObam preprocessing

proclist=''
while read -r f
	do 	fname=`echo $f | sed "s/.cram//g"`
		command="samtools view  -b -o $fname.bam -T $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa  $f"
		echo $command
		PROC1=$(submit_sbatch " " "${command}" | sed 's/Submitted batch job //g')
		#$command
		#proclist="${proclist}\|${PROC1}"
		sleep 0.5
	 done  < <(ls *.cram)
waitforcompletion "${proclist}"

#extraction
proclist=''
while read -r line
 do              fl=`echo $line | sed "s/.bam//g"`
                echo $fl
                 command1="samtools sort -n -T ${fl}_temp.bam -O BAM -o ${fl}_sorted.bam ${fl}.bam"
                 command2="bamToFastq -i ${fl}_sorted.bam -fq ${fl}.fq1 -fq2 ${fl}.fq2"
                 command3="pigz -9 -f ${fl}.fq1 ${fl}.fq2"
                 command4="rm ${fl}_sorted.bam"
                 PROC1=$(submit_sbatch " " "${command1}" | sed 's/Submitted batch job //g')
                 PROC2=$(submit_sbatch "--dependency=afterok:${PROC1}" "${command2}"  | sed 's/Submitted batch job //g')
                 PROC3=$(submit_sbatch "--dependency=afterok:${PROC2}" "${command3}"  | sed 's/Submitted batch job //g')
                 PROC4=$(submit_sbatch "--dependency=afterok:${PROC3}" "${command4}"  | sed 's/Submitted batch job //g')
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
