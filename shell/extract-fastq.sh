#!/bin/bash

submit_sbatch(){
sbatch --partition=LONG --output /dev/null --error /dev/null ${1} --wrap="${2}" 
}

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


proclist=''
while read -r line
 do              fl=`echo $line | sed "s/.bam//g"`
                echo $fl
                 command1="samtools sort -n -T ${fl}_temp.bam -O BAM -o ${fl}_sorted.bam ${fl}.bam"
                 command2="bamToFastq -i ${fl}_sorted.bam -fq ${fl}.fq1 -fq2 ${fl}.fq2"
                 command3="gzip -9 -f ${fl}.fq1 ${fl}.fq2"
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
