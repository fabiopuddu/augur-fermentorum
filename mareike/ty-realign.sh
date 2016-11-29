# /bin/sh

#sorting bam files based on read name
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

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
# cat bams_for_mpileup | while read line
# 	 do if [[ ! -a ${line}_sorted.bam ]]; 
# 	 		then sbatch  --wrap="samtools sort -n -T ${line}_temp -O BAM -o ${line}_sorted.bam ${line}";
# 	 		else echo "Skipping sort ${line}";
#      	fi
#      done	
# waitforcompletion
# cat bams_for_mpileup | while read line 	
# 	do 	percorso=`echo ${line} | sed "s|[A-Z_]*[0123456789b]*.[a-z]*.bam||g"`;
# 		name=`echo ${line} | tr '/' "\n" | grep "SC\|SD" | tail -n1`;
# 		echo "$name ... $percorso";
# 		if [[ ! -a ${line}.fq2 ]]; 
# 			then sbatch  --wrap="bamToFastq -i ${line}_sorted.bam -fq $percorso$name.fq1 -fq2 $percorso$name.fq2";
# 			else echo "Skipping bamToFastq ${line}";
# 		fi	
# 	done								
# waitforcompletion
#realigning
mkdir -p TR_BAMS
cat bams_for_mpileup | while read line
 	do 	percorso=`echo ${line} | sed "s|[A-Z_]*[0123456789b]*.[a-z]*.bam||g"`
		name=`echo ${line} | tr '/' "\n" | grep "SC\|SD" | tail -n1 | sed "s/\.bam//g"`		
		command1="bwa mem -R \"@RG\tID:$name\" $DIR/../mpileup_defaults/Ty_ref/Ty1-5.fa $percorso$name.fq1.gz $percorso$name.fq2.gz | samtools view -bS - | samtools sort -o TR_BAMS/$name.Ty.bam -O bam -T $name"
		command2="samtools index TR_BAMS/$name.Ty.bam"
		PROC=$(sbatch --wrap="${command1}" | sed 's/Submitted batch job //g') 
        sbatch --dependency=afterok:${PROC} --wrap="${command2}"
	done
waitforcompletion								   																												
