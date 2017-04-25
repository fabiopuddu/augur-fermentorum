submit_sbatch(){
sbatch --partition=LONG --output /dev/null --error /dev/null ${1} --wrap="${2}" 
}
ls *.bam | while read line
 do              fl=`echo $line | sed "s/.bam//g"`
                 command1="samtools sort -n -T ${fl}_temp.bam -O BAM -o ${fl}_sorted.bam ${fl}.bam"
                 command2="bamToFastq -i ${fl}_sorted.bam -fq ${fl}.fq1 -fq2 ${fl}.fq2"
                 command3="gzip -9 -f ${fl}.fq1 ${fl}.fq2"
                 command4="rm ${fl}_sorted.bam"
                 PROC1=$(submit_sbatch " " "${command1}" | sed 's/Submitted batch job //g') 
                 PROC2=$(submit_sbatch "--dependency=afterok:${PROC1}" "${command2}"  | sed 's/Submitted batch job //g')
                 PROC3=$(submit_sbatch "--dependency=afterok:${PROC2}" "${command3}"  | sed 's/Submitted batch job //g')
                 submit_sbatch "--dependency=afterok:${PROC3}" "${command4}"  
 done
