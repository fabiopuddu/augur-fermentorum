
#!/bin/bash


# Author:       Fabio Puddu  
# Maintainer:   Fabio Puddu
# Created:      Sep 2016
# Description:

#sorting bam files based on read name
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
if [[ ! -a bams_for_mpileup ]]
	then echo "Run the variant calling first!!!"
	exit 1
fi
mkdir -p TR_BAMS
cat bams_for_mpileup | while read line
 	do 	percorso=`echo ${line} | sed "s|/.*\.bam|/|g"`
		name=`echo ${line} | tr '/' "\n" | tail -n1 | sed "s/\.bam//g"`

		export SBATCH_CMD="bwa mem -R \"@RG\tID:$name\tSM:$name\" $DIR/../mpileup_defaults/repDNA_ref/repDNA.fa $percorso$name.fq1.gz $percorso$name.fq2.gz | samtools view -@ 8 -b -F 4 - | bamsormadup threads=16 SO=coordinate level=0 verbose=0 fixmate=1 adddupmarksupport=1 tmpfile=$name.tmp | samtools sort -@ 8 > TR_BAMS/$name.Ty.bam"
		PROC=$(sbatch submit_sbatch_ty_realign.sh | sed 's/Submitted batch job //g') 
		
		export SBATCH_CMD="samtools index TR_BAMS/$name.Ty.bam"
		PROC2=$(sbatch --dependency=afterok:${PROC} submit_sbatch_ty_realign.sh | sed 's/Submitted batch job //g') 
		
		export SBATCH_CMD=""

#		PROC=$(sbatch --partition=LONG --wrap="${command1}" | sed 's/Submitted batch job //g') 
#        	PROC2=$(sbatch --partition=LONG --dependency=afterok:${PROC} --wrap="${command2}" | sed 's/Submitted batch job //g')
#        	PROC3=$(sbatch --partition=LONG --dependency=afterok:${PROC2} --wrap="${command3}" | sed 's/Submitted batch job //g')
#        	PROC4=$(sbatch --partition=LONG --dependency=afterok:${PROC3} --wrap="${command4}" | sed 's/Submitted batch job //g')
#        	sbatch --partition=LONG --dependency=afterok:${PROC4} --wrap="${command5}"
	done
