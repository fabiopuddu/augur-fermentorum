#!/bin/bash

# written 22.09.2016 

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
		command1="bwa mem -R \"@RG\tID:$name\" $DIR/../mpileup_defaults/repDNA_ref/repDNA.fa $percorso$name.fq1.gz $percorso$name.fq2.gz | samtools view -bS - | samtools sort -o TR_BAMS/$name.Ty.bam -O bam -T $name-ty"
		command2="samtools index TR_BAMS/$name.Ty.bam"
		command3="samtools view -b -F 4 TR_BAMS/$name.Ty.bam > TR_BAMS/$name.Ty.bam.map.bam"
		command4="mv TR_BAMS/$name.Ty.bam.map.bam TR_BAMS/$name.Ty.bam"
		command5="samtools index TR_BAMS/$name.Ty.bam"
		PROC=$(sbatch --partition=LONG --wrap="${command1}" | sed 's/Submitted batch job //g') 
        PROC2=$(sbatch --partition=LONG --dependency=afterok:${PROC} --wrap="${command2}" | sed 's/Submitted batch job //g')
        PROC3=$(sbatch --partition=LONG --dependency=afterok:${PROC2} --wrap="${command3}" | sed 's/Submitted batch job //g')
        PROC4=$(sbatch --partition=LONG --dependency=afterok:${PROC3} --wrap="${command4}" | sed 's/Submitted batch job //g')
        sbatch --partition=LONG --dependency=afterok:${PROC4} --wrap="${command5}"
	done

mkdir -p TWOMICRON_BAMS
cat bams_for_mpileup | while read line
        do      percorso=`echo ${line} | sed "s|/.*\.bam|/|g"`
                name=`echo ${line} | tr '/' "\n" | tail -n1 | sed "s/\.bam//g"`
                command1="bwa mem -R \"@RG\tID:$name\" $DIR/../mpileup_defaults/2-micron_ref/2-micron.fa $percorso$name.fq1.gz $percorso$name.fq2.gz | samtools view -bS - | samtools sort -o TWOMICRON_BAMS/$name.2m.bam -O bam -T $name-2m"
                command2="samtools index TWOMICRON_BAMS/$name.2m.bam"
                command3="samtools view -b -F 4 TWOMICRON_BAMS/$name.2m.bam > TWOMICRON_BAMS/$name.2m.bam.map.bam"
                command4="mv TWOMICRON_BAMS/$name.2m.bam.map.bam TWOMICRON_BAMS/$name.2m.bam"
                command5="samtools index TWOMICRON_BAMS/$name.2m.bam"
                PROC=$(sbatch --partition=LONG --wrap="${command1}" | sed 's/Submitted batch job //g')
                PROC2=$(sbatch --partition=LONG --dependency=afterok:${PROC} --wrap="${command2}" | sed 's/Submitted batch job //g')
                PROC3=$(sbatch --partition=LONG --dependency=afterok:${PROC2} --wrap="${command3}" | sed 's/Submitted batch job //g')
                PROC4=$(sbatch --partition=LONG --dependency=afterok:${PROC3} --wrap="${command4}" | sed 's/Submitted batch job //g')
                sbatch --partition=LONG --dependency=afterok:${PROC4} --wrap="${command5}"
        done
