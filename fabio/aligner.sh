#!/bin/bash
submit_sbatch(){
sbatch --partition=LONG --output /dev/null --error %j.err ${1} --wrap="${2}" 
}
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) #get the directory this program is stored in
cat ../experiments/name\ conversion.tsv | cut -f5 | sort | uniq | while read line
	do 
	command="bwa mem -R \"@RG\tID:$ERSNO\tSM:$SNAME\" $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa $line.fq1.gz $line.fq2.gz | samtools view -bS - | samtools sort -o $line.bam  -O bam -T $line.temp"
#	echo $command
	submit_sbatch "" "$command"
	done
