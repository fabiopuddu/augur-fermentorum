# Author:       Fabio Puddu  
# Maintainer:   Fabio Puddu
# Created:      Jun 2018
# Description:



DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) 
mkdir -p BAM_repository
cd BAM_repository
cat "../name conversion.tsv" | while read line
	do ers=`echo "$line" | cut -f 6`
	   echo Getting $ers ...
	   enaDataGet $ers;
	done
for f in *.cram
	do out=`echo $f | sed 's/.cram//g'`
	   export SBATCH_CMD="samtools view -T $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa -b -o $out.bam $f"
           PROC=$(sbatch ${DIR}/submit_sbatch_cramtobam.sh | sed 's/Submitted batch job //g')
	   export SBATCH_CMD="samtools index $out.bam"
	   sbatch --dependency=afterok:${PROC} ${DIR}/submit_sbatch_cramtobam.sh	
           export SBATCH_CMD=""
	done
cd ..

