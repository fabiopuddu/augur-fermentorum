bwa mem -R "@RG\tID:sample name\tSM:SC_MFY/SD/ERS" /Applications/PF/mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa fastq/SRR800790_1.fastq.gz fastq/SRR800790_2.fastq.gz | samtools view -bS - | samtools sort -o SC_MFY000000050.bam  -O bam -T temp