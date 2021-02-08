##################################################
rule telomeres:                 #Calculate telomere length from fastq files
    version: "1.01+"+str(os.path.getmtime(installation_path + "perl/telomeres.pl"))
    input:
        "{DelFolder}/BAM/{sample}.cram"
    output:
        "{DelFolder}/repDNA/{sample}.tel"
    run:
        for line in shell(samtools+" view -F 2048 {input} | head -n 10000 | cut -f 10 | perl -ne 'chomp;print length($_)."+'"\\n"' + "' | sort | uniq -c || if [[ $? -eq 141 ]]; then true; else exit $?; fi", iterable=True):
            line=line.rstrip()
            result = [s for s in line.split(' ') if s]
            min_read_length=float(result[1])
        min_read_length=min_read_length-(min_read_length*0.10)
        shell(samtools+" view -F 2048 {input} | cut -f10 | " + installation_path + "perl/telomeres.pl 5 "+str(int(min_read_length))+" > {output}")
##################################################
rule repetitive_dna_estimation: #Estimate the length of various repetitive DNA elements
    version: "1.00+"+str(os.path.getmtime(installation_path + "perl/new-rDNA-cnv_estimate.pl"))
    input:
        bam="{DelFolder}/BAM/{sample}.cram",
        idx="{DelFolder}/BAM/{sample}.cram.crai",
        ref=genome_ref
    output:
        "{DelFolder}/repDNA/{sample}.txt"
    params:
        p=lambda wildcards: ploidy[wildcards.sample]
    shell:
        installation_path+"perl/new-rDNA-cnv_estimate.pl -p {params.p} -r {input.ref} -i {input.bam} > {output}"

##################################################
rule repetitive_table:          #Collect repDNA and telomere data in a file
    version: "1.00+"+str(os.path.getmtime(installation_path + "perl/get_repetitive_table.pl"))
    input:
        repDNA_file=lambda wildcards: expand("{{DelFolder}}/repDNA/{sample}.txt", sample=bcIDs[wildcards.DelFolder]),
        telomere_file=lambda wildcards:expand("{{DelFolder}}/repDNA/{sample}.tel", sample=bcIDs[wildcards.DelFolder]),
        nc=nameconversion
    output:
        repfile="{DelFolder}/repDNA/results.txt",
        reptable="{DelFolder}/reports/repDNA_table.txt"
    shell:
        installation_path+"perl/get_repetitive_table.pl -i {input.nc} -o {output.repfile}> {output.reptable} "
