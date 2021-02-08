##################################################
##CRAM FILE DOWNLOAD AND MANAGEMENT
#These rules are for dealing with cram files from external sources.
#They will download missing cram files from the specified link and store the data in the CRAM_repository
#When required the data will be linked into the secondary_cram_repository or, if more than one lane is available
#data from the different lanes will be merged in a single cram file
##################################################
ruleorder:
    download_crams_from_datastore>download_crams>bamTocram
##################################################
rule download_crams:  #Download crams from internet to scratch repository if possibl
  #version: "1.00"
  input:
      nc=lambda wildcards: ancient(nameconversion) if ".cram" in urls[wildcards.RLT] else ""
  params:
      u=lambda wildcards:urls[wildcards.RLT],
  output:
      CRAM_repository+"{run}/{RLT}.cram"
  shell:
      "curl ftp://anonymous@{params.u} > {output}"
#################################################
rule download_crams_from_datastore:
    resources:
      load=25
    input:
        CRAM_repository+"{run}/{RLT}.cram.datastore"
    output:
        temp(CRAM_repository+"{run}/{RLT}.cram")
    shell:
        "cp -L --preserve=all {input} {output}"
#################################################
rule link_crams:
  version: "1.00"
  input:
      crams=lambda wildcards: set(expand(CRAM_repository+"{run}/{RLT}.cram", zip, \
                      RLT=list(fln_for_merge[wildcards.sample]), \
                      run=map(lambda x:x.split("_")[0], list(fln_for_merge[wildcards.sample]))\
                      )) if len(list(fln_for_merge[wildcards.sample]))==1 else ""
  output:
      temp(mergedCRAM_repository+"{hash}/{sample}_merged.cram")
  shell:
      "ln -s {input} {output};"
#################################################
rule merge_crams:        #Merge BAMS with the same TubeBarcodeID, in scratch repository
  version: "1.00"
  input:
      crams=lambda wildcards: set(expand(CRAM_repository+"{run}/{RLT}.cram", zip, \
                      RLT=list(fln_for_merge[wildcards.sample]), \
                      run=map(lambda x:x.split("_")[0], list(fln_for_merge[wildcards.sample]))\
                      )) if len(list(fln_for_merge[wildcards.sample]))>1 else ""
  output:
      temp(mergedCRAM_repository+"{sample}_merged.temp.cram")
  shell:
      samtools+" merge -@ 12 --output-fmt cram {output} {input}"
#################################################
rule rehead_crams:        #Merge BAMS with the same TubeBarcodeID, in scratch repository
  version: "1.00"
  input:
      mergedCRAM_repository+"{sample}_merged.temp.cram"
  output:
      cr=mergedCRAM_repository+"{hash}/{sample}_merged.cram",
      hd=temp(mergedCRAM_repository+"{hash}/{sample}_merged.cram.header")
  shell:
      samtools+" view -@8 -H {input} > {output.hd} && "+
      samtools+" reheader {output.hd} {input} > {output.cr}"
#################################################
if ORGANISM == "yeast":
    rule extract_fastq:             #Extract fastq files from bam files
      version: "1.00"
      input:
          #the crams line is required because when merged cram files are a simple symlink to the main cram repo,
          # the files copied from datastore (that the links point to )should not be deleted before fastq have been extracted
          merg=mergedCRAM_repository+"{hash}/{sample}_merged.cram",
          crams=lambda wildcards: set(expand(CRAM_repository+"{run}/{RLT}.cram", zip, \
                          RLT=list(fln_for_merge[wildcards.sample]), \
                          run=map(lambda x:x.split("_")[0], list(fln_for_merge[wildcards.sample]))\
                          )) if len(list(fln_for_merge[wildcards.sample]))==1 else mergedCRAM_repository+"{hash}/{sample}_merged.cram"
      output:
          temp(workingCRAM_repository+"{hash}/{sample}.fq1.gz"),
          temp(workingCRAM_repository+"{hash}/{sample}.fq2.gz")
      params:
          fq1=workingCRAM_repository+"{hash}/{sample}.fq1",
          fq2=workingCRAM_repository+"{hash}/{sample}.fq2"
      shell:
          samtools+" collate -@ 8 -O {input.merg} | "+
          samtools+" fastq -1 {params.fq1} -2 {params.fq2} - &&"
          "pigz -9 -f {params.fq1} {params.fq2}"
    ##################################################
    rule realign_to_working_genome:
        version: "1.00"
        input:
            rep_dna_reference=genome_ref,
            fq1=workingCRAM_repository+"{hash}/{sample}.fq1.gz",
            fq2=workingCRAM_repository+"{hash}/{sample}.fq2.gz"
        output:
            workingCRAM_repository+"{hash}/{sample}.cram"
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}",
            tempname="{sample}",
            ref=genome_ref
        shell:
            "bwa mem -t 12 -T0 -R '{params.rg}' {input.rep_dna_reference} {input.fq1} {input.fq2} | " +
            samtools+" view -@ 8 -b - | "+
            samtools+" collate -O - tmp.{params.tempname} | "+
            "bamsormadup threads=16 SO=coordinate level=9 verbose=0 fixmate=1 tmpfile={params.tempname}.tmp |"+
            samtools+" view -@8 -O CRAM -T {input.rep_dna_reference} -o {output} -"
##################################################
else:
    rule link_crams_to_working_repo:
        version: "1.00"
        input:
            mergedCRAM_repository+"{hash}/{sample}_merged.cram"
        output:
            workingCRAM_repository+"{hash}/{sample}.cram"
        shell:
            "ln -s {input} {output}"

#################################################
rule link_bams:         #Link the BAMS to current analysis folder
    version: "1.00"
    input:
        lambda wildcards : workingCRAM_repository+wildcards.sample[0:7] +"/"+wildcards.sample+".cram"
    output:
        "{DelFolder}/BAM/{sample}.cram"
    shell:
        "ln -s {input} {output}"
##################################################
rule index_bams:        #Index the BAMS in current analysis folder
    version: "1.00"
    input:
        "{DelFolder}/BAM/{sample}.cram"
    output:
        "{DelFolder}/BAM/{sample}.cram.crai"
    shell:
        samtools+" index {input}"
