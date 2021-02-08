##################################################
##FASTQ FILE DOWNLOAD AND MANAGEMENT
#These rules are for dealing with fastq files from external sources
#They will download, trim, align to standard reference and mark duplicates
#And store the data in a compact file in the CRAM_repository
##################################################
ruleorder:
        rename_CI_fastq>download_fastq

rule rename_CI_fastq:           #Rename CI fastq to equivalent Sanger RLT in main cram repository folder (SangerLanes 1-8, CI Lanes 9-16)
    input:
        fq1=lambda wildcards : CRAM_repository+"../SLX-"+wildcards.RLT.split("_")[1].split("#")[0]+"."+SD[bcode[wildcards.RLT]].replace("-","_")+".r_1.fq.gz" if "SLX" in platLib[bcode[wildcards.RLT]] else "",
        fq2=lambda wildcards : CRAM_repository+"../SLX-"+wildcards.RLT.split("_")[1].split("#")[0]+"."+SD[bcode[wildcards.RLT]].replace("-","_")+".r_2.fq.gz" if "SLX" in platLib[bcode[wildcards.RLT]] else ""
    output:
        fq1=temp(CRAM_repository+"../"+"{RLT}xfq1_raw.gz"),
        fq2=temp(CRAM_repository+"../"+"{RLT}xfq2_raw.gz")
    shell:
        "mv {input.fq1} {output.fq1} && mv {input.fq2} {output.fq2}"
##################################################
rule download_fastq:            #Fetch fastq files from url
  version: "1.00"
  input:
      nc=lambda wildcards : ancient(nameconversion) if "None" not in urls[wildcards.RLT] else ""
  output:
      temp(CRAM_repository+"../"+"{RLT}x{fq}_raw.gz"),
  run:
      fq_url=dict()
      u=urls[wildcards.RLT].split("%3B")
      if 'fastq' in u[0]:
          try:
              fq_url['fq1']=u[0]
              fq_url['fq2']=u[1]
          except:
              fq_url['fq']=u[0]
      command="curl ftp://anonymous@"+fq_url[wildcards.fq]+" > {output}"
      shell(command)
##################################################
rule fastq_quality_control_pre: #optional pre-trimming fastqc
  version: "1.00"
  input:
      CRAM_repository+"../"+"{RLT}x{fq}_raw.gz"
  output:
      "pre_trimming_fastqc/{RLT}.{fq}_fastqc.html"
  shell:
      "mkdir -p pre_trimming_fastqc && "
      "fastqc --outdir ./pre_trimming_fastqc -t 8 {input}"
##################################################
rule trim:                      #Trim fastq files with trim-galore
  version: "1.00"
  input:
      CRAM_repository+"../"+"{RLT}xfq1_raw.gz",
      CRAM_repository+"../"+"{RLT}xfq2_raw.gz"
  output:
      fq1=temp(CRAM_repository+"../"+"{RLT}.fq1.gz"),
      fq2=temp(CRAM_repository+"../"+"{RLT}.fq2.gz"),
      report1="trim_reports/{RLT}.fq1.report.txt",
      report2="trim_reports/{RLT}.fq2.report.txt"
  shell:
      # "ENC=`zcat {input} | head -n 10000 |\
      # awk '{{if(NR%4==0) printf(\"%s\",$0);}}' |  od -A n -t u1 | \
      # awk 'BEGIN{{min=100;max=0;}} {{for(i=1;i<=NF;i++) \
      #         {{if($i>max) max=$i; if($i<min) min=$i;}} }}END \
      #         {{if(max<=74 && min<59) print \"33\"; \
      #         else if(max>73 && min>=64) print \"64\"; \
      #         else if(min>=59 && min<64 && max>73) print \"Solexa\"; }}'` &&"
      # " trim_galore -j 4 --phred${{ENC}} {input} &&"
      "trim_galore --paired -j 4 --phred33 {input} &&"
      " mv {wildcards.RLT}xfq1_raw.gz_val_1.fq.gz {output.fq1} &&"
      " mv {wildcards.RLT}xfq1_raw.gz_trimming_report.txt {output.report1} &&"
      " mv {wildcards.RLT}xfq2_raw.gz_val_2.fq.gz {output.fq2} &&"
      " mv {wildcards.RLT}xfq2_raw.gz_trimming_report.txt {output.report2}"
##################################################
rule fastq_quality_control_post:    #Optional fastqc post trimming
  version: "1.00"
  input:
      CRAM_repository+"../"+"{RLT}.{fq}.gz"
  output:
      "post_trimming_fastqc/{RLT}.{fq}_fastqc.html"
  shell:
      "mkdir -p post_trimming_fastqc && "
      "fastqc --outdir ./post_trimming_fastqc -t 8 {input}"
##################################################
rule align:                         #Align the fastq data to standard_genome_ref
  version: "1.00"
  input:
      fq1=CRAM_repository+"../"+"{RLT}.fq1.gz",
      fq2=CRAM_repository+"../"+"{RLT}.fq2.gz",
      ref=standard_genome_ref
  output:
      temp(CRAM_repository+"../"+"{RLT}.sam")
  params:
      rg=lambda wildcards : r"@RG\tID:"+wildcards.RLT+r"\tSM:"+bcode[wildcards.RLT],
  shell:
      "bwa mem -t 12 -T 0 -R '{params.rg}' {input.ref} {input.fq1} {input.fq2} > {output}"
#################################################
rule samTobam:                      #encode sam to bam
  version: "1.00"
  input:
      CRAM_repository+"../"+"{RLT}.sam"
  output:
      temp(CRAM_repository+"../"+"{RLT}.bam")
  params:
      rg=lambda wildcards : r"@RG\tID:"+bcode[wildcards.RLT]+r"\tSM:"+bcode[wildcards.RLT],
  shell:
      samtools+" view -@ 8 -b {input} > {output}"
#################################################
rule mark_duplicates_and_sort:      #MarkDuplicates and  sort (bamsormadup)
  version: "1.00"
  input:
      CRAM_repository+"../"+"{RLT}.bam"
  output:
      temp(CRAM_repository+"../"+"sort.mdup.{RLT}.bam")
  shell:
      samtools+" collate -O {input} {output}.collate.tmp |"
      "bamsormadup threads=12 SO=coordinate level=0 verbose=0 fixmate=1 adddupmarksupport=1 tmpfile={output}.tmp > {output}"
#################################################
rule bamTocram: #Otherwie create crams from aligned data and store in scratch repository
  version: "1.00"
  input:
      bam=lambda wildcards : CRAM_repository+"../"+"sort.mdup.{RLT}.bam" if "fastq" in urls[wildcards.RLT] or force_bam_conversion else "",
      ref=genome_ref,
  output:
      CRAM_repository+"{run}/{RLT}.cram"
  shell:
      samtools+" view -@8 -O CRAM -T {input.ref} -o {output} {input.bam}"
