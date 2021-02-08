rule basic_variant_caller:  # Call chromosomal variants chromosome by chromosome and store data in temporary files
  version: "1.00"
  input:
      ref=genome_ref,
      bam="{DelFolder}/BAM/{sample}.cram",
      idx="{DelFolder}/BAM/{sample}.cram.crai"
  output:
      temp("{DelFolder}/calling/{sample}/{chr}.{sample}.vcf.gz")
  params:
      chr=lambda wildcards:wildcards.chr
  shell:
      samtools+" mpileup -f {input.ref} -r {params.chr} -g -t DP,DV -C0 -p -m3 -F0.2 -d10000 {input.bam} |"+
      bcftools+" call -vm -f GQ | bgzip -@8 -f > {output}"
##################################################
rule freebayes_variant_caller:
  version: "1.00"
  input:
      ref=genome_ref,
      bam="{DelFolder}/BAM/{sample}.cram",
      idx="{DelFolder}/BAM/{sample}.cram.crai"
  output:
      temp("{DelFolder}/calling/freebayes/{sample}/{chr}.{sample}.vcf.gz")
  params:
      p=lambda wildcards: ploidy[wildcards.sample],
      mf=lambda wildcards: (2/int(ploidy[wildcards.sample]))/5,
      chr=lambda wildcards:wildcards.chr
  shell:
      "freebayes -f {input.ref} -p {params.p} -m0 -C3 -F {params.mf} -r {params.chr} --max-coverage 10000 -E-1 -J -= {input.bam} | bgzip -f > {output}"
##################################################
if variant_caller == "freebayes":
  version: "1.00"
  rule merge_chr_vcf:
      input:
          fby=expand("{{DelFolder}}/calling/freebayes/{{sample}}/{chr}.{{sample}}.vcf.gz", chr=chromosome_list),
      output:
          "{DelFolder}/calling/{sample}.vcf.gz"
      shell:
          bcftools+" concat {input.fby}  | bgzip > {output}"
else:
  version: "1.00"
  rule merge_chr_vcf:
      input:
          mpl=expand("{{DelFolder}}/calling/{{sample}}/{chr}.{{sample}}.vcf.gz", chr=chromosome_list),
      output:
          "{DelFolder}/calling/{sample}.vcf.gz"
      shell:
          bcftools+" concat {input.mpl}  | bgzip > {output}"
##################################################
rule rehead_vcf:
  version: "1.00"
  input:
      "{DelFolder}/calling/{sample}.vcf.gz",
  output:
      "{DelFolder}/calling/rehead.{sample}.vcf.gz",
  params:
      header=lambda wildcards: "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+wildcards.sample
  shell:
      "gunzip -dc {input} | sed \"s/^.*CHROM.*$/{params.header}/\" | sed '/##samtoolsCommand/d' | bgzip > {output}"
##################################################
#If we have a normal panel we need to subtract the mutations present in the normal panel from the mutations present in the samples
if normal_panel:
  rule subtract_normal_panel:
      version: "1.00"
      input:
          v="{DelFolder}/calling/rehead.{sample}.vcf.gz",
          n=normal_panel
      output:
          v="{DelFolder}/calling/np.rehead.{sample}.vcf.gz",
      shell:
          "tabix -f -p vcf {input.v} && "+
          bcftools+" isec -w1 -C {input.v} {input.n} | bgzip > {output.v} "
          # "tabix -f -p vcf {output.v}"
else:
  rule link_normal_panel:
      version: "1.00"
      input:
          v="{DelFolder}/calling/rehead.{sample}.vcf.gz",
      output:
          v="{DelFolder}/calling/np.rehead.{sample}.vcf.gz",
      shell:
          "ln -s ../../{input.v} {output.v}"
##################################################
rule variant_effect_predictor:
  version: "1.00"
  input:
      vcf="{DelFolder}/calling/np.rehead.{sample}.vcf.gz",
  output:
      vep="{DelFolder}/analysis/vep.{sample}.txt",
      vcf="{DelFolder}/analysis/{sample}.vcf.gz"
  params:
      #filename=lambda wildcards: wildcards.sample,
      sp=species
  shell:
      "vep --species {params.sp} -i {input.vcf} --format vcf -o {output.vep} --force_overwrite --offline &&"
      "vcf2consequences_vep -v {input.vcf} -i {output.vep} | bgzip -@8 > {output.vcf}"
##################################################
rule prepare_vcf:
  version: "1.00"
  input:
      vcf="{DelFolder}/analysis/{sample}.vcf.gz",
      ref=genome_ref,
  output:
      v="{DelFolder}/analysis/sort.norm.{sample}.vcf.gz",
      i="{DelFolder}/analysis/sort.norm.{sample}.vcf.gz.tbi",
  params:
      filename=lambda wildcards: wildcards.sample,
  shell:
      #no need anymore to remove mitochondrial mutations because we are calling mutations by chromosome and excluded mito
      bcftools+" norm -m-both --check-ref e -f {input.ref} {input.vcf} |"+    #normalise indels
      bcftools+" sort | bgzip -@8 > {output.v} &&"      #sort the vcf
      "tabix -f -p vcf {output.v}"
##################################################
###RULE TO CREATE OR COMBINE CONTROLS
##################################################
rule create_control:
  #Rule logic: if we have a control inside the DelName folder use it (=if there is one copy ; more than one merge them )
  #            if we do not have a control inside the DelName folder then look up a special folder called Del0000_control and use it; otherwise fail
  version: "1.01"
  input:
      v=lambda wildcards: expand("{{DelFolder}}/analysis/sort.norm.{c}.vcf.gz",    c=set(control[wildcards.DelFolder])) if wildcards.DelFolder in control.keys() and len(set(control[wildcards.DelFolder])) >= 1  \
                                                                                                                                        else expand("Del0000_WildtypeControl/analysis/sort.norm.{c}.vcf.gz",    c=set(control["Del0000_WildtypeControl"])),
      t=lambda wildcards: expand("{{DelFolder}}/analysis/sort.norm.{c}.vcf.gz.tbi",c=set(control[wildcards.DelFolder])) if wildcards.DelFolder in control.keys() and len(set(control[wildcards.DelFolder])) >= 1  \
                                                                                                                                        else expand("Del0000_WildtypeControl/analysis/sort.norm.{c}.vcf.gz.tbi",c=set(control["Del0000_WildtypeControl"]))
  output:
      v="{DelFolder}/analysis/sort.norm.control.vcf.gz",
      t="{DelFolder}/analysis/sort.norm.control.vcf.gz.tbi"
  run:
      print ("Delfolder",wildcards.DelFolder,"control keys",control.keys())
      if wildcards.DelFolder in control.keys() and len(set(control[wildcards.DelFolder])) == 1:
          shell("cp -v {input.v} {output.v} && cp {input.t} {output.t}")
      elif wildcards.DelFolder in control.keys() and len(set(control[wildcards.DelFolder])) > 1:
          shell(bcftools+" merge -m none {input.v} | bgzip > {output.v} && tabix -f -p vcf {output.v}") #merge control vcfs and split multiallelic sites on multiple lines
      elif wildcards.DelFolder not in control.keys():
          if "Del0000_WildtypeControl" in control.keys():
              if len(set(control["Del0000_WildtypeControl"])) == 1:
                  shell("cp -v {input.v} {output.v} && cp {input.t} {output.t}")
              elif len(set(control["Del0000_WildtypeControl"])) > 1:
                  shell(bcftools+" merge -m none {input.v} | bgzip > {output.v} && tabix -f -p vcf {output.v}") #merge control vcfs and split multiallelic sites on multiple lines

##################################################
rule pre_filtering:
  version: "1.00"
  input:
      s="{DelFolder}/analysis/sort.norm.{sample}.vcf.gz",
  output:
      s="{DelFolder}/analysis/filt.sort.norm.{sample}.vcf.gz",
      t="{DelFolder}/analysis/filt.sort.norm.{sample}.vcf.gz.tbi"
  shell:
      bcftools+" filter -m + -s 'minDP' -e 'INFO/DP<10' {input.s}| "+
      bcftools+" filter -m + -s 'minDV' -e 'FORMAT/DV<3' | "+
      bcftools+" filter -m + -s 'minQsnp' -e 'TYPE=\"snp\" & QUAL<100' | "+
      bcftools+" filter -m + -s 'minQind' -e 'TYPE=\"indel\" & QUAL<30' | "+
      bcftools+" filter -m + -s 'minGQ' -e 'FORMAT/GQ<40' | "+
      bcftools+" filter -m + -s 'SNPgap' -g 7 | bgzip -@8 > {output.s} &&"
      "tabix -f -p vcf {output.s}"
##################################################
rule filter_control:
  version: "1.00"
  input:
      c="{DelFolder}/analysis/sort.norm.control.vcf.gz"
  output:
      c="{DelFolder}/analysis/filt.sort.norm.control.vcf.gz",
      t="{DelFolder}/analysis/filt.sort.norm.control.vcf.gz.tbi"
  shell:
      bcftools+" filter -m + -s 'minDP' -e 'INFO/DP<10' {input.c}| "+
      bcftools+" filter -m + -s 'minDV' -e 'FORMAT/DV<3' | "+
      bcftools+" filter -m + -s 'minQsnp' -e 'TYPE=\"snp\" & QUAL<100' | "+
      bcftools+" filter -m + -s 'minQind' -e 'TYPE=\"indel\" & QUAL<30' | "+
      bcftools+" filter -m + -s 'minGQ' -e 'FORMAT/GQ<40' | "+
      bcftools+" filter -m + -s 'SNPgap' -g 7 | bgzip -@8 > {output.c} &&"
      "tabix -f -p vcf {output.c}"
##################################################
rule intersections:
  version: "1.00"
  input:
      s="{DelFolder}/analysis/filt.sort.norm.{sample}.vcf.gz",
      ts="{DelFolder}/analysis/filt.sort.norm.{sample}.vcf.gz.tbi",
      c="{DelFolder}/analysis/sort.norm.control.vcf.gz",
      tc="{DelFolder}/analysis/sort.norm.control.vcf.gz.tbi"
  output:
      v="{DelFolder}/analysis/isec.{sample}.vcf.gz",
      t="{DelFolder}/analysis/isec.{sample}.vcf.gz.tbi",
  shell:
      bcftools+" isec -w1 -C {input.s} {input.c} | bgzip -@8 > {output.v} && tabix -f -p vcf {output.v}"
##################################################
rule inverse_intersections:
  version: "1.00"
  input:
      s="{DelFolder}/analysis/sort.norm.{sample}.vcf.gz",
      c="{DelFolder}/analysis/filt.sort.norm.control.vcf.gz"
  output:
      v="{DelFolder}/analysis/inv_isec.{sample}.vcf.gz",
      t="{DelFolder}/analysis/inv_isec.{sample}.vcf.gz.tbi",
  shell:
      bcftools+" isec -w1 -C {input.c} {input.s} | bgzip -@8 > {output.v} && tabix -f -p vcf {output.v}"
      #"bedtools intersect -header -a {input.c} -b {input.s} -v | bgzip -@8 > {output}"
##################################################
rule merge_experiment_vcf:
  version: "1.00"
  input:
      data=lambda wildcards: expand("{{DelFolder}}/analysis/isec.{id}.vcf.gz", id=bcIDs[wildcards.DelFolder])if len(bcIDs[wildcards.DelFolder])>1 else "",
      idx=lambda wildcards: expand("{{DelFolder}}/analysis/isec.{id}.vcf.gz.tbi", id=bcIDs[wildcards.DelFolder]) if len(bcIDs[wildcards.DelFolder])>1 else ""
  output:
      "{DelFolder}/analysis/experiment_merge.vcf.gz"
  shell:
      bcftools+" merge -m none {input.data} | bgzip -@8 > {output}"
##################################################
rule link_experiment_vcf:  #in case there was only one sample
  version: "1.00"
  input:
      data=lambda wildcards: expand("{{DelFolder}}/analysis/isec.{id}.vcf.gz", id=bcIDs[wildcards.DelFolder]) if len(bcIDs[wildcards.DelFolder])==1 else "",
      idx=lambda wildcards: expand("{{DelFolder}}/analysis/isec.{id}.vcf.gz.tbi", id=bcIDs[wildcards.DelFolder]) if len(bcIDs[wildcards.DelFolder])==1 else ""
  output:
      "{DelFolder}/analysis/experiment_merge.vcf.gz"
  shell:
      "ln -s ../../{input.data} {output}"
##################################################
rule vcf_to_gene_list:      ##convert results in gene list format:
    version: "1.00+"+str(os.path.getmtime(installation_path + "perl/vcf_to_gene_list_het.pl"))
    input:
        "{DelFolder}/analysis/experiment_merge.vcf.gz"
    output:
        "{DelFolder}/analysis/het.table.file"
    shell:
        installation_path+"perl/vcf_to_gene_list_het.pl -i {input} > {output}"
