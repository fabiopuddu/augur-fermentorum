if ORGANISM == 'yeast':
  rule karyotype:             #Calculate chromosome copy number and rearrangements
      version: "1.00+"+str(os.path.getmtime(installation_path+"perl/CGH.pl"))
      input:
          bam="{DelFolder}/BAM/{sample}.cram",
          idx="{DelFolder}/BAM/{sample}.cram.crai"
      output:
          "{DelFolder}/ploidy/{sample}_ploidy_data.txt",
          "{DelFolder}/ploidy/{sample}_highlights.txt",
          "{DelFolder}/ploidy/{sample}_breakpoints.txt",
          "{DelFolder}/ploidy/{sample}_plstats.txt",
          "{DelFolder}/ploidy/{sample}_plstats_raw.txt",
      params:
          p=lambda wildcards: ploidy[wildcards.sample]
      shell:
          installation_path+"perl/CGH.pl -o "+ORGANISM+" -i {input.bam} -p {params.p} -f "
  ##################################################
  rule circos:
      version: "1.00+"+str(os.path.getmtime(circos_config))
      input:
         data="{DelFolder}/ploidy/{sample}_ploidy_data.txt",
         highl="{DelFolder}/ploidy/{sample}_highlights.txt",
         conf=circos_config
      output:
         png="{DelFolder}/ploidy/{sample}.png",
         web="{DelFolder}/ploidy/{sample}_web.jpg",
         svg="{DelFolder}/ploidy/{sample}.svg"
      params:
          filename=lambda wildcards: wildcards.sample,
          d=lambda wildcards: delname[wildcards.sample],
          s=lambda wildcards: wildcards.sample,
          e=lambda wildcards: aka[wildcards.sample]
      shell:
         "circos -silent -conf {input.conf} -param highlights/highlight/file={input.highl} -param plots/plot/file={input.data} -outputfile {wildcards.DelFolder}/ploidy/{params.filename} &&"
         "magick {output.png} -font DejaVu-Sans -weight 70  -gravity center -pointsize 60 -annotate 0 \"\n{params.d}\n\n \"  -pointsize 30 -annotate 0 \"\n{params.s}   {params.e}\" {wildcards.DelFolder}/ploidy/out.{params.filename}.png &&"
         "mv {wildcards.DelFolder}/ploidy/out.{params.filename}.png {output.png} &&"
         "magick {output.png} -quality 96 -resize 500x500  {output.web}"
else:
  rule karyotype:             #Calculate chromosome copy number and rearrangements
      version: "1.00+"+str(os.path.getmtime(installation_path+"perl/CGH.pl"))
      input:
          bam="{DelFolder}/BAM/{sample}.cram",
          idx="{DelFolder}/BAM/{sample}.cram.crai"
      output:
          "{DelFolder}/ploidy/{sample}_ploidy_data.txt"
      params:
          p=lambda wildcards: ploidy[wildcards.sample]
      shell:
          installation_path+"perl/CGH.pl -o "+ORGANISM+" -i {input.bam} -p {params.p}"
  ##################################################
  rule circos:
      version: "1.00+"+str(os.path.getmtime(circos_config))
      input:
         data="{DelFolder}/ploidy/{sample}_ploidy_data.txt",
         conf=circos_config
      output:
         png="{DelFolder}/ploidy/{sample}.png",
         web="{DelFolder}/ploidy/{sample}_web.jpg",
         svg="{DelFolder}/ploidy/{sample}.svg"
      params:
          filename=lambda wildcards: wildcards.sample,
          d=lambda wildcards: delname[wildcards.sample],
          s=lambda wildcards: wildcards.sample,
          e=lambda wildcards: aka[wildcards.sample]
      shell:
         "circos -silent -conf {input.conf} -param plots/plot/file={input.data} -outputfile {wildcards.DelFolder}/ploidy/{params.filename} &&"
         "magick {output.png} -font DejaVu-Sans -weight 70  -gravity center -pointsize 60 -annotate 0 \"\n{params.d}\n\n \"  -pointsize 30 -annotate 0 \"\n{params.s}   {params.e}\" {wildcards.DelFolder}/ploidy/out.{params.filename}.png &&"
         "mv {wildcards.DelFolder}/ploidy/out.{params.filename}.png {output.png} &&"
         "magick {output.png} -quality 96 -resize 500x500  {output.web}"
