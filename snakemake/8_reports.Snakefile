rule aneuploidy_report:     #produce an aneuploidy report
    version: "1.00"
    input:
        lambda wildcards: expand("{{DelFolder}}/ploidy/{sample}.png", sample=bcIDs[wildcards.DelFolder])
    output:
        "{DelFolder}/reports/aneuploidy_report.png"
    shell:
        "montage {input} -font DejaVu-Sans -limit memory 2GB -define registry:temporary-path=/tmp -geometry 1200x1200  {output}"
##################################################
#### CREATING RESULT REPORTS
##################################################
##################################################
rule mutation_stats:
    version: "1.00+"+str(os.path.getmtime(installation_path + "perl/vcf_stats_table_all.pl"))+"+"+str(os.path.getmtime(installation_path + "perl/vcf_stats_by_colony.pl"))
    input :
        "{DelFolder}/analysis/experiment_merge.vcf.gz"
    output:
        m="{DelFolder}/reports/mutation_stats_table.txt",
        s="{DelFolder}/reports/stats_table_col.txt"
    shell:
        installation_path+"perl/vcf_stats_table_all.pl -i {input} > {output.m} && "+
        installation_path+"perl/vcf_stats_by_colony.pl -i {input} > {output.s}"
##################################################
rule mutation_table:
    version: "1.01+"+str(os.path.getmtime(installation_path + "perl/get_mut_summary_table.pl"))
    input:
        data=lambda wildcards: expand("{{DelFolder}}/analysis/isec.{sample}.vcf.gz", sample=bcIDs[wildcards.DelFolder]),
        data_inv=lambda wildcards: expand("{{DelFolder}}/analysis/inv_isec.{sample}.vcf.gz", sample=bcIDs[wildcards.DelFolder]),
        nc=nameconversion,
    output:
        "{DelFolder}/reports/mutation_table.txt"
    params:
        p="1",
        c=lambda wildcards: control[wildcards.DelFolder] if wildcards.DelFolder in control.keys() else ""
    run:
        if wildcards.DelFolder in control.keys():
            shell(installation_path+"perl/get_mut_summary_table.pl -n \"{input.nc}\" -p {params.p} -c {params.c} -d {wildcards.DelFolder} > {output}")
        else:
            shell(installation_path+"perl/get_mut_summary_table.pl -n \"{input.nc}\" -p {params.p} -d {wildcards.DelFolder} > {output}")
##################################################
rule genotype:
    version: "1.00+"+str(os.path.getmtime(installation_path + "perl/get_genotype_table.pl"))
    input:
        het_table="{DelFolder}/analysis/het.table.file",
        nc=nameconversion
    output:
        gt="{DelFolder}/reports/genotype_table.txt",
        pg="{DelFolder}/analysis/predicted_genotypes.txt"
    shell:
        installation_path+"perl/get_genotype_table.pl -i {input.het_table} -n \"{input.nc}\" > {output.gt}"
##################################################
rule mutation_location_file:
    version: "1.00+"+str(os.path.getmtime(installation_path + "perl/mutationloc.pl"))
    input:
        "{DelFolder}/analysis/experiment_merge.vcf.gz"
    output:
        temp("{DelFolder}/reports/hits_genomemap.txt")
    shell:
        installation_path+"perl/mutationloc.pl {input} > {output}"
##################################################
rule mutation_location_plot:
    version: "1.00+"+str(os.path.getmtime(installation_path + "gnuplot/results_genomeplot.gpl"))
    input:
        "{DelFolder}/reports/hits_genomemap.txt"
    output:
        "{DelFolder}/reports/mut_genomemap.png"
    shell:
        "gnuplot -e \"inp='{input}'; outp='{output}'\" "+installation_path+"gnuplot/results_genomeplot.gpl"
##################################################
rule rDNA_plots:
    version: "1.00+"+str(os.path.getmtime(installation_path + "gnuplot/rep_hist.gpl"))
    input:
        "{DelFolder}/repDNA/results.txt"
    output:
        expand("{{DelFolder}}/reports/{id}.png", id=(2,3,4,5,6,7,8,9,10,11))
    shell:
        "gnuplot -e \"DELfolder='{wildcards.DelFolder}/'\" "+installation_path+"gnuplot/rep_hist.gpl"
##################################################
rule strand_bias:
    version: "1.00+"+str(os.path.getmtime(installation_path + "perl/strand_bias.pl"))
    input:
        "{DelFolder}/analysis/experiment_merge.vcf.gz",
        "{DelFolder}/reports/mutation_stats_table.txt"
    output:
        "{DelFolder}/analysis/trsc_strand.txt",
        "{DelFolder}/analysis/trsc_asym.txt",
        "{DelFolder}/reports/strand_bias.png",
        "{DelFolder}/analysis/rplc_loc.txt",
    shell:
        installation_path+"perl/strand_bias.pl {wildcards.DelFolder}"
##################################################
rule make_report:
    version: "1.00+"+str(os.path.getmtime(installation_path + "snakemake/make_report.py"))
    input:
        expand("{{DelFolder}}/reports/{id}.png", id=(2,3,4,5,6,7,8,9,10,11)),
        "{DelFolder}/repDNA/results.txt",
        "{DelFolder}/reports/aneuploidy_report.png",
        "{DelFolder}/analysis/predicted_genotypes.txt",
        "{DelFolder}/reports/strand_bias.png",
        "{DelFolder}/reports/genotype_table.txt",
        "{DelFolder}/reports/mutation_stats_table.txt",
        "{DelFolder}/reports/mutation_table.txt",
        "{DelFolder}/reports/repDNA_table.txt",
        "{DelFolder}/reports/deletion_report.txt",
        "{DelFolder}/reports/mut_genomemap.png",
    output:
        o="{DelFolder}/reports/WGS_report.pdf",
        l="all_reports/{DelFolder}_WGS_report.pdf"
    shell:
        installation_path+"snakemake/make_report.py {wildcards.DelFolder} &&"
        "ln -s ../{output.o} {output.l}"
