#################
##GLOBAL TARGETS
#################

rule global_all:
    version: "1.00"
    input:
        "GO",
        "Plots",
        "xref/kobayashi",
rule clean_working_crams:
    version: "1.00"
    params:
        cram=expand(workingCRAM_repository+"{hashsample}", hashsample=map(lambda x:x[0:7]+"/"+x+".cram",  set(delname.keys()))),
        fastq=expand(workingCRAM_repository+"{hashsample}", hashsample=map(lambda x:x[0:7]+"/"+x+".fq*",  set(delname.keys()))),
    run:
        for f in params.cram:
            command=("rm -vf "+f)
            shell (command)
        for f in params.fastq:
            command=("rm -vf "+f)
            shell (command)

rule clean_merged_crams:
    version: "1.00"
    params:
        cram=expand(mergedCRAM_repository+"{sample}", sample=map(lambda x:x[0:7]+"/"+x+"_merged.cram",  set(delname.keys())))
    run:
        for f in params.cram:
            command=("rm -vf "+f)
            shell (command)

rule standard_cram:
    version: "1.00"
    input:
        lambda wildcards : expand(CRAM_repository+"{run}/{RLT}.cram", \
                                                   RLT=[item for sublist in fln_for_merge.values() for item in sublist], \
                                                   run=map(lambda x:x.split("_")[0], [item for sublist in fln_for_merge.values() for item in sublist])\
                                                   )
rule fastqc:
    version: "1.00"
    input:
        lambda wildcards : expand("post_trimming_fastqc/{RLT}.{fq}_fastqc.html",\
                            RLT=[item for sublist in fln_for_merge.values() for item in sublist],\
                            fq=["fq1", "fq2"])
rule IGV:
    version: "1.00"
    input:
        expand("{DelFolder}/BAM/{sample}.cram.crai", zip, sample=map(lambda x:bcode[x], [item for sublist in sorted(fln_for_merge.values()) for item in sublist]), DelFolder=map(lambda x:delname[bcode[x]], [item for sublist in sorted(fln_for_merge.values()) for item in sublist]))

rule reassign:
    version: "1.00"
    input:
        expand("{DelFolder}/___deletion_checked___", DelFolder=uniq_delnames_unconfirmed)
rule snapshots:
    version: "1.00"
    input:
        expand("{DelFolder}/snapshots/{sample}_"+SNAP_GENE+".png",
            zip,
            sample=map(lambda x:bcode[x], [item for sublist in sorted(fln_for_merge.values()) for item in sublist]),
            DelFolder=map(lambda x:delname[bcode[x]], [item for sublist in sorted(fln_for_merge.values()) for item in sublist])
            ) if SNAP_GENE else map(lambda x: expand(
                                                "{DelFolder_sample}_{snap_gene}.png",
                                                DelFolder_sample=delname[bcode[x]]+"/snapshots/"+bcode[x],
                                                snap_gene=list(map(lambda y:y.split('-')[0], delname[bcode[x]].upper().split("_")[1].split("+")[0:]))
                                                ) ,
                                [item for sublist in sorted(fln_for_merge.values()) for item in sublist]
                                )

rule mutations:
    version: "1.00"
    input:
        expand("{DelFolder}/reports/genotype_table.txt", DelFolder=uniq_delnames),
        expand("{DelFolder}/reports/mutation_stats_table.txt", DelFolder=uniq_delnames),
        expand("{DelFolder}/reports/mutation_table.txt", DelFolder=uniq_delnames),
        expand("{DelFolder}/reports/mut_genomemap.png", DelFolder=uniq_delnames)
rule aneuploidy:
    version: "1.00"
    input:
        expand("{DelFolder}/reports/aneuploidy_report.png", DelFolder=uniq_delnames)
rule repDNA:
    version: "1.00"
    input:
        expand("{DelFolder}/reports/repDNA_table.txt", DelFolder=uniq_delnames),
rule mouse:
    version: "1.00"
    input:
        expand("{DelFolder}/reports/mutation_stats_table.txt", DelFolder=uniq_delnames),
        expand("{DelFolder}/reports/mutation_table.txt", DelFolder=uniq_delnames),
        expand("{DelFolder}/reports/aneuploidy_report.png", DelFolder=uniq_delnames)
rule man:
    version: "1.00"
    input:
        expand("{DelFolder}/reports/mutation_stats_table.txt", DelFolder=uniq_delnames),
        expand("{DelFolder}/reports/mutation_table.txt", DelFolder=uniq_delnames),
        expand("{DelFolder}/reports/aneuploidy_report.png", DelFolder=uniq_delnames),
        expand("{DelFolder}/analysis/isec.{sample}.snps.gz", zip, sample=map(lambda x:bcode[x], [item for sublist in sorted(fln_for_merge.values()) for item in sublist]), DelFolder=map(lambda x:delname[bcode[x]], [item for sublist in sorted(fln_for_merge.values()) for item in sublist]))

