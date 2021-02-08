#Snakemake file
#################
###FUNCTIONS AND MODULES
#################
from mysql_connector import get_metadata, read_nc, write_nc, update_del_database, update_repDNA_db, suppressor_choice, update_mutDB, update_repl_asymDB, update_vcf_db
import os
cwd = os.getcwd()
###############
### CHOICE OF ORGANISM AND EXPERIMENT
###############

try:
    NO_SUPP_CHOICE=config["no_supp_choice"]
except:
    NO_SUPP_CHOICE=0
try:
    ORGANISM=config["org"]
except:
    ORGANISM="yeast"
try:
    NORMPAN=config["n_p"]
except:
    NORMPAN=None
try:
    SNAP_GENE=config["snap"].upper()
except:
    SNAP_GENE=""
try:
    experiment_number = config["exper"]
except:
    try:
        experiment_number = int(cwd.split("/")[-1].split("_")[0])
    except ValueError:
        exit('Could not determine experiment number and one was not provided')
try:
    variant_caller = config["vc"]
except:
    variant_caller="mpileup"
try:
    config["force_bam_conversion"]
    force_bam_conversion=True
except:
    force_bam_conversion=False

###############
### PATHS
###############
CRAM_repository="/mnt/home1/jackson/fp305/scratch/cram_repositories/primary_cram_repository/"        #Where CRAM files aligned to standard references are stored, either from Sanger, ENA, or de-novo aligned from downloaded fastq
mergedCRAM_repository="/mnt/home1/jackson/fp305/scratch/cram_repositories/merged_cram_repository/"   #Where CRAM files are stored after data is merged (when 2 or more lanes have the same barcode ID)
workingCRAM_repository="/mnt/home1/jackson/fp305/scratch/cram_repositories/working_cram_repository/" #Where CRAM files are stored after they have been re-aligned to the custom working copy of the reference genome
installation_path=workflow.basedir+"/../"  #where augur-fermentorum is installed

#PROGRAMS
samtools="/mnt/home1/jackson/fp305/sw/samtools_1.9/samtools/samtools"
bcftools="/mnt/home1/jackson/fp305/sw/samtools_1.9/bcftools/bcftools"

IGV_path="/mnt/home1/jackson/fp305/scratch/IGV_Linux_2.5.0/"

nameconversion="name_conversion.tsv"                                                 #Standard name of the file containing sample information
chromosome_list=list()
genome_ref=str()
standard_genome_ref=str()
species=""
print("The organism chosen is:",ORGANISM,file=sys.stderr)
print("The variant caller chosen is:",variant_caller,file=sys.stderr)
if ORGANISM == "yeast":
    genome_ref = installation_path+"mpileup_defaults/new_reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm_MASKED+REPDNA.toplevel.fa"
    chromosome_list=["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
    standard_genome_ref = installation_path+"mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa"
    species = "saccharomyces_cerevisiae"
    normal_panel=""
    circos_config=installation_path+"defaults/circos_yeast.conf"                                                               #the circos config file                                                  #the reference for repetitive DNA
elif ORGANISM == "mouse":
    genome_ref =  "/mnt/home1/jackson/fp305/scratch/reference_genomes/mouse/Mus_musculus.GRCm38.68.dna.toplevel.fa"
    standard_genome_ref = "/mnt/home1/jackson/fp305/scratch/reference_genomes/mouse/Mus_musculus.GRCm38.68.dna.toplevel.fa"
    chromosome_list=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,"X", "Y"]
    species = "mus_musculus"
    normal_panel="/mnt/home1/jackson/fp305/scratch/reference_genomes/mouse_background_vcf/AN3-12.snps+indels.vcf.gz"
    circos_config=installation_path+"defaults/circos_mouse.conf"                                                               #the circos config file
elif ORGANISM == "man":
    np_db={
        "HAP1":"/mnt/home1/jackson/fp305/scratch/reference_genomes/1000Genomes_hs37d5/all/normal_panels/HAP1+dbSNP_panel.vcf.gz",
        "RPE1":"/mnt/home1/jackson/fp305/scratch/reference_genomes/1000Genomes_hs37d5/all/normal_panels/RPE1+dbSNP_panel.vcf.gz"
    }
    genome_ref =  "/mnt/home1/jackson/fp305/scratch/reference_genomes/1000Genomes_hs37d5/all/fasta/hs37d5.fa.gz"
    standard_genome_ref = "/mnt/home1/jackson/fp305/scratch/reference_genomes/1000Genomes_hs37d5/all/fasta/hs37d5.fa.gz"
    chromosome_list=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X", "Y"]
    species = "homo_sapiens"
    if NORMPAN != None:
        normal_panel=np_db[NORMPAN]
    else:
        exit('No normal panel was specified!')
    circos_config=installation_path+"defaults/circos_human.conf"                                                               #the circos config file
else:
    exit('ORGANISM NOT SUPPORTED')
print("The normal panel chosen is:",normal_panel,file=sys.stderr)
######################
#PRELIMINARY ACTIONS
######################
#cwd = cwd.split("/")[-1]
# print ("Working on: "+cwd,file=sys.stderr)
if not os.path.isfile(nameconversion):
    print ("Need to create nameconversion file")
    localrules: write_nameconverion_file
    rule n_c:
        version: "1.00"
        input:
            nameconversion
    rule write_nameconverion_file:
        version: "1.00"
        output:
            nameconversion
        run:
            write_nc(nameconversion, experiment_number)
else:
    print ("Proceeding with analysis",file=sys.stderr)
    print ("Force bam conversion?", force_bam_conversion ,file=sys.stderr)
    ploidy,delname,platLib,ers,bcode,urls,bcIDs,fln_for_merge,control,aka,SD,delcheck = read_nc(nameconversion)
    uniq_delnames=set(delname.values())
    if 'Del0000_WildtypeControl' in uniq_delnames:
    	uniq_delnames.remove('Del0000_WildtypeControl')
    uniq_delnames_unconfirmed=list()
    for k,v in delcheck.items():
        if v != "+" and v != "X":
            uniq_delnames_unconfirmed.append(delname[k])
    uniq_delnames_unconfirmed=set(uniq_delnames_unconfirmed)
    print ("Deletions to be confirmed:", uniq_delnames_unconfirmed,file=sys.stderr)
    #print (control,file=sys.stderr)
    #print (uniq_delnames,file=sys.stderr)
    #print (bcIDs,file=sys.stderr)
    #print (platLib,file=sys.stderr)
    if len(bcIDs) == 0:
        exit("Error")
    #################
    ##LOCAL RULES
    #################
    #These rules will be executed on the head node and not submitted to slurm
    if ORGANISM == "yeast":
        localrules:
            download_crams_from_datastore,\
            complete,\
            take_snapshot,\
            clean_merged_crams, \
            clean_working_crams, \
            download_crams,\
            download_fastq, \
            link_crams, \
            link_bams, \
            index_bams,\
            repetitive_table, \
            create_control, \
            link_experiment_vcf, \
            genotype, \
            mutation_location_file, \
            mutation_location_plot,\
            rDNA_plots, \
            strand_bias, \
            make_report, \
            attempt_reassign, \
            insert_mysql_results,  \
            aneuploidy_report,\
            mutation_stats,\
            mutation_table,\
#            prepare_vcf,\
#            pre_filtering,\
#            filter_control,\
#            intersections,\
#            inverse_intersections,\
            merge_experiment_vcf,\
            vcf_to_gene_list,
            global_replication_bias_plots,\
            rename_CI_fastq,\
            fastq_quality_control_pre,\
            fastq_quality_control_post,
            deletion_report,\
            rehead_vcf,\
#            merge_chr_vcf,\
            link_normal_panel
#            variant_effect_predictor
    else:
        localrules: download_crams_from_datastore, \
                    clean_merged_crams, \
                    clean_working_crams, \
                    download_crams,\
                    download_fastq, \
                    link_crams, \
                    link_bams,\
                    repetitive_table, \
                    create_control, \
                    link_experiment_vcf, \
                    genotype, \
                    mutation_location_file, \
                    mutation_location_plot,rDNA_plots, \
                    strand_bias, \
                    make_report, \
                    attempt_reassign, \
                    insert_mysql_results, \
                    complete, \
                    aneuploidy_report
    rule all:
        version: "1.00"
        input:
            "ALL_DONE"
    include:
        "0_targets.Snakefile"
    ##################################################
    ##WILDCARDS CONSTRAINTS
    ##################################################
    wildcard_constraints:
        sample="SC_MFY[0-9]+|[0-9]+STDY[0-9]+|FP[0-9]+",
        RLT="[0-9]+_[0-9]+#[0-9]+",
        DelFolder="Del[0-9]+_.+"
    ruleorder:
        link_crams > rehead_crams
    ruleorder:
        merge_experiment_vcf > link_experiment_vcf

    ##################################################
    include:
        "2_fastq_to_cram.Snakefile"
    include:
        "3_cram_management.Snakefile"
    ##################################################
    ####SAMPLE CHECK
    ##################################################
    if ORGANISM == "yeast":
        include:
            "4_deletion_check.Snakefile"
    ##################################################
    ##RULES EXTRACTING REPETITIVE DNA INFORMATION FROM DATA
    ##################################################
    include:
        "6_repDNA.Snakefile"
    ##################################################
    ##RULES EXTRACTING CHROMOSOMAL POINT/SMALL MUTATION INFORMATION FROM DATA
    ##################################################
    include:
        "7_small_variants.Snakefile"
    ##################################################
    ##RULES EXTRACTING MITOCHONDRIAL POINT/SMALL MUTATION INFORMATION FROM DATA
    ##################################################
    rule create_temp_bam:
        version: "1.00"
        input:
            bam="{DelFolder}/BAM/{sample}.cram",
            idx="{DelFolder}/BAM/{sample}.cram.crai"
        output:
            bam="{DelFolder}/BAM/mt.{sample}.bam",
            idx="{DelFolder}/BAM/mt.{sample}.bam.bai"
        shell:
            samtools+" view -h -@ 8 -O BAM -o {output.bam} {input.bam} Mito && " +
            samtools+" index {output.bam}"
    ##################################################
    rule lofreq_caller_mt:
        version: "1.00"
        input:
            ref=genome_ref,
            bam="{DelFolder}/BAM/mt.{sample}.bam",
            idx="{DelFolder}/BAM/mt.{sample}.bam.bai"
        output:
            "{DelFolder}/calling/mt.{sample}.vcf.gz"
        shell:
            "lofreq call -f {input.ref} -r Mito {input.bam} | bgzip > {output}"
    ##################################################
    rule filter_mt_mutations:
        version: "1.00"
        input:
            "{DelFolder}/calling/mt.{sample}.vcf.gz",
            "{DelFolder}/repDNA/{sample}.txt"
        output:
            "{DelFolder}/analysis/filt.mt.{sample}.vcf.gz"
        shell:
            "touch x"
    ##################################################
    ####RULE EXTRACTING INFORMATION ON PLOIDY FROM DATA
    ##################################################
    include:
        "5_karyotype.Snakefile"
    ##################################################
    include:
        "8_reports.Snakefile"
    ##################################################
    rule global_replication_bias_plots:
        version: "1.00+"+str(os.path.getmtime(installation_path + "gnuplot/replication_bias_global.gpl"))
        input:
            expand("{DelName}/reports/mutation_stats_table.txt",DelName=uniq_delnames),
            expand("{DelName}/analysis/rplc_loc.txt",DelName=uniq_delnames)
        output:
            "replication_strand_bias.png"
        shell:
            "gnuplot -e 'DelName=\""+" ".join(uniq_delnames)+"\"' "+installation_path+"gnuplot/replication_bias_global.gpl "
    ##################################################
    rule insert_mysql_results:
        version: "2.0"
        resources:
            load=50
        input:
            mut="{DelFolder}/analysis/predicted_genotypes.txt",
            mutstat="{DelFolder}/reports/stats_table_col.txt",
            repl_asym="{DelFolder}/analysis/rplc_loc.txt",
            vcf_f=lambda wildcards: expand("{{DelFolder}}/analysis/isec.{sample}.vcf.gz", sample=list(bcIDs[wildcards.DelFolder]) )
        output:
           "{DelFolder}/__SQL_REP_UPDATE__"
        run:
            import vcfpy
            supp_answer=None
            for file in input.vcf_f:
                smpl=file.split("/")[-1].split(".")[1]
                vcf_reader=vcfpy.Reader.from_path(file)
                sample=vcf_reader.header.samples.names[0]
                for Record in vcf_reader:
                    if 'PASS' in Record.FILTER:
                        update_vcf_db(experiment_number,smpl,\
                                  Record.CHROM, \
                                  Record.POS, \
                                  Record.ID, \
                                  Record.REF, \
                                  str(",".join([alt.type for alt in Record.ALT])),\
                                  str(",".join([alt.value for alt in Record.ALT])),\
                                  Record.QUAL, \
                                  str(",".join(Record.FILTER)), \
                                  Record.INFO['CSQ'][0],\
                                  Record.calls[0].data.get('GT'),\
                                  Record.calls[0].data.get('DP'),\
                                  Record.calls[0].data.get('DV'),\
                                  Record.calls[0].data.get('GQ'))
            with open(input.mutstat) as f:
                header=f.readline().rstrip().split("\t")
                d     =f.readline().rstrip().split("\t")
                bcids=d[-1].split(":")
                data={}
                for b in bcids:
                    data[b]={}
                    for n in range(0,len(header)-1):
                        if n==4:
                            continue
                        data[b][header[n]]=d[n].split(":")[bcids.index(b)]
                for b in bcids:
                    update_mutDB(b,data[b],experiment_number)
            if supp_answer == "NoToAll" or NO_SUPP_CHOICE:
                pass
            else:
                supp_answer=suppressor_choice(input.mut, "Screening")
                print(supp_answer)
            with open (input.repl_asym) as f:
                delname= file.split("/")[0]
                exper=experiment_number
                header=f.readline().rstrip().split("\t")
                decile=10
                for line in f:
                    line=line.rstrip()
                    riga=line.split("\t")
                    if decile <=100:
                        update_repl_asymDB(decile,exper,delname,header,riga)
                    decile=decile+10
            shell('touch {output}')
    ##################################################
    rule complete:
        version: "1.00"
        input:
            a=expand("all_reports/{DelFolder}_WGS_report.pdf", DelFolder=uniq_delnames),
            b=expand("{DelName}/__SQL_REP_UPDATE__", DelName=uniq_delnames),
            rep="rep.txt",
        output:
            "ALL_DONE"
        run:
            with open (input.rep) as f:
                header=list()
                for line in f:
                    line=line.rstrip()
                    riga=line.split("\t")
                    if "#" in riga[0]:
                        header=riga
                    else:
                        if header:
                            update_repDNA_db(header,riga,experiment_number)
            shell("rm -rf slurm_files && touch {output}")