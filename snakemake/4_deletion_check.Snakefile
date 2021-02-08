rule detect_deletion:
    version: "1.00+"+str(os.path.getmtime(installation_path+"perl/new_deletion_check.pl"))
    input:
        bam="{DelFolder}/BAM/{sample}.cram",
        idx="{DelFolder}/BAM/{sample}.cram.crai",
    output:
        d="{DelFolder}/del+bcd/{sample}.del",
    params:
        gene=lambda wildcards : delname[wildcards.sample].split("_")[1] #Del101_CDC34 > CDC34
    shell:
        installation_path+"perl/new_deletion_check.pl -i {input.bam} -d {wildcards.DelFolder} -o yeast > {output.d}"
##################################################
rule attempt_reassign:
    version: "1.00"
    input:
        bc=lambda wildcards: expand("{{DelFolder}}/del+bcd/{sample}.bcd", sample=bcIDs[wildcards.DelFolder]),
        de=lambda wildcards: expand("{{DelFolder}}/del+bcd/{sample}.del",sample=bcIDs[wildcards.DelFolder])
    output:
        "{DelFolder}/___deletion_checked___"
    run:
        best=dict()
        for f in bcIDs[wildcards.DelFolder]:
            with open(wildcards.DelFolder+"/del+bcd/"+f+".del") as gh:
                li=gh.readline().rstrip()
                answ=li.split("\t")[5].split(":")[1] if not "FAILED" in li else "?"
                x=gh.readline().split("\t")[2]
            if answ == "N":
                with open(wildcards.DelFolder+"/del+bcd/"+f+".bcd") as fh:
                    l=fh.readline().rstrip()
                    gene=l.split("\t")[0]
                    reads_supporting=0
                    try:
                        most_likely_candidate=l.split("\t")[2].split(";")[0].split(" ")[0]
                        reads_supporting=int(''.join([c for c in l.split("\t")[2].split(";")[0].split(" ")[2] if c in '1234567890']))
                    except IndexError:
                        most_likely_candidate="X"    #mark with X the samples that we cannot tell anything about and we will not analyse further
                    if most_likely_candidate.upper() == gene.upper() or reads_supporting <3:
                        most_likely_candidate="X"
                    print ("Reassigning sample", f, "from", x , "to", most_likely_candidate)
                    update_del_database(f,most_likely_candidate, experiment_number)
            elif answ == "Y":
                update_del_database(f,"+",experiment_number)
        shell("touch {output}")
##################################################
rule detect_tags:
    version: "1.00+"+str(os.path.getmtime(installation_path+"perl/find_deletion_tag.pl"))
    input:
        bam="{DelFolder}/BAM/{sample}.cram"
    output:
        "{DelFolder}/del+bcd/{sample}.bcd"
    shell:
        installation_path+"perl/find_deletion_tag.pl {wildcards.sample} {wildcards.DelFolder}"
##################################################
rule deletion_report:
    version: "1.00"
    input:
        d=lambda wildcards:expand("{{DelFolder}}/del+bcd/{sample}.del",sample=bcIDs[wildcards.DelFolder]),
        b=lambda wildcards:expand("{{DelFolder}}/del+bcd/{sample}.bcd",sample=bcIDs[wildcards.DelFolder])
    output:
        "{DelFolder}/reports/deletion_report.txt"
    shell:
        "echo =====Deletion Check Report======= > {output} &&"
        "cat {input.d} >> {output} &&"
        "echo =====Barcode Report============== >> {output} &&"
        "cat {input.b} >> {output}"
##################################################
rule take_snapshot:
    version: "1.00"
    input:
        fi="{DelFolder}/BAM/{sample}.cram",
        ind="{DelFolder}/BAM/{sample}.cram.crai"
    output:
        "{DelFolder}/snapshots/{sample}_{gene}.png"
    run:
        from telnetlib import Telnet
        gene_db=dict()
        sysname=dict()
        with open (installation_path+"defaults/all_yeast_genes.txt") as f:
            for line in f:
                line=line.rstrip()
                riga=line.split("\t")
                if '#' in riga[0]:
                    continue
                start=int(riga[2])-5000
                end=int(riga[3])+5000
                gene_db[riga[0]]="chr"+riga[1]+":"+str(start)+"-"+str(end)
                try:
                    sysname[riga[4]]=riga[0]
                except:
                    sysname[riga[0]]=riga[0]
        gene=wildcards.gene
        sn=sysname[gene]
        position=""
        try:
            position=gene_db[sn]
        except:
            pass
        if position:
            with Telnet('cb-head2.gurdon.private.cam.ac.uk', 60151) as tn:
                tn.write(b"new\n");
                tn.read_until(b"OK");
                tn.write(b"genome cerevisiae+repDNA\n");
                tn.read_until(b"OK");
                tn.write(b"snapshotDirectory "+bytes(cwd, encoding='utf-8')+b"/"+bytes(wildcards.DelFolder, encoding='utf-8')+b"/snapshots \n");
                # print("snapshotDirectory "+cwd+"/"+wildcards.DelFolder+"/snapshots \n");
                tn.read_until(b"OK");
                tn.write(b"load "+bytes(cwd, encoding='utf-8')+b"/"+bytes(input.fi, encoding='utf-8')+b"\n");
                tn.read_until(b"OK");
                tn.write(b"goto "+bytes(position, encoding='utf-8')+b"\n");
                tn.read_until(b"OK");
                tn.write(b"snapshot "+ bytes(wildcards.sample, encoding='utf-8')+b"_"+bytes(gene,encoding='utf-8')+b".png\n");
                # print ("snapshot "+wildcards.sample+"_"+gene+".png\n");
                tn.read_until(b"OK");
