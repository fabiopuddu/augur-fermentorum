## Augur Fermentorum
This repository contains a set of tools to analyse yeast whole genome sequencing data for repetitive DNA genome variation and aneuploidies.

The system is based on snakeake and is designed to interact with an SQL database to read data and store results. This guide does not (yet) cover setting up the database.

### Requirements
* A unix system
* HPC based on SLURM workload manager
* Perl 
* Python
* Samtools 1.9
* Bcftools 1.9
* Pigz 2.3.1
* Conda
* Snakemake
* IGV


### Build environment

Use the file snakemake/snakemake2.yaml to build the environment required for this tool to work

### Adapt workflow to your environment

Edit the file snakemake/1_main.Snakefile to ensure these pathways are properly set

* **CRAM_repository**: where CRAM files aligned to standard references are to be stored when downloaded from ENA, or de-novo aligned from downloaded fastq
* **mergedCRAM_repository**: where intermediate CRAM files are to be stored after data is merged (when the same library is split in 2 or more lanes)
* **workingCRAM_repository**: where CRAM files are to be stored after they have been re-aligned to the working copy of the reference genome
* **samtools**: path to samtools executable
* **bcftools**: path to bcftools executable
* **IGV_path**: path to IGV folder

### Usage
#### Data input requirements:

* A tab separated file named "name conversion.tsv" containing details of all the sequencing data to be analysed. One line per BAM file and the following columns:
   * Tube barcode identifier (unique key in your sequencing database)
   * Sample group name: this takes a standard format e.g. Del1_TEL1 (see note  below)
   * Plate ID
   * Human-friendly Sample name, ideally unique
   * Name of the BAM file
   * ERS Sample ID; can also be "None"
   * Ploidy of the sample (1,2)
   * URL to download the data (bam or cram file; or a pair of fastq files, semicolon separated); can also be "None"
   * Binary field indicating wether a bam file corresponds to a sample(0) or a control (1)
   * Sample Name Barcode (unique)
   * Control field: "+" to ensure the sample is processed

#### Running the Pipeline

Use the script snakemake/brew; here is a set of options that one can add in addition to the standard snakemake options

#####Targets:
 * **n\_c**: query the sql database and create a name conversion file
 * **reassign**: check expected deletion or mutation and reassign non-deleted samples based on their barcodes (yeast only). Remember to delete all the data and rerun brew n_c after reassigning
 * **snapshots**: take IGV snapshots of the genes in the DelFolder name(default behaviour)
 * **all**: run the pipeling producing the full report (yeast only)
 * **clean_all**:             remove all the data from the working directory
 * **clean\_working\_crams**: remove all the data of the current experiment from the repository of realigned data
 * **clean\_merged\_crams**:  remove all the data of the current experiment from the repository of merged data
 
#####Config options:
 * **org="organism"**  choose from yeast(Default), mouse, man; remember to use the appropriate target for each organims (yeast=all; mouse=mouse; man=man)
 * **snap="gene"**   take a snapshot from the selected gene instead of default
 * **exper="experiment"**    experiment number; must be used in conjunction with n_c, if it is not used the pipeline will attempt to guess it from the folder name e.g. 1_ThisIsExperimentOne
 * **vc="variant_caller"**   choose from mpileup(Default) or freebayes (under development)
 * **n_p="normal_panel"** before anything subtract all the mutations found in the indicated normal panel (currently available for man: HAP1(HAP1+dbSNP), RPE1(RPE1+dbSNP); normal panel for mouse is by default AN3-12)
 





