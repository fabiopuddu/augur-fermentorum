## Augur Fermentorum
This repository contains a set of tool to analyse yeast whole genome sequencing data and measure genomic stability

### Requirements
* A unix system
* HPC based on SLURM workload manager
* Perl 
* Python
* Samtools 1.3.1
* Pigz 2.3.1

### Usage
#### Data input requirements:
* A tab-separated file named "samples.tsv" containing on each line an isogenic group of samples to be analysed together. One line per sample and the following columns:
   * Sample number (e.g. "1")
   * Sample name(e.g. "TEL1")
* A tab separated file named "name conversion.tsv" containing details of all the sequencing data to be analysed. One line per BAM file and the following columns:
   * Sequencing Name (barcode)
   * Sample name, derived from samples.tsv and in the format "Del1_TEL1"
   * Sequencing Lane
   * Human-friendly Sample name, ideally unique
   * Name of the BAM file
   * ERS Accession number; has to be unique in each line if not available can be a custom number in the format "ERS000000"
* A folder containing all the BAM files to be analysed. 

#### Creation of the analysis folder structure:

#### Fastq extraction and realignment to custom genomes:

#### Variant Calling

#### Running the Pipeline

#### Recovering information on repetitive DNA and ploidy

#### Adjusting repetitive DNA estimates from samples across different sequencing lanes
