## Augur Fermentorum
This repository contains a set of tools to analyse yeast whole genome sequencing data for repetitive DNA genome variation and aneuploidies.

### Requirements
* A unix system
* HPC based on SLURM workload manager
* Perl 
* Python
* Samtools 1.3.1
* Pigz 2.3.1
*...


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

* A folder containing the BAM files to be analysed. 


#### Fastq extraction:

In the directory containing the BAM files run:

`shell/extract-fastq.sh`

This detects BAM or CRAM files in the current directory, converts them in the appropriate form (cram>bam>fastq) and sorts the files in the structure that is required by subsequent programs:

#### Creation of the analysis folder structure:

In the directory where the "samples.tsv" and "name conversion.tsv" files have been created run:

`shell/create_folder_structure.sh <absolute/path/to/folder/with/BAM/files>`

This will create symlinks to the appropriate BAM and fastq files.

#### Realignment to custom genomes for Ty and 2µ counting and mating type estimation:

In the directory where the "samples.tsv" and "name conversion.tsv" files have been created run:
`shell/experiment-ty-realigner.sh`

#### Variant Calling
In the directory where the "samples.tsv" and "name conversion.tsv" files have been created run:
`shell/experiment-caller.sh`

#### Running the Pipeline
In the directory where the "samples.tsv" and "name conversion.tsv" files have been created run:
`shell/experiment-analyser.sh`

#### Recovering information on repetitive DNA and ploidy
In the directory where the "samples.tsv" and "name conversion.tsv" files have been created run:
`perl perl/gather-rep-data.pl` and redirect the output to a file (e.g. `>repDNA.txt`)

#### Adjusting repetitive DNA estimates from samples across different sequencing lanes
In the directory where the "samples.tsv" and "name conversion.tsv" files have been created run:
`perl perl/lane-median-adjust.pl <repDNA file> <name conversion file>`, where <repDNA file> is the file containing the repetitive DNA data (e.g. "repDNA.txt") and <name conversion file> is "name conversion.tsv"


