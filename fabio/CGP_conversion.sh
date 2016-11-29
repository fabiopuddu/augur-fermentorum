#!/bin/sh
#This scripts imports the CGP bam files data into the Yeast folder, mocking the TK-DA pipeline data structure with symlinks
#this script has to be run inside the Esperimento_150 folder 

for strain in ../../BAM\ files/Yeast_CGP/SD*.bam
		do 	
			strain_name=$(echo $strain | sed 's/.bam//g'| sed 's|../../BAM files/Yeast_CGP/||g')
			sample_code=`cat name\ conversion.tsv | grep "$strain_name\t" | awk '{print $1}'`
			if [[ "$sample_code" == "" ]]
				then empty="1"
				else empty="0"
			fi			
			if [[ ! -a "../../BAM files/Yeast/$sample_code" && "$empty" == "0" ]]
				then mkdir "../../BAM files/Yeast/$sample_code"
					 ln -s "../../Yeast_CGP/$strain_name.bam" "../../BAM files/Yeast/$sample_code/$strain_name.bam" 
					 ln -s "../../Yeast_CGP/$strain_name.bai" "../../BAM files/Yeast/$sample_code/$strain_name.bai"
					 printf "processed $sample_code $strain_name\n"
			fi
		done
		