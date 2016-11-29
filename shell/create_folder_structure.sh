#!/bin/bash
#  vcftoresult.sh
#  Created by Fabio on 23/11/2014.
#
sort -nu samples.tsv > sorted_samples.tsv
mv sorted_samples.tsv samples.tsv
cat samples.tsv | while read line
	do 	t=1
			for tab in $line
				do	if [[ $t -eq 1 ]] 
						then n=$tab
					fi
					if [[ $t -eq 2 ]]
						then name=$tab
					fi
					t=$(($t+1))
				done
		folders=`cat name\ conversion.tsv | grep "Del$n"'_' | awk '{print $1}' | tr "\n" "\t"`
		fnames=`cat name\ conversion.tsv | grep "Del$n"'_' | awk '{print $5}' | tr "\n" "\t"`
		mkdir -p "Del$n"'_'"$name"
			cd "Del$n"'_'"$name"
			mkdir -p "BAM"
				cd BAM
				for fname in $fnames
					do  ln -s /mnt/scratch/jackson/fp305/from_CGP/CGP_BAM/$fname/$fname.bam $PWD/$fname.bam
						ln -s /mnt/scratch/jackson/fp305/from_CGP/CGP_BAM/$fname/$fname.bam.bai $PWD/$fname.bam.bai
						ln -s /mnt/scratch/jackson/fp305/from_CGP/CGP_BAM/$fname/$fname.fq1.gz $PWD/$fname.fq1.gz
						ln -s /mnt/scratch/jackson/fp305/from_CGP/CGP_BAM/$fname/$fname.fq2.gz $PWD/$fname.fq2.gz
					done
				cd ..	
			cd ..
	done
