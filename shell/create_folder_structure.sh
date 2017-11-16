#!/bin/bash
#  vcftoresult.sh
#  Created by Fabio on 23/11/2014.
#
sort -nu samples.tsv > sorted_samples.tsv
mv sorted_samples.tsv samples.tsv
#specify default source path
source_path='/mnt/scratch/jackson/fp305/from_CGP/CGP_BAM'
if [[ -z $1 ]]
	then source_path='/mnt/scratch/jackson/fp305/from_CGP/CGP_BAM'
	else source_path="../../$1"
fi
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
		folders=`cat name\ conversion.tsv | grep "Del$n"'_' | grep -oE '\bSD[[:alnum:]]*\b|\bSC_MFY[[:alnum:]]*\b' | tr "\n" "\t"`
		fnames=`cat name\ conversion.tsv | grep "Del$n"'_' | grep -oE '\bSD[[:alnum:]]*\b|\bSC_MFY[[:alnum:]]*\b' | tr "\n" "\t"`
		mkdir -p "Del$n"'_'"$name"
			cd "Del$n"'_'"$name"
			mkdir -p "BAM"
				cd BAM
				for fname in $fnames
					do  ln -s "${source_path}/$fname/$fname.bam" "$PWD/$fname.bam"
						ln -s "${source_path}/$fname/$fname.bam.bai" "$PWD/$fname.bam.bai"
						ln -s "${source_path}/$fname/$fname.fq1.gz" "$PWD/$fname.fq1.gz"
						ln -s "${source_path}/$fname/$fname.fq2.gz" "$PWD/$fname.fq2.gz"
					done
				cd ..	
			cd ..
	done
