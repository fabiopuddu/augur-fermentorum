#!/bin/bash

for folder in Del*
		do if [[ -a $folder/analysis/stat_table.txt ]]
				then cat $folder/analysis/stat_table.txt | grep Del >> mutation_statistics.txt		
			fi	
		done
		
