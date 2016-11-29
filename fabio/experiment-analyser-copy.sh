#!/bin/bash
#This program runs through all the folders matching the pattern Del123_Yfg1
#then if a proper calling has been done, runs the analyser with standard options
	for folder in Del*
		do 	
			echo $folder
			del_number=`echo $folder | tr '_' "\n" | head -n1` #Identify which deletion strain we are working on
			control_samples=`cat name\ conversion.tsv | grep "$del_number"'_' | grep C0[0123456789][0123456789] | grep -o ERS.* | tr "\n" ',' | sed 's/,*$//' | sed 's/ERS//g' `  #identify the ERS numbers of the control samples.
			control_samples=`cat name\ conversion.tsv | grep "$del_number"'_' | grep C0[0123456789][0123456789] | grep -o ERS.* | head -n1 | tr "\n" ',' | sed 's/,*$//' | sed 's/ERS//g' `  #identify the ERS numbers of the control samples.
			echo $control_samples
			cd $folder
			if [[ -a bams_for_mpileup ]]
				then if [[ $control_samples == '' ]] 
								then analyser-multi.sh -n2 -t -r  -F -x  >> results_nov.txt &#$folder
								else analyser-multi.sh -n2 -t -r -F -C $control_samples >> results_nov.txt &#$folder
						fi
			fi
			cd ..
			sleep 10
		done
