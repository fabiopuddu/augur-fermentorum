#!/bin/bash
#This program runs through all the folders matching the pattern Del123_Yfg1
#then if a proper calling has been done, runs the analyser with standard options
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) #get the directory this program is stored in
	for folder in Del*
		do 	
			echo $folder
			del_number=`echo $folder | tr '_' "\n" | head -n1` #Identify which deletion strain we are working on
		#	Uncoment the next line to use any strain marked as C in the plate field as control
		#	control_samples=`cat name\ conversion.tsv | grep "${del_number}_" | grep C[0123456789] | awk -F"\t" '$6 ~ ERS {print $6}' | sed -e 's/[[:space:]]*$//' | tr "\n" ',' | sed 's/,*$//' `  #identify the ERS numbers of the control samples.
		#	Uncomment the next line to use only the *first* strain  marked as C in the plate field as control
		#      control_samples=`cat name\ conversion.tsv | grep "${del_number}_" | grep C[0123456789] | awk -F"\t" '$6 ~ ERS {print $6}' | sed -e 's/[[:space:]]*$//' | head -n1 | tr "\n" ',' | sed 's/,*$//' ` #added head -n1 to only get one control #identify the ERS numbers of the control samples.
		#	Uncomment the next line to always use the same wild-type control for everything
			control_samples='ERS1042082'
			echo $control_samples
			cd $folder
			if [[ -a bams_for_mpileup ]]
				then	if [[ $control_samples == '' ]] 
								then analyser-multi.sh -a -r -n2 -x -F >> results.txt &#$folder
								else export SBATCH_CMD_CIAO="analyser-multi.sh -a -r -n2 -F -c $control_samples > results.txt 2>&1"
						fi
			fi
			sbatch ${DIR}/submit_sbatch.sh
			SBATCH_CMD_CIAO=""
			cd ..
			sleep 15
		done
