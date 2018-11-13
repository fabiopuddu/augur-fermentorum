#!/bin/bash
#This program runs through all the folders matching the pattern Del123_Yfg1
#then if a proper calling has been done, runs the analyser with standard options

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) #get the directory this program is stored in
OPTIND=1         # Reset in case getopts has been used previously in the shell.
control_defined=0
automatic_defined=0
multiple=0
function hlp {
   printf "############   HELP   ###############\nOPTIONS\n"
   printf "This program automates the running of the analyser pipeline,\n"
   printf "providing it with the reference to the proper control\n"
   printf "\t-h\tThis Help\n"    
   printf "\t-a\tAutomatic Detection of the control\n"
   printf "\t-m\tAllow multiple automatic controls\n"
   printf "\t-c\tDefine a standard wild-type control for all samples\n"  
   exit 1
}

while getopts "hac:m" opt
    do  case "$opt" in
                        h)  hlp
                        ;;
                        c)  control="$OPTARG"
                            control_defined=1
                        ;;
                        a)  automatic_defined=1
                        ;;
			m) multiple=1
	esac
    done


if  ls Del*_* 1> /dev/null 2>&1
	then	if [ $control_defined == 0 ] && [ $automatic_defined == 0 ]
    			then hlp
		fi
	for folder in Del*
		do
			del_number=`echo $folder | tr '_' "\n" | head -n1` #Identify which deletion strain we are working on
			if  [ $control_defined == 1 ]
				then control_samples=$control #WT control is ERS1076728
				elif [ $automatic_defined == 1 ]
						ploidy=`cat name\ conversion.tsv | grep "${del_number}_" | grep C[0123456789] | awk -F"\t" '{print $7}' | sed -e 's/[[:space:]]*$//' | head -n1 `
				  		then if [ $multiple == 1 ]
							then    control_samples=`cat name\ conversion.tsv | grep "${del_number}_" | grep C[0123456789] | awk -F"\t" '$6 ~ ERS {print $6}' | sed -e 's/[[:space:]]*$//' | tr "\n" ',' | sed 's/,*$//' `
							else	control_samples=`cat name\ conversion.tsv | grep "${del_number}_" | grep C[0123456789] | awk -F"\t" '$6 ~ ERS {print $6}' | sed -e 's/[[:space:]]*$//' | head -n1    | sed 's/,*$//'`
					    	fi
			fi

			echo "Analysing $folder"
			echo "Control(s): $control_samples"
		
			cd $folder
			if [[ -a bams_for_mpileup ]]
				then 
					if [[ $ploidy > 0 || $ploidy < 5 ]]
						then
							if [ $multiple == 1 ]
				 				then export SBATCH_CMD_CIAO="analyser-multi.sh -a -r -n${ploidy} -F -C $control_samples > results.txt 2>&1" 
								else export SBATCH_CMD_CIAO="analyser-multi.sh -a -r -n${ploidy} -F -c $control_samples > results.txt 2>&1"
				     			fi
						else
							if [ $multiple == 1 ]
                                                	        then export SBATCH_CMD_CIAO="analyser-multi.sh -a -r -n1 -F -C $control_samples > results.txt 2>&1"
                                                        	else export SBATCH_CMD_CIAO="analyser-multi.sh -a -r -n1 -F -c $control_samples > results.txt 2>&1"
                                               	 	fi
					fi
			fi
			sbatch ${DIR}/submit_sbatch.sh
			echo "Command: $SBATCH_CMD_CIAO"
			export SBATCH_CMD_CIAO=""
			cd ..
			sleep 2
		done

else printf "Please create a correct folder structure\n"
fi



 #       del_number=`echo $folder | tr '_' "\n" | head -n1` #Identify which deletion strain we are working on
                #       Uncoment the next line to use any strain marked as C in the plate field as control
                #       control_samples=`cat name\ conversion.tsv | grep "${del_number}_" | grep C[0123456789] | awk -F"\t" '$6 ~ ERS {print $6}' | sed -e 's/[[:space:]]*$//' | tr "\n" ',' | sed 's/,*$//' `  #identify the ERS numbers of the control samples.
                #       Uncomment the next line to use only the *first* strain  marked as C in the plate field as control
                #       control_samples=`cat name\ conversion.tsv | grep "${del_number}_" | grep C[0123456789] | awk -F"\t" '$6 ~ ERS {print $6}' | sed -e 's/[[:space:]]*$//' | head -n1 | tr "\n" ',' | sed 's/,*$//' ` #added head -n1 to only get one control #identify the ERS numbers of the control samples.
                #       Uncomment the next line to always use the same wild-type control for everything
                #       control_samples='ERS1076728'
