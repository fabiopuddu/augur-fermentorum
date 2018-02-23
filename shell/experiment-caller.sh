#!/bin/bash
OPTIND=1         # Reset in case getopts has been used previously in the shell.
caller=""

while getopts "SF" opt
    do  case "$opt" in
            S)
            caller='-S'
            ;;
            F)
            caller='-F'
            ;;
esac
done



for directory in Del*
    do  cd $directory
    	printf "Processing $directory..."
        wt_name=`echo $directory | grep "WT-"`
        echo $wt_name
        if [ -z $wt_name ]
            then
                printf "Mutant sample\n"
                sbatch --partition=LONG --wrap="variant_caller.sh -k $caller" #the k option enforce a check of the deletion before calling the mutations
            else
                printf "Wild type sample\n"
                sbatch --partition=LONG --wrap="variant_caller.sh $caller"      # No deletion check for wild type samples
        fi
	printf "Done\n"
        cd ..
	sleep 5
    done
#  experiment-caller.sh
#  
#
#  Created by Fabio on 13/02/2015.
#
