#!/bin/bash
OPTIND=1         # Reset in case getopts has been used previously in the shell.
caller=""
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
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
                export SBATCH_CMD_CIAO="variant_caller.sh -k $caller" #the k option enforce a check of the deletion before calling the mutations
            else
                printf "Wild type sample\n"
                export SBATCH_CMD_CIAO="variant_caller.sh $caller"      # No deletion check for wild type samples
        fi
	sbatch ${DIR}/submit_sbatch.sh
        echo "Command: $SBATCH_CMD_CIAO"
	printf "Done\n"
	export SBATCH_CMD_CIAO=""
        cd ..
	sleep 1
    done
#  experiment-caller.sh
#  
#
#  Created by Fabio on 13/02/2015.
#
