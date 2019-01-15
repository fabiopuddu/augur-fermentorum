#!/bin/bash
OPTIND=1         # Reset in case getopts has been used previously in the shell.
caller=""
delcheck=0
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
while getopts "SFhnk" opt
    do  case "$opt" in
	   h)  printf "############   HELP   ###############\nOPTIONS\n"
                            printf "\t-h\tThis Help\n"
                            printf "\t-k\tPerform gene KO check in folders\n"
                            printf "\t-n\tDo not perform gene KO check in folders\n"
                            printf "\t-S\tUse Samtools as variant caller\n"
			    printf "\t-S\tUse Freebayes as variant caller\n"
                            exit 0
	    ;;
	    k)
	    delcheck='yes';
	    ;;
	    n)
	    delcheck='no';
	    ;;
            S)
            caller='-S'
            ;;
            F)
            caller='-F'
            ;;
esac
done
if [ $delcheck == 0 ]
	then experiment-caller.sh -h  
	     exit 0
fi

for directory in Del*
    do  cd $directory
    	printf "Processing $directory..."
        wt_name=`echo $directory | grep "WT-"`
        echo $wt_name
        if [ -z $wt_name ]
            then
                printf "Mutant sample\n"
                if [ $delcheck == 'yes' ]
			then export SBATCH_CMD_CIAO="variant_caller.sh -k $caller" #the k option enforce a check of the deletion before calling the mutations
			else export SBATCH_CMD_CIAO="variant_caller.sh $caller"
		fi
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
