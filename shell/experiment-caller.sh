#!/bin/sh
for directory in Del*
    do  cd $directory
    	printf "Processing $directory..."
        wt_name=`echo $directory | grep "WT-"`
        echo $wt_name
        if [ -z $wt_name ]
            then
                printf "Mutant sample\n"
                variant_caller.sh -k  #the k option enforce a check of the deletion before calling the mutations
            else
                printf "Wild type sample\n"
                variant_caller.sh       # No deletion check for wild type samples
        fi
	printf "Done\n"
        cd ..
    done
#  experiment-caller.sh
#  
#
#  Created by Fabio on 13/02/2015.
#
