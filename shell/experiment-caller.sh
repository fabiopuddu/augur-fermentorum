#!/bin/sh
for directory in Del*
    do  cd $directory
    	printf "Processing $directory..."
          variant_caller.sh -k  #the k option enforce a check of the deletion before calling the mutations
	printf "Done\n"
        cd ..
    done
#  experiment-caller.sh
#  
#
#  Created by Fabio on 13/02/2015.
#
