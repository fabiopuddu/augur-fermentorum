#!/bin/bash

for directory in Del*
    do  cd $directory
        printf "Processing $directory..."
	sbatch --partition=LONG  --wrap="find_deletion_tag.pl"
	cd ..
    done
