#!/bin/sh
while getopts "hs:e:" opt
    do  case "$opt" in
                        h)  printf "############   HELP   ###############\nOPTIONS\n"
                            printf "\t-h\tThis Help\n\t-s\tstart: initial SC_MFY folder\n\t-e\tend: final SC_MFY folder"
                            exit 0
                        ;;
                        s) start=$OPTARG
                        ;; 
                        e) fine=$OPTARG
                        ;;                        
        esac
    done
num=$(($fine-$start))
echo $range
for x in $(seq 0 $num)
	do  	n=$(($start+$x))
		echo $n
		#ln -s "/Volumes/LuCia/BAM\ files/Yeast/SC_MFY$n/" "$PWD/SC_MFY$n "
		ln -s /Volumes/LuCia/BAM\ files/Yeast/SC_MFY$n/ $PWD/SC_MFY$n
	done

#  Script.sh
#
#
#  Created by Fabio on 25/01/2015.
#
