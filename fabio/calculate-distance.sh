#!/bin/sh
combine_first=0
while getopts "he:c" opt
    do  case "$opt" in
                        h)  printf "############   HELP   ###############\nOPTIONS\n"
                            printf "\t-h\tThis Help\n\t-f\tForce Rewrite\n\t-e\tName of the input experiment (if not specified it will be \"all\")"
                            exit 0
                        ;;
                        e) input_file=$OPTARG
                        ;; 
                        c) combine_first=1
                        ;;                        
        esac
    done
#collate all the sorted/intersected ERS files from the subfolders A and B corresponding to the two propagation lines    
mkdir $input_file
cp $input_file\_A/all/sort.ERS*.isec.vcf $input_file
cp $input_file\_B/all/sort.ERS*.isec.vcf $input_file
cd $input_file
strain=`echo $PWD | tr '/' '\t'| rev | cut -f1 | rev`
echo $n
rm *.dist
rm *.dat
#for every file, convert the relative coordinates to absolute coordinates
for file in sort.ERS*.isec.vcf
    do  n=`echo $file | tr -d '.isec.vcf'`
        n=`echo $n | tr -d 'sort.'`
        echo $n
        cat $file | grep '##' -v | grep "CHROM" -v | while read line
                    do  if [[ $line =~ '1/1' || $line =~ '0/1' ]]
                            then    chr=`echo $line | cut -d ' ' -f1`
                                    relpos=`echo $line | cut -d ' ' -f2`
                                    pos=0
                                    if [[ $chr == 'I' ]] #size 230,218
                                        then pos=$relpos
                                    elif [[ $chr == 'II' ]] #size 813,184
                                        then pos=$(($relpos+230218))
                                    elif [[ $chr == 'III' ]] #size 316,620
                                        then pos=$(($relpos+1043402))
                                    elif [[ $chr == 'IV' ]] #size 1,531,933
                                        then pos=$(($relpos+1360022))
                                    elif [[ $chr == 'V' ]] #size 576,874
                                        then pos=$(($relpos+2891955))
                                    elif [[ $chr == "VI" ]] #size 270,161
                                        then pos=$(($relpos+3468829))
                                    elif [[ $chr == "VII" ]] #size 1,090,940
                                        then pos=$(($relpos+3738990))
                                    elif [[ $chr == "VIII" ]] #size 562,643
                                        then pos=$(($relpos+4829930))
                                    elif [[ $chr == "IX" ]] #size 439,888
                                        then pos=$(($relpos+5392573))
                                    elif [[ $chr == "X" ]] #size 745,751
                                        then pos=$(($relpos+5832461))
                                    elif [[ $chr == "XI" ]] #size 666,816
                                        then pos=$(($relpos+6578212))
                                    elif [[ $chr == "XII" ]] #size 1,078,177
                                        then pos=$(($relpos+7245028))
                                    elif [[ $chr == "XIII" ]] #size 924,431
                                        then pos=$(($relpos+8323205))
                                    elif [[ $chr == "XIV" ]] #size 784,333
                                        then pos=$(($relpos+9247636))
                                    elif [[ $chr == "XV" ]] #size 1,091,291
                                        then pos=$(($relpos+10031969))
                                    elif [[ $chr == "XVI" ]] #size 948,066
                                        then pos=$(($relpos+11123260))
                                        else pos='error'
                                    fi
                                    if [[ $line =~ 'INDEL' ]] # if the line contains an indel, add a i to the line to mark it
                                        then pos="$pos i"
                                    fi
                                    echo $pos >> $n.dist
                        fi
            done
        echo $n.dist
        sort -n $n.dist > sorted.$n.dist #sort the numbers
        mv sorted.$n.dist $n.dist
        if [[ $combine_first == 0 ]] # only calculate the distances inbetween each file if combine_first option is not set
       		old_position=0
        	then cat $n.dist | while read line
            		do	if [[ ! $line =~ 'i' ]]
                    			then    pos=$line
                            			distance=$(($pos-$old_position))
                            			echo "$pos\t$distance\t0" >> $n.dist.dat
                    			else    pos=`echo $line | tr -d ' i'`
                            			distance=$(($pos-$old_position))
                            			echo "$pos\t0\t$distance" >> $n.dist.dat
                		fi
                		old_position=$pos
         			done
		fi
    done

if [[ $combine_first == 1 ]]
	then 	old_position=0
			cat *.dist > experiment.dist	#combine all the files with mutation position
			sort -n experiment.dist > sorted.experiment.dist #sort the numbers
        	mv sorted.experiment.dist experiment.dist
			cat experiment.dist | while read line
            		do	if [[ ! $line =~ 'i' ]]
                    			then    pos=$line
                            			distance=$(($pos-$old_position))
                            			echo "$pos\t$distance\t0" >> experiment.dist.dat
                    			else    pos=`echo $line | tr -d ' i'`
                            			distance=$(($pos-$old_position))
                            			echo "$pos\t0\t$distance" >> experiment.dist.dat
                		fi
                		old_position=$pos
                	done
mv experiment.dist.dat distances.dat
fi
if [[ $combine_first == 0 ]]
	then cat ERS*.dist.dat > distances.dat # combine all the distances files
fi
echo 'Plotting'
if [[ $combine_first == 0 ]]
	then gnuplot /Applications/premiata-forneria/fabio/raindrop-plot.gpl > "$strain.raindrop.png" # combine all the distances files
	else gnuplot /Applications/premiata-forneria/fabio/raindrop-plot.gpl > "$strain.combine_fist.raindrop.png"
fi
echo 'Done'
cd ..
 
#list=''
#for file in ERS*.dist.dat
#    do  list="$list\"$file\","
#    done
#echo $list
#sed "12s|.*|plot $list using 1:2 title \"SNV\"|" /Applications/WGS_software/fabio/raindrop-plot > raindrop-plot

#  Script.sh
#  
#
#  Created by Fabio on 03/02/2015.
#
