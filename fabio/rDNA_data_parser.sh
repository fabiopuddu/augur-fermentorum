#!/bin/bash
for folder in Del*
		do 	sum_length=0
			counter=0
			suffix=0
			if [[ -a $folder/rDNA/results.txt ]]
				then	cat $folder/rDNA/results.txt | { while read line	
							do	suffix=$(($suffix+1))
                    			declare "length_$suffix"=`echo $line | awk '{print $2}'`
							done
						until [[ $counter == $suffix ]]
							do 	counter=$(($counter + 1))
								l=length_$counter
								sum_length=`echo "scale=1; $sum_length + ${!l}" | bc`
							done	
						average=`echo "scale=1; $sum_length / $counter" | bc `
						counter=0
						sumdif=0
						until [[ $counter == $suffix ]]
							do 	counter=$(($counter + 1))
								l=length_$counter
								dif=`echo "scale=1; (${!l} - $average)^2" | bc`
								sumdif=`echo "scale=1; $sumdif + $dif" | bc`
							done
						
						stdev=`echo "scale=1; sqrt($sumdif / $counter)" | bc`
						printf "$folder\t$average\t$stdev\n"
						
			}
			fi		
		done
		
		
		
