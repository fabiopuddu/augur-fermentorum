#/bin/sh
#
cat Human_GeCKOv2_Library_*_09Mar2015.csv | parallel 'line={}
												seq=`echo "$line" | cut -f3`
												gene=`echo "$line" | cut -f1`
												cov=`LC_ALL=C fgrep $seq index3.fq | wc -l`
												printf "$gene\t$cov\t$seq\n" >> results_index_3_ter.txt
												'												
												