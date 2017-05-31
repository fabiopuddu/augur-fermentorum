outfile=$1
function gather_rep(){
cat name\ conversion.tsv | while read line
	do SD=`echo "$line" | cut -f5`
	   delname=`echo "$line" | cut -f2`
	   rep_group=`cat $delname/repDNA/$SD.txt | grep -v Sample`
	   telo=`cat $delname/repDNA/$SD.tel`
	   ERS=`echo "$line" | cut -f7`
	   printf "$rep_group\t$telo\t$ERS\t$delname\n"
	done 
}
gather_rep | grep "SD...." > $1
