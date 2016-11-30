function realign_rdna()
{
if [[ $v == "1" ]]; then echo "Processing... ${1}"; fi
name=`echo ${1} | grep -o "SC_MFY.......\|SD......"| sed "s|\.||g" | head -n1`
if [[ ! -a $name.rDNA ]]; then	rDNA_cov_extract.pl -i ../$line | sort -n -k1 >  $name.rDNA; fi
}

function rdna_extract()
{
name=`echo ${1} | grep -o "SC_MFY.......\|SD......"| sed 's|\.||g' | head -n1`
rDNA_repeat_estimate.pl -i $name.rDNA >> results.txt
ERSnum=`cat "../../name conversion.tsv" | grep -w $name | cut -f6`
Delname=`cat "../../name conversion.tsv" | grep -w $name | cut -f2`
printf "$ERSnum\t$Delname\n" >> results.txt
}