printf "Sample name\tchr01\tchr02\tchr03\tchr04\tchr05\tchr06\tchr07\tchr08\tchr09\tchr10\tchr11\tchr12\tchr13\tchr14\tchr15\tchr16\n"
for file in */ploidy_data/*_plstats.txt
	do
	fn=`echo $file | tr '/' "\n" | tail -n1 | sed 's/_plstats.txt//g'`
	printf "$fn\t"
	for chr in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 chr13 chr14 chr15 chr16
		do
		cn=`cat $file | grep $chr | cut -f2`
		printf -- "$cn\t"
		done
	printf "\n"
	done

