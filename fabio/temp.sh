mkdir repetitive
for dir in Del*
	do  gene_name=`echo $dir | tr '/' "\n" | grep Del | sed 's/Del[0123456789]*_//g' | tr '[:lower:]' '[:upper:]'` 
	      cat $dir/bams_for_mpileup | while read line
	        do
		name=`echo $line | grep -o "SD.*\." | tr -d '.'`
		#printf "Gene: $gene_name\tProcessing...$name\t"
		#detect_deletion_chr_region.pl $gene_name $line | tr "\t" "\n" | grep Deleted: | sed 's/Deleted://g'
		command="rDNA-cnv_estimate.pl -i $dir/$line -p 2 > repetitive/$name.txt"
		sbatch  --wrap="$command"
		printf "$command\n"
	        done
		sleep 5
	done
		
