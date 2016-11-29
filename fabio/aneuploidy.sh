#!/bin/sh
mkdir aneuploidy
cd aneuploidy

# "gnuplot -e \"as=$arrow_start[$i];ae=$arrow_end[$i];s=$start[$i];e=$end[$i];of=\'$name[$i].pdf\';\" plot.gpl" );


function aneup(){
	line="$1"
	ploidy="$2"
	strain=`echo $line | tr '/' "\n" | tail -n1 | tr '.' "\n" | head -n1`
	mkdir $strain
	cd $strain
	CGH.pl -i ../../$line -p $ploidy
	gnuplot -e "of='../plot_${strain}.png'" /Applications/PF/fabio/aneuploidy_plot.gpl
	cd ..
}
export -f aneup		

plo=$1

export -f aneup	
export plo	

cat ../bams_for_mpileup | parallel 'aneup {} $plo';
cd ..
	
