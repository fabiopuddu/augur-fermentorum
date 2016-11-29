#!/bin/bash

BEDTOOLSDIR="/software/vertres/bin-external"
IGVBED="$BEDTOOLSDIR/bedtools igv"

outdir=$1
indir=$2

echo "Generating IGV batch scripts"
for SAMPLE in `ls $indir`
do 
	SAMPLENAME=`basename $SAMPLE`
	$IGVBED -i $indir/$SAMPLE -path $outdir/$SAMPLENAME -sort position -slop 50 > $outdir/$SAMPLENAME.igvscript.txt
done