#!/bin/bash 
#  analyser-multi.sh
#  Created by Fabio on 23/11/2014.
#############################################
#  INITIALISE VARIABLES                        #
#############################################
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) #get the directory this program is stored in
control="0000"
ploidy=0
hetfilt=0
OPTIND=1         # Reset in case getopts has been used previously in the shell.
force_rewrite=0
low=0
multiple_ref=0
experiment_name=`echo $PWD | tr "/" "\n" | tail -n 1`
no_ctrl=0
v=0
show_syno=0
hard_force_rewrite=0
#############################################
#  FUNCTIONS                                #
#############################################
waitforcompletion(){
printf "Waiting for jobs $1  to complete"
list="ciao${1}"
finito=`squeue -u fp305 | grep "${list}" |wc -l | tr -d "\t"`
while [[ $finito != '0' ]]    
    do  finito=`squeue -u fp305 | grep "${list}" |wc -l | tr -d "\t"`
        printf ".${finito}"
        sleep 5
    done  
printf "\n"    
}
#############################################
#   GET   AND CHECK COMMAND LINE OPTIONS    #
#############################################
while getopts "sn:hfle:c:C:xvrtF" opt
    do  case "$opt" in
                        h)  printf "############   HELP   ###############\nOPTIONS\n"
                            printf "\t-h\tThis Help\n"    
                            printf "\t-f\tForce Rewrite\n"
                            printf "\t-x\tno control sample(an artificial control containing all the mutations shared between all the samples will be used\n"
                            printf "\t-c\tReference File\n"
                            printf "\t-C\tMultiple reference files: ERS numbers separated with a comma\n"
                            printf "\t-l\tmask only very very low quality variants\n"
                            printf "\t-n\tSet ploidy of samples (1 or 2)\n"
                            printf "\t-r\tanalyse rDNA in the samples\n"
                            printf "\t-t\ttransposons: analyse transposons in the samples\n"
                            printf "\t-v\tverbose: print out detailed analysis progression\n"
                            printf "\t-s\tsynonimous: print out detailed analysis progression\n"
                            exit 0
                        ;;
                        F)
                        hard_force_rewrite=1
                        ;;
                        s)
                        show_syno=1
                        ;;
                        c)
                        control="ERS$OPTARG"
                        ;;
                        n)
                        ploidy="$OPTARG"
                        ;;
                        l)
                        low=1
                        ;;
                        C)
                        control="$OPTARG"
                        multiple_ref=1
                        ;;
                        x)
                        no_ctrl=1
                        ;;
                        v)
                        v=1
                        ;;
                        r)
                        rDNA=1
                        ;;
                        t)
                        trnspn=1
                        ;;
                        
        esac
    done
if [[ $control == '0000' && $no_ctrl == '0' ]]
    then    echo "Please specify the ERS number of the reference sample with the flag -c or -C"
            exit 1
    else echo "Reference: $control"
fi
if [[ $ploidy == '0' ]]
    then    echo "Please specify the ploidy of the sample with the flag -n"
            exit 1
    else printf "Ploidy:   $ploidy"
         printf "n\n"
fi
if [[ ! -d calling ]]
    then printf "Run the variant calling first!"
         exit 1 
fi    

if [ ! "$(ls -A calling)" ]
    then printf "No vcf files\n"
         exit 1
fi         
# if [[ ${trnspn} == '1' && ! -d TR_BAMS ]]
#     then echo "Please run ty-realign.sh first"
#          exit 1   
# fi     
if [[ $multiple_ref == '1' ]]                                   ###in case of multiple reference files, 
    then    suffix=0                                                #cicle to the reference numbers and assign 
            control=$(echo $control | tr "," "\n")                    #each number to variables control_1, control_2, control_3, etc.

            for line in $control
                do  suffix=$(($suffix+1))
                    declare "control_$suffix"="ERS$line"
                done
fi
#
#################################################
#                                               #
#        PROGRAM STARTS                         #
#                                               #
#################################################
###################################
##CLEANUP IN CASE OF FORCE REWRITE
###################################
if [[ $hard_force_rewrite == '1' ]]
    then echo 'force rewrite' 
        rm -rf analysis
        rm -rf rDNA
        rm -rf transposons
fi
if [[ -d analysis ]] 
    then first_run=0
    else first_run=1
         mkdir analysis
fi    
if [[ -a analysis/csq.file ]]; then rm analysis/*.file; fi
######Before doing anything else, check that the sample names 
######in the vcf file are in the ERS format; if not change them appropriately
cd calling
if [[ $first_run == 1 ]]  
 then for f in ERS*.vcf.gz
        do headr=`zcat $f | grep "CHROM"  | cut -f10-`
            new_headr=''
            for tab in $headr
                do     if [[ ! $tab =~ "ERS" ]]
                        then     printf "Fixing headers...\n"
                        #name=`echo $tab | tr -d '..' 
                        ERSnumber=`grep -w $tab ../../name\ conversion.tsv | tr "\t" "\n" | tail -n 1`
                        new_headr=`printf "$new_headr$ERSnumber\t"`
                    else     new_headr=`printf "$new_headr$tab\t"`

                    fi
                done
            new_headr=`printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$new_headr"`
            gunzip -dc $f | sed "s/^.*CHROM.*$/$new_headr/" | sed '/##samtoolsCommand/d' | bgzip > new.$f  #the samtoolsCommand lines are removed                                                                                     
            rm $f                                                                                                #to allow bedtools intersect to work 
            mv new.$f $f                                                                                    #without complaining about mismatching columns
        done
fi        
###################################
######Remove mitochondrial mutation and move files to the analysis folder
###################################    
printf "Removing mitochondrial mutation: "
for f in ERS*.vcf.gz
        do if [[ ! -a ../analysis/mito.$f ]]
                then gunzip -dc $f | grep -v 'Mito' | bgzip > ../analysis/mito.$f 
                     printf    "\r\033[Kdone sample... $f"    
                else printf    "\r\033[Kskipping... $f "
            fi    
        done
        mito=0    
        printf "\n"        
cd ../analysis
####################################
###VARIANT EFFECT PREDICTION
####################################
proclist=''
ls mito.ERS*.vcf.gz | { while read line 
		do if [[ ! -a csq.$line ]]
				then command1="variant_effect_predictor.pl --species saccharomyces_cerevisiae -i $line  --format vcf -o vep.$line.txt --no_progress --force_overwrite --offline"
					 command2="vcf2consequences_vep -v $line -i vep.$line.txt 2>/dev/null | bgzip > csq.$line"
					 PROC1=$(sbatch --partition=LONG --wrap="${command1}" | sed 's/Submitted batch job //g') 
					 PROC2=$(sbatch --partition=LONG --dependency=afterok:${PROC1} --wrap="${command2}" | sed 's/Submitted batch job //g')
		   			 proclist="${proclist}\|${PROC2}"
		   fi
		done
waitforcompletion "${proclist}" 
}
######################################
##        ANALYSE RDNA LOCUS         ##
######################################
if [[ $rDNA == 1 ]] 
    then    printf "Analysing rDNA...\n"
            mkdir -p ../rDNA
        	cd ../rDNA
        	proclist=''
            cat ../bams_for_mpileup |{ while read line
            	do  name=`echo $line | grep -o "SC_MFY.......\|SD......"| sed "s|\.||g" | head -n1`
                    command1="rDNA_cov_extract.pl -i ../$line | sort -n -k1 >  $name.rDNA"
                    command2="rDNA_repeat_estimate.pl -i $name.rDNA >> $name.txt"
                    PROC1=$(sbatch --partition=LONG --wrap="${command1}" | sed 's/Submitted batch job //g') 
                    PROC2=$(sbatch --partition=LONG --dependency=afterok:${PROC1} --wrap="${command2}" | sed 's/Submitted batch job //g')
                    proclist="${proclist}\|${PROC2}"	 
                done
                waitforcompletion "${proclist}"
                }                
                ls *.txt | while read line
                       do 	result=`cat $line`
		                	name=`echo $line | sed 's|\.txt||g'`
                        	ERSnum=`cat "../../name conversion.tsv" | grep -w $name | cut -f6`
                        	Delname=`cat "../../name conversion.tsv" | grep -w $name | cut -f2`
                        	printf "$result\t$ERSnum\t$Delname\n" >> results.txt
                    	done
                cd ../analysis
fi
sleep 5
######################################
##        ANALYSE TRANSPOSONS         ##
######################################
##trnspn=1
if [[ $trnspn == 1 ]] 
    then	printf "Analysing Transposons...    "
        	cd ..
        	#ty-realign.sh
        	mkdir -p transposons
			cd transposons
            cat ../bams_for_mpileup | sed 's|.bam|.bam.Ty.bam|g' | while read line
                             do name=`echo ${line}| tr "/" "\n" | tail -n1`;
                                SC=`echo ${name} | sed 's/.bam.Ty.bam//g' | sed 's/.merged//g' `
                                ERSnum=`cat "../../name conversion.tsv" | grep -w $SC | cut -f6`;
                                command="Ty1-5Estimator.pl -i ../TR_BAMS/$name > ${ERSnum}.ty"
                                sbatch --partition=LONG -e slurm.%N.%j.err --wrap="${command}"
                             done    
			waitforcompletion
			sleep 15 
		 	ls ERS*.ty | while read line
             	do	ERSnum=$(echo ${line} | sed 's/\.ty//g');
                	ty=`cat $line`
                	printf "$ERSnum\t$ty\n" >> results.txt
             	done        
    		cd ../analysis
fi         
#####################################
#### START MUTATION ANALYSIS     #####
#####################################
###################################
####NORMALISATION AND PRELIMINARY FILTERS
####################################
#NORMALIZATION OF INDEL ANNOTATION
sleep 10
printf "Normalising indels.... \n"
for f in csq.*.vcf.gz 
        do     if [[ ! -a norm.$f ]]
                then printf "Normalising $f ..." 
                     n=`echo $f | sed 's|csq.mito.||g' | sed 's|.vcf.gz||g'`
                     bcftools norm -f $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa $f | bgzip > $n.vcf.gz
                else printf "\n Skipping $f"
            fi                     
        done         
printf 'Sorting mutations...'
isec_files=(*isec.vcf)
    if [[ ! -e "${isec_files[0]}" || "$force_rewrite" == "1" ]]
      then for x in ERS*.vcf.gz
        do  n=$(echo $x | sed 's/.vcf.gz//g')
            zcat $x | vcf-sort > sort.$n.vcf 2>/dev/null
            cp sort.$n.vcf sorted.$n.vcf
            bgzip -f sort.$n.vcf
            tabix -f -p vcf sort.$n.vcf.gz
        done
       printf "done\n"
       else printf "Files already sorted: skipping\n"
     fi
############################################################################################################################################################
#
#               COMBINE ALL THE CONTROLS IF THE NUMBER OF CONTROLS IS >1s   sort, bgzip and tabix the control files
#
#############################################################################################################################################################
if [[ $multiple_ref == '1' ]]
    then    if [[ $v == '1' ]] ; then echo '......Combining controls: ' ;fi
            counter=1
            file1=control_$counter
            if [[ -a sort.${!file1}.vcf.gz ]] #make sure all the control files exists before trying to merge them
                then     cp sort.${!file1}.vcf.gz sort.control.vcf.gz
                        until [[ $counter == $suffix ]] # cycle through each control sample (control_$suffix), combine the last combined sample (which is equivalent
                            do  counter=$(($counter+1)) #  to control_1 in the first cycle) with control_$counter. Append the intersection to the control_$counter sample, then sort, bgzip and tabix
                                file2=control_$counter
                                if [[ -a sort.${!file2}.vcf.gz ]] #make sure all the control files exists before trying to merge them
                                    then     bedtools intersect -a sort.control.vcf.gz -b sort.${!file2}.vcf.gz -v > isec.control.vcf
                                            cat sorted.${!file2}.vcf isec.control.vcf > combined.control.vcf
                                            cat combined.control.vcf | vcf-sort > sort.control.vcf
                                            cp sort.control.vcf control.vcf
                                            bgzip -f sort.control.vcf
                                            tabix -f -p vcf sort.control.vcf.gz
                                fi
                            done
            fi                
    #if multiple_ref is not 1 then either there is one control or no control was specified
    #so if no_ctrl = 0 then only one control was specified and its file should be renamed as sort.control.vcf.gz
    else     if [[  $no_ctrl == '0' ]] 
                then    cp sort.$control.vcf.gz sort.control.vcf.gz
                        cp sort.$control.vcf.gz.tbi sort.control.vcf.gz.tbi
            fi
fi
# if the process of combining multiple references failed because either the first control was missing or the first was present but the second was missing, 
# no sort.control.vcf.gz file will have been created, therefore if the file is not there the missing_control flag is set to 1 and the control variable is
# set to the same value that it has when no control has been specified
if [[ ! -a sort.control.vcf.gz ]] 
    then     missing_control='1'
            control="0000"
fi
############################################################################################################################################################
#
#              NO CONTROL SPECIFIED
#
#############################################################################################################################################################
# the fake control is produced by going through all the samples and intersecting every sample with the previous one
# so that at the end it will contain only the mutations shared by all samples
if [[ $no_ctrl == '1' || $missing_control == '1' ]] #if no control was specified or if the specified control files do not exist
    then    if [[ $v == '1' ]] ; then printf "no control\n"; fi
            pattern="sort.ERS*.vcf.gz"
            files=( $pattern )
            cp "${files[0]}" fake.control.vcf.gz
            for x in sort.ERS*.vcf.gz
                do  n=$(echo $x | sed 's/.vcf.gz//g' | sed 's|^sort.||g')
                    #vcf-isec -f -a -c
                    bedtools intersect -header -a fake.control.vcf.gz -b sort.$n.vcf.gz > fake.control.vcf
                    cat fake.control.vcf | vcf-sort > sort.fake.control.vcf
                    bgzip -f sort.fake.control.vcf
                    tabix -f -p vcf sort.fake.control.vcf.gz
                done               
            cp fake.control.vcf filt.$control.vcf
fi
# if 
# if [[ ! $multiple_ref == '1' || $missing_control == '1' ]]    
#     then     cp filt.$control.vcf control.vcf
#             cat control.vcf | vcf-sort > sort.control.vcf
#             bgzip -f sort.control.vcf
#             tabix -f -p vcf sort.control.vcf.gz
# fi
#############################################################################################################################################################
printf "Filtering mutations with mapping quality <40 and SnpGap=7"
if [[ ! -a already_filtered || $force_rewrite == '1' ]] 
    then     for x in sort.ERS*vcf.gz
                do  n=$(echo $x | sed 's/.gz//g')
                    zcat $x | vcf-annotate -H -f +/d=15/q=25/SnpGap=7 > filt.$x ######!!!!!!!! FILTER Q>30, SNPGAP<7
                    mv filt.$x $x
                    printf '.'
                done
            printf "Finished\n" 
            touch already_filtered
fi    
printf "All SNV in the control sample: "
zcat  sort.control.vcf.gz | grep '##' -v  | grep '#CHROM' -v | grep 'INDEL' -v | wc -l
############################################################################################################################################################
#
#               INTERSECTIONS
#
#############################################################################################################################################################
#DIRECT INTERSECTION
echo 'Intersecting samples with control: '$control
    for x in sort.ERS*vcf.gz
                do  n=$(echo $x | sed 's/.vcf.gz//g')
                    if [[ $v == '1' ]] ; then printf "\r\033[K                       Intersecting $n";fi
                    bedtools intersect -header -a $x -b sort.control.vcf.gz -v > $n.isec.vcf #actual intersect command
                    if [[ $low = 1 ]] #filtering mutations from non-control samples based on quality
                            then cat $n.isec.vcf  | vcf-annotate -f $DIR/../mareike/gt-filter-lax.pl > $n.isec.filt.vcf #lax filter
                            else cat $n.isec.vcf  | vcf-annotate -f $DIR/../mareike/gt-filter.pl > $n.isec.filt.vcf #strict filter
                    fi
                    mv $n.isec.filt.vcf $n.isec.vcf
                done
#INVERSE INTERSECTION
if [[ ( ! -a inverse_intersection || "$force_rewrite" == "1" ) ]]
    then    printf "\n Performing inverse intersections"
            mkdir inverse_intersection
            cd inverse_intersection
            for x in ../sort.ERS*.vcf.gz
                do  n=$(echo $x | sed 's/.vcf.gz//g')
                    n=$(echo $n | sed 's/..\///g')
                    if [[ $v == '1' ]] ; then printf "\r\033[K ........Inverse intersection: Intersecting control(s) with sample $n" ;fi
                    #vcf-isec -f -a -c ../sort.control.vcf.gz $x 2>/dev/null > $n.isec.vcf
                    bedtools intersect -header -a ../sort.control.vcf.gz -b $x -v > $n.isec.vcf #actual intersect command
                    if [[ $low = 1 ]] #filtering mutations from non-control samples based on quality
                        then cat $n.isec.vcf  | vcf-annotate -f $DIR/../mareike/gt-filter-lax.pl > $n.isec.filt.vcf #lax filter
                        else cat $n.isec.vcf  | vcf-annotate -f $DIR/../mareike/gt-filter.pl > $n.isec.filt.vcf #strict filter
                    fi
                    mv $n.isec.filt.vcf $n.isec.vcf
                done
            cd ..
fi
########
#MASKED INTERSECTIONS
#######
#MASKED DIRECT INTERSECTION
#create a control reference with masked hets: this allows for detection of LOH events (reference: 0/1, sample 1/1)
if [[ $ploidy == '2' ]]
    then    zcat sort.control.vcf.gz | vcf-annotate -f $DIR/../mareike/mask-hets.pl > combined.control.masked.vcf
            grep -Ev "\./\." combined.control.masked.vcf > combined.control.masked2.vcf #These two lines are needed to remove ./. from the masked control file, otherwise they will be intersected
            mv combined.control.masked2.vcf combined.control.masked.vcf                 # and the corresponding mutations in the sample files will be removed leading to  incorrect numbers of LOH TO ALT.
#vcf-subset -c $control combined.control.masked.vcf -e > control.masked.vcf.subset
#           mv control.masked.vcf.subset control.masked.vcf
            if [[ $v == '1' ]] ; then printf "Of which homozygous: " ;fi
            if [[ $v == '1' ]] ; then cat  combined.control.masked.vcf | grep '##' -v  | grep '#CHROM' -v | grep 'INDEL' -v | wc -l ;fi
fi
#intersection with masked het-masked reference
if [[ ( ! -a intersect_masked || "$force_rewrite" == "1" ) && $ploidy == "2" ]]
    then    printf "\n"
            mkdir intersect_masked
            cd intersect_masked
            cat ../combined.control.masked.vcf | vcf-sort > sort.control.masked.vcf
            bgzip -f sort.control.masked.vcf
            tabix -f -p vcf sort.control.masked.vcf.gz
            echo '......Intersecting samples with masked control: '$control
            for x in ../sort.ERS*vcf.gz
                do  n=$(echo $x | sed 's/.vcf.gz//g')
                    n=$(echo $n | sed 's/..\///g')
                    if [[ $v == '1' ]] ; then printf "\r\033[K                       Intersecting $n";fi
                    bedtools intersect -header -a $x -b sort.control.masked.vcf.gz -v > $n.isec.vcf
                    #vcf-isec -f -a -c $x sort.control.masked.vcf.gz 2>/dev/null > $n.isec.vcf #actual intersect command
                    if [[ $low = 1 ]] #filtering mutations from non-control samples based on quality
                            then cat $n.isec.vcf  | vcf-annotate -f $DIR/../mareike/gt-filter-lax.pl > $n.isec.filt.vcf #lax filter
                            else cat $n.isec.vcf  | vcf-annotate -f $DIR/../mareike/gt-filter.pl > $n.isec.filt.vcf #strict filter
                    fi
                    mv $n.isec.filt.vcf $n.isec.vcf
                done
                cd ..
fi
#TO BE USED WHEN INVERSE INTERSECTING WITH MASKED REF
#for x in ../ERS*.vcf
#do  cat $x | vcf-annotate -f $DIR/../mareike/mask-hets.pl > masked.$x #mask all the hets in the samples but not in the control
#cat masked.$x| vcf-sort > sort.masked.$x
#bgzip -f sort.masked.$x
#tabix -f -p vcf sort.masked.$x.gz
#done
#############################################################################################################################################################
##merge all intersections if merged file does not exist or if rewrite is forced.
if [[ ! -a experiment_merge.vcf || "$force_rewrite" == "1" ]]
     then   printf "\nMerging intersected files..."
            for x in sort.ERS*.isec.vcf
                do  n=$(echo $x | sed 's/.vcf//g')
                    cat $x | vcf-sort > sort.$n.vcf #sort
                    bgzip -f sort.$n.vcf #gzip
                    tabix -p vcf sort.$n.vcf.gz
                done
            list=''
            for x in sort.sort.ERS*.isec.vcf.gz #generate a list of all the files to be merged
               do  list=$list`echo "$x"`
                   list=$list' '
               done
            vcf-merge $list 2>/dev/null > experiment_merge.vcf #actual merge command
fi
printf "done\n"
#############################################################################################################################################################

##################################################
#                                                #
# CALCULATING AND DISPLAYING RESULTS             #
#                                                #
##################################################

##################################################
#                                                #
# MUTATION SUMMARY TABLE                         #
#                                                #
##################################################

#####
##HEADER
#####
number_of_samples=0
if [[ $v == '1' ]] ; then echo '...Finished';fi
if [[ $v == '1' ]] ; then     printf "\n\n\n" ;fi
    echo 'MUTATION SUMMARY'
    width=43
    headerdip2="\n %-10s %9s %12s %60s %7s %23s %4s %40s %12s %7s %3s %12s %10s "
    headerdip1="\n %-10s %10s %5s %4s \e[32m%5s \e[33m%9s \e[34m%7s \e[0m%12s %3s %7s %10s %7s %10s %10s %6s \e[32m%10s \e[0m%8s %10s %3s %5s %3s %8s %10s %6s %5s %6s %4s %5s\n"
    formatdip=" %-10s %8s %4s %4s \e[32m%5s \e[33m%8s \e[34m%9s \e[0m%9s %8s %4s %12s %8s %8s %10s %8s \e[32m%6s \e[0m%10s %8s %7s %9s %8s %5s %9s %10s %4s\n"
    headerhap="\n %-10s %4s %8s %6s \e[32m%5s \e[33m%9s \e[34m%7s \e[0m%12s %3s %4s %1s %6s %5s \e[32m%12s \e[0m%10s %12s %2s %6s %1s %2s %2s %5s %4s %4s %4s \n"
    formathap=" %-12s %8s %8s %7s \e[32m%5s \e[33m%8s \e[34m%9s \e[0m%9s %8s %4s %4s %5s %7s \e[32m%7s \e[0m%12s %11s %6s %4s %4s %5s %7s\n"
        if [[ $ploidy == '2' ]]
            then    printf "=========================================================================================================================================================================================================\n"
                    printf " $experiment_name \n"
                    printf "========================================================================================================================================================================================================="
                    printf "$headerdip2" "." "." "║" "SINGLE NUCLEOTIDE VARIANTS GAINED (OF WHICH HOMOZIGOUS)" "║" "LOSS OF HETEROZIGOSITY" "║" " INSERTIONS/DELETIONS GAINED" "║" "LOSS OF HETEROZIGOSITY" "║"
                    printf "$headerdip1" "ERS NO." "SAMPLE NAME" "REF" "║" "NONSENSE" "MISSENSE" "SENSE" "INTERGENIC" "│" "TO REF." "TOTAL SNV" "║" "TO ALT." "TO REF." "║" "FRAMESHIFT" "INFRAME" "INTERGENIC" "|" "TOT.INDEL(HOM)" "║" "TO ALT." "TO REF." "║"
                    printf "=========================================================================================================================================================================================================\n"
            else    printf "===================================================================================================================================================================\n"
                    printf " $experiment_name \n"
                    printf "==================================================================================================================================================================="
                    printf "$headerhap" "ERS NO." "SAMPLE NAME" "REF" "║" "NONSENSE" "MISSENSE" "SENSE" "INTERGENIC" "│" "TO REF." "|" "TOT.SNV" "║" "FRAMESHIFT" "INFRAME" "INTERGENIC" "|" "TO REF." "|" "TOT.INDEL" "║"
                    printf "===================================================================================================================================================================\n"

        fi
#########
##MUTATION SUMMARY TABLE
#########
        for x in sort.ERS*.isec.vcf
        do
            n=$(echo $x | sed 's/.isec.vcf//g' | sed 's/sort.//g')
            num=$(echo $n | sed 's/ERS//g')
            name=`cat ../../name\ conversion.tsv | grep $n | cut -f2`
            #count SNPs
            SNPtot=`grep '##' $x -v | grep '#CHROM' -v | grep 'INDEL' -v | grep '\./\.' -v | wc -l` #homo+hetero mutations from het-unmasked control
            SNPHOMREV=`grep '##' inverse_intersection/$x -v | grep '#CHROM' -v | grep 'INDEL' -v | grep "1/1\|2/2" | grep PASS | wc -l | tr -d ' '` #homozigous reversion to ref (i.e. control 1/1, sample 0/0)
            if [[ $ploidy == "2" ]]
                then     SNPhom=`grep '##' $x -v | grep '#CHROM' -v | grep 'INDEL' -v | grep "1/1\|2/2" | wc -l | tr -d ' '` #homo mutations from het-unmasked control
                        SNP="$SNPtot($SNPhom)"
                        SNPhomomask=`grep '##' intersect_masked/$x -v | grep '#CHROM' -v | grep 'INDEL' -v | grep "1/1\|2/2" | wc -l | tr -d ' '` #hetero mutations from het-masked control
                        SNPLOH1=$(( $SNPhomomask-$SNPhom )) #loss of heterozygosity towards 1/1
                        SNPLOH2=`grep '##' inverse_intersection/$x -v | grep '#CHROM' -v | grep 'INDEL' -v | grep "0/1" | wc -l | tr -d ' '` #loss of heterozygosity towards 0/0
                else    SNP="$(($SNPtot+$SNPHOMREV))"
            fi
        #count INDELs
            INDtot=`grep '##' $x -v | grep '#CHROM' -v | grep 'INDEL' | grep '\./\.' -v | wc -l` #homo+hetero indels from het-unmasked control
            INDHOMREV=`grep '##' inverse_intersection/$x -v | grep '#CHROM' -v | grep 'INDEL' | grep "1/1\|2/2" | grep PASS | wc -l | tr -d ' '` #homozigous reversion to ref (i.e. control 1/1, sample 0/0)
            if [[ $ploidy == "2" ]]
                then    INDhom=`grep '##' $x -v | grep '#CHROM' -v | grep 'INDEL' | grep "1/1\|2/2" | wc -l | tr -d ' '` #homo mutations from het-unmasked control
                        IND="$INDtot($INDhom)"
                        INDhomomask=`grep '##' intersect_masked/$x -v | grep '#CHROM' -v | grep 'INDEL' | grep "1/1\|2/2" | wc -l | tr -d ' '` #hetero indels from het-unmasked control
                        INDLOH1=$(( $INDhomomask-$INDhom )) #loss of heterozygosity towards 1/1
                        INDLOH2=`grep '##' inverse_intersection/$x -v | grep '#CHROM' -v | grep 'INDEL' | grep "0/1" | wc -l | tr -d ' '` #loss of heterozygosity towards 0/0
                else    IND="$INDtot"
            fi
    #count consequences
            perl $DIR/../mareike/counting_consequences.pl -i $x > csq.file
            STOP=`grep "stop_gained" csq.file | grep SNP | head -c 2 | tr -d '\t'`
            MISS=`grep "missense_variant" csq.file | grep SNP | head -c 2 | tr -d "\t"`
            FS=$((`grep "frameshift_variant" csq.file | grep INDEL | head -c 2 | tr -d '\t'` + `grep "stop_gained" csq.file | grep INDEL | head -c 2 | tr -d '\t'`))
            INFRAME=$((`grep "inframe_insertion" csq.file | grep INDEL | head -c 2 | tr -d '\t'` + `grep "inframe_deletion" csq.file | grep INDEL | head -c 2 | tr -d '\t'`))
            SENSE=`grep "synonymous_variant" csq.file | grep SNP | head -c 2 | tr -d '\t'`
            UP_DOWN_SNV=`grep '##' $x -v | grep '#CHROM' -v |grep 'INDEL' -v | grep stop_gained -v | grep missense -v | grep synonymous -v | grep "\./\." -v | wc -l`
            UP_DOWN_SNV_HOM=`grep '##' $x -v | grep '#CHROM' -v |grep 'INDEL' -v | grep stop_gained -v | grep missense -v | grep synonymous -v | grep "1/1\|2/2" | wc -l|tr -d ' '`
            UP_DOWN_SNV="$UP_DOWN_SNV($UP_DOWN_SNV_HOM)"    
            UP_DOWN_INDEL=`grep '##' $x -v | grep '#CHROM' -v  | grep 'INDEL' | grep inframe_deletion -v | grep inframe_insertion -v | grep frameshift -v | grep "\./\." -v | wc -l`
    #Display results
            if [[ "$control" =~ "$num"  ]]
                then
                    contr='+++'
                        if [[ $multiple_ref == '1' ]]
                            then    apex1='x'
                                    apex2='y'
                        fi
                else
                    contr='-'
                    number_of_samples=$(($number_of_samples+1))
                    apex1=''
                    apex2=''
            fi
            if [[ $ploidy == '1' ]]
                then printf "$formathap" \ $n $name $contr "║" $STOP $MISS $SENSE $UP_DOWN_SNV "│" $SNPHOMREV "|" $SNP "║" $FS $INFRAME $UP_DOWN_INDEL "|"  $INDHOMREV "|" $IND "║"
                else printf "$formatdip" \ $n $name $contr "║" $STOP $MISS $SENSE $UP_DOWN_SNV "│" $SNPHOMREV $SNP "║" $SNPLOH1$apex1 $SNPLOH2$apex2 "║" $FS $INFRAME $UP_DOWN_INDEL "|" $IND "║" $INDLOH1 $INDLOH2 "║" 
            fi
    done
if [[ $ploidy == '1' ]]
    then printf "===================================================================================================================================================================\n"
    else printf "=========================================================================================================================================================================================================\n"
fi
if [[ $multiple_ref == '1' && $ploidy == '2' ]]
    then    printf "x = number of 0/1 mutations in the combined control that are 1/1 in this control sample\n"
            printf "y = number of 0/1 mutations in the combined control that are 0/0 in this control sample\n"
fi

##convert results in gene list:
    ##You will need to amend this script. I will include a file called
    ##yeast_genelist_nameconv.tsv which contains the information for the more common names for genes.
    ##At "my $conversion_file=" you will need to specify the path to this file to make the script run well)
printf "\n"
perl $DIR/../mareike/vcf_to_gene_list.pl -i experiment_merge.vcf > hom.table.file
perl $DIR/../mareike/vcf_to_gene_list_het.pl -i experiment_merge.vcf > het.table.file

#############################################
#                                            #
#    REPETITIVE SEQUENCES TABLE                #
#                                            #
#############################################

##############
#HEADER
##############
number_of_samples=0
echo 'Repetitive Sequences'
width=43
headerhap="\n %-10s %8s %6s %5s %6s %5s %5s %5s %5s %5s %5s %5s \n"
formathap=" %-10s %12s %6s %5s %6.1f %5s %5s %5s %5s %5s %5s %5s\n"
printf "===============================================================================\n"
printf " $experiment_name \n"
printf "==============================================================================="
printf "$headerhap" "ERS NO." " SAMPLE NAME" "REF" "║" "rDNA" "║" "Ty1" "Ty2" "Ty3" "Ty4" "Ty5" "║" 
printf "===============================================================================\n"
    
#########
## TABLE
#########
        for x in sort*.isec.vcf
        do
            n=$(echo $x | sed 's/.isec.vcf//g' | sed 's/sort.//g')
            num=$(echo $n | sed 's/ERS//g')
            name=`cat ../../name\ conversion.tsv | grep $n | cut -f2`
            if [[ $rDNA == 1 ]] ; then RDNA=`cat ../rDNA/results.txt | grep "$n" | cut -f2`; else RDNA='-1'; fi; #if the rDNA option was set, go in the result file and grab the number
            if [[ $trnspn == 1 ]] 
                 then TY1=`cat ../transposons/results.txt | grep "$n" | cut -f3`; 
                       TY2=`cat ../transposons/results.txt | grep "$n" | cut -f4`
                       TY3=`cat ../transposons/results.txt | grep "$n" | cut -f5`
                       TY4=`cat ../transposons/results.txt | grep "$n" | cut -f6`
                       TY5=`cat ../transposons/results.txt | grep "$n" | cut -f7`
                 else TY1='-1'; fi; #if the transposon option was set, go in the result file and grab the number
            #Display results
            if [[ "$control" =~ "$num"  ]]
                then
                    contr='+++'
                        if [[ $multiple_ref == '1' ]]
                            then    apex1='x'
                                    apex2='y'
                        fi
                else
                    contr='-'
                    number_of_samples=$(($number_of_samples+1))
                    apex1=''
                    apex2=''
            fi
            printf "$formathap" \ $n $name $contr "║" $RDNA "║"  $TY1 $TY2 $TY3 $TY4 $TY5 "║" 
         done   
printf "===============================================================================\n"
#
#################################################
#                                               #
#  CALCULATING/DISPLAYING GENOTYPE TABLE        #
#                                               #
#################################################
#
##convert output into genotypes
#create genotype tables with both common names and systematic names; add 'NONAME' whenever a non standard gene name is found
#
#rad5=`zgrep "YLR032W:missense_variant:1603:535:G>R" sort.control.vcf.gz | wc -l | tr -d " "` #check for the presence of rad5-535
#mut_check='rad50S'
#mut=`zgrep "YNL250W:missense_variant:242:81:K>I" sort.control.vcf.gz | wc -l | tr -d " "` #check for the presence of the additional mutation
mut_check='sae2-F276S'
mut=`zgrep "YGL175C:missense_variant:827:276:F>S" sort.control.vcf.gz | wc -l | tr -d " "`




for x in ERS*.vcf.gz
    do  n=$(echo $x | sed 's/.vcf.gz//g')
        samp_name=`cat "../../name conversion.tsv" | grep -w $n | cut -f2`;    
        num=$(echo $n | sed 's/ERS//g')
        geno_common=`grep $n hom.table.file |grep Y.[LR]........... -o |  tr "\t" ";" | sed 's/Y[ABCDEFGHIJKLMNOP][LR][0123456789][0123456789][0123456789][WC]-*[ABCD]*;//g' | sed 's/;.*$//' | sed 's/""/NONAME/g' | sed 's/-[ABCDEF]""/NONAME/g' | tr '\n' '\t'`
        geno_systematic=`grep $n hom.table.file | grep Y.[LR].......... -o | tr "\t" ";" |  sed -e 's/(//g' |sed -e 's/[[:space:]][A-Z][A-Z][A-Z][0123456789]*//g' | sed 's/""//g' | sed 's/;.*$//g' | tr '\n' '\t'`
        geno_details=`cat hom.table.file | grep $n | grep -o £.*£ | sed 's/£//g' | tr '\n' '\t'`
        if [[ "$control" =~ "$num" ]]
            then    if [[ $rad5 != "0" ]] #print rad5-535 if the mutation is present
                        then add_gen="\e[31mrad5-535\e[0m"
                        else add_gen=''
                    fi
                    if [[ $mut != "0" ]] #print the mutation to be checked, if the mutation is present
                        then add_gen="$add_gen\t\e[96m$mut_check\e[0m"
                        else add_gen="$add_gen"
                    fi
                    printf "$n\t$samp_name\t++++\t$add_gen\n" >> common_name.file
            else    printf "$n\t$samp_name\t.\t$geno_common\n" >> common_name.file
        fi
        printf "$n\t$samp_name\t$geno_systematic\n" >>systematic_name.file
        printf "$n\t$samp_name\t$geno_details\n" >>mutation.file
    done
# identify the samples with known resistant genes
top1=`cat hom.table.file | grep TOP1 | grep -o ERS...... | tr '\n' '\t'`
pdr=`cat hom.table.file | grep PDR | grep -o ERS...... | tr '\n' '\t'` #pleiotropic drug resistance
#zds=`cat hom.table.file | grep ZDS | grep -o ERS...... | tr '\n' '\t'` #zillion different screens
#print headers
printf "\nDIFFERENTIAL GENOTYPE\n"
if [[ $ploidy == "2" ]]
        then printf "ERS NUMBER\tSAMPLE NAME\tREF\tHOMOZYGOUS GENOTYPE"
        else printf "ERS NUMBER\tSAMPLE NAME\tREF\tGENOTYPE"
fi
printf "\n=========================================================================================================================================================\n"
l=0
r=0
#parse common name genotype table and substitute the corresponding systematic name when 'NONAME' is found
cat common_name.file | while read line
    do
        l=$(($l+1))
        c=0
        if [[ "$l" != '1' ]]
            then printf "\n"
        fi
        for tab in $line
                do
                    c=$(($c+1))
                    column=$(($c-1))

                    if [[ $c -eq 1 ]]
                            then    if [[ $top1 =~ $tab || $pdr =~ $tab || $zds =~ $tab ]] #if the sample contains mutations in known suppressors (top1, pdr1, pdr3,zds1,zds2) write the sample name in grey
                                    then    printf "\e[31m$l $tab\t\e[0m" | tr "\n" "\t"
                                    else    printf "\e[0m$l $tab\t" | tr "\n" "\t"
                                fi
                    fi
                    if [[ $c -eq 2 ]]
                            then printf "\e[0m$tab\t" | tr "\n" "\t"
                    fi
                    if [[ $c -eq 3 ]]
                            then printf "\e[0m$tab\t" | tr "\n" "\t"
                    fi
                    if [[ $c -gt 3 ]]
                        then    gen="`sed -n "$l"p mutation.file | cut -f$column | tr "\n" "\t" | tr -d "\t"`"
                                if [[ "$tab" == 'NONAME' ]] #print in colors: red for frameshift and stop gained, yellow point mutations, gray synonymous
                                    then    sist=`sed -n "$l"p systematic_name.file | cut -f$column | tr "\n" "\t" | tr -d "\t"`
                                            if [[ $gen =~ 'Δ' ]] || [[ $gen =~ 'FS@' ]]
                                                then printf "\e[32m$sist-$gen\t"
                                            elif [[ $gen =~ ':' ]]
                                                then if [[ "$show_syno" == '1' ]]    
                                                        then printf "\e[34m$sist-$gen\t"
                                                     fi    
                                            elif [[ $gen =~ 'x' ]]
                                                then PMCID_NUM=`scraper.py ${sist} ${gen:1:-1} | wc -l` #${gen::-1} returns the genotipe minus the last character
                                                	 printf "\e[33m$sist-$gen(${PMCID_NUM})\t" | tr -d 'x' ######################
                                                else printf "\e[97m$sist$gen\t"
                                            fi
                                    else    if [[ $gen =~ 'Δ' ]] || [[ $gen =~ 'FS@' ]]
                                                then printf "\e[32m$tab-$gen\t" | tr "\n" "\t"
                                            elif [[ $gen =~ ':' ]]
                                                then if [[ "$show_syno" == '1' ]]    
                                                         then printf "\e[34m$tab-$gen\t"
                                                     fi
                                            elif [[ $gen =~ 'x' ]]
                                                then PMCID_NUM=`scraper.py ${tab} ${gen:1:-1} | wc -l` #${gen::-1} returns the genotipe minus the last character
                                                	 printf "\e[33m$tab-$gen(${PMCID_NUM})\t" | tr -d 'x'
                                                else printf "\e[97m$tab-$gen\t" | tr "\n" "\t"
                                            fi
                                fi
                    fi
                sleep 1
                done
    done
#################################################
#                                               #
#  CALCULATING/DISPLAYING OVERLAP TABLE         #
#                                               #
#################################################
#calculating number of shared and unique mutations
# c=1
# for x in sort.ERS*.isec.vcf
#     do     n=$(echo $x | sed 's/.isec.vcf//g' | sed 's/sort.//g')
#         num=$(echo $n | sed 's/ERS//g')
#         if [[ ! "$control" =~ "$num"  ]]
#             then     declare "file_$c"="$x"    
#                     c=$(($c+1))
#         fi
#         
#     done
# cat $file_1 | grep "\./\." -v > colonia1.vcf
# cat $file_2 | grep "\./\." -v > colonia2.vcf    
# bedtools intersect -a colonia1.vcf -b colonia2.vcf > shared_mutations.vcf



printf "\e[0m\n=======================================================================================================================================================================\n"
# printf "Shared SNV" 
# cat shared_mutations.vcf | grep 'INDEL' -v | wc -l
printf "Total Unique SNV" 
cat experiment_merge.vcf | grep -v '##' | grep "0/1\|1/1\|2/2" | grep -v 'CHROM' | grep -v 'INDEL' | wc -l
# printf "Shared INDELS" 
# cat shared_mutations.vcf | grep 'INDEL' | grep "0/1\|1/1" | wc -l
printf "Total Unique INDELS" 
cat experiment_merge.vcf | grep -v '##' | grep "0/1\|1/1\|2/2" | grep -v 'CHROM' | grep  'INDEL' | wc -l
printf "\n"
printf "\nSTATISTICS\n"
printf "=======================================================================================================================================================================\n"
perl $DIR/../mareike/vcf_stats_table_all.pl experiment_merge.vcf > stat_table.txt
echo $number_of_samples >>stat_table.txt
cat stat_table.txt 
printf "\nHETEROZYGOUS SNV MUTATION OVERLAP\n"
printf "=======================================================================================================================================================================\n"
line_counter=0
# cat experiment_merge.vcf | grep -v '##' | grep -v 'INDEL' | while read line
#                                     do  tab_counter=1
#                                         line_counter=$(($line_counter+1))
#                                         for tab in $line
#                                             do  if [[ $tab_counter == '1' || $tab_counter == '2' ]]
#                                                     then printf "$tab\t\t" >> overlap_table.txt
#                                                 elif [[ $tab_counter -gt 9 && $line_counter == '1' ]]
#                                                     then     printf "$tab\t" >> overlap_table.txt
#                                                 elif [[ $tab_counter -gt 9 && $line_counter -gt 1 ]]
#                                                     then    g=$(echo $tab | grep -o "./.")
#                                                             if [[ $g == '' ]]
#                                                                 then   printf ".\t\t" >> overlap_table.txt
#                                                                 else   printf "$g\t\t" >> overlap_table.txt
#                                                             fi
#                                                 fi
#                                                 tab_counter=$(($tab_counter+1))
#                                             done
#                                         printf "\n" >> overlap_table.txt
#                                     done
# cat overlap_table.txt | grep "0/1\|1/1\|2/2\|0/2\|CHROM"
#rm sort*.vcf.gz*
#rm csq.file
cd ..
printf '\e[32m\n********************************************************************************************************\n\e[0m'



