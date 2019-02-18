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
aneup=0;
#############################################
#  FUNCTIONS                                #
#############################################
waitforcompletion(){
if [[ ${v} -eq '1' ]]; then printf "Waiting for jobs $1  to complete"; else printf "Waiting for jobs to complete"; fi
list="ffsaknerogd${1}"
finito=`squeue | grep "${list}" |wc -l | tr -d "\t"`
while [[ $finito != '0' ]]    
    do  finito=`squeue  | grep "${list}" |wc -l | tr -d "\t"`
        if [[ ${v} -eq '1' ]]; then printf ".${finito}"; else printf "."; fi
        sleep $[ ( $RANDOM % 20 )  + 1 ]s
    done  
printf "\n"    
}
#############################################
#   GET   AND CHECK COMMAND LINE OPTIONS    #
#############################################
while getopts "sn:hfle:c:C:xvrtFa" opt
    do  case "$opt" in
                        h)  printf "############   HELP   ###############\nOPTIONS\n"
                            printf "\t-h\tThis Help\n"    
                            printf "\t-c\tReference File\n"
                            printf "\t-n\tSet ploidy of samples (1 or 2)\n"
                            printf "\t-C\tMultiple reference files: ERS numbers separated with a comma\n"
                            printf "\t-x\tno control sample(an artificial control containing all the mutations shared between all the samples will be used\n"
                            printf "\t-r\tanalyse repetitive regions in the samples\n"
                            printf "\t-a\taneuploidy: generate ploidy circos charts for every file analysed\n"
                            printf "\t-F\tForce Rewrite: clear all previous analysis before restarting\n"
                            printf "\t-v\tverbose: print out detailed analysis progression\n"
                            printf "\t-l\tmask only very very low quality variants\n"
                            printf "\t-s\tsynonimous: print out details on synonimous mutations\n"  
                            exit 0
                        ;;
                        F)
                        hard_force_rewrite=1
                        ;;
                        s)
                        show_syno=1
                        ;;
                        c)
                        control="$OPTARG"
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
                        a)
                        aneup=1
                        ;;
                        
        esac
    done
if [[ $control == '0000' && $no_ctrl == '0' ]]
    then    echo "Please specify the ERS number of the reference sample with the flag -c or -C"
            exit 1
fi
if [[ ! $control =~ "ERS"  ]]
    then    echo "Please specify the ERS number of the reference sample in the ERSXXXXXX format"
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
                    declare "control_$suffix"="$line"
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
echo "Verbosity ${v}"
if [[ $hard_force_rewrite == '1' ]]
    then echo 'force rewrite' 
        rm -rf analysis
        rm -rf repDNA
        rm -rf transposons
        rm -rf ploidy_data
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
        do headr=`zcat $f | grep "CHROM"  | cut -f10- | sed "s/\t*$//g"`
            new_headr=''
            for tab in "$headr"
                do     if [[ ! $tab =~ "ERS" ]]
                        then     printf "Fixing headers...\n"
                        #name=`echo $tab | tr -d '..' 
                        #ERSnumber=`grep -w $tab ../../name\ conversion.tsv | awk -F"\t" '$6 ~ ERS {print $6}'`
			ERSnumber=`echo $f | sed 's/.vcf.gz//g'`
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
while read -r line 
    do if [[ ! -a csq.$line ]]
           then export SBATCH_CMD="variant_effect_predictor.pl --species saccharomyces_cerevisiae -i $line  --format vcf -o vep.$line.txt --no_progress --force_overwrite --offline"
                PROC1=$(sbatch submit_sbatch_analyser.sh | sed 's/Submitted batch job //g') 
		export SBATCH_CMD="vcf2consequences_vep -v $line -i vep.$line.txt 2>/dev/null | bgzip > csq.$line"
                PROC2=$(sbatch --dependency=afterok:${PROC1} submit_sbatch_analyser.sh | sed 's/Submitted batch job //g')
                export SBATCH_CMD=""
		proclist="${proclist}\|${PROC2}"
       fi
    done < <(ls mito.ERS*.vcf.gz)
cd ..    
#########################################
##        ANALYSE repetitive DNA       ##
#########################################
#if [[  "$force_rewrite" == "1"  && $rDNA == 1 ]] 
if [[   $rDNA == 1 ]] 
    then    printf "Analysing repetitive DNA sequences ...\n"
            mkdir -p repDNA
        	  cd repDNA
            while read -r line ; 
            do
                check_pool=`echo $line | grep 'pool'`
                    if [ -z $check_pool ]
                        then sampleploidy=$ploidy
                        else sampleploidy=$(( $ploidy * 2 ))
                    fi
               name=`echo $line | tr "/" "\n" | tail -n1 | sed "s|\.bam||g" | head -n1`
	       export SBATCH_CMD="rDNA-cnv_estimate.pl -p ${sampleploidy} -i ../$line > $name.txt"
               if [[ ${v} -eq '1' ]]; then echo $command; fi
               PROC1=$(sbatch submit_sbatch_analyser.sh | sed 's/Submitted batch job //g') 
               proclist="${proclist}\|${PROC1}"
      	       export SBATCH_CMD="zcat ../BAM/${name}.fq1.gz ../BAM/${name}.fq2.gz | telomeres.pl 5  > ${name}.tel"
               PROC1=$(sbatch submit_sbatch_analyser.sh | sed 's/Submitted batch job //g') 
               proclist="${proclist}\|${PROC1}"  	 
		 export SBATCH_CMD=""
             done < <(cat ../bams_for_mpileup)
             cd ..
fi
################################
##        Analyse ploidy      ##
################################
if [[ $aneup == 1 ]]
    then printf "Analysing ploidy...\n"
         mkdir ploidy_data 
         cd ploidy_data
         while read -r line
            do
                check_pool=`echo $line | grep 'pool'`
                    if [ -z $check_pool ]
                        then sampleploidy=$ploidy
                        else sampleploidy=$(( $ploidy * 2 ))
                    fi
         	    code1=`echo $line | tr '/' "\n" | tail -n1 | sed "s|\.bam||g"` #SCMFY or #SD code
                    code2=`grep -w ${code1} ../../name\ conversion.tsv | awk -F"\t" '$6 ~ ERS {print $6}'` # Fetch ERS code (cut -f 7 gave issues depending on the tabulation of name\ conversion.tsv)
         	    name=`grep -w ${code1} ../../name\ conversion.tsv | awk -F"\t" '{print $2}'`
         	    export SBATCH_CMD="CGH.pl -i ../$line -p $sampleploidy -f -l \"${name}:${code1}:${code2}\""
         	    if [[ ${v} -eq '1' ]]; then echo ${command}; fi
         	    PROC1=$(sbatch submit_sbatch_analyser.sh | sed 's/Submitted batch job //g') 
         	    proclist="${proclist}\|${PROC1}"
		    export SBATCH_CMD=""
           done < <(cat ../bams_for_mpileup)
         cd ..	
fi
######## Wait for the parallel computing to be finished ########
waitforcompletion "${proclist}"      

########  ANEUPLOIDY DNA POSTPROCESSING  #############
if [[ $aneup == 1 ]]
	then cd ploidy_data
	        montage -geometry 1200x1200 png/*.png aneuploidy_report.png #mount all the images in a single file
                cd ..
fi

#####################################
#### START MUTATION ANALYSIS    #####
#####################################

###################################
####NORMALISATION AND PRELIMINARY FILTERS
####################################

#NORMALIZATION OF INDEL ANNOTATION
sleep 10
cd analysis
printf "Normalising indels.... \n"
for f in csq.*.vcf.gz 
   do  if [[ ! -a norm.$f ]]
         then if [[ ${v} -eq '1' ]]; then printf "Normalising $f ..."; fi
              n=`echo $f | sed 's|csq.mito.||g' | sed 's|.vcf.gz||g'`
	    if [[ ${v} -eq '1' ]]; 
	        then  bcftools norm -f $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa $f | bgzip > $n.vcf.gz
      	        else  bcftools norm -f $DIR/../mpileup_defaults/reference_genome/Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa $f 2>/dev/null| bgzip > $n.vcf.gz
              fi       
         else if [[ ${v} -eq '1' ]]; then printf "\n Skipping $f"; fi 
       fi                     
   done  
#INITIAL SORTING OF MUTATIONS IN VCF FILES, GZIPPING AND TABIX               
printf 'Sorting mutations...'
isec_files=(*isec.vcf)
    if [[ ! -e "${isec_files[0]}" ]]
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
                            then cat $n.isec.vcf  | vcf-annotate -f $DIR/../perl/gt-filter-lax.pl > $n.isec.filt.vcf #lax filter
                            else cat $n.isec.vcf  | vcf-annotate -f $DIR/../perl/gt-filter.pl > $n.isec.filt.vcf #strict filter
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
                        then cat $n.isec.vcf  | vcf-annotate -f $DIR/../perl/gt-filter-lax.pl > $n.isec.filt.vcf #lax filter
                        else cat $n.isec.vcf  | vcf-annotate -f $DIR/../perl/gt-filter.pl > $n.isec.filt.vcf #strict filter
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



# if [[ $ploidy == '2' ]]
#     then    zcat sort.control.vcf.gz | vcf-annotate -f $DIR/../perl/mask-hets.pl > combined.control.masked.vcf
#             grep -Ev "\./\." combined.control.masked.vcf > combined.control.masked2.vcf #These two lines are needed to remove ./. from the masked control file, otherwise they will be intersected
#             mv combined.control.masked2.vcf combined.control.masked.vcf                 # and the corresponding mutations in the sample files will be removed leading to  incorrect numbers of LOH TO ALT.




#vcf-subset -c $control combined.control.masked.vcf -e > control.masked.vcf.subset
#           mv control.masked.vcf.subset control.masked.vcf


# 
#             if [[ $v == '1' ]] ; then printf "Of which homozygous: " ;fi
#             if [[ $v == '1' ]] ; then cat  combined.control.masked.vcf | grep '##' -v  | grep '#CHROM' -v | grep 'INDEL' -v | wc -l ;fi
# fi




#intersection with masked het-masked reference



# if [[ ( ! -a intersect_masked || "$force_rewrite" == "1" ) && $ploidy == "2" ]]
#     then    printf "\n"
#             mkdir intersect_masked
#             cd intersect_masked
#             cat ../combined.control.masked.vcf | vcf-sort > sort.control.masked.vcf
#             bgzip -f sort.control.masked.vcf
#             tabix -f -p vcf sort.control.masked.vcf.gz
#             echo '......Intersecting samples with masked control: '$control
#             for x in ../sort.ERS*vcf.gz
#                 do  n=$(echo $x | sed 's/.vcf.gz//g')
#                     n=$(echo $n | sed 's/..\///g')
#                     if [[ $v == '1' ]] ; then printf "\r\033[K                       Intersecting $n";fi
#                     bedtools intersect -header -a $x -b sort.control.masked.vcf.gz -v > $n.isec.vcf
#                     #vcf-isec -f -a -c $x sort.control.masked.vcf.gz 2>/dev/null > $n.isec.vcf #actual intersect command
#                     if [[ $low = 1 ]] #filtering mutations from non-control samples based on quality
#                             then cat $n.isec.vcf  | vcf-annotate -f $DIR/../perl/gt-filter-lax.pl > $n.isec.filt.vcf #lax filter
#                             else cat $n.isec.vcf  | vcf-annotate -f $DIR/../perl/gt-filter.pl > $n.isec.filt.vcf #strict filter
#                     fi
#                     mv $n.isec.filt.vcf $n.isec.vcf
#                 done
#                 cd ..
# fi

#TO BE USED WHEN INVERSE INTERSECTING WITH MASKED REF
#for x in ../ERS*.vcf
#do  cat $x | vcf-annotate -f $DIR/../perl/mask-hets.pl > masked.$x #mask all the hets in the samples but not in the control
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

printf "\n"
formatcontrol=$(echo $control | tr " " ",") # Fixing for control format
get_mut_summary_table.pl -i ../../name\ conversion.tsv -c ${formatcontrol} -p $ploidy


##convert results in gene list:
    ##You will need to amend this script. I will include a file called
    ##yeast_genelist_nameconv.tsv which contains the information for the more common names for genes.
    ##At "my $conversion_file=" you will need to specify the path to this file to make the script run well)
printf "\n"
perl $DIR/../perl/vcf_to_gene_list.pl -i experiment_merge.vcf > hom.table.file
perl $DIR/../perl/vcf_to_gene_list_het.pl -i experiment_merge.vcf > het.table.file

#################################################
#                                               #
#  CALCULATING/DISPLAYING GENOTYPE TABLE        #
#                                               #
#################################################
printf "\n"
if [[ ${show_syno} == 1 ]]
	then get_genotype_table.pl -i hom.table.file -s
	else get_genotype_table.pl -i hom.table.file
fi
printf "\n"
#############################################
#                                           #
#    DISPLAY REPETITIVE DNA TABLE           #
#                                           #
#############################################
if [[ ${rDNA} ]]
	then
        cd ../repDNA
		get_repetitive_table.pl -i ../../name\ conversion.tsv
#            get_repetitive_table.pl -i ../../name\ conversion.tsv -c ${formatcontrol}
        cd ../analysis
fi
printf "\n"
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
#printf "Total Unique SNV"
#cat experiment_merge.vcf | grep -v '##' | grep "0/1\|1/1\|2/2" | grep -v 'CHROM' | grep -v 'INDEL' | wc -l
# printf "Shared INDELS"
# cat shared_mutations.vcf | grep 'INDEL' | grep "0/1\|1/1" | wc -l
#printf "Total Unique INDELS"
#cat experiment_merge.vcf | grep -v '##' | grep "0/1\|1/1\|2/2" | grep -v 'CHROM' | grep  'INDEL' | wc -l
printf "\n"
printf "\nSTATISTICS\n"
printf "=======================================================================================================================================================================\n"
perl $DIR/../perl/vcf_stats_table_all.pl experiment_merge.vcf > stat_table.txt
cat stat_table.txt
vcf_stats_by_colony.pl  > stats_table_col.txt
cat stats_table_col.txt

#printf "\nHETEROZYGOUS SNV MUTATION OVERLAP\n"
#printf "=======================================================================================================================================================================\n"
#line_counter=0
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

echo 'DONEOK'

