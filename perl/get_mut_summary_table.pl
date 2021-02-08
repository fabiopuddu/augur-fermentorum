#!/usr/bin/env perl
use strict;
use Cwd;
use warnings;
use Data::Dumper;
use Term::ANSIColor qw( colored );
use Getopt::Long;
use List::Util 1.33 'any';
use File::Basename;
####################################################
#                                                  #
#  CALCULATING/DISPLAYING MUT SUMMARY TABLE        #
#                                                  #
####################################################
# Arguments: ploidy, input (or just use directory), -i name\ conversion.tsv
my ($input); my ($control);my ($ploidy);my ($nameconv);my ($DELdir);
GetOptions
(
'n|name_conversion=s'       => \$nameconv,
'c|control=s'               => \$control,
'p|ploidy=s'                => \$ploidy,
'd|delfolder=s'             => \$DELdir
);
( $nameconv && -f $nameconv) or die qq[Usage: $0
                -n name conversion file
                -c control ERS number(s), comma separated if > 1
                -p ploidy of the organism
];
my $multiref;
if (defined $control){
    $multiref = $control =~ m/,/ ? 0:1 ;
}
else{
    $multiref=0;
}

open (my $nameconversion, '<', $nameconv);
my %nc;my %bcID;my %flnm;my %delnm;my %is_control;
while (my $line = <$nameconversion>){
    chomp($line);
    (my $barcodeID, my $delname, my $plate, my $aka, my $filename, my $ERSNO, my $sample_ploidy)=split("\t", $line);
    $nc{$barcodeID}=substr($aka,0,10);
    $bcID{$barcodeID}=$barcodeID;
    #if ($filename){                         # If variable is empty, use the barcodeID (for SC_MFY samples)
    #    $flnm{$barcodeID}=$filename;
    #} else{
        $flnm{$barcodeID}=$barcodeID;
    #}
    $delnm{$barcodeID}=$delname;
    if (defined $control and $control =~ $barcodeID){$is_control{$barcodeID}='+'}
    else{$is_control{$barcodeID}='-'}
}
close ($nameconversion);
######
###HEADER
######
my $dir = cwd(); # Directory where the script is executed
my $dirscript = dirname(__FILE__); # Directory where the script is stored
my @experiment_name = split('/', $dir);
my $experiment_name = $experiment_name[-1];
my $number_of_samples = 0;

my $headerdip2="\n %-10s %9s %12s %60s %7s %23s %4s %40s %12s %7s %3s %12s %5s ";
my $headerdip1="\n %-10s %10s %5s %4s \e[32m%5s \e[33m%9s \e[34m%7s \e[0m%12s %4s %8s %11s %7s %10s %10s %6s \e[32m%10s \e[0m%8s %10s %3s %5s %3s %12s %12s %6s \n";
my $formatdip=" %-10s %8s %4s %4s \e[32m%5s \e[33m%8s \e[34m%9s \e[0m%9s %8s %4s %12s %8s %8s %10s %8s \e[32m%6s \e[0m%10s %8s %7s %9s %8s %5s %9s %10s \n";

my $headerhap0="%-30s %-25s %-27s %-3s %-6s\n";
my $headerhap="\n %-14s %-10s %-3s %-2s \e[32m%5s \e[33m%9s \e[36m%7s \e[0m%12s %3s %4s %1s %9s %5s \e[32m%12s \e[0m%10s %12s %2s %6s %1s %7s %5s \n";
my $formathap=   "%-14s %-10s %-4s %-2s \e[32m%5s \e[33m%8s \e[36m%9s \e[0m%9s %8s %4s %4s %6s %8s \e[32m%7s \e[0m%12s %11s %6s %4s %4s %7s %7s \n";

if ($ploidy eq 2){
    print "=========================================================================================================================================================================================================\n";
    print "$experiment_name\n";
    print "=========================================================================================================================================================================================================";
    printf "$headerdip2", '.', '.', '║', 'SINGLE NUCLEOTIDE VARIANTS GAINED (OF WHICH HOMOZIGOUS)', '║', 'LOSS OF HETEROZIGOSITY', '║', ' INSERTIONS/DELETIONS GAINED', '║', 'LOSS', 'OF', 'HETEROZIGOSITY', '║';
    printf "$headerdip1", 'ERS NO.', 'SAMPLE NAME', 'REF', '║', colored('NONSENSE', 'green'), colored('MISSENSE', 'yellow'), colored('SENSE', 'blue'), 'INTERGENIC', '│', 'TO REF.', 'TOTAL SNV', '║', 'TO ALT.','TO REF.', '║', colored('FRAMESHIFT', 'green'), 'INFRAME', 'INTERGENIC', '|', 'TOT.INDEL(HOM)', '║', 'TO ALT.', 'TO REF.', '║';
    print "=========================================================================================================================================================================================================\n";
}
elsif ($ploidy eq 1){
    print "======================================================================================================================================================================\n";
    printf "$headerhap0", "$experiment_name", "║", "SNV TO ALT.", "|", "TO REF.";
    print "======================================================================================================================================================================";
    printf "$headerhap", 'bcID', 'NAME', 'REF', '║', 'NONSENSE', 'MISSENSE', 'SENSE', 'INTERGENIC', '│', 'TOT.SNV', '|', 'TOT.SNV.', '║', 'FRAMESHIFT', 'INFRAME', 'INTERGENIC', '|', 'TO REF.', '|', 'TOT.INDEL', '║';
    print "======================================================================================================================================================================\n";
}


opendir (DIR, "$dir/$DELdir/analysis") or die "Can't open $dir: $!";
my @files = grep {/^isec.*.vcf.gz$/} readdir(DIR);
@files = sort @files;
chomp @files;
closedir(DIR);

foreach my $cfile (@files) {
    my $n = $cfile;
    $n =~ s/.vcf.gz//g;
    $n =~ s/isec\.//g;
    my $num = $n;
    my $name = $nc{$n};
    my $sampleploidy = 0;
    if ($flnm{$n} =~ m/pool/){
        $sampleploidy = $ploidy * 2;
    } else {
        $sampleploidy = $ploidy;
    }
    my $line;
    ## Count SNPs
    my $SNPtot = 0;       #homo+hetero mutations from het-unmasked control
    my $SNPhom = 0;       #homo mutations from het-unmasked control
    my $INDtot = 0;       #homo+hetero indels from het-unmasked control
    my $INDhom = 0;       #homo indels from het-unmasked control
    my $SNPHOMREV = 0;    # homozigous reversion to ref (i.e. control 1/1, sample 0/0)
    my $SNPLOH2 = 0;      # loss of heterozygosity towards 0/0
    my $INDHOMREV = 0;    #homozigous reversion to ref (i.e. control 1/1, sample 0/0)
    my $INDLOH2 = 0;      #loss of heterozygosity towards 0/0
    my $STOP = 0;         #csq=stop_gained
    my $MISS = 0;         #csq=missense_variant
    my $FS = 0;           #csq=frameshift_variant
    my $INFRAME = 0;      #csq=inframe_insertion+inframe_deletion
    my $SENSE = 0;        #csq=synonymous_variant
    my $UP_DOWN_SNV = 0;        #intergenic snv
    my $UP_DOWN_INDEL = 0;      #intergenic indels
    my $UP_DOWN_SNV_HOM = 0;



    open(my $fhandle, qq[gunzip -c "$DELdir/analysis/$cfile"|]) or die;
    while ($line = <$fhandle>) {
            next unless $line =~ /PASS/ and $line !~ /\.\/\./;
            if ($line =~ /INDEL|ins|del/){
                  $INDtot++;
                  $INDhom++           if $sampleploidy eq 2 and $line =~ /([1-9]+)\/(\1)/;#match genotypes such as 1/1 or 2/2 but not 0/0 (\1 is a backreference)
                  $FS++               if $line =~ /frameshift_variant/;
                  $INFRAME++          if $line =~ /inframe_insertion\|inframe_deletion/;
                  $UP_DOWN_INDEL++    if $line !~ /frameshift_variant/ and $line !~ /inframe_insertion\|inframe_deletion/;
            }
            else{
                  $SNPtot++;
                  $SNPhom++           if $sampleploidy eq 2 and $line =~ /([1-9]+)\/(\1)/;#match genotypes such as 1/1 or 2/2 but not 0/0 (\1 is a backreference)
                  $STOP++             if $line =~ /stop_gained/;
                  $MISS++             if $line =~ /missense_variant/;
                  $SENSE++            if $line =~ /synonymous_variant/;
                  if ($line !~ /stop_gained/ and $line !~ /missense_variant/ and $line !~ /synonymous_variant/){
                                      $UP_DOWN_SNV++;
                                      $UP_DOWN_SNV_HOM++  if $sampleploidy eq 2 and $line =~ /([1-9]+)\/(\1)/;
                  }
            }
    }
    close($fhandle);
    open(my $fhandle2, qq[gunzip -c "$DELdir/analysis/inv_$cfile"|]) or die;
    while ($line = <$fhandle2>) {
            next unless $line =~ /PASS/ and $line !~ /\.\/\./;
            if ($line =~ /INDEL|ins|del/){
                  $INDHOMREV++ if $line =~ /([1-9]+)\/(\1)/;
                  $INDLOH2++ if $line =~ /([0-9]+)\/(?!\1)/;
            }
            else{
                  $SNPHOMREV++  if $line =~ /([1-9]+)\/(\1)/; #match genotypes such as 1/1 or 2/2 but not 0/0 (\1 is a backreference)
                  $SNPLOH2++    if $line =~ /([0-9]+)\/(?!\1)/;#match genotypes such as 0/1 or 1/2 but not 0/0 1/1 2/2 (?! is a negative lookahead)
            }
    }
    close($fhandle2);





    my $SNP;
    my $SNPLOH1;
    if ($ploidy eq 2) {
        my $SNPhomomask = 0; #hetero mutations from het-masked control
        open(my $fhandle3, '<', "$DELdir/analysis/intersect_masked/$cfile") ;#or die "Unable to open file, $!";
        while ($line = <$fhandle3>) {
            next if $line =~ /##|#CHROM|INDEL|ins|del/;
                if ($sampleploidy eq 2){
                    next unless $line =~ /1\/1|2\/2/;
                    $SNPhomomask++;
                }

        }
        close($fhandle3);

        $SNP="$SNPtot($SNPhom)";
        $SNPLOH1 = $SNPhomomask-$SNPhom; #loss of heterozygosity towards 1/1
    } else {
        $SNP="$SNPtot";
    }


    my $IND;
    my $INDLOH1;
    if ($ploidy eq 2) {
        my $INDhomomask = 0; #hetero indels from het-unmasked control
        open(my $fhandle6, '<', "$DELdir/analysis/intersect_masked/$cfile"); #or die "Unable to open file, $!";
        while ($line = <$fhandle6>) {
            next if $line =~ /##|#CHROM/;
                next unless $line =~ /INDEL|ins|del/;
            if ($sampleploidy eq 2){
                next unless $line =~ /1\/1|2\/2/;
                $INDhomomask++;
            }
            if ($sampleploidy eq 4){
                next unless ($line =~ /1\/1\/1\/1|2\/2\/2\/2/);
                $INDhomomask++;
            }
        }
        close($fhandle6);

        $IND="$INDtot($INDhom)";
        $INDLOH1 = $INDhomomask-$INDhom; #loss of heterozygosity towards 1/1
    } else {
        $IND="$INDtot";
    }

    if ($ploidy eq 2){
        $UP_DOWN_SNV="$UP_DOWN_SNV($UP_DOWN_SNV_HOM)";
    } else {
        $UP_DOWN_SNV="$UP_DOWN_SNV";
    }


    ## Display results
    my $contr;
    my $apex1;
    my $apex2;
    if (defined $control and $control =~ /$num/){
        $contr = '+++';
        if ($multiref eq '1'){
            $apex1 = 'x';
            $apex2 = 'y';
        }
    } else {
        $contr = '-';
        $number_of_samples++;
        $apex1 = '';
        $apex2 = '';
    }

    if ($ploidy eq '1'){
        printf "$formathap", "$n", "$name", "$contr", "║", "$STOP", "$MISS", "$SENSE", "$UP_DOWN_SNV", "│", "$SNP", "|", "$SNPHOMREV", "║", "$FS", "$INFRAME", "$UP_DOWN_INDEL", "|",  "$INDHOMREV", "|", "$IND", "║";
    } else {
        printf "$formatdip", "$n", "$name", "$contr", "║", "$STOP", "$MISS", "$SENSE", "$UP_DOWN_SNV", "│", "$SNPHOMREV", "$SNP", "║", "$SNPLOH1$apex1", "$SNPLOH2$apex2", "║", "$FS", "$INFRAME", "$UP_DOWN_INDEL", "|", "$IND", "║", "$INDLOH1", "$INDLOH2", "║";
    }
}

if ($ploidy eq '1'){
    print "======================================================================================================================================================================\n";
} else {
   print "=========================================================================================================================================================================================================\n"
}

if ($multiref eq '1' & $ploidy eq '2'){
    print "x = number of 0/1 mutations in the combined control that are 1/1 in this control sample\n";
    print "y = number of 0/1 mutations in the combined control that are 0/0 in this control sample\n";
}
