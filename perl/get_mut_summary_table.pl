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

my ($input);
my ($control);
my ($ploidy);

GetOptions
(
'i|input=s'      => \$input,
'c|control=s'      => \$control,
'p|ploidy=s'        => \$ploidy,
);
( $input && -f $input && $control) or die qq[Usage: $0 \n
-i name conversion file\n
-c control ERS number(s), comma separated if > 1
-p ploidy of the organism
];


open (my $nameconversion, '<', $input) or die "cannot open name conversion file";
chomp(my @NC = <$nameconversion>);
close ($nameconversion);
my %nc;
my %bcID;
my %flnm;
my %delnm;
my %is_control;

my $multiref=0;
if ($control =~ m/,/) {
    $multiref=1;
}

for my $line (@NC){
    (my $barcodeID, my $delname, my $plate, my $aka, my $filename, my $ERSNO)=split("\t", $line);
    # Removed 'undef' between my $filename and my $ERSNO, because the empty column was not being recognised
    $nc{$ERSNO}=substr($aka,0,10);
    $bcID{$ERSNO}=$barcodeID;
    if ($filename){ # If variable is empty, use the barcodeID (for SC_MFY samples)
        $flnm{$ERSNO}=$filename;
    } else{
        $flnm{$ERSNO}=$barcodeID;
    }
    $delnm{$ERSNO}=$delname;
    if ($control =~ $ERSNO){$is_control{$ERSNO}='+'}
    else{$is_control{$ERSNO}='-'}
}


######
###HEADER
######
my $number_of_samples = 0;

my $headerdip2="\n %-10s %9s %12s %60s %7s %23s %4s %40s %12s %7s %3s %12s %5s ";
my $headerdip1="\n %-10s %10s %5s %4s \e[32m%5s \e[33m%9s \e[34m%7s \e[0m%12s %4s %8s %11s %7s %10s %10s %6s \e[32m%10s \e[0m%8s %10s %3s %5s %3s %12s %12s %6s \n";
my $formatdip=" %-10s %8s %4s %4s \e[32m%5s \e[33m%8s \e[34m%9s \e[0m%9s %8s %4s %12s %8s %8s %10s %8s \e[32m%6s \e[0m%10s %8s %7s %9s %8s %5s %9s %10s \n";
my $headerhap="\n %-10s %4s %7s %6s \e[32m%5s \e[33m%9s \e[34m%7s \e[0m%12s %3s %4s %1s %9s %5s \e[32m%12s \e[0m%10s %12s %2s %6s %1s %7s %5s \n";
my $formathap=" %-12s %8s %8s %6s \e[32m%5s \e[33m%8s \e[34m%9s \e[0m%9s %8s %4s %4s %7s %7s \e[32m%7s \e[0m%12s %11s %6s %4s %4s %7s %7s \n";


########################################################################################################
#
###################################################
##                                                #
## MUTATION SUMMARY TABLE                         #
##                                                #
###################################################
#

my $dir = cwd(); # Directory where the script is executed
my $dirscript = dirname(__FILE__); # Directory where the script is stored

my @experiment_name = split('/', $dir);
my $experiment_name = $experiment_name[-2];

if ($ploidy eq 2){
    print "=========================================================================================================================================================================================================\n";
    print "$experiment_name\n";
    print "=========================================================================================================================================================================================================";
    printf "$headerdip2", '.', '.', '║', 'SINGLE NUCLEOTIDE VARIANTS GAINED (OF WHICH HOMOZIGOUS)', '║', 'LOSS OF HETEROZIGOSITY', '║', ' INSERTIONS/DELETIONS GAINED', '║', 'LOSS', 'OF', 'HETEROZIGOSITY', '║';
    printf "$headerdip1", 'ERS NO.', 'SAMPLE NAME', 'REF', '║', colored('NONSENSE', 'green'), colored('MISSENSE', 'yellow'), colored('SENSE', 'blue'), 'INTERGENIC', '│', 'TO REF.', 'TOTAL SNV', '║', 'TO ALT.','TO REF.', '║', colored('FRAMESHIFT', 'green'), 'INFRAME', 'INTERGENIC', '|', 'TOT.INDEL(HOM)', '║', 'TO ALT.', 'TO REF.', '║';
    print "=========================================================================================================================================================================================================\n";
}
else{
    print "======================================================================================================================================================================\n";
    print "$experiment_name\n";
    print "======================================================================================================================================================================";
    printf "$headerhap", 'ERS NO.', 'SAMPLE NAME', 'REF', '║', 'NONSENSE', 'MISSENSE', 'SENSE', 'INTERGENIC', '│', 'TO REF.', '|', 'TOT.SNV', '║', 'FRAMESHIFT', 'INFRAME', 'INTERGENIC', '|', 'TO REF.', '|', 'TOT.INDEL', '║';
    print "======================================================================================================================================================================\n";
}
    

##########
## MUTATION SUMMARY TABLE
#########
opendir (DIR, $dir) or die "Can't open $dir: $!";
my @files = grep {/^sort\.ERS.+\.isec\.vcf$/} readdir(DIR);
my @files_sorted = sort @files;
@files = @files_sorted;
#print "@files\n";
closedir(DIR);

foreach my $cfile (@files) {
    my $n = $cfile;
    $n =~ s/\.isec\.vcf$//g;
    $n =~ s/sort\.//g;
    my $num = $n;
    $num =~ s/ERS//g;
    my $name = $nc{$n};
    my $sampleploidy = 0;
    if ($flnm{$n} =~ m/pool/){
        $sampleploidy = $ploidy * 2;
    } else {
        $sampleploidy = $ploidy;
    }
    my $line;
    ## Count SNPs
    
    my $SNPtot = 0; #homo+hetero mutations from het-unmasked control
    my $SNPhom = 0; #homo mutations from het-unmasked control
    open(my $fhandle, '<', $cfile) or die "Unable to open file, $!";
    while ($line = <$fhandle>) {
        next if $line =~ /##|#CHROM|INDEL|ins|del|\.\/\./;
            $SNPtot++;
        if ($sampleploidy eq 2){
            next unless ($line =~ /1\/1|2\/2/);
            $SNPhom++;
        }
        if ($sampleploidy eq 4){
            next unless ($line =~ /1\/1\/1\/1|2\/2\/2\/2/);
            $SNPhom++;
        }
    }
    close($fhandle);
    
    my $SNPHOMREV = 0; #homozigous reversion to ref (i.e. control 1/1, sample 0/0)
    my $SNPLOH2 = 0; #loss of heterozygosity towards 0/0
    open(my $fhandle2, '<', "inverse_intersection/$cfile") or die "Unable to open file, $!";
    while ($line = <$fhandle2>) {
        next if $line =~ /##|#CHROM|INDEL|ins|del/;
            if ($line =~ /1\/1|2\/2/){
                if ($line =~ /PASS/){
                    $SNPHOMREV++;
                }
            }
        if ($line =~ /0\/1/){
            $SNPLOH2++;
        }
    }
    close($fhandle2);
    
    my $SNP;
    my $SNPLOH1;
    if ($ploidy eq 2) {
        my $SNPhomomask = 0; #hetero mutations from het-masked control
        open(my $fhandle3, '<', "intersect_masked/$cfile") or die "Unable to open file, $!";
        while ($line = <$fhandle3>) {
            next if $line =~ /##|#CHROM|INDEL|ins|del/;
                if ($sampleploidy eq 2){
                    next unless $line =~ /1\/1|2\/2/;
                    $SNPhomomask++;
                }
            if ($sampleploidy eq 4){
                next unless ($line =~ /1\/1\/1\/1|2\/2\/2\/2/);
                $SNPhomomask++;
            }
        }
        close($fhandle3);
    
        $SNP="$SNPtot($SNPhom)";
        $SNPLOH1 = $SNPhomomask-$SNPhom; #loss of heterozygosity towards 1/1
    } else {
        $SNP="$SNPtot";
    }
    
    ## Count INDELs
    
    my $INDtot = 0; #homo+hetero indels from het-unmasked control
    my $INDhom = 0; #homo indels from het-unmasked control
    open(my $fhandle4, '<', $cfile) or die "Unable to open file, $!";
    while ($line = <$fhandle4>) {
        next if $line =~ /##|#CHROM|\.\/\./;
            next unless $line =~ /INDEL|ins|del/;
        $INDtot++;
        if ($sampleploidy eq 2){
            next unless ($line =~ /1\/1|2\/2/);
            $INDhom++;
        }
        if ($sampleploidy eq 4){
            next unless ($line =~ /1\/1\/1\/1|2\/2\/2\/2/);
            $INDhom++;
        }
    }
    close($fhandle4);
    
    my $INDHOMREV = 0; #homozigous reversion to ref (i.e. control 1/1, sample 0/0)
    my $INDLOH2 = 0; #loss of heterozygosity towards 0/0
    open(my $fhandle5, '<', "inverse_intersection/$cfile") or die "Unable to open file, $!";
    while ($line = <$fhandle5>) {
        next if $line =~ /##|#CHROM/;
            next unless $line =~ /INDEL|ins|del/;
        if ($line =~ /1\/1|2\/2/){
            if ($line =~ /PASS/){
                $INDHOMREV++;
            }
        }
        if ($line =~ /0\/1/){
            $INDLOH2++;
        }
    }
    close($fhandle5);
    
    my $IND;
    my $INDLOH1;
    if ($ploidy eq 2) {
        my $INDhomomask = 0; #hetero indels from het-unmasked control
        open(my $fhandle6, '<', "intersect_masked/$cfile") or die "Unable to open file, $!";
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
    
    ## Count consequences

    system ("perl $dirscript/../mareike/counting_consequences.pl -i $cfile > csq.file");
    my $STOP = 0;
    my $MISS = 0;
    my $FS = 0;
    my $INFRAME = 0;
    my $SENSE = 0;
    my $UP_DOWN_SNV = 0;
    my $UP_DOWN_SNV_HOM = 0;
    my $UP_DOWN_INDEL = 0;
    my $foocount = 0;
    
    open(my $csqhandle, '<', "csq.file") or die "Unable to open file, $!";
    while ($line = <$csqhandle>) {
        if ($line =~ /SNP/){
            if ($line =~ /stop_gained/){
                ($foocount, undef, undef) = split("\t", $line);
                $STOP += $foocount;
            }
            if ($line =~ /missense_variant/){
                ($foocount, undef, undef) = split("\t", $line);
                $MISS += $foocount;
            }
            if ($line =~ /synonymous_variant/){
                ($foocount, undef, undef) = split("\t", $line);
                $SENSE += $foocount;
            }
        }
        if ($line =~ /INDEL/){
            if ($line =~ /frameshift_variant/){
                ($foocount, undef, undef) = split("\t", $line);
                $FS += $foocount;
            }
            if ($line =~ /stop_gained/){
                ($foocount, undef, undef) = split("\t", $line);
                $FS += $foocount;
            }
            if ($line =~ /inframe_insertion/){
                ($foocount, undef, undef) = split("\t", $line);
                $INFRAME += $foocount;
            }
            if ($line =~ /inframe_deletion/){
                ($foocount, undef, undef) = split("\t", $line);
                $INFRAME += $foocount;
            }
            
        }
    }
    close($csqhandle);
    
    open(my $fhandle7, '<', $cfile) or die "Unable to open file, $!";
    while ($line = <$fhandle7>) {
        next if $line =~ /##|#CHROM|\.\/\./;
        if ($line =~ /INDEL|ins|del/){
            next if $line =~ /inframe_deletion|inframe_insertion|frameshift/;
            $UP_DOWN_INDEL++;
        }
        next if $line =~ /INDEL|ins|del|stop_gained|missense|synonymous/;
        $UP_DOWN_SNV++;
        if ($sampleploidy eq 2){
            next unless $line =~ /1\/1|2\/2/;
            $UP_DOWN_SNV_HOM++;
        }
        if ($sampleploidy eq 4){
            next unless ($line =~ /1\/1\/1\/1|2\/2\/2\/2/);
            $UP_DOWN_SNV_HOM++;
        }
        
    }
    close($fhandle7);
    
    if ($ploidy eq 2){
        $UP_DOWN_SNV="$UP_DOWN_SNV($UP_DOWN_SNV_HOM)";
    } else {
        $UP_DOWN_SNV="$UP_DOWN_SNV";
    }
    
    
    ## Display results
    my $contr;
    my $apex1;
    my $apex2;
    if ($control =~ /$num/){
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
        printf "$formathap", "$n", "$name", "$contr", "║", "$STOP", "$MISS", "$SENSE", "$UP_DOWN_SNV", "│", "$SNPHOMREV", "|", "$SNP", "║", "$FS", "$INFRAME", "$UP_DOWN_INDEL", "|",  "$INDHOMREV", "|", "$IND", "║";
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


