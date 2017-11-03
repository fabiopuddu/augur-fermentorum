#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Term::ANSIColor;
use Getopt::Long;

#################################################
#                                               #
#  CALCULATING/DISPLAYING REP DNA TABLE	        #
#                                               #
#################################################

my ($input);
my ($control);
GetOptions
(
'i|input=s'	  => \$input,
'c|control=s'	  => \$control,
);
( $input && -f $input && $control) or die qq[Usage: $0 \n
					 	-i name conversion file\n
						-c control ERS number(s), comma separated if > 1
						];
						
						
open (my $nameconversion, '<', $input) or die "cannot open name conversion file";
chomp(my @NC = <$nameconversion>);
close ($nameconversion);
my %nc;
my %bcID;
my %flnm;
my %delnm;
my %is_control;

my $multiref;
if ($control =~ m/,/) {
    $multiref=1;
}
else{
    $multiref=0;
}

for my $line (@NC){
		(my $barcodeID, my $delname, my $plate, my $aka, my $filename, my $ERSNO)=split("\t", $line);
    # Removed 'undef' between my $filename and my $ERSNO, because the empty column was not being recognised
		$nc{$ERSNO}=substr($aka,0,10);
		$bcID{$ERSNO}=$barcodeID;
		$flnm{$ERSNO}=$filename;
		$delnm{$ERSNO}=$delname;
		if ($control =~ $ERSNO){$is_control{$ERSNO}='+'}
		else{$is_control{$ERSNO}='-'}
		# if ($plate =~ /C/){$is_control{$ERSNO}='+'}
# 		elsif ($plate =~ /S/){$is_control{$ERSNO}='-'}
# 		else{$is_control{$ERSNO}='?'}
}

my %rDNA;
my %CUP1;
my %mito;
my %ty1;
my %ty2;
my %ty3;
my %ty4;
my %ty5; 
my %gwm;
my %telo;
my @SAMPLES;
my @LINE;
foreach my $ERSNO (keys %bcID) {
	if (-e "$flnm{$ERSNO}".'.txt'){
		open (my $res_file, '<', $flnm{$ERSNO}.'.txt');
		chomp(my @LNR = <$res_file>);	
		close($res_file);
		my @LINE=split("\t", $LNR[1]);
		$rDNA{$ERSNO}= (defined $LINE[1] and $LINE[1] ne "-1" and $LINE[1] ne "") ? sprintf "%.1f",$LINE[1] : "N/A" ;
		$CUP1{$ERSNO}=(defined $LINE[2] and $LINE[2] ne "-1" and $LINE[2] ne "") ? sprintf "%.1f",$LINE[2] : "N/A" ;
		$mito{$ERSNO}=(defined $LINE[3] and $LINE[3] ne "-1" and $LINE[3] ne "") ? sprintf "%.1f",$LINE[3] : "N/A" ;
		$ty1{$ERSNO}=(defined $LINE[4] and $LINE[4] ne "-1" and $LINE[4] ne "") ? sprintf "%.1f",$LINE[4] : "N/A" ;
		$ty2{$ERSNO}=(defined $LINE[5] and $LINE[5] ne "-1" and $LINE[5] ne "") ? sprintf "%.1f",$LINE[5] : "N/A" ;;
		$ty3{$ERSNO}=(defined $LINE[6] and $LINE[6] ne "-1" and $LINE[6] ne "") ? sprintf "%.1f",$LINE[6] : "N/A" ;;
		$ty4{$ERSNO}=(defined $LINE[7] and $LINE[7] ne "-1" and $LINE[7] ne "") ? sprintf "%.1f",$LINE[7] : "N/A" ;;
		$ty5{$ERSNO}=(defined $LINE[8] and $LINE[8] ne "-1" and $LINE[8] ne "") ? sprintf "%.1f",$LINE[8] : "N/A" ;;
		$gwm{$ERSNO}=(defined $LINE[9] and $LINE[9] ne "-1" and $LINE[9] ne "") ? sprintf "%.1f",$LINE[9]  : "N/A" ;;
		push @SAMPLES, $ERSNO, ;
	}
	if (-e "$flnm{$ERSNO}".'.tel'){
		open (my $tel_file, '<', $flnm{$ERSNO}.'.tel');
		chomp(my @LNT = <$tel_file>);	
		close($tel_file);
		@LINE=split("\t", $LNT[0]);
		$telo{$ERSNO}=(defined $LINE[0] and $LINE[0] ne "-1" and $LINE[0] ne "") ? sprintf "%.1f", $LINE[0] : "N/A" ;;
	}
}

print "REPETITIVE DNA QUANTIFICATION";
my   $header="%-13s%10s%7s%6s%11s%12s%7s%6s%6s%6s%6s%6s%6s%20s%7s%14s%8s\n";
my $format="\n%-13s%10s%8s%7s%13s%9s%8s%6s%6s%6s%6s%6s%5s%18s%10s%12s%10s";
print "\n===============================================================================================================================================\n";
printf $header, ("ERS NO."," SAMPLE NAME", "REF", "║", "rDNA", "CUP1", "║", "Ty1", "Ty2", "Ty3", "Ty4", "Ty5", "║", "TELOMERES (rpm)", "║", "Mitochondria", "║");
print "===============================================================================================================================================";

my $refrDNA=0;
my $refmito=0;
my $reftelo=0;

if ($multiref){
    my @allrefs;
    my $ncontrols;
    my $i=0;
    @allrefs=split(',', $control);
    foreach my $ncontrol (sort @allrefs){
        $refrDNA+=$rDNA{$ncontrol};
        $refmito+=$mito{$ncontrol};
        $reftelo+=$telo{$ncontrol};
        $i++;
    }
    $refrDNA/=$i;
    $refmito/=$i;
    $reftelo/=$i;
}
else{
    $refrDNA=$rDNA{$control};
    $refmito=$mito{$control};
    $reftelo=$telo{$control};
}

foreach my $ERSNO (sort @SAMPLES){
    my $rDNA_var=sprintf "%.2f", $rDNA{$ERSNO}/$refrDNA;
    my $mito_var=sprintf "%.2f", $mito{$ERSNO}/$refmito;
    my $tel_var=sprintf "%.2f", $telo{$ERSNO}/$reftelo;
    if ($control =~ $ERSNO){
        printf $format, ("$ERSNO","$nc{$ERSNO}","$is_control{$ERSNO}","║","$rDNA{$ERSNO} (-)","$CUP1{$ERSNO}","║","$ty1{$ERSNO}","$ty2{$ERSNO}","$ty3{$ERSNO}","$ty4{$ERSNO}","$ty5{$ERSNO}","║","$telo{$ERSNO} (-)","║","$mito{$ERSNO} (-)", "║");#print the key (current sample ERS number)
    }
    else{
        printf $format, ("$ERSNO","$nc{$ERSNO}","$is_control{$ERSNO}","║","$rDNA{$ERSNO} ($rDNA_var)","$CUP1{$ERSNO}","║","$ty1{$ERSNO}","$ty2{$ERSNO}","$ty3{$ERSNO}","$ty4{$ERSNO}","$ty5{$ERSNO}","║","$telo{$ERSNO} ($tel_var)","║","$mito{$ERSNO} ($mito_var)", "║");#print the key (current sample ERS number)

    }
}
print "\n===============================================================================================================================================\n";

open (my $outfile, '>', 'results.txt');
#print $outfile, "REPETITIVE DNA QUANTIFICATION";
#print $outfile, "\n=========================================================================================================================================================================\n";
#print $outfile, "ERSNumber\tSample name\trDNA\tCUP1\tTy1\tTy2\tTy3\tTy4\tTy5\tTelomeres\tMitochondrial DNA\n";
#print $outfile, "=========================================================================================================================================================================\n";
#run through the has sorted by keys name
foreach my $ERSNO (sort @SAMPLES){
	print $outfile "$flnm{$ERSNO}\t$rDNA{$ERSNO}\t$CUP1{$ERSNO}\t$mito{$ERSNO}\t$ty1{$ERSNO}\t$ty2{$ERSNO}\t$ty3{$ERSNO}\t$ty4{$ERSNO}\t$ty5{$ERSNO}\t$gwm{$ERSNO}\t$telo{$ERSNO}\t$ERSNO\t$delnm{$ERSNO}\n";
}
close ($outfile);
