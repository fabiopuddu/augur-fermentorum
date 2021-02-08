#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Term::ANSIColor;
use Getopt::Long;
use Sort::Naturally;
use Cwd;
#################################################
#                                               #
#  CALCULATING/DISPLAYING REP DNA TABLE	        #
#                                               #
#################################################
my ($input);
my ($output);
my $control='';
GetOptions
(
'i|input=s'	  => \$input,
'c|control=s'	  => \$control,
'o|input=s'	  => \$output,
);
( $input && -f $input) or die qq[Usage: $0 \n
					 	-i name conversion file\n
						-c control ERS number(s), comma separated if > 1
						];
my $DELdir=(split "/", $output)[0];
open (my $nameconversion, '<', $input) or die "cannot open name conversion file";
chomp(my @NCfile = <$nameconversion>);
close ($nameconversion);
my %nc;my %NC;
my %bcID;
my %flnm;
my %delnm;
my %is_control;
my @CONTROLS;
for my $line (@NCfile){
	my ($barcodeID,$delname,$plate,$aka,$filename,$ERSNO,$plo,$url,$ctrl)=split("\t", $line);
	next unless $delname eq $DELdir;
	$nc{$barcodeID}=substr($aka,0,10);
	$NC{$barcodeID}=$aka;
	$delnm{$barcodeID}=$delname;
	$is_control{$barcodeID} = ($ctrl == 1) ? '+' : '-';
	push @CONTROLS, $barcodeID if $ctrl;
}
my %rDNA;my %CUP1;my %mito;my %twom;
my %ty1;my %ty2;my %ty3;my %ty4;my %ty5;
my %gwm;my %telo;my %sex;
my @SAMPLES;my @LINE;
foreach my $bcID (keys %nc) {
	if (-e "$DELdir/repDNA/$bcID.txt"){
		open (my $res_file, '<', "$delnm{$bcID}/repDNA/$bcID.txt");
		chomp(my @LNR = <$res_file>);
		close($res_file);
		my @LINE=split("\t", $LNR[1]);
		$rDNA{$bcID}=(defined $LINE[1]  and $LINE[1]  ne "-1" and $LINE[1]  ne "") ? sprintf "%.1f",$LINE[1]  : "N/A" ;
		$CUP1{$bcID}=(defined $LINE[2]  and $LINE[2]  ne "-1" and $LINE[2]  ne "") ? sprintf "%.1f",$LINE[2]  : "N/A" ;
		$mito{$bcID}=(defined $LINE[3]  and $LINE[3]  ne "-1" and $LINE[3]  ne "") ? sprintf "%.1f",$LINE[3]  : "N/A" ;
		$twom{$bcID}=(defined $LINE[4]  and $LINE[4]  ne "-1" and $LINE[4]  ne "") ? sprintf "%.1f",$LINE[4]  : "N/A" ;
		$ty1{$bcID}= (defined $LINE[5]  and $LINE[5]  ne "-1" and $LINE[4]  ne "") ? sprintf "%.1f",$LINE[5]  : "N/A" ;
		$ty2{$bcID}= (defined $LINE[6]  and $LINE[6]  ne "-1" and $LINE[5]  ne "") ? sprintf "%.1f",$LINE[6]  : "N/A" ;;
		$ty3{$bcID}= (defined $LINE[7]  and $LINE[7]  ne "-1" and $LINE[6]  ne "") ? sprintf "%.1f",$LINE[7]  : "N/A" ;;
		$ty4{$bcID}= (defined $LINE[8]  and $LINE[8]  ne "-1" and $LINE[7]  ne "") ? sprintf "%.1f",$LINE[8]  : "N/A" ;;
		$ty5{$bcID}= (defined $LINE[9]  and $LINE[9]  ne "-1" and $LINE[8]  ne "") ? sprintf "%.1f",$LINE[9]  : "N/A" ;;
		$gwm{$bcID}= (defined $LINE[10] and $LINE[10] ne "-1" and $LINE[9]  ne "") ? sprintf "%.1f",$LINE[10] : "N/A" ;;
		$sex{$bcID}= (defined $LINE[11] and $LINE[11] ne "-1" and $LINE[11] ne "") ? 		    $LINE[11] : "N/A" ;;
		push @SAMPLES, $bcID, ;
	}
	if (-e "$DELdir/repDNA/$bcID.tel"){
		open (my $tel_file, '<', "$DELdir/repDNA/$bcID.tel");
		chomp(my @LNT = <$tel_file>);
		close($tel_file);
		@LINE=split("\t", $LNT[0]);
		$telo{$bcID}=(defined $LINE[0] and $LINE[0] ne "-1" and $LINE[0] ne "") ? sprintf "%.1f", $LINE[0] : "N/A" ;;
	}
}

print "REPETITIVE DNA QUANTIFICATION";
my   $header="%-16s %-14s %-4s %-6s %-14s %-5s %-4s %-6s %-6s %-6s %-6s %-6s %-6s %-16s %-4s %-16s %-13s %-1s %-2s %-2s\n";
my $format=  "\n%-16s %-14s %-4s %-4s %-25s %-5s %-4s %-6s %-6s %-6s %-6s %-6s %-6s %-25s %-4s %-24s %-23s %-5s %-9s %-2s";
print "\n============================================================================================================================================================================\n";
printf $header, ("ERS NO.","SAMPLE NAME", "REF", "║", "rDNA", "CUP1", "║", "Ty1", "Ty2", "Ty3", "Ty4", "Ty5", "║", "TELOMERES (rpm)", "║", "Mitochondria", "2-micron", "║", "Mating Type", "║");
print "============================================================================================================================================================================";

my $refrDNA=0;
my $refmito=0;
my $reftelo=0;
my $reftwom=0;

if (scalar @CONTROLS>0){
	foreach my $control (@CONTROLS){
		$refrDNA+=$rDNA{$control};
		$refmito+=$mito{$control};
		$reftelo+=$telo{$control};
		$reftwom+=$twom{$control};
	}
	$refrDNA/=scalar @CONTROLS;
	$refmito/=scalar @CONTROLS;
	$reftelo/=scalar @CONTROLS;
	$reftwom/=scalar @CONTROLS;
}
else{
	$refrDNA='N/A';
	$refmito='N/A';
	$reftelo='N/A';
	$reftwom='N/A';
}

foreach my $bcID (sort {ncmp($NC{$a},$NC{$b})} @SAMPLES){
	my $rDNA_var='N/A'; my $mito_var='N/A'; my $tel_var='N/A'; my $twom_var='N/A';
		if ($refrDNA ne 'N/A' and $rDNA{$bcID} ne 'N/A'){
			$rDNA_var=sprintf "%.2f", $rDNA{$bcID}/$refrDNA;
			$rDNA_var= ($rDNA_var >1.5 || $rDNA_var <0.5) ? sprintf ("%4s", colored($rDNA_var, 'red') ): sprintf ("%10s", colored($rDNA_var, 'white') );
		}
	if ($refmito ne 'N/A' and $mito{$bcID} ne 'N/A'){
			if ($refmito == 0 and $mito{$bcID} ==0){
				$mito_var = 0;
			}
			else{
				$mito_var= $refmito > 0 ? sprintf "%.2f", $mito{$bcID}/$refmito : "inf";
			}
			$mito_var= ($mito_var >1.5 || $mito_var <0.5) ? sprintf ("%4s", colored($mito_var, 'red') ): sprintf ("%10s", colored($mito_var, 'white') );
		}
	if ($reftelo ne 'N/A' and $telo{$bcID} ne 'N/A'){
			$tel_var =sprintf "%.2f", $telo{$bcID}/$reftelo;
			$tel_var= ($tel_var >1.2 || $tel_var <0.8) ? sprintf ("%4s", colored($tel_var, 'red') ): sprintf ("%10s", colored($tel_var, 'white') );
		}
		if ($reftwom ne 'N/A' and $twom{$bcID} ne 'N/A'){
			if ($reftwom == 0 and $twom{$bcID} ==0){
				$twom_var = 1;
			}
			else{
				$twom_var= $reftwom>0 ? sprintf "%.2f", $twom{$bcID}/$reftwom : "inf";
			}
			$twom_var= ($twom_var >1.2 || $twom_var <0.8) ? sprintf ("%4s", colored($twom_var, 'red') ): sprintf ("%10s", colored($twom_var, 'white') );
		}
	if ($control =~ $bcID){
			printf $format, ("$bcID","$NC{$bcID}","$is_control{$bcID}","║","$rDNA{$bcID} (-)","$CUP1{$bcID}","║","$ty1{$bcID}","$ty2{$bcID}","$ty3{$bcID}","$ty4{$bcID}","$ty5{$bcID}","║","$telo{$bcID} (-)","║","$mito{$bcID} (-)",  "$twom{$bcID} (-)", "║", $sex{$bcID},"║");#print the key (current sample ERS number)
		}
		else{
			printf $format, ("$bcID","$NC{$bcID}","$is_control{$bcID}","║","$rDNA{$bcID} ($rDNA_var)","$CUP1{$bcID}","║","$ty1{$bcID}","$ty2{$bcID}","$ty3{$bcID}","$ty4{$bcID}","$ty5{$bcID}","║","$telo{$bcID} ($tel_var)","║","$mito{$bcID} ($mito_var)", "$twom{$bcID} ($twom_var)", "║", $sex{$bcID},"║");#print the key (current sample ERS number)

		}
}
print "\n============================================================================================================================================================================\n";
open (my $outfile, '>', $output);
#run through the has sorted by keys name
foreach my $bcID (sort {ncmp($NC{$a},$NC{$b})} @SAMPLES){
    if ($bcID ~~ @CONTROLS){
			print $outfile "*$bcID\t$rDNA{$bcID}\t$CUP1{$bcID}\t$mito{$bcID}\t$twom{$bcID}\t$ty1{$bcID}\t$ty2{$bcID}\t$ty3{$bcID}\t$ty4{$bcID}\t$ty5{$bcID}\t$telo{$bcID}\t$sex{$bcID}\t$gwm{$bcID}\t$NC{$bcID}\n";
		}
		else{
			print $outfile "$bcID\t$rDNA{$bcID}\t$CUP1{$bcID}\t$mito{$bcID}\t$twom{$bcID}\t$ty1{$bcID}\t$ty2{$bcID}\t$ty3{$bcID}\t$ty4{$bcID}\t$ty5{$bcID}\t$telo{$bcID}\t$sex{$bcID}\t$gwm{$bcID}\t$NC{$bcID}\n";

		}
}
close ($outfile);
