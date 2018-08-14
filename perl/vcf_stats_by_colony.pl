#!/usr/bin/env perl

use strict; 
use warnings;
use Cwd;


#Declare variables
 
my $pwd = (split /\//, cwd() )[-2];

my @indel_count; my @snp_count;
my @transitions; my @transversions;
my @A_C;my @A_G;my @A_T;my @C_A;my @C_G;my @C_T;my @G_A;my @G_C;my @G_T;my @T_A;my @T_C;my @T_G;
my @ERS;
#INDELs
my @one;my @two;my @three;my @more;
my @minone;my @mintwo;my @minthree;my @less;

my @files=`ls filt.isec.sort.norm.csq.genomic.ERS*.vcf`;
chomp @files;

for my $f (@files){ 
	my $command = "cat $f". '| grep "[0-9]/[0-9]\|[0-9]/[0-9]\|#" | grep "PASS\|#" |  bcftools stats';
	my @stats_out =  readpipe("$command");
	my $MM=0; my $LL=0;
	my $M1=0; my $M2=0; my $M3=0;
	my $L1=0; my $L2=0; my $L3=0;
	my $ers;
	if ($f =~ /(ERS[0-9]+)/) {$ers=$1}
	if  ($pwd ne 'Del4733_WT-1'){
		next if $ers eq 'ERS1076728';
	}
	push @ERS, $ers;
	foreach my $l (@stats_out){
	        	my $k = $l;
   		my @s;
		chomp( $l );
		next if ($l =~ /#/);
        		if ($l =~ /SN/){
            		if ($l =~ /number of indels/){
                			@s = split("\t", $l); push @indel_count, $s[3];
            		}
            		if ($l =~ /number of SNP/){
                			@s = split("\t", $l); push @snp_count, $s[3];
            		}
    		}

    		if ($l =~ /TSTV/){
                		@s = split("\t", $l); push @transitions, $s[2]; push @transversions, $s[3];
    		}	
		if ($l =~ /ST/){
			if ($l =~ /A>C/){
				@s = split("\t", $l); push @A_C, $s[3];
		        	}
		        	if ($l =~ /A>G/){
				@s = split("\t", $l); push @A_G,$s[3];
		        	}
		        	if ($l =~ /A>T/){
				@s = split("\t", $l); push @A_T,$s[3];
		        	}
		        	if ($l =~ /C>A/){
				@s = split("\t", $l); push @C_A,$s[3];
		        	}
		        	if ($l =~ /C>G/){
				@s = split("\t", $l); push @C_G,$s[3];
		        	}
		        	if ($l =~ /C>T/){
				@s = split("\t", $l); push @C_T,$s[3];
		        	}
		        	if ($l =~ /G>A/){
				@s = split("\t", $l); push @G_A,$s[3];
		        	}
		        	if ($l =~ /G>C/){
				@s = split("\t", $l); push @G_C,$s[3];
		        	}
		        	if ($l =~ /G>T/){
				@s = split("\t", $l); push @G_T,$s[3];
		        	}
		        	if ($l =~ /T>A/){
				@s = split("\t", $l); push @T_A,$s[3];
		        	}
		        	if ($l =~ /T>C/){
				@s = split("\t", $l); push @T_C,$s[3];
		        	}
		        	if ($l =~ /T>G/){
				@s = split("\t", $l); push @T_G,$s[3];
		        	}
	    	}

    		if ($l =~ /IDD/){
    			@s = split("\t", $l);
	       	   	my $num = $s[2];
  	      		if ($num == 1){$M1+=$s[3];}
  	      		if ($num == 2){$M2+=$s[3];}	
  	      		if ($num == 3){$M3+=$s[3];}
  	      		if ($num == -1){$L1+=$s[3];}
 	      		if ($num == -2){$L2+=$s[3];}
  	      		if ($num == -3){$L3+=$s[3];}
  	      		if ($num > 3) {$MM = $MM + $s[3];}
 	      		if ($num < -3) {$LL = $LL + $s[3];}
   		}
	}
	push @more, $MM;
	push @less, $LL;
	push @one , $M1;
	push @two , $M2;
	push @three,$M3;
	push @minone,$L1;
	push @mintwo,$L2;
	push @minthree,$L3;
# close
}
foreach (@files){
		if ($_ =~ /(ERS[0-9]+)/) {$_=$1}
		else{$_='ERROR'}
}

my @output=( (join ':', @indel_count) ,
	   (join ':', @snp_count) ,
	   (join ':', @transitions) ,
	   (join ':', @transversions) ,
	   $pwd,
	   (join ':', @C_T),
	   (join ':', @A_G), 
	   (join ':', @A_T),
	   (join ':', @C_G),
	   (join ':', @G_T), 
	   (join ':', @A_C),
	   (join ':', @less),
	   (join ':', @minthree) ,
	   (join ':', @mintwo) ,
	   (join ':', @minone) ,
	   (join ':', @one) ,
	   (join ':', @two) ,
	   (join ':', @three) ,
	   (join ':', @more), 
	   scalar @ERS,
	   (join ':', @ERS)   );


print("INDEL\tSNP\tTs\tTv\tSample    \tC>T\tA>G\tA>T\tC>G\tG>T\tA>C\t<-3\t-3\t-2\t-1\t1\t2\t3\t>3\tNumber of samples\n");
print join "\t", @output;
print "\n";

