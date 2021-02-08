#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use List::Util qw(sum);
use POSIX qw(ceil floor);

my @six_class=qw(C>T A>G A>T C>G G>T A>C);
my @twelve_class=qw(C>T G>A A>G T>C A>T T>A C>G G>C G>T C>A A>C T>G);
#my @twelve_class=qw(C>A C>G C>T A>C A>G A>T G>A G>C G>T T>A T>G T>C);
my $DEL_folder=shift;
my %mut;
my $id=0;
#read positional info from experiment_merge
my $fh;
open ($fh, "<", "$DEL_folder/analysis/experiment_merge.vcf") or open ($fh, qq[gunzip -c "$DEL_folder/analysis/experiment_merge.vcf.gz" |]) or die;
while (my $line = <$fh>){
    next unless $line =~ /PASS/ and $line =~ /[0-9]\/[0-9]/ and $line !~ /INDEL/;
    my ($chr, $pos, undef, $ref, $alt, undef, undef, $desc, @rest)=split "\t", $line;
    foreach (split ",", $alt){
      $id++;
      $mut{$id}{'chr'}="chr".$chr;
      $mut{$id}{'pos'}=$pos;
      $mut{$id}{'ref'}=$ref;
      $mut{$id}{'alt'}=$_;
    }
}
close ($fh);
#read positional info info from list of RNAP POL2 orf
my $dirscript = dirname(__FILE__); # Directory where the script is stored]
my %trsc_db;
open ($fh, "<", $dirscript."/../defaults/RNApII_ORFs.tsv") or die;
while (my $line = <$fh>){
    my ($chr, $start, $end, $strand, $sysname, @rest)=split "\t", $line;
    $trsc_db{$chr}{$start."-".$end}=$strand;
}
close ($fh);

my %ori_db;
my %ori_db2;
open ($fh, "<", $dirscript."/../defaults/ORI_loc.txt") or die;
my $prev_chr="I",
my $prev_ars=0;
while (my $line = <$fh>){
    my ($c, $start, $end, $name, $othernames, $status)=split "\t", $line;
    next if $line =~ /^#/;
    #print $status;
    next unless $status =~ /Confirmed/i;
    my $chr=get_as_rom($c);
    #print "$c\t$chr\n";
    my $ars=($start+$end)/2;
    my $max_span=1000000000;
    if ( $chr eq $prev_chr){
          my $midpoint = ($ars + $prev_ars)/2;
          my $lead_span=$midpoint-$prev_ars;
          my $leading= $lead_span <$max_span ? $prev_ars.'-'.$midpoint : $prev_ars.'-'.($prev_ars+$max_span);
          my $lag_span=$ars-$midpoint;
          my $lagging= $lag_span < $max_span? $midpoint.'-'.$ars : ($ars-$max_span).'-'.$ars;
          $ori_db{$chr}{$leading}='1';
          $ori_db{$chr}{$lagging}='-1';
          push @{$ori_db2{$chr}}, $prev_ars.'-'.$ars;
          $prev_ars=$ars;
    }
    else{
        $prev_ars=0;
    }
    $prev_chr=$chr;
}
close ($fh);

#crossreference and Count
my %trsc_counts=%{get_counts(\%mut, \%trsc_db)};
# print Dumper \%counts;
my %replc_counts=%{get_counts(\%mut, \%ori_db)};

my %replc_locations=%{get_location(\%mut, \%ori_db2)};
my %bin_replc_locations=%{hist_locations(\%replc_locations)};

print Dumper \%bin_replc_locations;
my %trsc_frqs=%{get_frq(\%trsc_counts)};
my %rplc_frqs=%{get_frq(\%replc_counts)};

#print Dumper \%rplc_frqs;

my %trsc_asym=%{get_asym(\%trsc_frqs)};

my %rplc_asym=%{get_asym(\%rplc_frqs)};

my $trsc_s=sum (values %trsc_counts) || 0;
my $rplc_s=sum (values %replc_counts) || 0;


my @out;
foreach (@twelve_class){
  push @out, $trsc_frqs{$_};
}

open ($fh, ">", $DEL_folder."/analysis/trsc_strand.txt");
print $fh (join "\t", @twelve_class)."\ttot\tTOT\n";
print $fh (join "\t", @out)."\t$trsc_s\t$id\n";
close ($fh);

@out=();
foreach (@six_class){
  push @out, $trsc_asym{$_};
}
open ($fh, ">", $DEL_folder."/analysis/trsc_asym.txt");
print $fh (join "\t", @six_class)."\n";
print $fh (join "\t", @out)."\n";
close ($fh);

@out=();
foreach (@twelve_class){
  push @out, $rplc_frqs{$_};
}
open ($fh, ">", $DEL_folder."/analysis/rplc_strand.txt");
print $fh (join "\t", @twelve_class)."\ttot\tTOT\n";
print $fh (join "\t", @out)."\t$rplc_s\t$id\n";
close ($fh);


open ($fh, ">", $DEL_folder."/analysis/rplc_loc.txt");
# my @header = @twelve_class;
# print $fh join ("\t", @header), "\n";
# #print by columns
# while ( map {@$_} values %replc_locations ) {
#    my @row;
#    push( @row, shift @{ $replc_locations{$_} } // '' ) for @header;
#    print $fh join ("\t", @row ), "\n";
# }
print $fh join ("\t", @twelve_class), "\n";
for my $row (0,10,20,30,40,50,60,70,80,90,100){
        for my $col (@twelve_class){
                my $o=$bin_replc_locations{$col}{$row} || '0';
                print $fh "$o\t";
        }
        print $fh "\n";
}
close ($fh);


@out=();
foreach (@six_class){
  push @out, $rplc_asym{$_};
}
open ($fh, ">", $DEL_folder."/analysis/rplc_asym.txt");
print $fh (join "\t", @six_class)."\n";
print $fh (join "\t", @out)."\n";
close ($fh);

#calculate angles
my $t_a=0;
my $r_a=0;
if ($id>0){
	my $t_a = (360 * $trsc_s )/ $id;
	my $r_a = (360 * $rplc_s )/ $id;
}
#print "gnuplot -e \"a=$a; b=$id; df=$DEL_folder;\" $dirscript/../gnuplot/plot_bias.gpl";
print ("gnuplot -e \"a=$t_a; b=$id; c=$r_a; df='$DEL_folder';\" $dirscript/../gnuplot/plot_bias.gpl");
      `gnuplot -e "a=$t_a; b=$id; c=$r_a; df='$DEL_folder';" $dirscript/../gnuplot/plot_bias.gpl`;

sub complement{
  my $base = shift;
  return "T" if $base eq "A";
  return "A" if $base eq "T";
  return "C" if $base eq "G";
  return "G" if $base eq "C";
}

sub get_as_rom{
	my $chrom=shift;
	my $cn;
	if ($chrom eq '1'){$cn='chrI'}
	elsif ($chrom eq '2'){ $cn='chrII'}
	elsif ($chrom eq '3'){ $cn='chrIII'}
	elsif ($chrom eq '4'){ $cn='chrIV'}
	elsif ($chrom eq '5'){ $cn='chrV'}
	elsif ($chrom eq '6'){ $cn='chrVI'}
	elsif ($chrom eq '7'){ $cn='chrVII'}
	elsif ($chrom eq '8'){ $cn='chrVIII'}
	elsif ($chrom eq '9'){ $cn='chrIX'}
	elsif ($chrom eq '10'){ $cn='chrX'}
	elsif ($chrom eq '11'){ $cn='chrXI'}
	elsif ($chrom eq '12'){ $cn='chrXII'}
	elsif ($chrom eq '13'){ $cn='chrXIII'}
	elsif ($chrom eq '14'){ $cn='chrXIV'}
	elsif ($chrom eq '15'){ $cn='chrXV'}
	elsif ($chrom eq '16'){ $cn='chrXVI'}
        else {die};
        return $cn;

}

sub get_counts{
        my $m=shift;
        #print "$m\n";
        my %mut=%{$m};
        #print Dumper \%mut;
        my $db=shift;
        my %database=%{$db};
        my %counts;
        foreach my $idf (keys %mut){
          my $strand;
          my $chr_of_mutation = $mut{$idf}{'chr'};
          foreach my $coord (keys %{$database{$chr_of_mutation}}){
            my ($start, $end) = split "-", $coord;
            next unless $mut{$idf}{'pos'}>=$start and $mut{$idf}{'pos'}<=$end;
            $strand = $database{$chr_of_mutation}{$coord};
            last;
          }
          if ($strand){
            if ($strand == "1"){
              $counts{$mut{$idf}{'ref'}.">".$mut{$idf}{'alt'}}++;
            }
            elsif ($strand == "-1"){
              $counts{complement($mut{$idf}{'ref'}).">".complement($mut{$idf}{'alt'})}++;
            }
            #print "$idf\t$strand\t$mut{$idf}{'chr'}\t$mut{$idf}{'pos'}\t$mut{$idf}{'ref'}\t$mut{$idf}{'alt'}\n"
          }
        }
return \%counts;
}

sub get_location{
        my $m=shift;
        #print "$m\n";
        my %mut=%{$m};
        #print Dumper \%mut;
        my $db=shift;
        my %output;
        my %database=%{$db};
        foreach my $idf (keys %mut){
                my $position=$mut{$idf}{'pos'};
                my $chr_of_mutation = $mut{$idf}{'chr'};
                foreach my $coord (@{$database{$chr_of_mutation}}){
                        my ($start, $end) = split "-", $coord;
                        next unless $position>=$start and $position<=$end;
                        my $percentage= ($position-$start)/($end-$start)*100;
                        push @{$output{$mut{$idf}{'ref'}.'>'.$mut{$idf}{'alt'}}}, $percentage;
                        last;
                }
        }
return \%output;
}

sub hist_locations{
        my $m=shift;
        my %mut=%{$m};
        my %output;
        foreach my $mm(keys %mut){
                my @list=@{$mut{$mm}};
                my $bin_width=10;
                my %histogram;
                $histogram{(ceil((floor($_)+1)/$bin_width)-1)*$bin_width}++ for @list;
                $histogram{$_}=$histogram{$_}/scalar @list foreach keys %histogram;
                $output{$mm}=\%histogram;
        }
        return \%output;
}


sub get_frq{
        my $c=shift;
        my %counts=%{$c};
        my %frqs;
        my $s=sum (values %counts) || 0;
        foreach (@twelve_class){
          my $c=$counts{$_} || 0;
          if ($s>0){
            $frqs{$_} = $c/$s;
          }
          else{
            $frqs{$_} = 0;
          }
        }
return \%frqs;
}





sub get_asym{
        my $f=shift;
        my %frqs=%{$f};
        my %asym;
        foreach my $class (@six_class){
          my $a = $frqs{$class};
          my @INV=split ">", $class;
          @INV = map {complement($_)} @INV;
          my $inv=join ">", @INV;
          my $b = $frqs{$inv};
          if ($a>0 and $b > 0){
            $asym{$class} = log($a/$b)/log(2);
          }
          else {
            $asym{$class} = 0;
          }
        }
        return \%asym;
}
