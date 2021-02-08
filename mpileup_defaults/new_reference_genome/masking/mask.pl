use strict;
use warnings;
use Data::Dumper;
my %genome;
open (my $fh, "<", 'Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa');
my $chr;
my @sequence;
my @chromosome_list;
while (my $line = <$fh>){
        chomp $line;
        if ($line =~ /^>/){
            if (defined $chr){
                $genome{$chr}=join "", @sequence;
            }
            $chr=substr($line,1);
            push @chromosome_list, $chr;
            @sequence=();
        }
        else{
            push @sequence, $line;
        }

}
$genome{$chr}=join "", @sequence;

close $fh;
open ($fh, "<", 'region-filter-slim.txt');
my %filter;
while (my $line = <$fh>){
        chomp $line;
        my @riga=split "\t", $line;
        $filter{$riga[3]}{'ch'}=$riga[0];
        $filter{$riga[3]}{'st'}=$riga[1];
        $filter{$riga[3]}{'en'}=$riga[2];
}
close $fh;
#VIII	212000	216000	#CUP1
#IV	527000	539000	#ENA
print "Pre-filtering\n";
genome_stats(\%genome);

foreach my $c (keys %genome){
        my $chrom=(split " ", $c)[0];
        foreach my $region (keys %filter){
                next unless $filter{$region}{'ch'} eq $chrom;
                my $offset=$filter{$region}{'st'};
                my $length=$filter{$region}{'en'}-$filter{$region}{'st'}+1;
                my $new_sequence=("n") x $length;
                substr($genome{$c},$offset,$length)=$new_sequence;
        }
}

print "\nPost-filtering\n";
genome_stats(\%genome);


sub genome_stats{
        my %input=%{$_[0]};
        foreach (sort keys %input){
                my $l=length($input{$_});
                print "$_\t$l\n";
        }
}
open (my $outfh, ">", "Saccharomyces_cerevisiae.EF4.69.dna_sm_MASKED.toplevel.fa");
foreach (@chromosome_list){
        print $outfh ">$_\n";
        print $outfh "$genome{$_}\n";
}
