#!/usr/bin/perl

# These two examples show that the VCF output line can be edited. (Thanks to Shane McCarthy)
{
    tag      => 'FORMAT/DP',
    name     => 'MinSampleDP',
    desc     => 'Genotypes set to . for samples with DP < 5',
    apply_to => 'all',
    test     => sub {
        my $i = 8;
        for my $dp (@$MATCH)
        {
            $i++;
            next unless ($dp<5);
            my @format = split(/:/,$$RECORD[$i]);
            $format[0] = $format[0] =~ /\// ? "./." : ".";
            $$RECORD[$i] = join(":",@format);
        }
        return $PASS;
    },
},
{
    tag      => 'FORMAT/GQ',
    name     => 'MinSampleGQ',
    desc     => 'Genotypes set to . for samples with GQ < 10', #FABIO EDITED, was 15
    apply_to => 'all',
    test     => sub {
        my $i = 8;
        for my $gq (@$MATCH)
        {
            $i++;
            next unless ($gq<10); #FABIO EDITED, was 30
            my @format = split(/:/,$$RECORD[$i]);
            $format[0] = $format[0] =~ /\// ? "./." : ".";
            $$RECORD[$i] = join(":",@format);
        }
        return $PASS;
    },
},
# In this example, a minimum value of MQ>20 is required
{
    tag  => 'INFO/MQ',                       # The VCF tag to apply this filter on
    name => 'MinMQ',                          # The filter ID
    desc => 'Minimum MQ [20]',             # Description for the VCF header
    test => sub { return $MATCH < 20 ? $FAIL : $PASS },
},

