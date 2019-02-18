#!/usr/bin/perl
#Require a depth >=5 to call a valid mutation 
{
    tag      => 'FORMAT/DP',						# The VCF tag to apply this filter on
    name     => 'MinSampleDP',						# The filter ID
    desc     => 'Genotypes set to . for samples with DP < 5',		# Description for the VCF header
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
#Require a GQ (GenotypeQuality)>= 10 to call a valid mutation (90% accuracy) 
{
    tag      => 'FORMAT/GQ',
    name     => 'MinSampleGQ',
    desc     => 'Genotypes set to . for samples with GQ < 10',
    apply_to => 'all',
    test     => sub {
        my $i = 8;
        for my $gq (@$MATCH)
        {
            $i++;
            next unless ($gq<10);
            my @format = split(/:/,$$RECORD[$i]);
            $format[0] = $format[0] =~ /\// ? "./." : ".";
            $$RECORD[$i] = join(":",@format);
        }
        return $PASS;
    },
},
#Require a MQ (MappingQuality)>= 20 to call a valid mutation (99% accuracy) 
{
    tag  => 'INFO/MQ',                 
    name => 'MinMQ',                
    desc => 'Minimum MQ [20]',           
    test => sub { return $MATCH < 20 ? $FAIL : $PASS },
},

