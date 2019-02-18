#!/usr/bin/perl
{
    tag      => 'FORMAT/GT',
    name     => 'MaskHets',
    desc     => 'Genotypes set to . for samples that are het',
    apply_to => 'all',
    test     => sub {
        my $i = 8;
        for my $gt (@$MATCH)
        {
            $i++;
            my $gtype = $gt;
            my @g = split( /\//, $gtype);
            next if ($g[0] eq $g[1]);
            my @format = split(/:/,$$RECORD[$i]);
            $format[0] = $format[0] =~ /\// ? "./." : ".";
            $$RECORD[$i] = join(":",@format);
        }
        return $PASS;
    },
},

