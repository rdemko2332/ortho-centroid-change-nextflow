#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

=pod

=head1 Description

Flatten a combined groups file ("GROUPID: seq1 seq2 ...", core/peripheral/residual
intermixed) into a flat protein -> group lookup, one row per member.

=head1 Input Parameters

=over 4

=item groupFile

Combined groups file (e.g. fullGroupFile.txt).

=item outputFile

Output tsv: seqId\tgroupId, one row per member.

=back

=cut

my ($groupFile, $outputFile);

&GetOptions("groupFile=s" => \$groupFile,
            "outputFile=s" => \$outputFile);

open(my $in, '<', $groupFile) || die "Could not open file $groupFile: $!";
open(my $out, '>', $outputFile) || die "Could not open file $outputFile for writing: $!";

while (my $line = <$in>) {
    chomp $line;
    next unless length($line);

    if ($line =~ /^(\S+):\s(.+)/) {
        my ($groupId, $seqString) = ($1, $2);
        foreach my $seqId (split(/\s+/, $seqString)) {
            print $out "$seqId\t$groupId\n";
        }
    }
    else {
        die "Improper group file format: $line";
    }
}

close($in);
close($out);
