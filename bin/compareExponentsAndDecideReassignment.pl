#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

=pod

=head1 Description

Decide which flagged proteins actually get reassigned: a flagged protein
moves to its single best-matching group (from the global flagged-vs-all-
centroids search) iff that group differs from its current one AND the match
is at least --exponentImprovement orders of magnitude better (e.g. 2 =
100-fold) than its own recorded distance to its current group's centroid.

Scored as exponent = -log10(max(evalue, 1e-300)), which is mathematically
equivalent to "100-fold better" for any evalue representation (including
diamond's literal "0"/"0.0" and the "no_hit" sentinel used when step 2/3
found nothing at all for that protein, treated as the worst possible score,
exponent 0) without any evalue-string parsing.

=head1 Input Parameters

=over 4

=item flaggedProteins

protein_id\tgroup_id\tevalue_or_no_hit (flagProteinsByDistance.pl output --
evalue/marker here is the protein's distance to its CURRENT group's centroid).

=item flaggedHits

protein_id\tbestMatchGroupId\tevalue (globalFlaggedSearchDiamond output --
each flagged protein's single best hit across every group's centroid, old
and new alike).

=item exponentImprovement

Minimum -log10(new) - -log10(old) to qualify for reassignment.

=item outputFile

protein_id\told_group_id\tnew_group_id\told_evalue\tnew_evalue\texponent_improvement,
one row per protein that actually qualifies to move.

=back

=cut

my ($flaggedProteinsFile, $flaggedHitsFile, $exponentImprovement, $outputFile);

&GetOptions("flaggedProteins=s" => \$flaggedProteinsFile,
            "flaggedHits=s" => \$flaggedHitsFile,
            "exponentImprovement=f" => \$exponentImprovement,
            "outputFile=s" => \$outputFile);

sub exponentOf {
    my ($evalue) = @_;
    return 0 if (!defined($evalue) || $evalue eq 'no_hit');
    my $numeric = 0 + $evalue;
    $numeric = 1e-300 if ($numeric < 1e-300);
    return -(log($numeric) / log(10));
}

open(my $fp, '<', $flaggedProteinsFile) || die "Could not open file $flaggedProteinsFile: $!";
my %oldGroupOf;
my %oldEvalueOf;
while (my $line = <$fp>) {
    chomp $line;
    next unless length($line);
    my ($proteinId, $groupId, $evalue) = split(/\t/, $line);
    $oldGroupOf{$proteinId} = $groupId;
    $oldEvalueOf{$proteinId} = $evalue;
}
close($fp);

open(my $fh, '<', $flaggedHitsFile) || die "Could not open file $flaggedHitsFile: $!";
open(my $out, '>', $outputFile) || die "Could not open file $outputFile for writing: $!";

while (my $line = <$fh>) {
    chomp $line;
    next unless length($line);
    my ($proteinId, $newGroupId, $newEvalue) = split(/\t/, $line);

    my $oldGroupId = $oldGroupOf{$proteinId};
    next unless defined($oldGroupId);        # not actually one of our flagged proteins; ignore
    next if $newGroupId eq $oldGroupId;       # best global match is still its own group; nothing to do

    my $oldEvalue = $oldEvalueOf{$proteinId};
    my $improvement = exponentOf($newEvalue) - exponentOf($oldEvalue);

    if ($improvement >= $exponentImprovement) {
        my $oldEvalueDisplay = defined($oldEvalue) ? $oldEvalue : 'no_hit';
        print $out "$proteinId\t$oldGroupId\t$newGroupId\t$oldEvalueDisplay\t$newEvalue\t$improvement\n";
    }
}

close($fh);
close($out);
