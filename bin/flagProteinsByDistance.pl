#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

=pod

=head1 Description

Flag every protein whose distance to its own (new) group centroid is worse
than --threshold, or for which no distance was recorded at all (no diamond
hit found). Walks every group's full membership from --groupsFile, not just
rows present in --blastFile, so "no hit at all" members are caught. A
group's own centroid is never flagged (its distance to itself is
structurally ~0 and is never even recorded, since diamond excludes self-hits).

=head1 Input Parameters

=over 4

=item bestRepFile

group_id\tcentroid_id, one per group.

=item groupsFile

Combined groups file, "GROUPID: seq1 seq2 ...".

=item blastFile

Collected batchDistances.tsv (group_id\tqseq\tsseq\tevalue) across all batches.

=item threshold

Numeric evalue threshold; qseq flagged iff recorded evalue > threshold, or no
evalue was recorded at all.

=item outputFile

protein_id\tgroup_id\tevalue (evalue is the literal string "no_hit" when
nothing was recorded).

=back

=cut

my ($bestRepFile, $groupsFile, $blastFile, $threshold, $outputFile);

&GetOptions("bestRepFile=s" => \$bestRepFile,
            "groupsFile=s" => \$groupsFile,
            "blastFile=s" => \$blastFile,
            "threshold=s" => \$threshold,
            "outputFile=s" => \$outputFile);

my $numericThreshold = 0 + $threshold;

open(my $best, '<', $bestRepFile) || die "Could not open file $bestRepFile: $!";
my %centroidOf;
while (my $line = <$best>) {
    chomp $line;
    next unless length($line);
    my ($groupId, $repSeqId) = split(/\t/, $line);
    $centroidOf{$groupId} = $repSeqId;
}
close($best);

open(my $blast, '<', $blastFile) || die "Could not open file $blastFile: $!";
my %evalueOf;
while (my $line = <$blast>) {
    chomp $line;
    next unless length($line);
    my ($groupId, $qseq, $sseq, $evalue) = split(/\t/, $line);
    $evalueOf{$qseq} = $evalue;
}
close($blast);

open(my $groups, '<', $groupsFile) || die "Could not open file $groupsFile: $!";
open(my $out, '>', $outputFile) || die "Could not open file $outputFile for writing: $!";

while (my $line = <$groups>) {
    chomp $line;
    next unless length($line);
    if ($line =~ /^(\S+):\s(.+)/) {
        my ($groupId, $seqString) = ($1, $2);
        my $centroid = $centroidOf{$groupId};

        foreach my $seqId (split(/\s+/, $seqString)) {
            next if defined($centroid) && $seqId eq $centroid;

            my $evalue = $evalueOf{$seqId};
            if (!defined($evalue)) {
                print $out "$seqId\t$groupId\tno_hit\n";
            }
            elsif ((0 + $evalue) > $numericThreshold) {
                print $out "$seqId\t$groupId\t$evalue\n";
            }
        }
    }
    else {
        die "Improper group file format: $line";
    }
}

close($groups);
close($out);
