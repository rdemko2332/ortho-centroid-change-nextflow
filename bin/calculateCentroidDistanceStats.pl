#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Statistics::Descriptive;

=pod

=head1 Description

Per-group distance-to-centroid summary stats (min, 25th pct, median, mean,
75th pct, max, count), generalizing calculateGroupStatistics.pl /
calculateResidualGroupStatistics.pl (which this codebase already near-
duplicates for core vs. residual) to be OG/OGR-agnostic and to add a mean
column. Run once on the OLD centroid's recorded distances
(intraGroupBlastFile.tsv + intraResidualGroupBlastFile.tsv, no new diamond
work needed -- those files already hold every member's distance to the
then-current rep) and once on the NEW centroid's (the freshly-computed
batchDistances.tsv), so the two runs are directly comparable.

Evalues are clamped to [1e-200, 1e-5] (matches calculateGroupStatistics.pl
exactly) and any group member with no recorded distance to its rep (not
present in the blast file at all) is backfilled at the 1e-5 non-significant
floor, same convention. A group with ZERO recorded distances at all
(--missingGroups) gets a degenerate all-zero row, also matching the
existing convention.

=head1 Input Parameters

=over 4

=item blastFile

4-column group_id\tqseq\tsseq\tevalue (sseq is unused here -- already
guaranteed equal to the group's rep by construction).

=item groupsFile

Combined groups file, "GROUPID: seq1 seq2 ...".

=item missingGroups

One group ID per line -- groups with no blastFile rows at all.

=item outputFile

group_id\tmin\tq25\tmedian\tmean\tq75\tmax\tcount

=back

=cut

my ($blastFile, $groupsFile, $missingGroups, $outputFile);

&GetOptions("blastFile=s" => \$blastFile,
            "groupsFile=s" => \$groupsFile,
            "missingGroups=s" => \$missingGroups,
            "outputFile=s" => \$outputFile);

open(my $missing, '<', $missingGroups) || die "Could not open file $missingGroups: $!";
my %missingGroupsHash;
while (my $line = <$missing>) {
    chomp $line;
    next unless length($line);
    $missingGroupsHash{$line} = 1;
}
close($missing);

open(my $groups, '<', $groupsFile) || die "Could not open file $groupsFile: $!";
my %groupToGroupSize;
while (my $line = <$groups>) {
    chomp $line;
    next unless length($line);
    if ($line =~ /^(\S+):\s(.+)/) {
        my ($groupId, $seqString) = ($1, $2);
        $groupToGroupSize{$groupId} = scalar(split(/\s+/, $seqString));
    }
    else {
        die "Improper group file format: $line";
    }
}
close($groups);

open(my $blast, '<', $blastFile) || die "Could not open file $blastFile: $!";
my %groupEvalues;
while (my $line = <$blast>) {
    chomp $line;
    next unless length($line);
    my ($groupId, $qseq, $sseq, $evalue) = split(/\t/, $line);
    next if $qseq eq $sseq;  # defensive; diamond runs here already exclude self-hits

    my $numericEvalue = 0 + $evalue;
    if ($numericEvalue < 1.0e-200) {
        push(@{$groupEvalues{$groupId}}, 1.0e-200);
    }
    elsif ($numericEvalue > 1.0e-5) {
        push(@{$groupEvalues{$groupId}}, 1.0e-5);
    }
    else {
        push(@{$groupEvalues{$groupId}}, $numericEvalue);
    }
}
close($blast);

open(my $out, '>', $outputFile) || die "Could not open file $outputFile for writing: $!";

foreach my $groupId (sort keys %groupToGroupSize) {
    if ($missingGroupsHash{$groupId}) {
        print $out "$groupId\t0\t0\t0\t0\t0\t0\t1\n";
        next;
    }

    my @evalues = @{$groupEvalues{$groupId} || []};
    my $groupSize = $groupToGroupSize{$groupId};
    # Same accounting convention as calculateGroupStatistics.pl: backfill the
    # shortfall between recorded hits and group size at the non-significant
    # floor, rather than group size minus one (excluding the centroid itself).
    my $nonSignificantCount = $groupSize - scalar(@evalues);
    push(@evalues, 1.0e-5) foreach (1..$nonSignificantCount);

    if (!@evalues) {
        print $out "$groupId\t0\t0\t0\t0\t0\t0\t1\n";
        next;
    }

    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@evalues);

    my $min = $stat->min();
    my $q25 = $stat->quantile(1);
    my $median = $stat->median();
    my $mean = $stat->mean();
    my $q75 = $stat->quantile(3);
    my $max = $stat->max();
    my $count = scalar(@evalues);

    print $out "$groupId\t$min\t$q25\t$median\t$mean\t$q75\t$max\t$count\n";
}

close($out);
