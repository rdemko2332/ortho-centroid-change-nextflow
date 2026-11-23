#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

=pod

=head1 Description

Build the final per-group centroid comparison report (deliverable a):
group_id, old_centroid_id, new_centroid_id, old_avg_distance,
new_mean_distance, method, n_proteins, changed.

=head1 Input Parameters

=over 4

=item oldBestReps / newBestReps

group_id\tcentroid_id (combined core+peripheral+residual, old and new).

=item oldStats / newStats

group_id\tmin\tq25\tmedian\tmean\tq75\tmax\tcount (calculateCentroidDistanceStats.pl output).

=item centroidsTsv

group_id\tn_proteins\tmethod\tcentroid_id\t... (step 1's computeGroupCentroids
output -- only covers groups actually recomputed there).

=item groupsFile

Combined groups file, for member counts of groups absent from centroidsTsv
(untouched <=2-member groups, and all residual groups).

=item outputFile

=back

=cut

my ($oldBestRepsFile, $newBestRepsFile, $oldStatsFile, $newStatsFile,
    $centroidsTsvFile, $groupsFile, $outputFile);

&GetOptions("oldBestReps=s" => \$oldBestRepsFile,
            "newBestReps=s" => \$newBestRepsFile,
            "oldStats=s" => \$oldStatsFile,
            "newStats=s" => \$newStatsFile,
            "centroidsTsv=s" => \$centroidsTsvFile,
            "groupsFile=s" => \$groupsFile,
            "outputFile=s" => \$outputFile);

sub loadBestReps {
    my ($file) = @_;
    my %map;
    open(my $fh, '<', $file) || die "Could not open file $file: $!";
    while (my $line = <$fh>) {
        chomp $line;
        next unless length($line);
        my ($groupId, $repId) = split(/\t/, $line);
        $map{$groupId} = $repId;
    }
    close($fh);
    return %map;
}

sub loadStats {
    my ($file) = @_;
    my %map;
    open(my $fh, '<', $file) || die "Could not open file $file: $!";
    while (my $line = <$fh>) {
        chomp $line;
        next unless length($line);
        my ($groupId, $min, $q25, $median, $mean, $q75, $max, $count) = split(/\t/, $line);
        $map{$groupId} = $mean;
    }
    close($fh);
    return %map;
}

my %oldCentroid = loadBestReps($oldBestRepsFile);
my %newCentroid = loadBestReps($newBestRepsFile);
my %oldMean = loadStats($oldStatsFile);
my %newMean = loadStats($newStatsFile);

open(my $ct, '<', $centroidsTsvFile) || die "Could not open file $centroidsTsvFile: $!";
my %method;
my %nProteinsRecomputed;
while (my $line = <$ct>) {
    chomp $line;
    next unless length($line);
    my ($groupId, $nProteins, $recomputeMethod) = split(/\t/, $line);
    $method{$groupId} = $recomputeMethod;
    $nProteinsRecomputed{$groupId} = $nProteins;
}
close($ct);

open(my $groups, '<', $groupsFile) || die "Could not open file $groupsFile: $!";
open(my $out, '>', $outputFile) || die "Could not open file $outputFile for writing: $!";

print $out "group_id\told_centroid_id\tnew_centroid_id\told_avg_distance\tnew_mean_distance\tmethod\tn_proteins\tchanged\n";

while (my $line = <$groups>) {
    chomp $line;
    next unless length($line);
    if ($line =~ /^(\S+):\s(.+)/) {
        my ($groupId, $seqString) = ($1, $2);

        my $recomputeMethod = $method{$groupId} || 'unchanged';
        my $nProteins = $nProteinsRecomputed{$groupId};
        if (!defined($nProteins)) {
            $nProteins = scalar(split(/\s+/, $seqString));
        }

        my $oldRep = $oldCentroid{$groupId} || '';
        my $newRep = $newCentroid{$groupId} || '';
        my $changed = ($oldRep ne '' && $newRep ne '' && $oldRep ne $newRep) ? 'true' : 'false';

        print $out "$groupId\t$oldRep\t$newRep\t" .
            (defined($oldMean{$groupId}) ? $oldMean{$groupId} : '') . "\t" .
            (defined($newMean{$groupId}) ? $newMean{$groupId} : '') . "\t" .
            "$recomputeMethod\t$nProteins\t$changed\n";
    }
    else {
        die "Improper group file format: $line";
    }
}

close($groups);
close($out);
