#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

=pod

=head1 Description

Rewrite a combined groups file with reassignments applied: each reassigned
protein is removed from its old group's member list and added to its new
group's. Unlike mergeAssignedIntoStableGroups.pl (which only ever adds a
previously-unassigned sequence), this both removes and adds, since every
reassigned protein already belonged to a group.

A group's own centroid is never reassigned (see compareExponentsAndDecide
Reassignment.pl -- it's never even a candidate), so a group can shrink but
never empty out entirely. A group reduced to exactly 1 member (always that
lone group's centroid) is dissolved: written to --dissolvedSingletons
instead of staying in the output group file.

=head1 Input Parameters

=over 4

=item groupFile

Combined groups file, "GROUPID: seq1 seq2 ...".

=item reassignments

protein_id\told_group_id\tnew_group_id\t... (compareExponentsAndDecide
Reassignment.pl output; only the first 3 columns are used).

=item outputFile

Rewritten groups file, same format.

=item dissolvedSingletons

protein_id\told_group_id for every group dissolved down to 1 member.

=back

=cut

my ($groupFile, $reassignmentsFile, $outputFile, $dissolvedSingletonsFile);

&GetOptions("groupFile=s" => \$groupFile,
            "reassignments=s" => \$reassignmentsFile,
            "outputFile=s" => \$outputFile,
            "dissolvedSingletons=s" => \$dissolvedSingletonsFile);

open(my $r, '<', $reassignmentsFile) || die "Could not open file $reassignmentsFile: $!";
my %removeFrom;   # groupId -> { proteinId -> 1 }
my %addTo;        # groupId -> [ proteinId, ... ]
while (my $line = <$r>) {
    chomp $line;
    next unless length($line);
    my ($proteinId, $oldGroupId, $newGroupId) = split(/\t/, $line);
    $removeFrom{$oldGroupId}{$proteinId} = 1;
    push(@{$addTo{$newGroupId}}, $proteinId);
}
close($r);

open(my $in, '<', $groupFile) || die "Could not open file $groupFile: $!";
open(my $out, '>', $outputFile) || die "Could not open file $outputFile for writing: $!";
open(my $dis, '>', $dissolvedSingletonsFile) || die "Could not open file $dissolvedSingletonsFile for writing: $!";

while (my $line = <$in>) {
    chomp $line;
    next unless length($line);

    if ($line =~ /^(\S+):\s(.+)/) {
        my ($groupId, $seqString) = ($1, $2);
        my @members = split(/\s+/, $seqString);

        if ($removeFrom{$groupId}) {
            @members = grep { !$removeFrom{$groupId}{$_} } @members;
        }
        if ($addTo{$groupId}) {
            push(@members, @{$addTo{$groupId}});
        }

        if (scalar(@members) == 1) {
            print $dis "$members[0]\t$groupId\n";
        }
        elsif (scalar(@members) >= 2) {
            print $out "$groupId: " . join(' ', @members) . "\n";
        }
        # scalar(@members) == 0 should be structurally impossible (a group's
        # centroid is never reassigned away), so there's nothing to do for it.
    }
    else {
        die "Improper group file format: $line";
    }
}

close($in);
close($out);
close($dis);
