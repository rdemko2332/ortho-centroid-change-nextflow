#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

=pod

=head1 Description

Filter one batch's diamond output (query = pooled group-member fastas, subject
database = that same batch's centroids, headered by group id) down to, for
each query protein, the one row whose subject is specifically that protein's
OWN group's centroid -- discarding any hit to a different group's centroid
that happened to also land in the batch's small reference db. Output shape
matches intraGroupBlastFile.tsv (group_id\tqseq\tsseq\tevalue) so old and new
distance-to-centroid data can be processed identically downstream.

=head1 Input Parameters

=over 4

=item diamondOut

Diamond blastp output, -f 6 qseqid sseqid evalue.

=item seqToGroup

protein -> group lookup (makeSeqToGroupMap.pl output).

=item bestReps

group -> centroid seq id (group_id\trepSeqId, one per group).

=item outputFile

Filtered output: group_id\tqseq\trepSeqId\tevalue.

=back

=cut

my ($diamondOut, $seqToGroupFile, $bestRepsFile, $outputFile);

&GetOptions("diamondOut=s" => \$diamondOut,
            "seqToGroup=s" => \$seqToGroupFile,
            "bestReps=s" => \$bestRepsFile,
            "outputFile=s" => \$outputFile);

open(my $s2g, '<', $seqToGroupFile) || die "Could not open file $seqToGroupFile: $!";
my %seqToGroup;
while (my $line = <$s2g>) {
    chomp $line;
    next unless length($line);
    my ($seqId, $groupId) = split(/\t/, $line);
    $seqToGroup{$seqId} = $groupId;
}
close($s2g);

open(my $br, '<', $bestRepsFile) || die "Could not open file $bestRepsFile: $!";
my %bestReps;
while (my $line = <$br>) {
    chomp $line;
    next unless length($line);
    my ($groupId, $repSeqId) = split(/\t/, $line);
    $bestReps{$groupId} = $repSeqId;
}
close($br);

open(my $in, '<', $diamondOut) || die "Could not open file $diamondOut: $!";
open(my $out, '>', $outputFile) || die "Could not open file $outputFile for writing: $!";

while (my $line = <$in>) {
    chomp $line;
    next unless length($line);
    my ($qseq, $sseq, $evalue) = split(/\t/, $line);

    # sseq is a group id (the batch's reference db is headered by group id);
    # keep this hit only if it's the query's own group.
    my $ownGroup = $seqToGroup{$qseq};
    next unless defined($ownGroup) && $ownGroup eq $sseq;

    print $out "$sseq\t$qseq\t$bestReps{$sseq}\t$evalue\n";
}

close($in);
close($out);
