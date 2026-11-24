include { computeGroupCentroids as computeTouchedGroupCentroids } from './centroidRecompute.nf'
include { mergeByGroupId as mergeTouchedCoreBestReps } from './shared.nf'
include { mergeByGroupId as mergeTouchedResidualBestReps } from './shared.nf'

process extractFlaggedFasta {
  container 'veupathdb/orthofinder:1.9.97'

  input:
    path flaggedProteins
    path fullProteome

  output:
    path 'flagged.fasta'

  script:
    template 'extractFlaggedFasta.bash'
}

/**
* Near-verbatim reuse of the sibling's bestRepsSelfDiamond shape: every
* flagged protein, batched, queried against the WHOLE new centroid db
* (core+peripheral+residual). -k 1 is sufficient here (unlike steps 2-3's
* own-centroid search) -- the single best hit across every group's centroid
* is automatically at least as good as the hit to the protein's own
* centroid, so there's no need to separately guarantee the own-group row
* survives a target-count cutoff.
*/
process globalFlaggedSearchDiamond {
  container 'veupathdb/orthofinder:1.9.97'

  input:
    path flaggedSubset
    path newBestRepsFull
    val globalSearchEvalue
    val diamondOutputFields

  output:
    path 'flaggedHits.tsv'

  script:
    template 'globalFlaggedSearchDiamond.bash'
}

process compareExponentsAndDecideReassignment {
  container 'veupathdb/orthofinder:1.9.97'

  publishDir "$params.outputDir", mode: 'copy'

  input:
    path flaggedProteins
    path flaggedHits
    val exponentImprovement

  output:
    path 'reassignments.tsv'

  script:
    template 'compareExponentsAndDecideReassignment.bash'
}

process applyReassignments {
  container 'veupathdb/orthofinder:1.9.97'

  publishDir "$params.outputDir", mode: 'copy'

  input:
    path groupFile
    path reassignments

  output:
    path 'newGroupFile.txt', emit: groupFile
    path 'dissolvedSingletons.tsv', emit: dissolvedSingletons

  script:
    template 'applyReassignments.bash'
}

process splitTouchedGroupFastas {
  container 'veupathdb/orthofinder:1.9.97'

  input:
    path groupFile
    path touchedGroups
    path fullProteome

  output:
    path 'touchedGroupFastas/*.fasta', optional: true

  script:
    template 'splitTouchedGroupFastas.bash'
}

process splitTouchedBestRepsByType {
  container 'veupathdb/orthofinder:1.9.97'

  input:
    path touchedBestReps

  output:
    path 'touchedCoreBestReps.txt', emit: core
    path 'touchedResidualBestReps.txt', emit: residual

  script:
    template 'splitTouchedBestRepsByType.bash'
}

process buildCentroidComparisonReport {
  container 'veupathdb/orthofinder:1.9.97'

  publishDir "$params.outputDir", mode: 'copy'

  input:
    path oldBestReps
    path newBestReps
    path oldStats
    path newStats
    path centroidsTsv
    path groupsFile

  output:
    path 'centroidComparison.tsv'

  script:
    template 'buildCentroidComparisonReport.bash'
}

workflow reassignmentWorkflow {
  take:
    flaggedProteins       // step 3 output
    newBestRepsFullFasta  // steps 2-3 output (core+peripheral NEW + residual unchanged)
    newBestRepsFullTxt    // steps 2-3 output, group_id\tcentroid_id for every group
    centroidsTsv          // step 1 output (touched core/peripheral groups only)
    oldStats              // steps 2-3 output
    newStats              // steps 2-3 output

  main:
    flagged = extractFlaggedFasta(flaggedProteins, file(params.fullProteome))

    // Same collectFile-merge-order concern as centroidRecompute.nf/
    // distanceToCentroid.nf: sort actual line content before merging so
    // this stays a stable cache input for compareExponentsAndDecideReassignment.
    flaggedHits = globalFlaggedSearchDiamond(
        flagged.splitFasta(by: 1000, file: true),
        newBestRepsFullFasta.collect(),
        params.globalSearchEvalue,
        params.distanceDiamondOutputFields)
        .splitText()
        .toSortedList()
        .flatMap { it }
        .collectFile(name: 'flaggedHits.tsv', storeDir: params.outputDir)

    reassignments = compareExponentsAndDecideReassignment(
        flaggedProteins, flaggedHits, params.exponentImprovement)

    (newGroupFile, dissolvedSingletons) = applyReassignments(file(params.fullGroupFile), reassignments)

    touchedGroups = reassignments
        .splitCsv(sep: '\t')
        .flatMap { row -> [row[1], row[2]] }
        .unique()
        .collectFile(name: 'touchedGroups.txt', newLine: true)

    touchedFastas = splitTouchedGroupFastas(newGroupFile, touchedGroups, file(params.fullProteome))

    // Same determinism concern as centroidRecompute.nf's batching -- sort
    // by filename before collating so a later -resume doesn't reshuffle
    // batch composition and turn cache hits into misses.
    sortedTouchedFastas = touchedFastas.flatten().toSortedList { a, b -> a.name <=> b.name }.flatMap { it }

    touchedResults = computeTouchedGroupCentroids(
        sortedTouchedFastas.collate(50),
        file("${projectDir}/bin/approx_centroid.py"),
        params.exactCentroidLimit, params.centroidDiamondArgs)

    touchedCentroidsTsv = touchedResults.tsv
        .splitText()
        .toSortedList()
        .flatMap { it }
        .collectFile(name: 'touchedCentroids.tsv', storeDir: params.outputDir)

    touchedBestReps = touchedCentroidsTsv
        .splitCsv(sep: '\t')
        .map { row -> "${row[0]}\t${row[3]}" }
        .collectFile(name: 'touchedBestReps.txt', newLine: true)

    touchedByType = splitTouchedBestRepsByType(touchedBestReps)

    // Groups that didn't gain/lose a member from any reassignment never show
    // up in touchedGroups.txt at all, so splitTouchedBestRepsByType's output
    // can legitimately be empty for one or both types -- mergeByGroupId
    // handles a missing/empty "fresh" file fine (nothing to drop or add).
    finalCoreBestReps = mergeTouchedCoreBestReps(
        newBestRepsFullTxt, touchedGroups, touchedByType.core)
        .collectFile(name: 'finalCoreBestReps.txt', storeDir: params.outputDir)

    finalResidualBestReps = mergeTouchedResidualBestReps(
        file(params.residualBestReps), touchedGroups, touchedByType.residual)
        .collectFile(name: 'finalResidualBestReps.txt', storeDir: params.outputDir)

    finalBestReps = finalCoreBestReps.concat(finalResidualBestReps)
        .collectFile(name: 'finalBestRepsFull.txt', storeDir: params.outputDir)

    oldBestReps = Channel.fromPath(params.coreBestReps)
        .concat(Channel.fromPath(params.residualBestReps))
        .collectFile(name: 'oldBestRepsFull.txt', storeDir: params.outputDir)

    buildCentroidComparisonReport(
        oldBestReps, finalBestReps, oldStats, newStats, centroidsTsv, newGroupFile)

  emit:
    newGroupFile
    dissolvedSingletons
    reassignments
    finalCoreBestReps
    finalResidualBestReps
}
