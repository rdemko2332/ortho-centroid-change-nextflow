include { calculateCentroidDistanceStats as calculateOldCentroidDistanceStats } from './shared.nf'
include { calculateCentroidDistanceStats as calculateNewCentroidDistanceStats } from './shared.nf'

/**
* Combine the step-1 output (new core/peripheral centroids, already merged
* with untouched groups) with the unchanged residual centroids into one
* global group->rep map and one group-ID-headered fasta covering every
* group (core+peripheral+residual).
*/
process mergeCoreAndResidualCentroids {
  container 'veupathdb/orthofinder:1.9.97'

  publishDir "$params.outputDir", mode: 'copy'

  input:
    path newCoreBestReps
    path residualBestReps
    path fullProteome

  output:
    path 'newBestRepsFull.txt', emit: txt
    path 'newBestRepsFull.fasta', emit: fasta

  script:
    template 'mergeCoreAndResidualCentroids.bash'
}

process makeSeqToGroupMap {
  container 'veupathdb/orthofinder:1.9.97'

  input:
    path groupFile

  output:
    path 'seqToGroup.tsv'

  script:
    template 'makeSeqToGroupMap.bash'
}

/**
* For a batch of groups (core+peripheral+residual mixed), build a small
* per-batch diamond db from just that batch's centroids and blastp that
* batch's pooled member fastas against it -- one diamond call per batch
* rather than one per group.
*/
process distanceToOwnCentroidDiamond {
  container 'veupathdb/orthofinder:1.9.97'

  input:
    path fastas
    path newBestRepsFull
    val distanceEvalue
    val diamondOutputFields

  output:
    path 'batchHits.tsv'

  script:
    template 'distanceToOwnCentroidDiamond.bash'
}

process filterDistanceToOwnCentroidHits {
  container 'veupathdb/orthofinder:1.9.97'

  input:
    path batchHits
    path seqToGroup
    path newBestRepsFullTxt

  output:
    path 'batchDistances.tsv'

  script:
    template 'filterDistanceToOwnCentroidHits.bash'
}

process flagProteinsByDistance {
  container 'veupathdb/orthofinder:1.9.97'

  publishDir "$params.outputDir", mode: 'copy'

  input:
    path newBestRepsFullTxt
    path groupsFile
    path blastFile
    val threshold

  output:
    path 'flaggedProteins.tsv'

  script:
    template 'flagProteinsByDistance.bash'
}

workflow distanceToCentroidWorkflow {
  take:
    newCoreBestReps       // step 1 output: group_id\tcentroid_id, all core/peripheral groups
    coreGroupFastas       // channel of core+peripheral per-group fasta paths
    residualGroupFastas   // channel of residual per-group fasta paths

  main:
    distanceBatchSize = params.distanceBatchSize as Integer

    merged = mergeCoreAndResidualCentroids(
        newCoreBestReps, file(params.residualBestReps), file(params.fullProteome))

    seqToGroup = makeSeqToGroupMap(file(params.fullGroupFile))

    allGroupFastas = coreGroupFastas.mix(residualGroupFastas)

    batchHits = distanceToOwnCentroidDiamond(
        allGroupFastas.collate(distanceBatchSize),
        merged.fasta.collect(),
        params.distanceEvalue,
        params.distanceDiamondOutputFields)

    batchDistances = filterDistanceToOwnCentroidHits(
        batchHits, seqToGroup.collect(), merged.txt.collect())

    newDistances = batchDistances.collectFile(name: 'newDistances.tsv', storeDir: params.outputDir)

    // Old distance-to-(then-current)-centroid data already exists -- it's
    // exactly what intraGroupBlastFile.tsv/intraResidualGroupBlastFile.tsv
    // are, by construction (filterSimilaritiesByBestRepresentative.pl only
    // ever keeps rows whose subject is the group's mapped rep). No new
    // diamond work needed for the "old" side.
    oldDistances = Channel.fromPath(params.intraGroupBlastFile)
        .concat(Channel.fromPath(params.intraResidualGroupBlastFile))
        .collectFile(name: 'oldDistances.tsv', storeDir: params.outputDir)

    oldStats = calculateOldCentroidDistanceStats(oldDistances, file(params.fullGroupFile))
        .collectFile(name: 'oldDistanceStats.txt', storeDir: params.outputDir)

    newStats = calculateNewCentroidDistanceStats(newDistances, file(params.fullGroupFile))
        .collectFile(name: 'newDistanceStats.txt', storeDir: params.outputDir)

    flagged = flagProteinsByDistance(
        merged.txt, file(params.fullGroupFile), newDistances, params.distanceThreshold)

  emit:
    newBestRepsFullTxt = merged.txt
    newBestRepsFullFasta = merged.fasta
    seqToGroup
    newDistances
    oldStats
    newStats
    flagged
}
