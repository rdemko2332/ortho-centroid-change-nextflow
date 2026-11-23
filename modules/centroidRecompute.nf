include { mergeByGroupId } from './shared.nf'

/**
* Bucket a batch of per-group fasta files by member count:
*   <=2 members          -- skipped entirely (step 1 only recomputes >2-member groups)
*   3..exactCentroidLimit -- exact all-vs-all bucket
*   >exactCentroidLimit   -- two-round sampling bucket
*/
process classifyGroupsForRecompute {
  container 'veupathdb/orthofinder:1.9.97'

  input:
    path fastas
    val exactCentroidLimit

  output:
    path 'exact/*.fasta', optional: true, emit: exact
    path 'large/*.fasta', optional: true, emit: large

  script:
    template 'classifyGroupsForRecompute.bash'
}

/**
* Compute (approximate) centroids for a batch of group fastas via approx_centroid.py,
* appending one TSV row and one new-centroid fasta record per group.
*/
process computeGroupCentroids {
  container 'veupathdb/orthofinder:1.9.97'

  input:
    path fastas
    path approxCentroidScript
    val exactCentroidLimit
    val diamondArgs

  output:
    path 'centroids_batch.tsv', emit: tsv
    path 'newCentroids.fasta', emit: fasta

  script:
    template 'computeGroupCentroids.bash'
}

workflow centroidRecomputeWorkflow {
  take:
    coreGroupFastas     // channel of per-group fasta file paths (core+peripheral only, not residual)
    cachedCoreBestReps  // path to current coreBestReps.txt

  main:
    approxCentroidScript = file("${projectDir}/bin/approx_centroid.py")

    // CLI/params-file overrides of these arrive as Strings even when declared as
    // integers in nextflow.config; .collate() only resolves against a real int.
    classifyBatchSize = params.classifyBatchSize as Integer
    exactBatchSize = params.exactBatchSize as Integer
    largeBatchSize = params.largeBatchSize as Integer

    (exactFastas, largeFastas) = classifyGroupsForRecompute(
        coreGroupFastas.collate(classifyBatchSize),
        params.exactCentroidLimit)

    // Both buckets batched to their own size, then mixed into one channel so
    // computeGroupCentroids is only called (DSL2 "process already used" --
    // each process can only be invoked once per workflow scope) once, while
    // still preserving a different batch size per bucket.
    allBatches = exactFastas.flatten().collate(exactBatchSize)
        .mix(largeFastas.flatten().collate(largeBatchSize))

    results = computeGroupCentroids(
        allBatches, approxCentroidScript, params.exactCentroidLimit, params.centroidDiamondArgs)

    centroidsTsv = results.tsv
        .collectFile(name: 'centroids.tsv',
                     keepHeader: false,
                     storeDir: params.outputDir)

    newCentroidsFasta = results.fasta
        .collectFile(name: 'newCentroids.fasta', storeDir: params.outputDir)

    touchedGroups = centroidsTsv
        .splitCsv(sep: '\t')
        .map { row -> row[0] }
        .collectFile(name: 'touchedGroups.txt', newLine: true)

    // group_id\tcentroid_id, same shape as coreBestReps.txt
    freshBestReps = centroidsTsv
        .splitCsv(sep: '\t')
        .map { row -> "${row[0]}\t${row[3]}" }
        .collectFile(name: 'freshBestReps.txt', newLine: true)

    newCoreBestReps = mergeByGroupId(cachedCoreBestReps, touchedGroups, freshBestReps)
        .collectFile(name: 'newCoreBestReps.txt', storeDir: params.outputDir)

  emit:
    centroidsTsv       // group_id n_proteins method centroid_id mean_sim n_candidates n_refs n_alignments seed
    newCentroidsFasta  // fasta, headers = group_id, only for the >2-member groups actually recomputed
    newCoreBestReps    // group_id\tcentroid_id for ALL core/peripheral groups (touched + untouched merged)
}
