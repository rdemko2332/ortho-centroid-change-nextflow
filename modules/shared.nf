// Processes shared across more than one step. Copied/adapted from
// orthoFinder-nextflow's modules/shared.nf.

process mergeByGroupId {
  container 'veupathdb/orthofinder:1.9.97'

  input:
    path cached, stageAs: 'cached.txt'
    path touchedGroups
    path fresh, stageAs: 'fresh.txt'

  output:
    path 'merged.txt'

  script:
    template 'mergeByGroupId.bash'
}

/**
* Per-group distance-to-centroid summary stats. Defined once here and
* `include`'d twice (aliased) by distanceToCentroid.nf, since a DSL2 process
* can only be invoked once per importing workflow scope -- this needs to run
* on both the old-centroid distances and the new-centroid distances.
*/
process calculateCentroidDistanceStats {
  container 'veupathdb/orthofinder:1.9.97'

  input:
    path blastFile
    path groupsFile

  output:
    path 'distance_stats.txt'

  script:
    template 'calculateCentroidDistanceStats.bash'
}
