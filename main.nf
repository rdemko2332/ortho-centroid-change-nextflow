#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//---------------------------------------------------------------
// Including workflows
//---------------------------------------------------------------

include { centroidRecomputeWorkflow } from './modules/centroidRecompute.nf'
include { distanceToCentroidWorkflow } from './modules/distanceToCentroid.nf'
include { reassignmentWorkflow } from './modules/reassignment.nf'

//---------------------------------------------------------------
// Step 1: recompute centroids for core/peripheral groups with
// more than 2 members. Residual groups are untouched here --
// their centroids only change later, if a reassignment (steps
// 4-5) adds a member to one of them.
//
// Steps 2+3: distance to own (new) centroid for every protein in
// every group (core+peripheral+residual), plus old-vs-new stats
// and flagging.
//---------------------------------------------------------------

workflow {
  if (!params.coreGroupFastasDir) {
    throw new Exception("Missing params.coreGroupFastasDir")
  }
  if (!params.coreBestReps) {
    throw new Exception("Missing params.coreBestReps")
  }
  if (!params.residualGroupFastasDir) {
    throw new Exception("Missing params.residualGroupFastasDir")
  }

  coreGroupFastas = Channel.fromPath("${params.coreGroupFastasDir}/*.fasta", checkIfExists: true)
  residualGroupFastas = Channel.fromPath("${params.residualGroupFastasDir}/*.fasta", checkIfExists: true)

  centroidRecomputeWorkflow(coreGroupFastas, file(params.coreBestReps))

  distanceToCentroidWorkflow(
      centroidRecomputeWorkflow.out.newCoreBestReps,
      coreGroupFastas,
      residualGroupFastas)

  reassignmentWorkflow(
      distanceToCentroidWorkflow.out.flagged,
      distanceToCentroidWorkflow.out.newBestRepsFullFasta,
      distanceToCentroidWorkflow.out.newBestRepsFullTxt,
      centroidRecomputeWorkflow.out.centroidsTsv,
      distanceToCentroidWorkflow.out.oldStats,
      distanceToCentroidWorkflow.out.newStats)
}
