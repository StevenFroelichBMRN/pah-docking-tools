#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.mapping = "${projectDir}/assets/compound_mapping.csv"
params.outdir  = "${launchDir}/results/agg_test"

include { AGGREGATE } from './modules/aggregate/main.nf'

workflow {
    summaries = Channel.fromPath("${projectDir}/assets/agg_test/summary_batch_000.csv").collect()
    poses     = Channel.fromPath("${projectDir}/assets/agg_test/BDB_50156156.sdf").collect()
    mapping_ch = file(params.mapping)
    AGGREGATE(summaries, poses, mapping_ch)
}
