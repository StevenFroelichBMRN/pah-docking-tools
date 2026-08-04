#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Second-model structural cross-check: Boltz-2 co-folding (+ affinity head)
// of the 20-compound DiffDock-L consensus shortlist against the PAH
// Fe/BH4 active-site receptor sequence, followed by an aggregate step that
// joins results against the DiffDock consensus table and pushes them back
// to this repo via the GitHub Contents API.
//
// params.batches_dir : directory containing batch_*.csv shards
//                       (columns: safe_id,smiles)
// params.receptor     : receptor FASTA (single-chain PAH sequence, shared
//                        across all shards)
// params.mapping      : safe_id -> compound_id/smiles/category/DiffDock
//                        scores join table

params.batch_subdir   = 'boltz_batches'
params.receptor_name  = 'receptor_boltz.fasta'
params.mapping        = "${projectDir}/assets/boltz_compound_mapping.csv"
params.outdir         = "${launchDir}/results/boltz2_batch"

params.batches_dir = "${projectDir}/assets/${params.batch_subdir}"
params.receptor    = "${projectDir}/assets/${params.receptor_name}"

include { BOLTZ2_BATCH }     from './modules/boltz2_batch/main.nf'
include { BOLTZ2_AGGREGATE } from './modules/boltz2_aggregate/main.nf'

workflow {
    receptor_ch = file(params.receptor)
    batches_ch  = Channel.fromPath("${params.batches_dir}/batch_*.csv", checkIfExists: true)
                         .map { csv -> tuple(csv, receptor_ch) }

    BOLTZ2_BATCH(batches_ch)

    summaries_ch = BOLTZ2_BATCH.out.summary.collect()
    pdbs_ch      = BOLTZ2_BATCH.out.pdbs.collect()
    mapping_ch   = file(params.mapping)

    BOLTZ2_AGGREGATE(summaries_ch, pdbs_ch, mapping_ch)
}
