#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Second-model structural cross-check: Boltz-2 co-folding (+ affinity head)
// of the 20-compound DiffDock-L consensus shortlist against the PAH
// Fe/BH4 active-site receptor sequence.
//
// params.batches_dir : directory containing batch_*.csv shards
//                       (columns: safe_id,smiles)
// params.receptor     : receptor FASTA (single-chain PAH sequence, shared
//                        across all shards)

params.batch_subdir   = 'boltz_batches'
params.receptor_name  = 'receptor_boltz.fasta'
params.outdir         = "${launchDir}/results/boltz2_batch"

params.batches_dir = "${projectDir}/assets/${params.batch_subdir}"
params.receptor    = "${projectDir}/assets/${params.receptor_name}"

include { BOLTZ2_BATCH } from './modules/boltz2_batch/main.nf'

workflow {
    receptor_ch = file(params.receptor)
    batches_ch  = Channel.fromPath("${params.batches_dir}/batch_*.csv", checkIfExists: true)
                         .map { csv -> tuple(csv, receptor_ch) }

    BOLTZ2_BATCH(batches_ch)
}
