#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PAH pharmacological-chaperone virtual screen: DiffDock-L docking of a
// sharded compound library against the PAH Fe/BH4 active-site pocket.
//
// params.batches_dir : directory containing batch_*.csv shards
//                       (columns: complex_name,protein_path,ligand_description,protein_sequence)
// params.receptor     : receptor PDB (shared across all shards)
// params.samples       : DiffDock samples_per_complex (default 8)
// params.outdir        : publish directory

params.batches_dir = "${projectDir}/assets/batches"
params.receptor    = "${projectDir}/assets/receptor.pdb"
params.samples     = 8
params.outdir      = "${launchDir}/results/diffdock_batch"

include { DIFFDOCK_BATCH } from './modules/diffdock_batch/main.nf'

workflow {
    receptor_ch = file(params.receptor)
    batches_ch  = Channel.fromPath("${params.batches_dir}/batch_*.csv")
                         .map { csv -> tuple(csv, receptor_ch) }

    DIFFDOCK_BATCH(batches_ch)
}
