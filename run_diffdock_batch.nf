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

// batch_subdir/receptor_name are simple relative names so a param override
// (e.g. via params-text on Seqera) still resolves against projectDir instead
// of the launch/work directory.
params.batch_subdir  = 'batches'
params.receptor_name = 'receptor.pdb'
params.samples       = 8
params.outdir        = "${launchDir}/results/diffdock_batch"

params.batches_dir = "${projectDir}/assets/${params.batch_subdir}"
params.receptor    = "${projectDir}/assets/${params.receptor_name}"

params.mapping        = "${projectDir}/assets/compound_mapping.csv"
params.github_repo    = 'StevenFroelichBMRN/pah-docking-tools'
params.top_n_poses    = 50
params.results_subdir = 'results'

include { DIFFDOCK_BATCH } from './modules/diffdock_batch/main.nf'
include { AGGREGATE }      from './modules/aggregate/main.nf'

workflow {
    receptor_ch = file(params.receptor)
    mapping_ch  = file(params.mapping)
    batches_ch  = Channel.fromPath("${params.batches_dir}/batch_*.csv", checkIfExists: true)
                         .map { csv -> tuple(csv, receptor_ch) }

    DIFFDOCK_BATCH(batches_ch)

    AGGREGATE(
        DIFFDOCK_BATCH.out.summary.collect(),
        DIFFDOCK_BATCH.out.poses_named.collect(),
        mapping_ch
    )
}
