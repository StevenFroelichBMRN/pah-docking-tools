# PAH docking tools

Nextflow (DSL2) pipeline wrapping DiffDock-L (Corso et al., 2023/2024,
github.com/gcorso/DiffDock) for a virtual screening campaign against the
human phenylalanine hydroxylase (PAH) Fe/BH4 active-site pocket, as part of
a pharmacological-chaperone discovery program for the R408W folding-defect
variant.

## Layout

- `modules/diffdock_batch/main.nf` — the `DIFFDOCK_BATCH` process. One task
  docks one shard (CSV) of compounds against the shared receptor, using the
  upstream prebuilt `rbgcsail/diffdock:latest` image
  (`micromamba run -n diffdock python -m inference`). Emits a per-shard
  `summary_<shard>.csv` (complex_name, diffdock_confidence, pose_path,
  status) plus the rank-1 SDF poses (renamed to `<complex_name>.sdf` to
  avoid collisions across shards when published).
- `run_diffdock_batch.nf` — root launcher; `include`s the module (Wave's
  Dockerfile/container detection needs `mainScript` to point at a root file
  that includes the module, not at the module file directly).
- `assets/receptor.pdb` — prepared apo PAH monomer (PDB 6HYC chain B,
  BH4 removed).
- `assets/batches/batch_*.csv` — sharded compound library
  (`complex_name,protein_path,ligand_description,protein_sequence`).
- `assets/compound_mapping.csv` — `safe_id -> compound_id, smiles, category,
  mw, logp, source` join table (DiffDock complex names are filesystem-safe
  slugs of the original compound ids).
- `nextflow.config` — Wave/container settings and `errorStrategy = 'ignore'`
  so one failed shard doesn't abort the run.

## Running

```bash
nextflow run run_diffdock_batch.nf \
  --batches_dir assets/batches \
  --receptor assets/receptor.pdb \
  --samples 8 \
  --outdir results/diffdock_batch
```

On Seqera Platform / AWS Batch, each shard becomes an independent Batch job
scheduled onto the GPU-enabled autoscaling compute environment, giving
natural per-shard parallelism across the library.
