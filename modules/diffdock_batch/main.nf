#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// One task = one shard of a compound library docked against one receptor.
// Uses the upstream prebuilt rbgcsail/diffdock image (micromamba env "diffdock").

params.samples = 8

process DIFFDOCK_BATCH {
    tag "${batch_csv.baseName}"
    container 'rbgcsail/diffdock:latest'
    accelerator 1
    cpus 8
    memory '56 GB'
    errorStrategy 'ignore'
    maxRetries 0
    publishDir "${params.outdir}/${batch_csv.baseName}", mode: 'copy'

    input:
    tuple path(batch_csv), path(receptor)

    output:
    path "summary_${batch_csv.baseName}.csv", emit: summary
    path 'poses_out/*.sdf', emit: poses_named, optional: true

    script:
    """
    set -e
    export HOME=/home/appuser
    WORK=\$PWD
    DDIR=/home/appuser/DiffDock

    # rewrite protein_path column to point at the receptor staged in this task dir
    micromamba run -n diffdock python - <<'PY'
import csv
rows = list(csv.DictReader(open("${batch_csv}")))
for r in rows:
    r["protein_path"] = "${receptor}"
with open("batch_local.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["complex_name","protein_path","ligand_description","protein_sequence"])
    w.writeheader()
    w.writerows(rows)
PY

    # copy + edit the inference config for our sampling depth
    micromamba run -n diffdock python - <<PY
import yaml
cfg = yaml.safe_load(open('\$DDIR/default_inference_args.yaml'))
cfg['samples_per_complex'] = int(${params.samples})
yaml.safe_dump(cfg, open('inf.yaml', 'w'))
PY

    cd \$DDIR
    micromamba run -n diffdock python -m inference \\
        --config \$WORK/inf.yaml \\
        --protein_ligand_csv \$WORK/batch_local.csv \\
        --out_dir \$WORK/out || true
    cd \$WORK

    # extract rank1 confidence per complex, copy pose under a globally-unique name,
    # and write a per-shard summary csv. Tolerates missing/failed complexes.
    # The output dir always exists (mkdir -p) even if zero poses succeeded,
    # so this output declaration never needs `optional: true`.
    mkdir -p "poses_${batch_csv.baseName}"
    micromamba run -n diffdock python - <<PY
import csv, glob, os, re, shutil

rows = list(csv.DictReader(open("batch_local.csv")))
out_rows = []
for r in rows:
    name = r["complex_name"]
    d = os.path.join("out", name)
    conf = ""
    pose_rel = ""
    status = "failed"
    if os.path.isdir(d):
        cands = glob.glob(os.path.join(d, "rank1_confidence*.sdf"))
        if not cands:
            cands = glob.glob(os.path.join(d, "rank1.sdf"))
        if cands:
            pose = cands[0]
            m = re.search(r"confidence([\\-0-9.]+)", os.path.basename(pose))
            if m:
                conf = m.group(1)
            dest = os.path.join("poses_${batch_csv.baseName}", name + ".sdf")
            shutil.copy(pose, dest)
            pose_rel = dest
            status = "ok"
    out_rows.append({
        "complex_name": name,
        "diffdock_confidence": conf,
        "pose_path": pose_rel,
        "status": status,
    })

with open("summary_${batch_csv.baseName}.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["complex_name","diffdock_confidence","pose_path","status"])
    w.writeheader()
    w.writerows(out_rows)
PY
    """
}
