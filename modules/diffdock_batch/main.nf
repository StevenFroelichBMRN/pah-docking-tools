#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// One task = one shard of a compound library docked against one receptor.
// Uses the upstream prebuilt rbgcsail/diffdock image (micromamba env "diffdock").
//
// Results retrieval: rather than pushing files back to GitHub or reading S3
// directly (both unavailable/unreliable from the driving session), each
// task prints its per-shard summary CSV to its OWN stdout wrapped in
// ===SHARD_CSV_START:<shard>=== / ===SHARD_CSV_END:<shard>=== markers.
// `debug true` streams that stdout into the Nextflow console log, which IS
// retrievable via the Seqera Platform workflow /log API. All the noisy
// DiffDock inference output is redirected to a file so it never reaches
// that stream and drowns out the markers.
//
// params.dump_poses_b64 (bool): if true, also base64-dump each successful
// rank1 pose inline (wrapped in ===POSE_START:<id>===/===POSE_END:<id>===
// markers) -- only used for a small targeted top-N follow-up shard, not
// the full-library run, to keep console log volume bounded.

params.samples        = 8
params.dump_poses_b64 = false

process DIFFDOCK_BATCH {
    tag "${batch_csv.baseName}"
    container 'rbgcsail/diffdock:latest'
    accelerator 1
    cpus 8
    memory '56 GB'
    errorStrategy 'ignore'
    maxRetries 0
    debug true

    input:
    tuple path(batch_csv), path(receptor)

    output:
    stdout emit: log_out

    script:
    """
    set -e
    export HOME=/home/appuser
    WORK=\$PWD
    DDIR=/home/appuser/DiffDock

    # rewrite protein_path column to point at the receptor staged in this task dir
    micromamba run -n diffdock python - <<'PY' > noisy.log 2>&1
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
    micromamba run -n diffdock python - <<PY >> noisy.log 2>&1
import yaml
cfg = yaml.safe_load(open('\$DDIR/default_inference_args.yaml'))
cfg['samples_per_complex'] = int(${params.samples})
yaml.safe_dump(cfg, open('inf.yaml', 'w'))
PY

    cd \$DDIR
    micromamba run -n diffdock python -m inference \\
        --config \$WORK/inf.yaml \\
        --protein_ligand_csv \$WORK/batch_local.csv \\
        --out_dir \$WORK/out >> \$WORK/noisy.log 2>&1 || true
    cd \$WORK

    # extract rank1 confidence per complex, tolerating missing/failed complexes,
    # then print the per-shard summary CSV (and optionally base64 poses) to
    # the REAL (unredirected) stdout, which `debug true` streams to the
    # Seqera-fetchable console log.
    mkdir -p poses_out
    micromamba run -n diffdock python - <<PY
import base64, csv, glob, os, re, shutil

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
            dest = os.path.join("poses_out", name + ".sdf")
            shutil.copy(pose, dest)
            pose_rel = dest
            status = "ok"
    out_rows.append({
        "complex_name": name,
        "diffdock_confidence": conf,
        "pose_path": pose_rel,
        "status": status,
    })

shard = "${batch_csv.baseName}"
with open(f"summary_{shard}.csv", "w", newline="") as f:
    fieldnames = ["complex_name","diffdock_confidence","pose_path","status"]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(out_rows)

print(f"===SHARD_CSV_START:{shard}===")
with open(f"summary_{shard}.csv") as f:
    print(f.read(), end="")
print(f"===SHARD_CSV_END:{shard}===")

# if everything in this shard failed, dump the tail of the redirected
# DiffDock log so the driving session can diagnose without S3 access
n_ok = sum(1 for r in out_rows if r["status"] == "ok")
if n_ok == 0 and os.path.isfile("noisy.log"):
    with open("noisy.log") as f:
        lines = f.readlines()
    tail = "".join(lines[-60:])
    print(f"===NOISY_TAIL_START:{shard}===")
    print(tail)
    print(f"===NOISY_TAIL_END:{shard}===")

if ${params.dump_poses_b64 ? 'True' : 'False'}:
    for r in out_rows:
        if r["status"] == "ok" and r["pose_path"]:
            with open(r["pose_path"], "rb") as pf:
                b64 = base64.b64encode(pf.read()).decode()
            print(f"===POSE_START:{r['complex_name']}===")
            print(b64)
            print(f"===POSE_END:{r['complex_name']}===")
PY
    """
}
