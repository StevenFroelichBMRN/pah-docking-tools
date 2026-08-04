#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Collects all per-shard summaries + pose files, joins against the compound
// mapping table, ranks by DiffDock confidence, and pushes the results back
// to the pipeline's own GitHub repo (no direct S3/AWS credentials needed to
// retrieve results — this task has network egress and a GITHUB_TOKEN
// pipeline secret instead).

params.github_repo    = 'StevenFroelichBMRN/pah-docking-tools'
params.top_n_poses    = 50
params.results_subdir = 'results'

process AGGREGATE {
    container 'python:3.11-slim'
    secret 'GITHUB_TOKEN'
    cpus 2
    memory '4 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path summaries
    path poses
    path mapping

    output:
    path 'diffdock_results.csv'
    path 'top_poses.tar.gz'
    path 'push_log.txt'

    script:
    """
    set -e
    apt-get update -qq && apt-get install -y -qq git tar gzip > /dev/null

    python3 - <<'PY'
import csv, glob, os, tarfile

# 1. load compound mapping (safe_id -> compound_id, smiles, category, ...)
mapping = {}
with open("${mapping}") as f:
    for row in csv.DictReader(f):
        mapping[row["safe_id"]] = row

# 2. concat all per-shard summaries
all_rows = []
for path in sorted(glob.glob("summary_batch_*.csv")):
    with open(path) as f:
        for row in csv.DictReader(f):
            all_rows.append(row)

# 3. join + rank
joined = []
for row in all_rows:
    cid = row["complex_name"]
    m = mapping.get(cid, {})
    conf = row.get("diffdock_confidence", "")
    try:
        conf_f = float(conf) if conf != "" else float("-inf")
    except ValueError:
        conf_f = float("-inf")
    joined.append({
        "compound_id": m.get("compound_id", cid),
        "smiles": m.get("smiles", ""),
        "category": m.get("category", ""),
        "diffdock_confidence": conf,
        "status": row.get("status", "failed"),
        "pose_file": row.get("pose_path", ""),
        "_conf_sort": conf_f,
    })

joined.sort(key=lambda r: r["_conf_sort"], reverse=True)

with open("diffdock_results.csv", "w", newline="") as f:
    fieldnames = ["compound_id", "smiles", "category", "diffdock_confidence", "status", "pose_file"]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in joined:
        w.writerow({k: r[k] for k in fieldnames})

n_ok = sum(1 for r in joined if r["status"] == "ok")
n_total = len(joined)
print(f"docked_ok={n_ok} total_attempted={n_total}")

# 4. tar up the top-N rank1 poses (by confidence, status==ok only)
top_n = int(${params.top_n_poses})
ok_ranked = [r for r in joined if r["status"] == "ok"]
top = ok_ranked[:top_n]
with tarfile.open("top_poses.tar.gz", "w:gz") as tar:
    for r in top:
        pf = r["pose_file"]
        if pf and os.path.isfile(pf):
            tar.add(pf, arcname=f"{r['compound_id']}.sdf")

print(f"tarred {len(top)} top poses")
PY

    # push results + top pose tarball back to the pipeline's own repo
    git clone -q "https://x-access-token:\${GITHUB_TOKEN}@github.com/${params.github_repo}.git" repo > push_log.txt 2>&1 || true
    if [ -d repo ]; then
        mkdir -p "repo/${params.results_subdir}"
        cp diffdock_results.csv "repo/${params.results_subdir}/"
        cp top_poses.tar.gz "repo/${params.results_subdir}/"
        cd repo
        git config user.email "diffdock-bot@example.com"
        git config user.name "DiffDock Campaign Bot"
        git add "${params.results_subdir}"
        git commit -m "Docking campaign results (\$(date -u +%FT%TZ))" >> ../push_log.txt 2>&1 || echo "nothing to commit" >> ../push_log.txt
        git push origin main >> ../push_log.txt 2>&1 || echo "PUSH FAILED" >> ../push_log.txt
        cd ..
    else
        echo "CLONE FAILED" >> push_log.txt
    fi
    cat push_log.txt
    """
}
