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
    // GITHUB_TOKEN is NOT declared via Nextflow's `secret` directive: that
    // requires the compute env's AWS execution role to have
    // secretsmanager:CreateSecret, which this account's TowerForge role
    // does not grant. Instead the token is passed as a plain process env
    // var set only for this specific launch (via launch-time configText,
    // never committed to the repo) -- see launch calls in the driving code.
    //
    // Results are pushed back via the GitHub REST "contents" API using
    // only the Python standard library (urllib) -- no `git`/apt-get
    // install needed inside the container, which is what made the
    // slim-image git-based push flaky.
    container 'python:3.11-slim'
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
    python3 - <<'PY'
import base64, csv, glob, json, os, tarfile, urllib.request, urllib.error

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

# 5. push diffdock_results.csv + top_poses.tar.gz back to the pipeline's
#    own GitHub repo via the Contents API (stdlib urllib only).
token = os.environ.get("GITHUB_TOKEN", "")
repo = "${params.github_repo}"
subdir = "${params.results_subdir}"
log_lines = []

def gh_put_file(path_in_repo, local_path):
    api = f"https://api.github.com/repos/{repo}/contents/{path_in_repo}"
    req = urllib.request.Request(api, headers={
        "Authorization": f"token {token}",
        "Accept": "application/vnd.github+json",
    })
    sha = None
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            sha = json.load(resp).get("sha")
    except urllib.error.HTTPError as e:
        if e.code != 404:
            log_lines.append(f"GET {path_in_repo} failed: {e.code} {e.read()[:300]}")

    with open(local_path, "rb") as f:
        content_b64 = base64.b64encode(f.read()).decode()

    payload = {
        "message": f"Docking campaign results: {path_in_repo}",
        "content": content_b64,
        "branch": "main",
    }
    if sha:
        payload["sha"] = sha

    body = json.dumps(payload).encode()
    put_req = urllib.request.Request(api, data=body, method="PUT", headers={
        "Authorization": f"token {token}",
        "Accept": "application/vnd.github+json",
        "Content-Type": "application/json",
    })
    try:
        with urllib.request.urlopen(put_req, timeout=60) as resp:
            resp.read()
        log_lines.append(f"PUSHED {path_in_repo}")
    except urllib.error.HTTPError as e:
        log_lines.append(f"PUSH FAILED {path_in_repo}: {e.code} {e.read()[:500]}")

if not token:
    log_lines.append("NO GITHUB_TOKEN SET -- skipping push")
else:
    gh_put_file(f"{subdir}/diffdock_results.csv", "diffdock_results.csv")
    gh_put_file(f"{subdir}/top_poses.tar.gz", "top_poses.tar.gz")

with open("push_log.txt", "w") as f:
    f.write("\\n".join(log_lines) + "\\n")
print("\\n".join(log_lines))
PY
    """
}
