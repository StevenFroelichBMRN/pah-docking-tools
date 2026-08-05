#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Collects all Boltz-2 per-shard summary CSVs + co-folded PDB structures
// (real Nextflow file outputs from BOLTZ2_BATCH, not console-log markers),
// joins against the compound mapping table (which also carries the prior
// DiffDock-L consensus_score / diffdock_confidence for each compound), and
// pushes the joined results table + a tarball of co-folded PDBs back to
// the pipeline's own GitHub repo via the Contents API (stdlib urllib only)
// -- same strategy validated for the DiffDock AGGREGATE step, avoiding the
// need for S3/AWS credentials in the driving session to retrieve results.

params.github_repo    = 'StevenFroelichBMRN/pah-docking-tools'
params.results_subdir = 'results'

process BOLTZ2_AGGREGATE {
    // GITHUB_TOKEN passed as a plain process env var set only for this
    // specific launch (via launch-time configText), never committed to
    // the repo -- see the launch call in the driving code. IMPORTANT: a
    // `process { withName: 'X' { env.GITHUB_TOKEN = ... } }` configText
    // block alone does NOT reliably set the var -- Nextflow silently drops
    // a withName-scoped `env` directive unless an (even empty) top-level
    // `env {}` scope is ALSO present in the merged config
    // (nextflow-io/nextflow#1187) -- the driving code's configText must
    // include both scopes.
    //
    // Container: 'python:3.11-slim' lacks `ps`, and on this Fusion/AWS
    // Batch setup Fusion polls task metrics via `ps` *during* execution --
    // installing procps mid-script (apt-get) was too late and the task was
    // still marked FAILED by Nextflow even though the Python logic itself
    // exited 0 ("Command 'ps' required by nextflow to collect task metrics
    // cannot be found" -- nextflow-io/nextflow#3080/#1289). The
    // pytorch/pytorch image used for BOLTZ2_BATCH is Ubuntu/NVIDIA-base and
    // already ships procps (that process completed cleanly), so it's
    // reused here purely as a CPU-only container to sidestep the issue.
    container 'pytorch/pytorch:2.4.1-cuda12.1-cudnn9-runtime'
    cpus 2
    memory '4 GB'
    publishDir params.outdir, mode: 'copy'
    errorStrategy 'terminate'

    input:
    path summaries
    path pdbs
    path mapping

    output:
    path 'boltz2_results.csv'
    path 'boltz2_pdbs.tar.gz'
    path 'push_log.txt'
    path 'diag.txt'

    script:
    """
    set -e
    ls -la > diag.txt

    python3 - <<'PY' 2>&1 | tee -a diag.txt
import base64, csv, glob, json, os, sys, tarfile, traceback, urllib.request, urllib.error

def main():
    mapping = {}
    with open("${mapping}") as f:
        for row in csv.DictReader(f):
            mapping[row["safe_id"]] = row

    all_rows = []
    for path in sorted(glob.glob("summary_batch_*.csv")):
        with open(path) as f:
            for row in csv.DictReader(f):
                all_rows.append(row)

    joined = []
    for row in all_rows:
        sid = row["compound_safe_id"]
        m = mapping.get(sid, {})
        joined.append({
            "compound_id": m.get("compound_id", sid),
            "safe_id": sid,
            "category": m.get("category", ""),
            "smiles": m.get("smiles", ""),
            "diffdock_confidence": m.get("diffdock_confidence", ""),
            "consensus_score": m.get("consensus_score", ""),
            "qed": m.get("qed", ""),
            "mw": m.get("mw", ""),
            "logp": m.get("logp", ""),
            "boltz2_status": row.get("status", "failed"),
            "boltz2_iptm": row.get("iptm", ""),
            "boltz2_ptm": row.get("ptm", ""),
            "boltz2_complex_plddt": row.get("complex_plddt", ""),
            "boltz2_confidence_score": row.get("confidence_score", ""),
            "boltz2_affinity_pred_value": row.get("affinity_pred_value", ""),
            "boltz2_affinity_probability_binary": row.get("affinity_probability_binary", ""),
            "boltz2_pocket_dist": row.get("boltz_pocket_dist", ""),
        })

    fieldnames = ["compound_id","safe_id","category","smiles","diffdock_confidence",
                  "consensus_score","qed","mw","logp","boltz2_status","boltz2_iptm",
                  "boltz2_ptm","boltz2_complex_plddt","boltz2_confidence_score",
                  "boltz2_affinity_pred_value","boltz2_affinity_probability_binary",
                  "boltz2_pocket_dist"]
    with open("boltz2_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in joined:
            w.writerow(r)

    n_ok = sum(1 for r in joined if r["boltz2_status"] == "ok")
    n_total = len(joined)
    print(f"boltz2_ok={n_ok} total_attempted={n_total}")

    # Nextflow/Fusion stages the collected BOLTZ2_BATCH.out.pdbs channel
    # files into this task's work dir as SYMLINKS into the Fusion S3
    # virtual filesystem (only resolvable from inside a task container).
    # tarfile.add() does NOT dereference symlinks by default, so a naive
    # tar here silently produces dangling-symlink archive members with
    # zero real bytes once extracted anywhere else -- dereference=True
    # copies the actual referenced file content into the tar instead.
    pdb_files = glob.glob("*.pdb")
    with tarfile.open("boltz2_pdbs.tar.gz", "w:gz", dereference=True) as tar:
        for pdb in pdb_files:
            real_size = os.stat(pdb).st_size  # stat() follows symlinks -> real remote size
            tar.add(pdb, arcname=os.path.basename(pdb))
    print("tar contents (with real sizes):", [(p, os.stat(p).st_size) for p in pdb_files])

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
        except Exception as e:
            log_lines.append(f"GET {path_in_repo} exception: {type(e).__name__}: {e}")

        with open(local_path, "rb") as f:
            content_b64 = base64.b64encode(f.read()).decode()

        payload = {
            "message": f"Boltz-2 second-model consensus results: {path_in_repo}",
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
        except Exception as e:
            log_lines.append(f"PUSH EXCEPTION {path_in_repo}: {type(e).__name__}: {e}")

    if not token:
        log_lines.append("NO GITHUB_TOKEN SET -- skipping push")
    else:
        gh_put_file(f"{subdir}/boltz2_results.csv", "boltz2_results.csv")
        gh_put_file(f"{subdir}/boltz2_pdbs.tar.gz", "boltz2_pdbs.tar.gz")

    with open("push_log.txt", "w") as f:
        f.write("\\n".join(log_lines) + "\\n")
    print("\\n".join(log_lines))

try:
    main()
except Exception:
    traceback.print_exc()
    # still emit required output files (even if empty) so Nextflow doesn't
    # error on missing declared outputs -- exit 0 so the diagnostic and any
    # partial outputs actually get published for inspection.
    for fn, content in [("boltz2_results.csv", ""), ("boltz2_pdbs.tar.gz", b""), ("push_log.txt", "CRASHED, see diag.txt\\n")]:
        if not os.path.isfile(fn):
            mode = "wb" if isinstance(content, bytes) else "w"
            with open(fn, mode) as f:
                f.write(content)
PY
    """
}
