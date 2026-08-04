#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// One task = one shard of compounds co-folded with the PAH receptor using
// Boltz-2 (Passaro & Wohlwend et al. 2025, github.com/jwohlwend/boltz),
// requesting the affinity head alongside the co-folded structure. Serves as
// the "second structural model" consensus cross-check against the prior
// DiffDock-L virtual screen (see modules/diffdock_batch).
//
// params.his285_pos / his290_pos / glu330_pos: 1-based position of the
// active-site landmark residues (His285, His290, Glu330) WITHIN the
// receptor sequence as submitted to Boltz -- used to recompute the
// ligand-centroid-to-active-site-landmark distance directly in each
// Boltz-predicted structure's own coordinate frame.
//
// Retrieval: results are emitted as REAL Nextflow file outputs (a
// per-shard summary_<shard>.csv plus a pdb_out/ directory of co-folded
// structures), NOT via the console-log marker strategy used for
// DIFFDOCK_BATCH -- Seqera's live console log keeps only a bounded tail
// window, which silently drops earlier shards' marker blocks once later
// shards print. Downstream BOLTZ2_AGGREGATE reads these as proper
// Nextflow channel inputs and pushes the joined table back to GitHub.

params.his285_pos = 259
params.his290_pos = 264
params.glu330_pos = 304
params.recycling_steps  = 3
params.diffusion_samples = 1

process BOLTZ2_BATCH {
    tag "${batch_csv.baseName}"
    container 'pytorch/pytorch:2.4.1-cuda12.1-cudnn9-runtime'
    // Default Docker/AWS Batch shm is 64 MB, which crashed Boltz's PyTorch
    // DataLoader workers with "Bus error ... out of shared memory" on the
    // first launch attempt (20/20 failed identically). AWS Batch's
    // shm-size container option wants a plain integer of bytes with NO
    // '=' (nextflow-io/nextflow#5190/#5176).
    containerOptions '--shm-size 8000000000'
    accelerator 1
    cpus 8
    memory '30 GB'
    errorStrategy 'ignore'
    maxRetries 0

    input:
    tuple path(batch_csv), path(receptor_fasta)

    output:
    path "summary_${batch_csv.baseName}.csv", emit: summary
    path 'pdb_out/*.pdb', emit: pdbs, optional: true

    script:
    """
    set -e
    export HOME=/root
    WORK=\$PWD

    pip install --quiet --no-cache-dir boltz > pip.log 2>&1 || (cat pip.log && exit 1)

    python3 - <<'PY' > prep.log 2>&1
import csv, os

with open("${receptor_fasta}") as f:
    lines = [l.strip() for l in f if not l.startswith(">")]
receptor_seq = "".join(lines)

rows = list(csv.DictReader(open("${batch_csv}")))
os.makedirs("yamls", exist_ok=True)
for r in rows:
    sid = r["safe_id"]
    smiles = r["smiles"]
    yaml_text = f'''version: 1
sequences:
  - protein:
      id: A
      sequence: {receptor_seq}
  - ligand:
      id: B
      smiles: '{smiles}'
properties:
  - affinity:
      binder: B
'''
    with open(f"yamls/{sid}.yaml", "w") as yf:
        yf.write(yaml_text)

print(f"wrote {len(rows)} yaml inputs")
PY

    boltz predict yamls/ \\
        --use_msa_server \\
        --out_dir out \\
        --output_format pdb \\
        --recycling_steps ${params.recycling_steps} \\
        --diffusion_samples ${params.diffusion_samples} \\
        --num_workers 0 \\
        > noisy.log 2>&1 || true

    mkdir -p pdb_out
    python3 - <<PY
import csv, glob, json, math, os

his285_pos = ${params.his285_pos}
his290_pos = ${params.his290_pos}
glu330_pos = ${params.glu330_pos}

rows = list(csv.DictReader(open("${batch_csv}")))

def parse_pdb_atoms(path):
    by_chain = {}
    with open(path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            chain_id = line[21].strip()
            try:
                res_seq = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            by_chain.setdefault(chain_id, []).append((res_seq, atom_name, x, y, z))
    return by_chain

out_rows = []
for r in rows:
    sid = r["safe_id"]
    status = "failed"
    iptm = ptm = plddt = conf_score = ""
    affinity_pred_value = affinity_prob_binary = ""
    pocket_dist = ""
    pdb_rel = ""

    pred_dir = os.path.join("out", "boltz_results_yamls", "predictions", sid)
    conf_json = glob.glob(os.path.join(pred_dir, f"confidence_{sid}_model_0.json"))
    aff_json = glob.glob(os.path.join(pred_dir, f"affinity_{sid}.json"))
    pdb_file = glob.glob(os.path.join(pred_dir, f"{sid}_model_0.pdb"))

    if conf_json:
        try:
            cj = json.load(open(conf_json[0]))
            iptm = cj.get("iptm", "")
            ptm = cj.get("ptm", "")
            plddt = cj.get("complex_plddt", "")
            conf_score = cj.get("confidence_score", "")
        except Exception:
            pass

    if aff_json:
        try:
            aj = json.load(open(aff_json[0]))
            affinity_pred_value = aj.get("affinity_pred_value", "")
            affinity_prob_binary = aj.get("affinity_probability_binary", "")
        except Exception:
            pass

    if pdb_file:
        status = "ok"
        dest = os.path.join("pdb_out", sid + ".pdb")
        with open(pdb_file[0]) as fsrc, open(dest, "w") as fdst:
            fdst.write(fsrc.read())
        pdb_rel = dest

        by_chain = parse_pdb_atoms(pdb_file[0])
        protein_atoms = by_chain.get("A", [])
        ligand_atoms = []
        for cid, atoms in by_chain.items():
            if cid != "A":
                ligand_atoms.extend(atoms)

        def ca_coord(res_seq):
            for rs, aname, x, y, z in protein_atoms:
                if rs == res_seq and aname == "CA":
                    return (x, y, z)
            return None

        c285 = ca_coord(his285_pos)
        c290 = ca_coord(his290_pos)
        c330 = ca_coord(glu330_pos)

        if c285 and c290 and c330 and ligand_atoms:
            lx = sum(a[2] for a in ligand_atoms) / len(ligand_atoms)
            ly = sum(a[3] for a in ligand_atoms) / len(ligand_atoms)
            lz = sum(a[4] for a in ligand_atoms) / len(ligand_atoms)
            cx = (c285[0] + c290[0] + c330[0]) / 3.0
            cy = (c285[1] + c290[1] + c330[1]) / 3.0
            cz = (c285[2] + c290[2] + c330[2]) / 3.0
            pocket_dist = math.sqrt((lx - cx) ** 2 + (ly - cy) ** 2 + (lz - cz) ** 2)

    out_rows.append({
        "compound_safe_id": sid,
        "status": status,
        "iptm": iptm,
        "ptm": ptm,
        "complex_plddt": plddt,
        "confidence_score": conf_score,
        "affinity_pred_value": affinity_pred_value,
        "affinity_probability_binary": affinity_prob_binary,
        "boltz_pocket_dist": pocket_dist,
        "pdb_path": pdb_rel,
    })

shard = "${batch_csv.baseName}"
fieldnames = ["compound_safe_id","status","iptm","ptm","complex_plddt","confidence_score",
              "affinity_pred_value","affinity_probability_binary","boltz_pocket_dist","pdb_path"]
with open(f"summary_{shard}.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(out_rows)
PY

    # ensure the optional pdb glob output has at least a placeholder so
    # Nextflow doesn't error on a totally empty glob when every compound
    # in this shard failed (should not happen post shm-size fix, but
    # keep the shard robust to zero-success without aborting the run)
    if [ -z "\$(ls -A pdb_out 2>/dev/null)" ]; then
        touch pdb_out/.empty
    fi
    """
}
