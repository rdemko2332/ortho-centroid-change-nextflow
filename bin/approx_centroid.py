#!/usr/bin/env python3
"""
approx_centroid.py - approximate the centroid (medoid) protein of an orthogroup.

The centroid is the member with the highest mean normalised similarity to the
rest of the group. Large groups use two rounds of random sub-sampling instead
of a full all-vs-all.

  Groups <= --exact-limit proteins:
      DIAMOND all-vs-all; exact medoid.

  Larger groups (all sets drawn at random, without replacement, mutually disjoint):
      R1 = --r1 reference proteins                         (default 50)
      C  = --cand-fraction of the group as candidates      (default 10%, capped by --max-cand)
      Round 1: DIAMOND  C  vs R1 -> keep the --top-k best candidates (default 50)
      Round 2: DIAMOND top vs R2 = --r2 more references    (default 400)
      Centroid = finalist with the highest mean similarity over R1 + R2.

Similarity(a, b) = bitscore(a, b) / max(selfbits(a), selfbits(b)); pairs with no
DIAMOND hit count as 0. Self bitscores come from one DIAMOND run of each protein
against itself, using the same DIAMOND settings.

Requirements: Python >= 3.8 (standard library only) and `diamond` on PATH.

Usage:
  approx_centroid.py OG401.fasta --threads 8 > OG401.centroid.tsv
  approx_centroid.py OG401.fasta --group-id OG6_100401 --details OG401.candidates.tsv

Output (stdout, one row, with header):
  group_id  n_proteins  method  centroid_id  mean_sim  n_candidates  n_refs  n_alignments  seed

The random seed is derived from the group id, so reruns give the same answer;
override with --seed.
"""
import argparse
import hashlib
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


# ----------------------------------------------------------------------------- I/O
def read_fasta(path):
    """Return list of (id, sequence); id = first word of the header."""
    recs, pid, chunks = [], None, []
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = lambda p: gzip.open(p, "rt")
    with opener(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if pid is not None:
                    recs.append((pid, "".join(chunks)))
                pid, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if pid is not None:
        recs.append((pid, "".join(chunks)))
    ids = [r[0] for r in recs]
    if len(ids) != len(set(ids)):
        sys.exit(f"ERROR: duplicate sequence ids in {path}")
    return recs


def write_fasta(path, ids, seqs):
    with open(path, "w") as fh:
        for i in ids:
            fh.write(f">{i}\n{seqs[i]}\n")


def log(msg):
    print(f"[approx_centroid] {msg}", file=sys.stderr, flush=True)


# ----------------------------------------------------------------------------- DIAMOND
class Diamond:
    def __init__(self, workdir, threads, extra_args):
        self.workdir, self.threads = workdir, threads
        self.extra = extra_args.split()
        self.n_runs = 0

    def search(self, query_ids, ref_ids, seqs, tag, max_targets=0):
        """Return {(query, ref): best bitscore} for query_ids vs ref_ids."""
        q = os.path.join(self.workdir, f"{tag}.query.fa")
        r = os.path.join(self.workdir, f"{tag}.ref.fa")
        db = os.path.join(self.workdir, f"{tag}.ref")
        out = os.path.join(self.workdir, f"{tag}.hits.tsv")
        write_fasta(q, query_ids, seqs)
        write_fasta(r, ref_ids, seqs)
        self._run(["makedb", "--in", r, "-d", db])
        self._run(["blastp", "-q", q, "-d", db, "-o", out,
                   "-k", str(max_targets), "--max-hsps", "1", "-p", str(self.threads),
                   *self.extra, "--outfmt", "6", "qseqid", "sseqid", "bitscore"])
        hits = {}
        with open(out) as fh:
            for line in fh:
                a, b, s = line.rstrip("\n").split("\t")
                s = float(s)
                if s > hits.get((a, b), 0.0):
                    hits[(a, b)] = s
        self.n_runs += 1
        return hits

    def _run(self, args):
        cmd = ["diamond", *args, "--quiet"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            sys.exit(f"ERROR: DIAMOND failed: {' '.join(cmd)}\n{res.stderr}")


def self_bitscores(dmd, ids, seqs, tag):
    """Self-alignment bitscore for each id. Each protein is searched against a
    database containing it; the self hit is kept. If DIAMOND reports no self hit
    (e.g. a fully low-complexity-masked sequence), fall back to the protein's best
    hit to anything, or 1.0 if it has none (its similarities are then all 0)."""
    hits = dmd.search(ids, ids, seqs, tag, max_targets=1)   # top hit = self (or an identical copy)
    best_any = defaultdict(float)
    for (a, _), s in hits.items():
        best_any[a] = max(best_any[a], s)
    return {i: hits.get((i, i)) or best_any.get(i) or 1.0 for i in ids}


def norm_sim(bits, sa, sb):
    return bits / max(sa, sb)


# ----------------------------------------------------------------------------- methods
def exact_medoid(ids, seqs, dmd):
    hits = dmd.search(ids, ids, seqs, "exact")
    selfb = {i: hits.get((i, i), 0.0) for i in ids}
    best_any = defaultdict(float)
    for (a, b), s in hits.items():
        if a != b:
            best_any[a] = max(best_any[a], s)
    for i in ids:
        if selfb[i] == 0.0:
            selfb[i] = best_any.get(i) or 1.0
    # symmetrise: best score in either direction
    pair = {}
    for (a, b), s in hits.items():
        if a == b:
            continue
        k = (a, b) if a < b else (b, a)
        if s > pair.get(k, 0.0):
            pair[k] = s
    total = defaultdict(float)
    for (a, b), s in pair.items():
        v = norm_sim(s, selfb[a], selfb[b])
        total[a] += v
        total[b] += v
    n = len(ids)
    scores = {i: total[i] / (n - 1) for i in ids}
    best = max(ids, key=lambda i: (scores[i], i))
    return best, scores, n * n


def two_round_medoid(ids, seqs, dmd, a, rng):
    n = len(ids)
    order = ids[:]
    rng.shuffle(order)

    n_r1 = min(a.r1, n // 4)
    n_c = min(max(math.ceil(a.cand_fraction * n), a.top_k), a.max_cand, n - n_r1)
    r1 = order[:n_r1]
    cand = order[n_r1:n_r1 + n_c]
    rest = order[n_r1 + n_c:]
    r2 = rest[:a.r2]
    if not r2:
        log("WARNING: no proteins left for round-2 references; using round 1 only")

    selfb = self_bitscores(dmd, r1 + cand + r2, seqs, "self")

    # Round 1: all candidates vs R1
    h1 = dmd.search(cand, r1, seqs, "round1")
    s1 = defaultdict(float)
    for (c, r), bits in h1.items():
        s1[c] += norm_sim(bits, selfb[c], selfb[r])
    top = sorted(cand, key=lambda c: (-s1[c], c))[:a.top_k]

    # Round 2: finalists vs R2
    s2 = defaultdict(float)
    if r2:
        h2 = dmd.search(top, r2, seqs, "round2")
        for (c, r), bits in h2.items():
            s2[c] += norm_sim(bits, selfb[c], selfb[r])

    n_refs = len(r1) + len(r2)
    final = {c: (s1[c] + s2[c]) / n_refs for c in top}
    best = max(top, key=lambda c: (final[c], c))
    details = {c: (s1[c] / len(r1), final.get(c), c in top) for c in cand}
    n_aln = len(cand) * len(r1) + len(top) * len(r2)
    return best, final[best], len(cand), n_refs, n_aln, details


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="FASTA of one orthogroup (may be .gz)")
    ap.add_argument("--group-id", help="default: FASTA file name without extension")
    ap.add_argument("--r1", type=int, default=50, help="round-1 reference set size")
    ap.add_argument("--cand-fraction", type=float, default=0.10,
                    help="fraction of the group tested as candidates")
    ap.add_argument("--max-cand", type=int, default=3000,
                    help="cap on the number of candidates (0 = no cap)")
    ap.add_argument("--top-k", type=int, default=50, help="candidates kept for round 2")
    ap.add_argument("--r2", type=int, default=400, help="round-2 reference set size")
    ap.add_argument("--exact-limit", type=int, default=2000,
                    help="groups up to this size use the exact all-vs-all")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--diamond-args", default="--very-sensitive -e 1e-3",
                    help="extra DIAMOND blastp arguments (match production)")
    ap.add_argument("--seed", type=int, help="default: derived from the group id")
    ap.add_argument("--details", help="optional TSV of per-candidate scores")
    ap.add_argument("--workdir", help="temp directory (default: system temp, removed after)")
    ap.add_argument("--keep-workdir", action="store_true")
    ap.add_argument("--emit-centroid-seq",
                    help="append the chosen centroid's sequence to this FASTA, headered by group_id "
                         "(so a batch of groups can accumulate one combined new-centroids FASTA)")
    a = ap.parse_args()
    if a.max_cand <= 0:
        a.max_cand = sys.maxsize

    if shutil.which("diamond") is None:
        sys.exit("ERROR: diamond not found on PATH")

    recs = read_fasta(a.fasta)
    if not recs:
        sys.exit(f"ERROR: no sequences in {a.fasta}")
    ids = [r[0] for r in recs]
    seqs = dict(recs)
    gid = a.group_id or os.path.basename(a.fasta).split(".")[0]
    seed = a.seed if a.seed is not None else \
        int(hashlib.sha256(gid.encode()).hexdigest()[:8], 16)
    n = len(ids)

    workdir = a.workdir or tempfile.mkdtemp(prefix=f"centroid_{gid}_")
    os.makedirs(workdir, exist_ok=True)
    dmd = Diamond(workdir, a.threads, a.diamond_args)
    try:
        if n == 1:
            method, best, score, n_c, n_refs, n_aln, details = "single", ids[0], 1.0, 1, 0, 0, None
        elif n <= a.exact_limit:
            log(f"{gid}: {n} proteins -> exact all-vs-all")
            best, scores, n_aln = exact_medoid(ids, seqs, dmd)
            method, score, n_c, n_refs = "exact", scores[best], n, n - 1
            details = {i: (None, scores[i], True) for i in ids}
        else:
            log(f"{gid}: {n} proteins -> two-round sampling (seed {seed})")
            best, score, n_c, n_refs, n_aln, details = two_round_medoid(
                ids, seqs, dmd, a, random.Random(seed))
            method = "two_round"
    finally:
        if not a.keep_workdir and not a.workdir:
            shutil.rmtree(workdir, ignore_errors=True)

    print("group_id\tn_proteins\tmethod\tcentroid_id\tmean_sim\tn_candidates\tn_refs\tn_alignments\tseed")
    print(f"{gid}\t{n}\t{method}\t{best}\t{score:.4f}\t{n_c}\t{n_refs}\t{n_aln}\t{seed}")

    if a.emit_centroid_seq:
        with open(a.emit_centroid_seq, "a") as fh:
            fh.write(f">{gid}\n{seqs[best]}\n")

    if a.details and details:
        with open(a.details, "w") as fh:
            fh.write("protein\tround1_mean_sim\tfinal_mean_sim\tfinalist\n")
            for p, (r1, fin, is_top) in sorted(details.items(),
                                                key=lambda kv: -(kv[1][1] or kv[1][0] or 0)):
                fh.write(f"{p}\t{'' if r1 is None else f'{r1:.4f}'}\t"
                         f"{'' if fin is None else f'{fin:.4f}'}\t{int(is_top)}\n")
    log(f"{gid}: centroid {best} (mean sim {score:.4f}, {method})")


if __name__ == "__main__":
    main()
