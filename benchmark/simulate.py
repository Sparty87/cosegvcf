#!/usr/bin/env python3
"""Simulation benchmark for CoSegVCF.

For each scenario (pedigree x inheritance model) the script
  1. draws a pool of background sites with population allele frequencies,
  2. assigns genotypes to founders under Hardy-Weinberg equilibrium and transmits
     alleles to offspring following Mendelian rules (X chromosome hemizygous in males),
  3. spikes one causal event consistent with the pedigree phenotypes,
  4. optionally adds genotype errors and low-depth calls,
  5. writes one variant-only VCF per individual and a PED file,
  6. runs CoSegVCF and records whether the causal variant(s) are retained and how many
     candidate variants are reported.

Usage:
    python benchmark/simulate.py --reps 50 --out benchmark/results
    # or in parallel, then merge:
    python benchmark/simulate.py --reps 50 --only 0,2,4,6 &
    python benchmark/simulate.py --reps 50 --only 1,3,5,7 &
    wait; python benchmark/simulate.py --merge
Outputs: results*.tsv (one row per replicate) and summary.tsv (one row per scenario).
"""

import argparse
import csv
import os
import random
import shutil
import statistics
import sys
import tempfile
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from cosegvcf.core import load_samples, read_ped, run_analysis  # noqa: E402

AUTOSOMES = ["chr%d" % i for i in range(1, 23)]

# pedigrees: name -> list of (id, sex, status, father, mother)
PEDIGREES = {
    "trio": [
        ("F", "M", "unaffected", "0", "0"),
        ("M", "F", "unaffected", "0", "0"),
        ("P", "M", "affected", "F", "M"),
    ],
    "quartet": [
        ("F", "M", "unaffected", "0", "0"),
        ("M", "F", "unaffected", "0", "0"),
        ("P", "M", "affected", "F", "M"),
        ("S", "F", "unaffected", "F", "M"),
    ],
    "x_family": [
        ("F", "M", "unaffected", "0", "0"),
        ("M", "F", "unaffected", "0", "0"),
        ("P", "M", "affected", "F", "M"),
        ("B", "M", "unaffected", "F", "M"),
        ("S", "F", "unaffected", "F", "M"),
    ],
    "three_gen": [
        ("GF", "M", "affected", "0", "0"),
        ("GM", "F", "unaffected", "0", "0"),
        ("F", "M", "affected", "GF", "GM"),
        ("U", "M", "unaffected", "GF", "GM"),
        ("M", "F", "unaffected", "0", "0"),
        ("P", "M", "affected", "F", "M"),
        ("S", "F", "unaffected", "F", "M"),
    ],
}

SCENARIOS = [  # (label, pedigree, model)
    ("AD, 3 generations (7)", "three_gen", "AD"),
    ("AR hom, trio", "trio", "AR"),
    ("AR hom, quartet", "quartet", "AR"),
    ("Comp. het, trio", "trio", "CH"),
    ("Comp. het, quartet", "quartet", "CH"),
    ("De novo, trio", "trio", "DN"),
    ("De novo, quartet", "quartet", "DN"),
    ("XLR, 5 members", "x_family", "XLR"),
]


def draw_af(rng, rare_fraction):
    """Mixture spectrum: rare sites log-uniform in 1e-5..1e-2, common sites uniform 0.01..0.5."""
    if rng.random() < rare_fraction:
        return 10 ** rng.uniform(-5, -2)
    return rng.uniform(0.01, 0.5)


def make_sites(rng, n_sites, rare_fraction, x_fraction=0.04, sites_per_gene=4):
    sites = []
    for i in range(n_sites):
        chrom = "chrX" if rng.random() < x_fraction else rng.choice(AUTOSOMES)
        # keep X sites outside the GRCh38 PAR regions
        pos = rng.randint(3_000_000, 150_000_000)
        sites.append([chrom, pos, draw_af(rng, rare_fraction)])
    sites.sort(key=lambda s: (s[0], s[1]))
    uniq, seen = [], set()
    for s in sites:
        if (s[0], s[1]) not in seen:
            seen.add((s[0], s[1]))
            uniq.append(s)
    for i, s in enumerate(uniq):
        s.append("GENE%05d" % (i // sites_per_gene))
    return uniq  # [chrom, pos, af, gene]


def haplotypes(rng, ped, sites):
    """Return {id: [(a1, a2) per site]} with a2 = None for male X."""
    hap = {}
    for pid, s, _, fa, mo in ped:  # founders are listed before their children
        h = []
        for chrom, _, af, _ in sites:
            hemi = chrom == "chrX" and s == "M"
            if fa == "0":
                a1 = int(rng.random() < af)
                a2 = None if hemi else int(rng.random() < af)
            else:
                ma = hap[mo][len(h)]
                a_mat = ma[rng.randrange(2)] if ma[1] is not None else ma[0]
                if chrom == "chrX":
                    if hemi:
                        a1, a2 = a_mat, None
                    else:  # daughters receive the father's only X
                        a1, a2 = a_mat, hap[fa][len(h)][0]
                else:
                    pa = hap[fa][len(h)]
                    a1, a2 = a_mat, pa[rng.randrange(2)]
            h.append((a1, a2))
        hap[pid] = h
    return hap


def dosage(h):
    return h[0] + (h[1] or 0)


def spike(rng, model, ped, sites, hap):
    """Add a causal event consistent with the pedigree. Returns causal site keys."""
    status = {p[0]: p[2] for p in ped}
    ids = [p[0] for p in ped]

    def new_site(chrom, gene):
        pos = rng.randint(3_000_000, 150_000_000)
        sites.append([chrom, pos, 10 ** rng.uniform(-5, -3.5), gene])
        for pid in ids:
            hap[pid].append((0, None if chrom == "chrX" and
                             [p for p in ped if p[0] == pid][0][1] == "M" else 0))
        return len(sites) - 1

    if model == "AD":
        i = new_site(rng.choice(AUTOSOMES), "CAUSAL_AD")
        for pid in ids:
            if status[pid] == "affected":
                hap[pid][i] = (1, 0)
        return [i]
    if model == "AR":
        i = new_site(rng.choice(AUTOSOMES), "CAUSAL_AR")
        for pid in ids:
            if status[pid] == "affected":
                hap[pid][i] = (1, 1)
            elif pid in ("F", "M"):
                hap[pid][i] = (1, 0)
            else:  # unaffected sibling: carrier with p = 2/3, otherwise reference
                hap[pid][i] = (1, 0) if rng.random() < 2 / 3 else (0, 0)
        return [i]
    if model == "CH":
        chrom = rng.choice(AUTOSOMES)
        i, j = new_site(chrom, "CAUSAL_CH"), new_site(chrom, "CAUSAL_CH")
        sib = rng.choice([(1, 0), (0, 1), (0, 0)])  # never both (unaffected)
        for pid in ids:
            if status[pid] == "affected":
                hap[pid][i], hap[pid][j] = (1, 0), (0, 1)
            elif pid == "F":
                hap[pid][i] = (1, 0)
            elif pid == "M":
                hap[pid][j] = (1, 0)
            else:
                hap[pid][i], hap[pid][j] = (sib[0], 0), (sib[1], 0)
        return [i, j]
    if model == "DN":
        i = new_site(rng.choice(AUTOSOMES), "CAUSAL_DN")
        sites[i][2] = 0.0
        hap["P"][i] = (1, 0)
        return [i]
    if model == "XLR":
        i = new_site("chrX", "CAUSAL_XLR")
        hap["M"][i] = (1, 0)
        hap["P"][i] = (1, None)
        if "S" in hap:
            hap["S"][i] = (rng.randrange(2), 0)
        return [i]
    raise ValueError(model)


def write_vcfs(rng, directory, ped, sites, hap, err, miss):
    header = ("##fileformat=VCFv4.2\n"
              '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="Population AF">\n'
              '##INFO=<ID=ANN,Number=.,Type=String,Description="Annotation">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
              '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n'
              + "".join("##contig=<ID=%s>\n" % c for c in AUTOSOMES + ["chrX"]))
    order = sorted(range(len(sites)), key=lambda k: (sites[k][0], sites[k][1]))
    paths, rare = [], []
    for pid, sex, *_ in ped:
        n_rare = 0
        p = os.path.join(directory, pid + ".vcf")
        with open(p, "w") as fh:
            fh.write(header)
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % pid)
            for k in order:
                chrom, pos, af, gene = sites[k]
                h = hap[pid][k]
                d = dosage(h)
                hemi = h[1] is None
                if err and rng.random() < err:  # genotype error
                    if hemi:
                        d = 1 - d
                    else:
                        d = {0: 1, 1: rng.choice((0, 2)), 2: 1}[d]
                if d == 0:
                    continue
                gt = "1" if hemi else {1: "0/1", 2: "1/1"}[d]
                dp = rng.randint(3, 7) if (miss and rng.random() < miss) else rng.randint(15, 80)
                gq = 99 if dp >= 15 else 10
                n_rare += af <= 0.01
                fh.write("%s\t%d\t.\tA\tG\t60\tPASS\tgnomAD_AF=%.6g;ANN=G|missense_variant|"
                         "MODERATE|%s|x\tGT:DP:GQ\t%s:%d:%d\n" % (chrom, pos, af, gene, gt, dp, gq))
        paths.append(p)
        rare.append(n_rare)
    ped_path = os.path.join(directory, "family.ped")
    with open(ped_path, "w") as fh:
        for pid, sex, status, fa, mo in ped:
            fh.write("F\t%s\t%s\t%s\t%s\t%s\n" % (pid, fa, mo, "1" if sex == "M" else "2",
                                                  "2" if status == "affected" else "1"))
    return paths, ped_path, sum(rare) / len(rare)


def run_once(rng, ped_name, model, n_sites, rare_fraction, err, miss, ignore_missing, max_af, tmp):
    ped = PEDIGREES[ped_name]
    sites = make_sites(rng, n_sites, rare_fraction)
    hap = haplotypes(rng, ped, sites)
    causal = spike(rng, model, ped, sites, hap)
    causal_keys = {(sites[i][0], sites[i][1]) for i in causal}
    paths, ped_path, rare = write_vcfs(rng, tmp, ped, sites, hap, err, miss)
    samples = load_samples(paths)
    read_ped(ped_path, samples)
    out = os.path.join(tmp, "out.vcf")
    cfg = {"model": model, "max_af": max_af, "ignore_missing": ignore_missing,
           "af_fields": "gnomAD_AF", "output": out}
    t0 = time.perf_counter()
    run_analysis(samples, cfg, log=lambda m: None)
    elapsed = time.perf_counter() - t0
    found = set()
    n_out = 0
    with open(out) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n_out += 1
            f = line.split("\t", 2)
            found.add((f[0], int(f[1])))
    n_input = sum(1 for p in paths for l in open(p) if not l.startswith("#"))
    return {
        "recovered": int(causal_keys <= found),
        "candidates": n_out,
        "background": n_out - len(causal_keys & found),
        "records_per_sample": n_input / len(paths),
        "rare_per_sample": rare,
        "seconds": elapsed,
    }


KEYS = ["scenario", "model", "pedigree", "condition", "rep", "recovered", "candidates",
        "background", "records_per_sample", "rare_per_sample", "seconds"]

CONDITIONS = [  # label, genotype error rate, low-depth rate, ignore_missing
    ("error-free", 0.0, 0.0, False),
    ("noisy", 0.002, 0.02, False),
    ("noisy, missing ignored", 0.002, 0.02, True),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--sites", type=int, default=200000,
                    help="background sites in the population pool per replicate")
    ap.add_argument("--rare-fraction", type=float, default=0.78,
                    help="fraction of pool sites with allele frequency below 1%%")
    ap.add_argument("--max-af", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=20260930)
    ap.add_argument("--out", default="benchmark/results")
    ap.add_argument("--only", default="",
                    help="comma-separated scenario indices (0-%d), to run in parallel" % (len(SCENARIOS) - 1))
    ap.add_argument("--merge", action="store_true",
                    help="only merge results_*.tsv in --out into summary.tsv")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    only = [int(x) for x in a.only.split(",")] if a.only else range(len(SCENARIOS))
    if not a.merge:
        rows = []
        for idx in only:
            label, ped_name, model = SCENARIOS[idx]
            rng = random.Random(a.seed + idx)  # one stream per scenario: order-independent
            for cond, err, miss, ign in CONDITIONS:
                for rep in range(a.reps):
                    tmp = tempfile.mkdtemp(prefix="coseg_")
                    try:
                        r = run_once(rng, ped_name, model, a.sites, a.rare_fraction, err, miss,
                                     ign, a.max_af, tmp)
                    finally:
                        shutil.rmtree(tmp)
                    r.update(scenario=label, model=model, pedigree=ped_name, condition=cond,
                             rep=rep)
                    rows.append(r)
                sub = [r for r in rows if r["scenario"] == label and r["condition"] == cond]
                print("%-24s %-24s recall %.2f  median background %s" % (
                    label, cond, statistics.mean(r["recovered"] for r in sub),
                    statistics.median(r["background"] for r in sub)), flush=True)
        name = "results_%s.tsv" % "-".join(map(str, only)) if a.only else "results.tsv"
        with open(os.path.join(a.out, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, KEYS, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        if a.only:
            return
    rows = []
    for fn in sorted(os.listdir(a.out)):
        if fn.startswith("results") and fn.endswith(".tsv") and (a.merge or fn == "results.tsv"):
            with open(os.path.join(a.out, fn)) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    for k in ("recovered", "background", "rep"):
                        r[k] = int(r[k])
                    for k in ("records_per_sample", "rare_per_sample", "seconds"):
                        r[k] = float(r[k])
                    rows.append(r)
    write_summary(rows, a.out)


def write_summary(rows, out):
    with open(os.path.join(out, "summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["scenario", "condition", "n", "sensitivity", "median_background",
                    "q1_background", "q3_background", "mean_records_per_sample",
                    "mean_rare_per_sample", "median_seconds"])
        for label, _, _ in SCENARIOS:
            for cond, *_ in CONDITIONS:
                sub = [r for r in rows if r["scenario"] == label and r["condition"] == cond]
                if not sub:
                    continue
                bg = sorted(r["background"] for r in sub)
                q = statistics.quantiles(bg, n=4) if len(bg) > 1 else [bg[0]] * 3
                w.writerow([label, cond, len(sub),
                            "%.3f" % statistics.mean(r["recovered"] for r in sub),
                            statistics.median(bg), q[0], q[2],
                            "%.0f" % statistics.mean(r["records_per_sample"] for r in sub),
                            "%.0f" % statistics.mean(r["rare_per_sample"] for r in sub),
                            "%.2f" % statistics.median(r["seconds"] for r in sub)])

if __name__ == "__main__":
    main()
