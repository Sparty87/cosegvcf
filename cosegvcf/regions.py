"""Region-level outputs: BED files of genomic intervals shared by the affected individuals.

Two independent exports are provided:

* callable regions - intervals confidently genotyped (reference blocks or variant calls
  passing the DP/GQ thresholds) in *every* affected individual. This needs gVCF input for
  all affected samples, because only reference blocks tell where a sample is 0/0 with
  evidence rather than simply not covered.

* shared segments - runs of consecutive informative sites at which the genotypes of the
  affected individuals are compatible with a shared haplotype (identity by state):
    ibs1  no pair of affected individuals is homozygous for opposite alleles
          (at least one allele shared; dominant models)
    ibs2  all affected individuals have the same genotype (recessive models)
    roh   all affected individuals are homozygous for the same allele
          (autozygosity; works with a single affected individual)
"""

from collections import defaultdict

from .core import PAR_REGIONS, chrom_class, natural_key

SHARED_MODES = ("auto", "ibs1", "ibs2", "roh")
RECESSIVE_MODELS = ("AR", "CH", "AR_ANY", "XLR")
# a discordant site is tolerated as a probable genotype error when it is preceded by at
# least this many concordant sites since the previous discordant one
ISOLATED_GAP = 10


def resolve_mode(mode, model, n_affected):
    """Translate 'auto' into a concrete sharing mode; return None if not applicable."""
    if mode == "auto":
        if model in RECESSIVE_MODELS:
            mode = "ibs2" if n_affected > 1 else "roh"
        else:
            mode = "ibs1"
    if mode in ("ibs1", "ibs2") and n_affected < 2:
        return None
    return mode


def site_state(dosages, mode):
    """True (concordant), False (discordant) or None (not informative) for one site.

    `dosages` are alternate-allele dosages of the affected individuals, with hemizygous
    male chrX calls already recoded as homozygous.
    """
    if any(d is None for d in dosages):
        return None
    if mode == "ibs1":
        return not (0 in dosages and 2 in dosages)
    if mode == "ibs2":
        return len(set(dosages)) == 1
    if mode == "roh":
        return len(set(dosages)) == 1 and dosages[0] != 1
    raise ValueError(mode)


def segments(sites, min_sites=20, min_bp=100000, tolerate_isolated=True):
    """Collapse per-site states into segments.

    `sites` is a position-sorted list of (pos, state) for one chromosome, with state True
    or False. A segment spans the first to the last concordant site of a run. A single
    discordant site surrounded by concordant ones is tolerated and counted; two discordant
    sites within ISOLATED_GAP concordant sites end the segment.

    Returns a list of (start, end, n_concordant, n_discordant), 1-based inclusive.
    """
    out = []
    start = last_ok = None
    n = disc = since = 0

    def close():
        nonlocal start
        if start is not None:
            d = disc - (1 if disc and since == 0 else 0)  # trailing discordant is the border
            if n >= min_sites and last_ok - start + 1 >= min_bp:
                out.append((start, last_ok, n, d))
        start = None

    for pos, ok in sites:
        if ok:
            if start is None:
                start, n, disc, since = pos, 0, 0, ISOLATED_GAP
            n += 1
            since += 1
            last_ok = pos
        elif start is not None:
            if tolerate_isolated and since >= ISOLATED_GAP:
                disc += 1
                since = 0
            else:
                close()
    close()
    return out


def shared_segments(engine, mode, min_sites, min_bp):
    """Shared segments over all variant sites of the store. Returns {chrom: [segments]}."""
    by_chrom = defaultdict(list)
    aff = engine.aff
    for key, v in engine.store.vars.items():
        chrom, pos = key[0], key[1]
        cc = chrom_class(chrom, pos, engine.par)
        if cc in ("Y", "MT"):
            continue
        dos = []
        for s in aff:
            d = engine.genotype(v, key, s.name)
            if cc == "X" and d is not None:
                if s.sex == "M":
                    d = 2 if d >= 1 else 0
                elif s.sex != "F":
                    d = None  # unknown sex on chrX: not interpretable
            dos.append(d)
        st = site_state(dos, mode)
        if st is not None:
            by_chrom[chrom].append((pos, st))
    res = {}
    for chrom, lst in by_chrom.items():
        lst.sort()
        # several alleles at one position: the site is concordant only if all are
        merged = []
        for pos, st in lst:
            if merged and merged[-1][0] == pos:
                merged[-1] = (pos, merged[-1][1] and st)
            else:
                merged.append((pos, st))
        segs = segments(merged, min_sites, min_bp)
        if segs:
            res[chrom] = segs
    return res


def _merge(intervals):
    intervals.sort()
    out = []
    for s, e in intervals:
        if out and s <= out[-1][1] + 1:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return out


def _intersect(a, b):
    out = []
    i = j = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s <= e:
            out.append([s, e])
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def callable_regions(engine):
    """Intervals confidently genotyped in every affected individual.

    Returns ({chrom: [[start, end], ...]}, missing) where `missing` lists the affected
    samples without gVCF reference blocks (in which case the dict is empty).
    """
    store = engine.store
    missing = [s.name for s in engine.aff if s.name not in store.refblocks]
    if missing or not engine.aff:
        return {}, missing
    per_sample = {s.name: defaultdict(list) for s in engine.aff}
    for s in engine.aff:
        for chrom, (_, blocks) in store.refblocks[s.name].items():
            for start, end, dp, gq in blocks:
                if engine._qual_ok(dp, gq):
                    per_sample[s.name][chrom].append((start, end))
    for key, v in store.vars.items():
        chrom, pos, ref = key[0], key[1], key[2]
        for s in engine.aff:
            if s.name in v["gts"] and engine.genotype(v, key, s.name) is not None:
                per_sample[s.name][chrom].append((pos, pos + len(ref) - 1))
    result = {}
    names = [s.name for s in engine.aff]
    chroms = set(per_sample[names[0]])
    for n in names[1:]:
        chroms &= set(per_sample[n])
    for chrom in chroms:
        cur = _merge(per_sample[names[0]][chrom])
        for n in names[1:]:
            cur = _intersect(cur, _merge(per_sample[n][chrom]))
            if not cur:
                break
        if cur:
            result[chrom] = cur
    return result, missing


def _chrom_order(store, chroms):
    order = {c: i for i, c in enumerate(store.contigs)}
    return sorted(chroms, key=lambda c: (order.get(c, 10 ** 6), natural_key(c)))


def write_region_beds(engine, cfg, base, log=lambda m: None):
    """Write the BED files requested in `cfg`. Returns the list of paths written."""
    written = []
    store = engine.store
    aff_names = ",".join(s.name for s in engine.aff)

    if cfg.get("bed_shared"):
        mode = resolve_mode(cfg.get("bed_shared_mode", "auto"), cfg["model"], len(engine.aff))
        if mode is None:
            log("WARNING: shared-segment BED skipped: IBS sharing needs at least two "
                "affected individuals (use mode 'roh' for a single affected individual).")
        else:
            if not any(s.name in store.refblocks for s in engine.aff):
                log("NOTE: no gVCF among the affected samples; shared segments rely on "
                    "variant sites only, with absent calls read as set in the options.")
            min_sites = int(cfg.get("bed_min_sites", 20))
            min_bp = int(float(cfg.get("bed_min_kb", 100)) * 1000)
            segs = shared_segments(engine, mode, min_sites, min_bp)
            path = "%s.shared_%s.bed" % (base, mode)
            total = n_seg = 0
            with open(path, "w") as fh:
                fh.write("#CoSegVCF shared segments; mode=%s; affected=%s; min_sites=%d; "
                         "min_kb=%g\n" % (mode, aff_names, min_sites, min_bp / 1000.0))
                fh.write("#chrom\tstart\tend\tname\tconcordant_sites\tdiscordant_sites\t"
                         "length_bp\n")
                for chrom in _chrom_order(store, segs):
                    for start, end, n, d in segs[chrom]:
                        fh.write("%s\t%d\t%d\t%s\t%d\t%d\t%d\n"
                                 % (chrom, start - 1, end, mode.upper(), n, d, end - start + 1))
                        total += end - start + 1
                        n_seg += 1
            log("Shared segments (%s): %d segments, %.2f Mb -> %s"
                % (mode, n_seg, total / 1e6, path))
            written.append(path)

    if cfg.get("bed_callable"):
        regions, missing = callable_regions(engine)
        if missing:
            log("WARNING: callable-region BED skipped: no gVCF reference blocks for "
                "affected sample(s) %s." % ", ".join(missing))
        else:
            path = "%s.callable_affected.bed" % base
            total = n_int = 0
            with open(path, "w") as fh:
                fh.write("#CoSegVCF regions callable in all affected; affected=%s; "
                         "min_dp=%d; min_gq=%d\n" % (aff_names, cfg["min_dp"], cfg["min_gq"]))
                for chrom in _chrom_order(store, regions):
                    for start, end in regions[chrom]:
                        fh.write("%s\t%d\t%d\n" % (chrom, start - 1, end))
                        total += end - start + 1
                        n_int += 1
            log("Callable in all affected: %d intervals, %.2f Mb -> %s"
                % (n_int, total / 1e6, path))
            written.append(path)
    return written


__all__ = ["SHARED_MODES", "PAR_REGIONS", "resolve_mode", "site_state", "segments",
           "shared_segments", "callable_regions", "write_region_beds"]
