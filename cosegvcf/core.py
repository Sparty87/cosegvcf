"""Core engine of CoSegVCF: VCF/gVCF parsing, pedigree checks and inheritance-model filters.

This module has no third-party dependencies and no GUI code, so it can be imported
by the command-line interface, by the Tkinter front end and by the test suite.
"""

import bisect
import csv
import datetime
import gzip
import os
import re
from collections import OrderedDict, defaultdict

from . import __version__

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------
PAR_REGIONS = OrderedDict([
    ("GRCh38", [(10001, 2781479), (155701383, 156030895)]),
    ("GRCh37", [(60001, 2699520), (154931044, 155260560)]),
    ("none", []),
])

MODELS = OrderedDict([
    ("AD", "Autosomal dominant"),
    ("AR", "Autosomal recessive - homozygous"),
    ("CH", "Autosomal recessive - compound heterozygous"),
    ("AR_ANY", "Autosomal recessive - homozygous + compound heterozygous"),
    ("XLR", "X-linked recessive"),
    ("XLD", "X-linked dominant"),
    ("DN", "De novo"),
    ("MT", "Mitochondrial (maternal)"),
    ("SHARED", "No model: shared by affected, absent in unaffected"),
])

MODEL_HELP = {
    "AD": "Autosomes and PAR. Affected: carriers (0/1 or 1/1). Unaffected: 0/0 "
          "(unless incomplete penetrance is allowed).",
    "AR": "Autosomes and PAR. Affected: 1/1. Unaffected: not 1/1. Unaffected parents of "
          "affected individuals present in the dataset: obligate carriers (if enabled).",
    "CH": "Pairs of 0/1 variants in the same gene carried by every affected individual. "
          "When parents are available, the two variants must come from different parents "
          "(in trans). No unaffected individual may carry both. Requires gene annotation "
          "(ANN, CSQ, Gene.refGene, GENEINFO...) or a gene BED file.",
    "AR_ANY": "Union of the homozygous and compound heterozygous recessive models.",
    "XLR": "chrX outside PAR only. Affected males: hemizygous; affected females: 1/1. "
           "Unaffected males: reference; unaffected females: at most 0/1. Mothers of "
           "affected males: obligate carriers. Sex is required for samples with known status.",
    "XLD": "chrX outside PAR only. Affected: carriers. Unaffected: reference (unless "
           "incomplete penetrance is allowed). Sex is required for samples with known status.",
    "DN": "Needs at least one affected individual with both parents in the dataset. "
          "Affected: carrier; both parents: called 0/0 passing quality thresholds (always "
          "enforced); other unaffected: reference.",
    "MT": "chrM/MT only. Affected: carriers; mothers of affected (if present): carriers. "
          "Unaffected individuals are not constrained (variable heteroplasmy).",
    "SHARED": "Any chromosome. Variants carried by all affected and absent in all unaffected "
              "(unless incomplete penetrance is allowed), with no model assumption.",
}

ROLES = ["Proband", "Father", "Mother", "Brother", "Sister", "Son", "Daughter",
         "Partner", "Paternal grandfather", "Paternal grandmother", "Maternal grandfather",
         "Maternal grandmother", "Paternal uncle", "Paternal aunt", "Maternal uncle",
         "Maternal aunt", "Cousin", "Other"]
ROLE_SEX = {"Father": "M", "Mother": "F", "Brother": "M", "Sister": "F", "Son": "M",
            "Daughter": "F", "Paternal grandfather": "M", "Paternal grandmother": "F",
            "Maternal grandfather": "M", "Maternal grandmother": "F", "Paternal uncle": "M",
            "Paternal aunt": "F", "Maternal uncle": "M", "Maternal aunt": "F"}
SEXES = ["M", "F", "?"]
STATUSES = ["affected", "unaffected", "unknown"]

DEFAULT_AF_FIELDS = ("gnomAD_AF,gnomADe_AF,gnomADg_AF,gnomAD_exome_ALL,gnomAD_genome_ALL,"
                     "AF_popmax,MAX_AF,ExAC_ALL")

DEFAULT_CONFIG = {
    "model": "AD",
    "build": "GRCh38",
    "min_dp": 8,
    "min_gq": 20,
    "pass_only": True,
    "absent_as_ref": True,
    "ignore_missing": False,
    "incomplete_penetrance": False,
    "obligate_carriers": True,
    "max_af": None,
    "af_fields": DEFAULT_AF_FIELDS,
    "gene_bed": "",
    "output": "cosegvcf_output.vcf",
    "bed_shared": False,        # BED of segments shared (IBS) by the affected individuals
    "bed_shared_mode": "auto",  # auto, ibs1, ibs2 or roh
    "bed_min_sites": 20,
    "bed_min_kb": 100,
    "bed_callable": False,      # BED of regions callable in all affected (gVCF only)
}

SKIP_ALTS = {"<NON_REF>", "<*>", "*", "."}
META_ID_RE = re.compile(r"^##(INFO|FORMAT|FILTER|ALT|contig)=<ID=([^,>]+)")
NUMBER_RE = re.compile(r"Number=([^,>]+)")
CSQ_FMT_RE = re.compile(r"Format: ([^\">]+)")
END_RE = re.compile(r"(?:^|;)END=(\d+)")
GT_SPLIT_RE = re.compile(r"[/|]")
VCF_SUFFIX_RE = re.compile(r"(\.g)?\.vcf(\.b?gz)?$|\.gvcf(\.b?gz)?$")


class CoSegError(Exception):
    """Raised for user-facing input or configuration errors."""


# ----------------------------------------------------------------------------
# Samples
# ----------------------------------------------------------------------------
class Sample:
    """One sequenced individual. `name` is unique and acts as identifier."""

    def __init__(self, name, path, col, vcf_name=None):
        self.name = name
        self.path = path
        self.col = col
        self.vcf_name = vcf_name or name
        self.role = "Other"
        self.sex = "?"
        self.father = "0"
        self.mother = "0"
        self.status = "unknown"

    @property
    def affected(self):
        """True, False or None (unknown)."""
        return {"affected": True, "unaffected": False}.get(self.status)

    def __repr__(self):
        return "Sample(%s, sex=%s, status=%s, father=%s, mother=%s)" % (
            self.name, self.sex, self.status, self.father, self.mother)


def open_text(path, mode="rt"):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, mode)
    return open(path, mode[0])


def vcf_samples(path):
    """Return the sample names in the #CHROM line of a VCF."""
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                names = line.rstrip("\r\n").split("\t")[9:]
                if names:
                    return names
                break
            if not line.startswith("#"):
                break
    raise CoSegError("No #CHROM header line or no sample columns in %s" % path)


def file_stem(path):
    return VCF_SUFFIX_RE.sub("", os.path.basename(path))


def load_samples(paths, existing=None):
    """Create Sample objects for every sample column of every file.

    Duplicate names across files are made unique by appending the file stem.
    """
    samples = OrderedDict((s.name, s) for s in (existing or []))
    added = []
    for p in paths:
        stem = file_stem(p)
        for col, n in enumerate(vcf_samples(p)):
            name = n
            if name in samples:
                name = "%s_%s" % (n, stem)
            k = 2
            while name in samples:
                name = "%s_%s_%d" % (n, stem, k)
                k += 1
            s = Sample(name, p, col, n)
            samples[name] = s
            added.append(s)
    return added


def read_ped(path, samples):
    """Apply a 6-column PED file to Sample objects matched by individual ID.

    Returns the number of samples updated.
    """
    by_name = {s.name: s for s in samples}
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 6:
                continue
            s = by_name.get(f[1])
            if s is None:
                continue
            s.father = f[2] if f[2] in by_name else "0"
            s.mother = f[3] if f[3] in by_name else "0"
            s.sex = {"1": "M", "2": "F"}.get(f[4], "?")
            s.status = {"2": "affected", "1": "unaffected"}.get(f[5], "unknown")
            n += 1
    return n


def write_ped(path, samples, family="FAM1"):
    with open(path, "w") as fp:
        for s in samples:
            fp.write("\t".join([family, s.name, s.father, s.mother,
                                {"M": "1", "F": "2"}.get(s.sex, "0"),
                                {True: "2", False: "1"}.get(s.affected, "0")]) + "\n")


def deduce_relationships(samples):
    """Fill empty father/mother/sex fields from the declared role of each sample.

    Roles are relative to the single sample labelled 'Proband'. Fields that are already
    set are never overwritten. Returns a list of warnings.
    """
    warns = []
    by_role = defaultdict(list)
    for s in samples:
        by_role[s.role].append(s)
        if s.sex == "?" and s.role in ROLE_SEX:
            s.sex = ROLE_SEX[s.role]

    def one(role):
        lst = by_role.get(role, [])
        if len(lst) > 1:
            warns.append("More than one sample labelled '%s': not used for deduction." % role)
            return None
        return lst[0] if lst else None

    proband = one("Proband")
    if proband is None:
        raise CoSegError("Exactly one sample must be labelled 'Proband'.")

    def setpar(child, fa, mo):
        if child is None:
            return
        if fa is not None and child.father == "0" and fa is not child:
            child.father = fa.name
        if mo is not None and child.mother == "0" and mo is not child:
            child.mother = mo.name

    father, mother, partner = one("Father"), one("Mother"), one("Partner")
    if partner is not None and partner.sex == "?" and proband.sex in ("M", "F"):
        partner.sex = "F" if proband.sex == "M" else "M"
    setpar(proband, father, mother)
    for s in by_role["Brother"] + by_role["Sister"]:
        setpar(s, father, mother)
    for s in by_role["Son"] + by_role["Daughter"]:
        if proband.sex == "M":
            setpar(s, proband, partner)
        elif proband.sex == "F":
            setpar(s, partner, proband)
    pgf, pgm = one("Paternal grandfather"), one("Paternal grandmother")
    mgf, mgm = one("Maternal grandfather"), one("Maternal grandmother")
    setpar(father, pgf, pgm)
    for s in by_role["Paternal uncle"] + by_role["Paternal aunt"]:
        setpar(s, pgf, pgm)
    setpar(mother, mgf, mgm)
    for s in by_role["Maternal uncle"] + by_role["Maternal aunt"]:
        setpar(s, mgf, mgm)
    return warns


# ----------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------
def to_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def parse_info(info):
    d = {}
    if info in (".", ""):
        return d
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        elif item:
            d[item] = True
    return d


def chrom_class(chrom, pos, par):
    """Classify a position as autosomal (A), X non-PAR (X), X PAR (XPAR), Y or MT."""
    c = chrom.lower()
    if c.startswith("chr"):
        c = c[3:]
    if c in ("x", "23"):
        for s, e in par:
            if s <= pos <= e:
                return "XPAR"
        return "X"
    if c in ("y", "24"):
        return "Y"
    if c in ("m", "mt", "25", "26"):
        return "MT"
    return "A"


def natural_key(c):
    c2 = c.lower()
    if c2.startswith("chr"):
        c2 = c2[3:]
    special = {"x": 23, "y": 24, "m": 25, "mt": 25}
    if c2.isdigit():
        return (0, int(c2), "")
    if c2 in special:
        return (0, special[c2], "")
    return (1, 0, c2)


# ----------------------------------------------------------------------------
# VCF reading and merging
# ----------------------------------------------------------------------------
class VariantStore:
    """In-memory union of all input files, split into biallelic records."""

    def __init__(self):
        self.vars = {}             # (chrom, pos, ref, alt) -> record dict
        self.refblocks = {}        # sample -> chrom -> (starts, blocks)
        self.meta = []             # ## lines carried to the output
        self.meta_ids = set()
        self.contigs = []
        self.contig_set = set()
        self.info_number = {}
        self.csq_fields = None
        self.chr_prefix_styles = set()
        self._have_fileformat = False

    def _meta(self, line, first):
        m = META_ID_RE.match(line)
        if m:
            typ, mid = m.group(1), m.group(2)
            if typ == "contig" and mid not in self.contig_set:
                self.contigs.append(mid)
                self.contig_set.add(mid)
            if typ == "INFO":
                n = NUMBER_RE.search(line)
                if n:
                    self.info_number.setdefault(mid, n.group(1))
                if mid == "CSQ" and self.csq_fields is None:
                    cm = CSQ_FMT_RE.search(line)
                    if cm:
                        self.csq_fields = cm.group(1).strip().split("|")
            if (typ, mid) in self.meta_ids:
                return
            self.meta_ids.add((typ, mid))
            self.meta.append(line)
            return
        if line.startswith("##fileformat"):
            if not self._have_fileformat:
                self.meta.insert(0, line)
                self._have_fileformat = True
            return
        if first:
            self.meta.append(line)

    def load(self, path, samples, log=lambda m: None):
        first = not self.meta
        n = 0
        with open_text(path) as fh:
            for line in fh:
                if line.startswith("##"):
                    self._meta(line.rstrip("\r\n"), first)
                    continue
                if line.startswith("#"):
                    continue
                f = line.rstrip("\r\n").split("\t")
                if len(f) < 10:
                    continue
                n += 1
                self._record(f, samples)
        log("  %s: %d records (%s)" % (os.path.basename(path), n,
                                       ", ".join(s.name for s in samples)))
        return n

    @staticmethod
    def _fmt(vals, i):
        if i is None or i >= len(vals):
            return None
        return to_int(vals[i])

    def _record(self, f, samples):
        chrom, pos, vid, ref, alt_s, qual, filt, info, fmt = f[:9]
        pos = int(pos)
        self.chr_prefix_styles.add(chrom.startswith("chr"))
        if chrom not in self.contig_set:
            self.contigs.append(chrom)
            self.contig_set.add(chrom)
        alts = alt_s.split(",")
        fk = fmt.split(":")
        if "GT" not in fk:
            return
        gi = fk.index("GT")
        di = fk.index("DP") if "DP" in fk else None
        mdi = fk.index("MIN_DP") if "MIN_DP" in fk else None
        qi = fk.index("GQ") if "GQ" in fk else None
        ai = fk.index("AD") if "AD" in fk else None
        fti = fk.index("FT") if "FT" in fk else None
        real = [(i + 1, a) for i, a in enumerate(alts) if a not in SKIP_ALTS]
        end = None
        if not real:
            m = END_RE.search(info)
            end = int(m.group(1)) if m else pos + len(ref) - 1

        for s in samples:
            vals = f[9 + s.col].split(":")
            gt = vals[gi] if gi < len(vals) else "."
            sep = "|" if "|" in gt else "/"
            al = GT_SPLIT_RE.split(gt)
            dp = self._fmt(vals, di)
            if dp is None:
                dp = self._fmt(vals, mdi)
            gq = self._fmt(vals, qi)
            if not real:  # gVCF reference block or reference-only site
                if al and all(a == "0" for a in al):
                    self.refblocks.setdefault(s.name, {}).setdefault(chrom, []).append(
                        (pos, end, dp, gq))
                continue
            ft = vals[fti] if fti is not None and fti < len(vals) else "."
            sfilt = filt if ft in (".", "PASS", "") else ft
            ad = vals[ai].split(",") if ai is not None and ai < len(vals) else None
            missing = gt in ("", ".") or "." in al
            for idx, alt in real:
                key = (chrom, pos, ref, alt)
                v = self.vars.get(key)
                if v is None:
                    v = {"id": vid, "qual": qual, "filter": filt, "info": info,
                         "alt_idx": idx, "n_alts": len(alts), "gts": {}}
                    self.vars[key] = v
                if missing:
                    dos, gts = None, ("./." if len(al) > 1 else ".")
                else:
                    sidx = str(idx)
                    dos = sum(1 for a in al if a == sidx)
                    gts = sep.join("1" if a == sidx else "0" for a in al)
                ad_s = "."
                if ad and len(ad) > idx and ad[0] != ".":
                    ad_s = "%s,%s" % (ad[0], ad[idx])
                v["gts"][s.name] = (dos, gts, dp, gq, ad_s, sfilt)

    def finalize(self):
        for d in self.refblocks.values():
            for chrom, lst in list(d.items()):
                lst.sort()
                d[chrom] = ([b[0] for b in lst], lst)

    def ref_cover(self, sname, chrom, pos):
        """None: sample has no gVCF blocks; False: not covered; (dp, gq): covered as 0/0."""
        d = self.refblocks.get(sname)
        if d is None:
            return None
        t = d.get(chrom)
        if not t:
            return False
        starts, blocks = t
        i = bisect.bisect_right(starts, pos) - 1
        if i >= 0 and blocks[i][1] >= pos:
            return blocks[i][2], blocks[i][3]
        return False


# ----------------------------------------------------------------------------
# Optional gene BED
# ----------------------------------------------------------------------------
class GeneBed:
    def __init__(self, path):
        tmp = defaultdict(list)
        with open_text(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                p = line.rstrip("\r\n").split("\t")
                if len(p) < 4:
                    continue
                tmp[p[0]].append((int(p[1]), int(p[2]), p[3]))
        self.data = {}
        for c, lst in tmp.items():
            lst.sort()
            maxlen = max(e - s for s, e, _ in lst)
            self.data[c] = ([s for s, _, _ in lst], lst, maxlen)

    def lookup(self, chrom, pos):
        t = self.data.get(chrom)
        if not t:
            alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            t = self.data.get(alt)
        if not t:
            return set()
        starts, lst, maxlen = t
        out = set()
        i = bisect.bisect_right(starts, pos - 1) - 1
        while i >= 0 and lst[i][0] >= pos - 1 - maxlen:
            s, e, name = lst[i]
            if s <= pos - 1 < e:
                out.add(name)
            i -= 1
        return out


# ----------------------------------------------------------------------------
# Segregation engine
# ----------------------------------------------------------------------------
class SegregationEngine:
    def __init__(self, store, samples, cfg):
        self.store = store
        self.cfg = cfg
        self.samples = samples
        self.by_name = {s.name: s for s in samples}
        self.aff = [s for s in samples if s.affected is True]
        self.unaff = [s for s in samples if s.affected is False]
        self.par = PAR_REGIONS[cfg["build"]]
        self.bed = GeneBed(cfg["gene_bed"]) if cfg.get("gene_bed") else None
        self.af_fields = [x.strip() for x in cfg.get("af_fields", "").split(",") if x.strip()]
        self.csq_idx = {}
        if store.csq_fields:
            self.csq_idx = {n: i for i, n in enumerate(store.csq_fields)}

    # --- genotype access ----------------------------------------------------
    def parents(self, s):
        return self.by_name.get(s.father), self.by_name.get(s.mother)

    def _qual_ok(self, dp, gq):
        if dp is not None and dp < self.cfg["min_dp"]:
            return False
        if gq is not None and gq < self.cfg["min_gq"]:
            return False
        return True

    def genotype(self, v, key, sname):
        """Alternate-allele dosage (0, 1, 2), or None if missing or below thresholds."""
        g = v["gts"].get(sname)
        if g is None:
            cov = self.store.ref_cover(sname, key[0], key[1])
            if cov is None:
                return 0 if self.cfg["absent_as_ref"] else None
            if cov is False:
                return None
            return 0 if self._qual_ok(*cov) else None
        dos, _, dp, gq, _, filt = g
        if dos is None:
            return None
        if self.cfg["pass_only"] and filt not in ("PASS", "."):
            return None
        if not self._qual_ok(dp, gq):
            return None
        return dos

    def req(self, d, cond):
        if d is None:
            return self.cfg["ignore_missing"]
        return cond(d)

    def any_aff_carrier(self, G):
        return any(G[s.name] is not None and G[s.name] >= 1 for s in self.aff)

    # --- annotation -----------------------------------------------------------
    def max_af(self, v, info):
        vals = []
        idx, n = v["alt_idx"], v["n_alts"]
        for fld in self.af_fields:
            val = info.get(fld)
            if isinstance(val, str):
                parts = val.split(",")
                cand = [parts[idx - 1]] if len(parts) == n and n > 1 else parts
                for p in cand:
                    try:
                        vals.append(float(p))
                    except ValueError:
                        pass
        csq = info.get("CSQ")
        if isinstance(csq, str) and self.csq_idx:
            idxs = [self.csq_idx[f] for f in self.af_fields if f in self.csq_idx]
            for entry in csq.split(","):
                p = entry.split("|")
                for i in idxs:
                    if i < len(p):
                        for x in p[i].split("&"):
                            try:
                                vals.append(float(x))
                            except ValueError:
                                pass
        return max(vals) if vals else None

    def genes_for(self, key, info):
        if self.bed:
            return self.bed.lookup(key[0], key[1])
        g = set()
        ann = info.get("ANN")
        if isinstance(ann, str):
            for e in ann.split(","):
                p = e.split("|")
                if len(p) > 3 and p[3] and "intergenic" not in p[1]:
                    g.add(p[3])
        csq = info.get("CSQ")
        if isinstance(csq, str) and self.csq_idx:
            i = self.csq_idx.get("SYMBOL", self.csq_idx.get("Gene"))
            if i is not None:
                for e in csq.split(","):
                    p = e.split("|")
                    if i < len(p) and p[i]:
                        g.add(p[i])
        for k in ("Gene.refGene", "Gene.refgene", "Gene.ensGene", "GENE", "SYMBOL",
                  "GENEINFO", "Gene_Name"):
            val = info.get(k)
            if isinstance(val, str):
                for t in re.split(r"\\x3b|[,;|&]", val):
                    t = t.split(":")[0].strip()
                    if t and t not in (".", "NONE"):
                        g.add(t)
        return g

    # --- inheritance models ---------------------------------------------------
    def m_ad(self, G):
        for s in self.aff:
            if not self.req(G[s.name], lambda d: d >= 1):
                return False
        if not self.cfg["incomplete_penetrance"]:
            for s in self.unaff:
                if not self.req(G[s.name], lambda d: d == 0):
                    return False
        return self.any_aff_carrier(G)

    def m_ar(self, G):
        for s in self.aff:
            if not self.req(G[s.name], lambda d: d == 2):
                return False
        for s in self.unaff:
            if not self.req(G[s.name], lambda d: d <= 1):
                return False
        if self.cfg["obligate_carriers"]:
            for s in self.aff:
                for p in self.parents(s):
                    if p is not None and p.affected is not True:
                        if not self.req(G[p.name], lambda d: d >= 1):
                            return False
        return any(G[s.name] == 2 for s in self.aff)

    def m_xlr(self, G):
        for s in self.aff:
            need = (lambda d: d >= 1) if s.sex == "M" else (lambda d: d == 2)
            if not self.req(G[s.name], need):
                return False
        for s in self.unaff:
            need = (lambda d: d == 0) if s.sex == "M" else (lambda d: d <= 1)
            if not self.req(G[s.name], need):
                return False
        if self.cfg["obligate_carriers"]:
            for s in self.aff:
                if s.sex == "M":
                    _, mo = self.parents(s)
                    if mo is not None and not self.req(G[mo.name], lambda d: d >= 1):
                        return False
        return self.any_aff_carrier(G)

    def m_xld(self, G):
        return self.m_ad(G)

    def m_dn(self, G):
        trios = 0
        for s in self.aff:
            if not self.req(G[s.name], lambda d: d >= 1):
                return False
            fa, mo = self.parents(s)
            if fa is not None and mo is not None:
                if G[s.name] is None:
                    return False
                if G[fa.name] != 0 or G[mo.name] != 0:  # always enforced
                    return False
                trios += 1
        if trios == 0:
            return False
        if not self.cfg["incomplete_penetrance"]:
            for s in self.unaff:
                if not self.req(G[s.name], lambda d: d == 0):
                    return False
        return True

    def m_mt(self, G):
        for s in self.aff:
            if not self.req(G[s.name], lambda d: d >= 1):
                return False
            _, mo = self.parents(s)
            if mo is not None and not self.req(G[mo.name], lambda d: d >= 1):
                return False
        return self.any_aff_carrier(G)

    def m_shared(self, G):
        return self.m_ad(G)

    def ch_candidate(self, G):
        for s in self.aff:
            if not self.req(G[s.name], lambda d: d == 1):
                return False
        for s in self.unaff:
            if not self.req(G[s.name], lambda d: d <= 1):
                return False
        return any(G[s.name] == 1 for s in self.aff)

    def ch_pair(self, G1, G2):
        """Test a pair of candidate variants in one gene.

        Returns None if the pair is incompatible, otherwise "trans" when at least one
        affected individual is phased through genotyped parents, or "unphased".
        """
        im = self.cfg["ignore_missing"]
        phased = False
        for s in self.aff:
            if G1[s.name] is None or G2[s.name] is None:
                if not im:
                    return None
                continue
            fa, mo = self.parents(s)
            if fa is not None and mo is not None:
                f1, f2, m1, m2 = G1[fa.name], G2[fa.name], G1[mo.name], G2[mo.name]
                if None in (f1, f2, m1, m2):
                    if not im:
                        return None
                    continue
                trans = ((f1 >= 1 and f2 == 0 and m1 == 0 and m2 >= 1) or
                         (f1 == 0 and f2 >= 1 and m1 >= 1 and m2 == 0))
                if not trans:
                    return None
                phased = True
            elif fa is not None or mo is not None:
                p = fa if fa is not None else mo
                d1, d2 = G1[p.name], G2[p.name]
                if d1 is None or d2 is None:
                    if not im:
                        return None
                    continue
                if (d1 >= 1) == (d2 >= 1):  # parent carries both or neither: not in trans
                    return None
                phased = True
        for s in self.unaff:
            d1, d2 = G1[s.name], G2[s.name]
            if d1 is not None and d2 is not None and d1 >= 1 and d2 >= 1:
                return None
        return "trans" if phased else "unphased"

    # --- run --------------------------------------------------------------------
    def run(self, log=lambda m: None):
        model = self.cfg["model"]
        af_thr = self.cfg.get("max_af")
        want_ch = model in ("CH", "AR_ANY")
        selected = {}
        ch_cand = defaultdict(list)
        n_af = n_nogene = n_chcand = 0
        names = [s.name for s in self.samples]

        def mark(key, tag):
            selected.setdefault(key, {"models": set(), "partners": {}, "genes": set()})
            selected[key]["models"].add(tag)

        for key, v in self.store.vars.items():
            cc = chrom_class(key[0], key[1], self.par)
            info = parse_info(v["info"])
            if af_thr is not None:
                af = self.max_af(v, info)
                if af is not None and af > af_thr:
                    n_af += 1
                    continue
            G = {n: self.genotype(v, key, n) for n in names}
            auto = cc in ("A", "XPAR")
            if model == "AD" and auto and self.m_ad(G):
                mark(key, "AD")
            elif model in ("AR", "AR_ANY") and auto and self.m_ar(G):
                mark(key, "AR_HOM")
            elif model == "XLR" and cc == "X" and self.m_xlr(G):
                mark(key, "XLR")
            elif model == "XLD" and cc == "X" and self.m_xld(G):
                mark(key, "XLD")
            elif model == "DN" and cc in ("A", "XPAR", "X") and self.m_dn(G):
                mark(key, "DE_NOVO")
            elif model == "MT" and cc == "MT" and self.m_mt(G):
                mark(key, "MT")
            elif model == "SHARED" and self.m_shared(G):
                mark(key, "SHARED")
            if want_ch and auto and self.ch_candidate(G):
                n_chcand += 1
                genes = self.genes_for(key, info)
                if not genes:
                    n_nogene += 1
                for g in genes:
                    ch_cand[g].append((key, G))

        if af_thr is not None:
            log("Variants removed by population frequency > %g: %d" % (af_thr, n_af))
        if want_ch:
            log("Compound heterozygous candidates: %d (without gene: %d)"
                % (n_chcand, n_nogene))
            if n_chcand and n_nogene == n_chcand:
                log("WARNING: no gene annotation found. Annotate the VCFs "
                    "(VEP/SnpEff/ANNOVAR) or provide a gene BED file.")
            n_pairs = 0
            for gene, lst in ch_cand.items():
                for i in range(len(lst)):
                    k1, G1 = lst[i]
                    for j in range(i + 1, len(lst)):
                        k2, G2 = lst[j]
                        phase = self.ch_pair(G1, G2)
                        if phase:
                            n_pairs += 1
                            for a, b in ((k1, k2), (k2, k1)):
                                mark(a, "COMP_HET")
                                selected[a]["partners"]["%s:%d:%s:%s" % b] = phase
                                selected[a]["genes"].add(gene)
            log("Compatible compound heterozygous pairs: %d" % n_pairs)

        for key, rec in selected.items():
            if not rec["genes"]:
                rec["genes"] = self.genes_for(key, parse_info(self.store.vars[key]["info"]))
        return selected


# ----------------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------------
def validate(samples, cfg):
    """Return (errors, warnings) for a pedigree and configuration."""
    errs, warns = [], []
    if cfg["model"] not in MODELS:
        errs.append("Unknown model '%s'." % cfg["model"])
        return errs, warns
    if cfg["build"] not in PAR_REGIONS:
        errs.append("Unknown build '%s'." % cfg["build"])
    names = {s.name: s for s in samples}
    aff = [s for s in samples if s.affected is True]
    if not aff:
        errs.append("At least one sample must be marked as affected.")
    if not any(s.affected is False for s in samples) and cfg["model"] != "DN":
        warns.append("No unaffected sample: filtering will be weakly selective.")
    for s in samples:
        for rel, exp in (("father", "M"), ("mother", "F")):
            p = getattr(s, rel)
            if p == "0":
                continue
            if p == s.name:
                errs.append("%s cannot be its own parent." % s.name)
            elif p not in names:
                errs.append("%s: parent '%s' is not among the loaded samples." % (s.name, p))
            elif names[p].sex not in (exp, "?"):
                warns.append("%s: %s is given as %s but has sex %s."
                             % (s.name, p, rel, names[p].sex))
        if s.father != "0" and s.father == s.mother:
            errs.append("%s: father and mother are the same sample." % s.name)
    if cfg["model"] in ("XLR", "XLD"):
        unk = [s.name for s in samples if s.affected is not None and s.sex == "?"]
        if unk:
            errs.append("X-linked models need the sex of: " + ", ".join(unk))
    if cfg["model"] == "DN":
        if not any(s.father in names and s.mother in names for s in aff):
            errs.append("The de novo model needs at least one affected sample with both "
                        "parents loaded.")
    if cfg["model"] in ("CH", "AR_ANY") and not cfg.get("gene_bed"):
        warns.append("Compound heterozygosity relies on gene annotation in the VCFs "
                     "(ANN/CSQ/Gene.refGene/GENEINFO); without it no pairs are found.")
    if cfg.get("bed_shared_mode", "auto") not in ("auto", "ibs1", "ibs2", "roh"):
        errs.append("Unknown shared-segment mode '%s'." % cfg["bed_shared_mode"])
    if not cfg["absent_as_ref"]:
        warns.append("Variants absent from a VCF without gVCF blocks will be treated as "
                      "missing genotypes.")
    return errs, warns


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------
def split_info(info, idx, n_alts, number_map):
    """Reduce Number=A/R INFO fields to the allele kept after splitting."""
    out = []
    if info in (".", ""):
        return out
    for item in info.split(";"):
        if not item or item.startswith("COSEG_"):
            continue
        if "=" in item:
            k, val = item.split("=", 1)
            num = number_map.get(k)
            parts = val.split(",")
            if n_alts > 1:
                if num == "A" and len(parts) == n_alts:
                    val = parts[idx - 1]
                elif num == "R" and len(parts) == n_alts + 1:
                    val = parts[0] + "," + parts[idx]
                elif num == "G":
                    continue
            out.append("%s=%s" % (k, val))
        else:
            out.append(item)
    return out


def output_paths(out_vcf):
    base = re.sub(r"\.vcf(\.gz)?$", "", out_vcf)
    return base + ".tsv", base + ".ped"


def write_outputs(store, samples, selected, cfg, log=lambda m: None):
    out_vcf = cfg["output"]
    out_tsv, out_ped = output_paths(out_vcf)
    header_order = {c: i for i, c in enumerate(store.contigs)}

    def skey(k):
        return (header_order.get(k[0], 10 ** 6), natural_key(k[0]), k[1], k[2], k[3])

    keys = sorted(selected, key=skey)

    meta = list(store.meta)
    if not meta or not meta[0].startswith("##fileformat"):
        meta.insert(0, "##fileformat=VCFv4.2")
    fmt_defs = {
        "GT": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "AD": '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        "DP": '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        "GQ": '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    }
    for fid, line in fmt_defs.items():
        if ("FORMAT", fid) not in store.meta_ids:
            meta.append(line)
    for c in store.contigs:
        if ("contig", c) not in store.meta_ids:
            meta.append("##contig=<ID=%s>" % c)
    meta += [
        '##INFO=<ID=COSEG_MODEL,Number=.,Type=String,Description="Compatible inheritance '
        'model(s): AD,AR_HOM,COMP_HET,XLR,XLD,DE_NOVO,MT,SHARED">',
        '##INFO=<ID=COSEG_GENE,Number=.,Type=String,Description="Gene(s) assigned to the variant">',
        '##INFO=<ID=COSEG_CH_PARTNER,Number=.,Type=String,Description="Partner variant(s) '
        '(CHROM:POS:REF:ALT) of a compound heterozygous pair">',
        '##INFO=<ID=COSEG_CH_PHASE,Number=.,Type=String,Description="Phase of each pair in '
        'COSEG_CH_PARTNER order: trans (inferred from parental genotypes) or unphased">',
        '##INFO=<ID=COSEG_AFF_CARRIERS,Number=1,Type=Integer,Description="Affected carriers">',
        '##INFO=<ID=COSEG_UNAFF_CARRIERS,Number=1,Type=Integer,Description="Unaffected carriers">',
    ]
    params = ("Model=%s,Build=%s,MinDP=%d,MinGQ=%d,PassOnly=%s,AbsentAsRef=%s,"
              "IgnoreMissing=%s,IncompletePenetrance=%s,ObligateCarriers=%s,MaxAF=%s"
              % (cfg["model"], cfg["build"], cfg["min_dp"], cfg["min_gq"],
                 cfg["pass_only"], cfg["absent_as_ref"], cfg["ignore_missing"],
                 cfg["incomplete_penetrance"], cfg["obligate_carriers"], cfg.get("max_af")))
    meta.append("##CoSegVCF=<Version=%s,Date=%s,%s>"
                % (__version__, datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), params))
    for s in samples:
        meta.append('##SAMPLE=<ID=%s,Role="%s",Sex=%s,Status=%s,Father=%s,Mother=%s,'
                    'Source="%s">' % (s.name, s.role, s.sex, s.status, s.father, s.mother,
                                      os.path.basename(s.path)))

    def sample_field(v, key, s):
        g = v["gts"].get(s.name)
        if g is not None:
            _, gts, dp, gq, ad, _ = g
            return gts, "%s:%s:%s:%s" % (gts, ad, "." if dp is None else dp,
                                         "." if gq is None else gq)
        cov = store.ref_cover(s.name, key[0], key[1])
        if cov is None:
            gts = "0/0" if cfg["absent_as_ref"] else "./."
            return gts, gts + ":.:.:."
        if cov is False:
            return "./.", "./.:.:.:."
        dp, gq = cov
        return "0/0", "0/0:.:%s:%s" % ("." if dp is None else dp, "." if gq is None else gq)

    opener = gzip.open if out_vcf.endswith(".gz") else open
    with opener(out_vcf, "wt") as fo, open(out_tsv, "w", newline="") as ft:
        tw = csv.writer(ft, delimiter="\t")
        tw.writerow(["CHROM", "POS", "REF", "ALT", "GENE", "MODEL", "CH_PARTNER", "CH_PHASE", "QUAL",
                     "FILTER"] + ["%s (%s)" % (s.name, s.status) for s in samples])
        fo.write("\n".join(meta) + "\n")
        fo.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                            "FORMAT"] + [s.name for s in samples]) + "\n")
        for key in keys:
            v = store.vars[key]
            rec = selected[key]
            fields = [sample_field(v, key, s) for s in samples]
            n_aff = sum(1 for s, (g, _) in zip(samples, fields)
                        if s.affected is True and "1" in g)
            n_un = sum(1 for s, (g, _) in zip(samples, fields)
                       if s.affected is False and "1" in g)
            info = split_info(v["info"], v["alt_idx"], v["n_alts"], store.info_number)
            info.append("COSEG_MODEL=" + ",".join(sorted(rec["models"])))
            if rec["genes"]:
                info.append("COSEG_GENE=" + ",".join(sorted(rec["genes"])))
            partners = sorted(rec["partners"])
            if partners:
                info.append("COSEG_CH_PARTNER=" + ",".join(partners))
                info.append("COSEG_CH_PHASE=" + ",".join(rec["partners"][p] for p in partners))
            info.append("COSEG_AFF_CARRIERS=%d" % n_aff)
            info.append("COSEG_UNAFF_CARRIERS=%d" % n_un)
            fo.write("\t".join([key[0], str(key[1]), v["id"], key[2], key[3], v["qual"],
                                v["filter"], ";".join(info), "GT:AD:DP:GQ"]
                               + [x for _, x in fields]) + "\n")
            tw.writerow([key[0], key[1], key[2], key[3], ",".join(sorted(rec["genes"])),
                         ",".join(sorted(rec["models"])), ",".join(partners),
                         ",".join(rec["partners"][p] for p in partners), v["qual"], v["filter"]]
                        + ["%s (DP=%s)" % (x.split(":")[0], x.split(":")[2]) for _, x in fields])

    write_ped(out_ped, samples)
    log("Written: %s\n         %s\n         %s" % (out_vcf, out_tsv, out_ped))
    return len(keys)


def run_analysis(samples, cfg, log=print):
    """Load all files, apply the selected model and write VCF, TSV and PED outputs.

    Returns the number of variants written.
    """
    full = dict(DEFAULT_CONFIG)
    full.update(cfg)
    cfg = full
    errs, _ = validate(samples, cfg)
    if errs:
        raise CoSegError("\n".join(errs))
    t0 = datetime.datetime.now()
    log("=== CoSegVCF %s - model: %s ===" % (__version__, MODELS[cfg["model"]]))
    store = VariantStore()
    by_file = OrderedDict()
    for s in samples:
        by_file.setdefault(s.path, []).append(s)
    log("Reading %d file(s), %d sample(s)..." % (len(by_file), len(samples)))
    for path, ss in by_file.items():
        store.load(path, ss, log)
    store.finalize()
    log("Unique biallelic variants in the union: %d" % len(store.vars))
    if len(store.chr_prefix_styles) > 1:
        log("WARNING: mixed chromosome naming ('chr1' and '1'); variants will not be "
            "matched across these files.")
    for s in samples:
        if s.name in store.refblocks:
            log("  %s: gVCF reference blocks found (coverage checked)" % s.name)
    engine = SegregationEngine(store, samples, cfg)
    selected = engine.run(log)
    counts = defaultdict(int)
    for rec in selected.values():
        for m in rec["models"]:
            counts[m] += 1
    for m, c in sorted(counts.items()):
        log("  %s: %d variants" % (m, c))
    n = write_outputs(store, samples, selected, cfg, log)
    if cfg.get("bed_shared") or cfg.get("bed_callable"):
        from .regions import write_region_beds
        write_region_beds(engine, cfg, re.sub(r"\.vcf(\.gz)?$", "", cfg["output"]), log)
    log("Compatible variants written: %d (elapsed %s)"
        % (n, str(datetime.datetime.now() - t0).split(".")[0]))
    return n
