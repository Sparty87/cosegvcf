"""Command-line interface: ``cosegvcf --vcf a.vcf b.vcf --ped family.ped --model AR -o out.vcf``.

Without arguments the graphical interface is started.
"""

import argparse
import sys

from . import __version__
from .core import (DEFAULT_AF_FIELDS, MODELS, PAR_REGIONS, CoSegError, load_samples,
                   read_ped, run_analysis, validate)


def build_parser():
    p = argparse.ArgumentParser(
        prog="cosegvcf",
        description="Filter variants from one or more VCF/gVCF files according to a "
                    "suspected inheritance model and a pedigree. Run without arguments "
                    "to open the graphical interface.",
        epilog="Models: " + "; ".join("%s = %s" % kv for kv in MODELS.items()))
    p.add_argument("--vcf", nargs="+", metavar="FILE", required=True,
                   help="input VCF/gVCF files (plain or .gz), single- or multi-sample")
    p.add_argument("--ped", required=True, metavar="FILE",
                   help="6-column PED file; individual IDs must match VCF sample names")
    p.add_argument("--model", required=True, choices=list(MODELS))
    p.add_argument("-o", "--output", required=True, metavar="VCF",
                   help="output VCF (.vcf or .vcf.gz); a .tsv and a .ped are written alongside")
    p.add_argument("--build", default="GRCh38", choices=list(PAR_REGIONS),
                   help="genome build used to locate the X pseudoautosomal regions")
    p.add_argument("--min-dp", type=int, default=8, help="minimum genotype depth (default 8)")
    p.add_argument("--min-gq", type=int, default=20,
                   help="minimum genotype quality (default 20)")
    p.add_argument("--keep-filtered", action="store_true",
                   help="also use calls whose FILTER/FT is not PASS")
    p.add_argument("--absent-as-missing", action="store_true",
                   help="treat variants absent from a VCF without gVCF blocks as missing "
                        "instead of 0/0")
    p.add_argument("--ignore-missing", action="store_true",
                   help="missing or low-quality genotypes do not reject a variant")
    p.add_argument("--incomplete-penetrance", action="store_true",
                   help="unaffected individuals may carry dominant/de novo variants")
    p.add_argument("--no-obligate-carriers", action="store_true",
                   help="do not require carrier status in parents (AR) or mothers (XLR)")
    p.add_argument("--max-af", type=float, default=None,
                   help="discard variants with population frequency above this value")
    p.add_argument("--af-fields", default=DEFAULT_AF_FIELDS,
                   help="comma-separated INFO/CSQ fields holding population frequencies")
    p.add_argument("--gene-bed", default="",
                   help="BED file (chrom, start, end, gene) used for compound heterozygosity")
    p.add_argument("--shared-bed", action="store_true",
                   help="also write a BED of segments shared (identical by state) by the "
                        "affected individuals")
    p.add_argument("--shared-mode", default="auto", choices=["auto", "ibs1", "ibs2", "roh"],
                   help="ibs1: at least one allele shared; ibs2: identical genotypes; roh: "
                        "shared homozygosity; auto: ibs1 for dominant-type models, ibs2 "
                        "(or roh with one affected) for recessive ones")
    p.add_argument("--shared-min-sites", type=int, default=20,
                   help="minimum concordant sites per shared segment (default 20)")
    p.add_argument("--shared-min-kb", type=float, default=100,
                   help="minimum length of a shared segment in kb (default 100)")
    p.add_argument("--callable-bed", action="store_true",
                   help="also write a BED of regions confidently genotyped in all affected "
                        "individuals (requires gVCF input for every affected sample)")
    p.add_argument("-q", "--quiet", action="store_true", help="print errors only")
    p.add_argument("--version", action="version", version="%(prog)s " + __version__)
    return p


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        from .gui import main as gui_main
        return gui_main()
    args = build_parser().parse_args(argv)
    log = (lambda m: None) if args.quiet else (lambda m: print(m, file=sys.stderr))
    try:
        samples = load_samples(args.vcf)
        n = read_ped(args.ped, samples)
        if n == 0:
            raise CoSegError("No PED individual matches a VCF sample. VCF samples: "
                             + ", ".join(s.name for s in samples))
        cfg = {
            "model": args.model, "build": args.build, "min_dp": args.min_dp,
            "min_gq": args.min_gq, "pass_only": not args.keep_filtered,
            "absent_as_ref": not args.absent_as_missing,
            "ignore_missing": args.ignore_missing,
            "incomplete_penetrance": args.incomplete_penetrance,
            "obligate_carriers": not args.no_obligate_carriers, "max_af": args.max_af,
            "af_fields": args.af_fields, "gene_bed": args.gene_bed, "output": args.output,
            "bed_shared": args.shared_bed, "bed_shared_mode": args.shared_mode,
            "bed_min_sites": args.shared_min_sites, "bed_min_kb": args.shared_min_kb,
            "bed_callable": args.callable_bed,
        }
        errs, warns = validate(samples, cfg)
        for w in warns:
            log("WARNING: " + w)
        if errs:
            raise CoSegError("\n".join(errs))
        run_analysis(samples, cfg, log=log)
    except (CoSegError, OSError) as e:
        print("cosegvcf: error: %s" % e, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
