import gzip

import pytest

from cosegvcf.cli import main as cli_main
from cosegvcf.core import (CoSegError, Sample, SegregationEngine, VariantStore, deduce_relationships,
                           load_samples, read_ped, run_analysis, split_info, validate)


def records(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]


def positions(path):
    return [(r[0], int(r[1]), r[4]) for r in records(path)]


def run(tmp, samples, model, **kw):
    cfg = {"model": model, "max_af": 0.01, "output": str(tmp / ("out_%s.vcf" % model))}
    cfg.update(kw)
    run_analysis(samples, cfg, log=lambda m: None)
    return positions(cfg["output"])


EXPECTED = {
    "AD": [("chr1", 400, "G")],
    "AR": [("chr1", 100, "T")],
    "CH": [("chr1", 200, "A"), ("chr1", 300, "C")],
    "AR_ANY": [("chr1", 100, "T"), ("chr1", 200, "A"), ("chr1", 300, "C")],
    "XLR": [("chrX", 5000000, "T")],
    "XLD": [],
    "DN": [("chr1", 400, "G")],
    "MT": [("chrM", 150, "T")],
    "SHARED": [("chr1", 400, "G")],
}


@pytest.mark.parametrize("model", sorted(EXPECTED))
def test_models(family, model):
    tmp, samples, _, _ = family
    assert run(tmp, samples, model) == EXPECTED[model]


def test_compound_het_partners_and_gene(family):
    tmp, samples, _, _ = family
    run(tmp, samples, "CH")
    info = {int(r[1]): r[7] for r in records(str(tmp / "out_CH.vcf"))}
    assert "COSEG_CH_PARTNER=chr1:300:G:C" in info[200]
    assert "COSEG_GENE=GENEB" in info[300]
    assert "COSEG_CH_PHASE=trans" in info[300]


def test_compound_het_phasing_rules(family):
    _, samples, _, _ = family
    cfg = {"model": "CH", "build": "GRCh38", "ignore_missing": False, "gene_bed": ""}
    eng = SegregationEngine(VariantStore(), samples, cfg)
    g = lambda fa, mo, pro, sis: {"FATHER": fa, "MOTHER": mo, "PROBAND": pro, "SISTER": sis}
    # in trans: one variant from each parent
    assert eng.ch_pair(g(1, 0, 1, 0), g(0, 1, 1, 0)) == "trans"
    # in cis: both variants from the mother
    assert eng.ch_pair(g(0, 1, 1, 0), g(0, 1, 1, 0)) is None
    # unaffected sister carries both
    assert eng.ch_pair(g(1, 0, 1, 1), g(0, 1, 1, 1)) is None
    # a missing parental genotype rejects the pair unless missing data are ignored
    assert eng.ch_pair(g(None, 0, 1, 0), g(0, 1, 1, 0)) is None
    eng.cfg["ignore_missing"] = True
    assert eng.ch_pair(g(None, 0, 1, 0), g(0, 1, 1, 0)) == "unphased"


def test_compound_het_unphased_without_parents(family):
    tmp, samples, _, _ = family
    for s in samples:
        s.father = s.mother = "0"
    run(tmp, samples, "CH")
    info = {int(r[1]): r[7] for r in records(str(tmp / "out_CH.vcf"))}
    assert "COSEG_CH_PHASE=unphased" in info[200]


def test_multiallelic_split_and_af(family):
    tmp, samples, _, _ = family
    # without AF filter chr1:500 T is inherited from the unaffected father -> not AD
    out = run(tmp, samples, "AD", max_af=None)
    assert ("chr1", 500, "T") not in out
    store = VariantStore()
    store.info_number = {"AC": "A", "gnomAD_AF": "A"}
    kept = split_info("AC=3,1;gnomAD_AF=0.3,0.0001", 2, 2, store.info_number)
    assert kept == ["AC=1", "gnomAD_AF=0.0001"]


def test_incomplete_penetrance(family):
    tmp, samples, _, _ = family
    out = run(tmp, samples, "AD", incomplete_penetrance=True)
    assert ("chr1", 100, "T") in out and ("chr1", 400, "G") in out


def test_gvcf_blocks_control_absent_calls(family, tmp_path):
    tmp, samples, _, _ = family
    gvcf = tmp_path / "FATHER.g.vcf"
    gvcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFATHER\n"
        "chr1\t1\t.\tN\t<NON_REF>\t.\t.\tEND=350\tGT:DP:GQ:MIN_DP\t0/0:30:99:25\n"
        "chr1\t351\t.\tN\t<NON_REF>\t.\t.\tEND=450\tGT:DP:GQ:MIN_DP\t0/0:3:5:2\n")
    for s in samples:
        if s.name == "FATHER":
            s.path = str(gvcf)
    out = run(tmp, samples, "DN")
    # chr1:400 falls in a low-quality reference block of the father -> not called de novo
    assert ("chr1", 400, "G") not in out
    # chr1:200 is now covered as confident 0/0 in the father -> de novo candidate
    assert ("chr1", 200, "A") in out


def test_low_quality_genotype_is_missing(family):
    tmp, samples, _, _ = family
    assert run(tmp, samples, "DN", min_dp=40) == []
    assert run(tmp, samples, "AD", min_dp=40, ignore_missing=True) == []


def test_validation_errors(family):
    _, samples, _, _ = family
    for s in samples:
        s.status = "unaffected"
    errs, _ = validate(samples, {"model": "AD", "build": "GRCh38", "absent_as_ref": True})
    assert any("affected" in e for e in errs)
    with pytest.raises(CoSegError):
        run_analysis(samples, {"model": "AD"}, log=lambda m: None)


def test_dn_needs_trio(family):
    _, samples, _, _ = family
    for s in samples:
        s.father = s.mother = "0"
    errs, _ = validate(samples, {"model": "DN", "build": "GRCh38", "absent_as_ref": True})
    assert errs


def test_ped_and_duplicate_names(family, tmp_path):
    _, _, paths, ped = family
    samples = load_samples([paths["FATHER"], paths["FATHER"]])
    assert [s.name for s in samples] == ["FATHER", "FATHER_FATHER"]
    samples = load_samples(list(paths.values()))
    assert read_ped(ped, samples) == 4
    pro = [s for s in samples if s.name == "PROBAND"][0]
    assert (pro.father, pro.mother, pro.sex, pro.affected) == ("FATHER", "MOTHER", "M", True)


def test_deduce_relationships():
    names = {"P": "Proband", "F": "Father", "M": "Mother", "S": "Sister", "G": "Maternal grandmother"}
    samples = []
    for n, role in names.items():
        s = Sample(n, "x.vcf", 0)
        s.role = role
        samples.append(s)
    deduce_relationships(samples)
    by = {s.name: s for s in samples}
    assert (by["P"].father, by["P"].mother) == ("F", "M")
    assert (by["S"].father, by["S"].mother, by["S"].sex) == ("F", "M", "F")
    assert by["M"].mother == "G"


def test_cli_gz_output(family):
    tmp, _, paths, ped = family
    out = str(tmp / "cli.vcf.gz")
    rc = cli_main(["--vcf"] + list(paths.values()) + ["--ped", ped, "--model", "AR",
                                                       "-o", out, "-q"])
    assert rc == 0
    assert positions(out) == EXPECTED["AR"]
    assert (tmp / "cli.tsv").exists() and (tmp / "cli.ped").exists()


def test_cli_reports_errors(family, capsys):
    tmp, _, paths, ped = family
    rc = cli_main(["--vcf", paths["FATHER"], "--ped", ped, "--model", "AD",
                   "-o", str(tmp / "x.vcf"), "-q"])
    assert rc == 1
    assert "affected" in capsys.readouterr().err
