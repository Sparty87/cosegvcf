from cosegvcf.cli import main as cli_main
from cosegvcf.core import Sample, run_analysis
from cosegvcf.regions import _intersect, _merge, resolve_mode, segments, site_state

GVCF_HEADER = ("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##contig=<ID=chr2>\n"
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n")


def gvcf(path, name, records):
    """records: (chrom, pos, end_or_None, alt, gt, dp, gq); end set -> reference block."""
    with open(path, "w") as fh:
        fh.write(GVCF_HEADER % name)
        for chrom, pos, end, alt, gt, dp, gq in records:
            if end is not None:
                fh.write("%s\t%d\t.\tA\t<NON_REF>\t.\t.\tEND=%d\tGT:DP:GQ\t0/0:%d:%d\n"
                         % (chrom, pos, end, dp, gq))
            else:
                fh.write("%s\t%d\t.\tA\t%s,<NON_REF>\t50\tPASS\t.\tGT:DP:GQ\t%s:%d:%d\n"
                         % (chrom, pos, alt, gt, dp, gq))
    return str(path)


def bed(path):
    with open(path) as fh:
        return [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]


def test_site_state():
    assert site_state([1, 2], "ibs1") is True
    assert site_state([0, 2], "ibs1") is False
    assert site_state([1, 1], "ibs2") is True
    assert site_state([1, 2], "ibs2") is False
    assert site_state([2, 2], "roh") is True
    assert site_state([0], "roh") is True
    assert site_state([1], "roh") is False
    assert site_state([1, None], "ibs1") is None


def test_resolve_mode():
    assert resolve_mode("auto", "AD", 2) == "ibs1"
    assert resolve_mode("auto", "AR", 2) == "ibs2"
    assert resolve_mode("auto", "AR", 1) == "roh"
    assert resolve_mode("auto", "AD", 1) is None
    assert resolve_mode("roh", "AD", 1) == "roh"


def test_segments_break_and_tolerance():
    ok = lambda a, b: [(p * 1000, True) for p in range(a, b)]
    # one isolated discordant site inside a long run is tolerated and counted
    sites = ok(1, 31) + [(31000, False)] + ok(32, 61)
    assert segments(sites, min_sites=5, min_bp=1000) == [(1000, 60000, 59, 1)]
    # with tolerance off the run is split in two
    assert segments(sites, 5, 1000, tolerate_isolated=False) == [
        (1000, 30000, 30, 0), (32000, 60000, 29, 0)]
    # two close discordant sites end the segment; the first one is its border
    sites = ok(1, 31) + [(31000, False), (32000, False)] + ok(33, 50)
    assert segments(sites, 5, 1000) == [(1000, 30000, 30, 0), (33000, 49000, 17, 0)]
    # thresholds on number of sites and on length
    assert segments(ok(1, 5), min_sites=5, min_bp=0) == []
    assert segments(ok(1, 11), min_sites=5, min_bp=50000) == []


def test_interval_helpers():
    assert _merge([(5, 9), (1, 3), (4, 4), (20, 30)]) == [[1, 9], [20, 30]]
    assert _intersect([[1, 10], [20, 30]], [[5, 25]]) == [[5, 10], [20, 25]]


def two_sibs(tmp_path):
    # chr1: 30 sites where both sibs are 0/1 (shared), then both discordant (0/0 vs 1/1)
    a, b = [("chr1", 1, 50000, None, None, 30, 99)], [("chr1", 1, 50000, None, None, 30, 99)]
    for i in range(30):
        pos = 60000 + i * 10000
        a.append(("chr1", pos, None, "G", "0/1", 30, 99))
        b.append(("chr1", pos, None, "G", "0/1", 30, 99))
    for i in range(30):
        pos = 400000 + i * 10000
        a.append(("chr1", pos, None, "G", "1/1", 30, 99))
        b.append(("chr1", pos, pos, None, None, 30, 99))       # confident 0/0 block
    # chr2: sib A covered 1-1000 well, sib B only 1-400 well and 401-1000 at low depth
    a.append(("chr2", 1, 1000, None, None, 30, 99))
    b.append(("chr2", 1, 400, None, None, 30, 99))
    b.append(("chr2", 401, 1000, None, None, 3, 5))
    pa = gvcf(tmp_path / "A.g.vcf", "A", a)
    pb = gvcf(tmp_path / "B.g.vcf", "B", b)
    samples = []
    for n, p in (("A", pa), ("B", pb)):
        s = Sample(n, p, 0)
        s.sex, s.status = "F", "affected"
        samples.append(s)
    return samples


def test_shared_and_callable_beds(tmp_path):
    samples = two_sibs(tmp_path)
    out = str(tmp_path / "out.vcf")
    run_analysis(samples, {"model": "AD", "output": out, "bed_shared": True,
                           "bed_callable": True, "bed_min_sites": 10, "bed_min_kb": 50},
                 log=lambda m: None)
    shared = bed(str(tmp_path / "out.shared_ibs1.bed"))
    assert shared == [["chr1", "59999", "350000", "IBS1", "30", "0", "290001"]]
    call = bed(str(tmp_path / "out.callable_affected.bed"))
    assert ["chr1", "0", "50000"] in call
    assert ["chr2", "0", "400"] in call                 # low-depth block of B excluded
    assert not any(c[0] == "chr2" and int(c[2]) > 400 for c in call)
    assert ["chr1", "399999", "400000"] in call         # variant in A, confident 0/0 in B


def test_roh_mode_and_cli(tmp_path):
    samples = two_sibs(tmp_path)
    out = str(tmp_path / "cli.vcf")
    ped = tmp_path / "f.ped"
    ped.write_text("F\tA\t0\t0\t2\t2\nF\tB\t0\t0\t2\t1\n")
    rc = cli_main(["--vcf", samples[0].path, samples[1].path, "--ped", str(ped),
                   "--model", "AR", "-o", out, "-q", "--shared-bed",
                   "--shared-min-sites", "10", "--shared-min-kb", "50", "--callable-bed"])
    assert rc == 0
    roh = bed(str(tmp_path / "cli.shared_roh.bed"))      # single affected -> auto = roh
    assert roh == [["chr1", "399999", "690000", "ROH", "30", "0", "290001"]]


def test_callable_needs_gvcf(family):
    tmp, samples, _, _ = family
    logs = []
    run_analysis(samples, {"model": "AD", "output": str(tmp / "o.vcf"), "bed_callable": True,
                           "bed_shared": True}, log=logs.append)
    assert any("callable-region BED skipped" in m for m in logs)
    assert any("shared-segment BED skipped" in m for m in logs)   # one affected, AD
    assert not (tmp / "o.callable_affected.bed").exists()
