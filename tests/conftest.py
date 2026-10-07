"""Shared fixtures: a four-member family (father, mother, affected son, unaffected daughter)
with one variant per inheritance pattern, written as variant-only single-sample VCFs."""

import pytest

from cosegvcf.core import Sample

HEADER = """##fileformat=VCFv4.2
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotation">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##contig=<ID=chr1>
##contig=<ID=chrX>
##contig=<ID=chrM>
"""


def ann(gene):
    return "ANN=A|missense_variant|MODERATE|%s|x" % gene


# (chrom, pos, ref, alt, info, genotypes) - expected model in the comment
VARIANTS = [
    ("chr1", 100, "C", "T", ann("GENEA") + ";gnomAD_AF=0.001",
     {"FATHER": "0/1", "MOTHER": "0/1", "PROBAND": "1/1", "SISTER": "0/1"}),   # AR hom
    ("chr1", 200, "G", "A", ann("GENEB"),
     {"FATHER": "0/1", "MOTHER": "0/0", "PROBAND": "0/1", "SISTER": "0/0"}),   # CH, paternal
    ("chr1", 300, "G", "C", ann("GENEB"),
     {"FATHER": "0/0", "MOTHER": "0/1", "PROBAND": "0/1", "SISTER": "0/1"}),   # CH, maternal
    ("chr1", 400, "A", "G", ann("GENEC"),
     {"FATHER": "0/0", "MOTHER": "0/0", "PROBAND": "0/1", "SISTER": "0/0"}),   # de novo / AD
    ("chr1", 500, "A", "G,T", ann("GENED") + ";AC=3,1;gnomAD_AF=0.3,0.0001",
     {"FATHER": "1/2", "MOTHER": "0/0", "PROBAND": "0/2", "SISTER": "0/0"}),   # T inherited
    ("chrX", 5000000, "C", "T", ann("GENEX"),
     {"FATHER": "0", "MOTHER": "0/1", "PROBAND": "1", "SISTER": "0/0"}),       # XLR
    ("chrM", 150, "C", "T", ".",
     {"FATHER": "0", "MOTHER": "1", "PROBAND": "1", "SISTER": "1"}),           # MT
]

PEDIGREE = [  # name, sex, status, father, mother
    ("FATHER", "M", "unaffected", "0", "0"),
    ("MOTHER", "F", "unaffected", "0", "0"),
    ("PROBAND", "M", "affected", "FATHER", "MOTHER"),
    ("SISTER", "F", "unaffected", "FATHER", "MOTHER"),
]


def write_family(directory):
    paths = {}
    for name, *_ in PEDIGREE:
        p = directory / ("%s.vcf" % name)
        with open(p, "w") as fh:
            fh.write(HEADER)
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % name)
            for chrom, pos, ref, alt, info, gts in VARIANTS:
                gt = gts[name]
                if set(gt.replace("/", "")) == {"0"}:
                    continue  # variant-only VCF: reference calls are not written
                fh.write("%s\t%d\t.\t%s\t%s\t50\tPASS\t%s\tGT:DP:GQ\t%s:30:60\n"
                         % (chrom, pos, ref, alt, info, gt))
        paths[name] = str(p)
    ped = directory / "family.ped"
    with open(ped, "w") as fh:
        for name, sex, status, fa, mo in PEDIGREE:
            fh.write("FAM1\t%s\t%s\t%s\t%s\t%s\n" % (
                name, fa, mo, {"M": "1", "F": "2"}[sex], {"affected": "2", "unaffected": "1"}[status]))
    return paths, str(ped)


@pytest.fixture
def family(tmp_path):
    paths, ped = write_family(tmp_path)
    samples = []
    for name, sex, status, fa, mo in PEDIGREE:
        s = Sample(name, paths[name], 0)
        s.sex, s.status, s.father, s.mother = sex, status, fa, mo
        samples.append(s)
    return tmp_path, samples, paths, ped
