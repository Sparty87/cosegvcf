# CoSegVCF

**Pedigree-aware inheritance-model filtering of VCF and gVCF files, with a graphical interface.**

CoSegVCF reads the VCF or gVCF files of the members of a family, lets you describe who is
related to whom and who is affected, and writes a multi-sample VCF that keeps only the
variants compatible with the inheritance model you suspect: autosomal dominant,
autosomal recessive (homozygous and/or compound heterozygous), X-linked recessive or
dominant, *de novo*, mitochondrial, or simply "shared by the affected and absent in the
unaffected".

It is meant for diagnostic laboratories that hold **one VCF per individual** (often from
different runs or pipelines) and want the person who knows the family to run the
segregation analysis without writing filter expressions.

![CI](https://github.com/OWNER/cosegvcf/actions/workflows/tests.yml/badge.svg)

## Features

- Any number of plain or gzipped VCF/gVCF files, single- or multi-sample, mixed freely.
- gVCF reference blocks (`<NON_REF>`, `END=`) are used to tell a confident 0/0 from a
  position that was not covered; with variant-only VCFs you choose how absence is read.
- Pedigree table in the GUI: relationship to the proband, sex, father, mother, status
  (affected / unaffected / unknown). Parents and sex can be filled in from the declared
  relationships; PED import and export.
- Nine models, selectable from a drop-down menu (see below).
- Compound heterozygous pairs grouped by gene (SnpEff `ANN`, VEP `CSQ`, ANNOVAR
  `Gene.refGene`, `GENEINFO`, or a BED of genes) and phased through parental genotypes.
- Genotype quality thresholds (DP, GQ), FILTER/FT = PASS, obligate-carrier checks,
  incomplete penetrance, optional population allele-frequency cut-off on any INFO or CSQ
  field.
- Optional BED export of the genomic regions shared by the affected individuals:
  identical-by-state segments (IBS1, IBS2 or runs of shared homozygosity) and, with gVCF
  input, the regions confidently genotyped in all affected individuals.
- Multi-allelic records split into biallelic ones; Number=A/R INFO fields reduced to the
  kept allele.
- Output: sorted multi-sample VCF (`.vcf` or `.vcf.gz`) with added `COSEG_*` INFO tags,
  the full parameter set and pedigree in the header, plus a TSV summary and a PED file.
- Pure Python standard library: no dependencies to install.

## Installation

Python 3.8 or later is required. The GUI uses Tkinter, which ships with the Python
installers for Windows and macOS; on Linux install it with your package manager, e.g.
`sudo apt install python3-tk`.

```bash
git clone https://github.com/OWNER/cosegvcf.git
cd cosegvcf
pip install .
```

You can also run it without installing: `python run_cosegvcf.py` (or `python -m cosegvcf`)
from the repository folder.

## Usage

### Graphical interface

```bash
cosegvcf            # or: cosegvcf-gui, or: python -m cosegvcf
```

1. **Add VCF...** — select the files of all family members.
2. Fill the pedigree table (relationship, sex, father, mother, status), or use
   **Fill parents and sex from relationships**, or **Import PED...**.
3. Choose the suspected model from the drop-down menu; a short description of the rules
   applied appears underneath. Adjust filters if needed.
4. Choose the output file and press **Run**. Progress and counts appear in the log.

### Command line

```bash
cosegvcf --vcf father.vcf.gz mother.vcf.gz proband.vcf.gz sister.vcf.gz \
         --ped family.ped --model AR_ANY --max-af 0.01 -o candidates.vcf.gz
```

PED individual IDs must match the VCF sample names. Run `cosegvcf --help` for all
options. A worked example with a four-member family is in [`examples/`](examples/).

## Inheritance models

| Code | Model | Affected | Unaffected | Extra checks |
| --- | --- | --- | --- | --- |
| `AD` | Autosomal dominant | carrier | 0/0 (carrier allowed with `--incomplete-penetrance`) | — |
| `AR` | Recessive, homozygous | 1/1 | not 1/1 | unaffected parents of affected are carriers |
| `CH` | Recessive, compound heterozygous | 0/1 at both variants of a gene | not carrying both | variants from different parents when parents are present |
| `AR_ANY` | Recessive, combined | `AR` or `CH` | | |
| `XLR` | X-linked recessive (non-PAR) | male carrier, female 1/1 | male 0, female ≤ 0/1 | mothers of affected males are carriers |
| `XLD` | X-linked dominant (non-PAR) | carrier | reference | — |
| `DN` | *De novo* | carrier | reference | both parents called 0/0 and passing thresholds |
| `MT` | Mitochondrial | carrier | not constrained | mothers of affected are carriers |
| `SHARED` | No model | carrier | reference | — |

Samples with status *unknown* are written to the output but do not constrain filtering.
Missing or low-quality genotypes reject a variant unless `--ignore-missing` is set; for
`DN` the parental 0/0 calls are always required.

## Output

Added INFO fields:

| Field | Meaning |
| --- | --- |
| `COSEG_MODEL` | compatible model(s): `AD`, `AR_HOM`, `COMP_HET`, `XLR`, `XLD`, `DE_NOVO`, `MT`, `SHARED` |
| `COSEG_GENE` | gene(s) assigned to the variant |
| `COSEG_CH_PARTNER` | partner variant(s) of a compound heterozygous pair |
| `COSEG_CH_PHASE` | `trans` (phased through parents) or `unphased`, one per partner |
| `COSEG_AFF_CARRIERS` / `COSEG_UNAFF_CARRIERS` | number of carriers by status |

`FORMAT` is `GT:AD:DP:GQ`. A position confirmed by a gVCF reference block is written as
`0/0` with that block's depth and quality; a position without evidence is written as `./.`.

## Region BED export

Two optional BED files describe *where* the affected individuals share their genome,
independently of the variant-level filter. Tick the boxes under **Region BED export** in
the GUI, or use `--shared-bed` / `--callable-bed` on the command line.

**`<output>.shared_<mode>.bed` — segments shared by the affected.** Runs of consecutive
informative sites at which the genotypes of all affected individuals are compatible with
a shared haplotype:

| Mode | A site is concordant when | Typical use |
| --- | --- | --- |
| `ibs1` | no two affected are homozygous for opposite alleles | dominant models, distant relatives |
| `ibs2` | all affected have the same genotype | affected siblings, recessive models |
| `roh` | all affected are homozygous for the same allele | autozygosity mapping, consanguinity; works with one affected |
| `auto` | `ibs1` for AD/XLD/DN/MT/SHARED; `ibs2` for AR/CH/XLR (`roh` if one affected) | default |

Sites with a missing or low-quality genotype in any affected individual are skipped. One
isolated discordant site is tolerated as a probable genotyping error and counted; two
discordant sites within 10 informative sites end the segment. Segments shorter than
`--shared-min-kb` (default 100 kb) or with fewer than `--shared-min-sites` (default 20)
concordant sites are discarded. Columns: chrom, start, end, mode, concordant sites,
discordant sites, length. Hemizygous chrX calls in males are treated as homozygous; chrY
and chrM are not analysed.

**`<output>.callable_affected.bed` — regions callable in all affected.** Intervals in
which every affected individual has a confident genotype (a gVCF reference block or a
variant call passing the DP/GQ thresholds). It requires gVCF input for all affected
samples and tells you where the absence of a shared variant is informative.

Notes: sharing is by state, not by descent; with exome data, sites are sparse and segment
boundaries are approximate, so loosen or tighten the thresholds to the design. gVCF input
is recommended, because reference blocks let homozygous-reference genotypes count as
evidence; with variant-only VCFs the result depends on how absent calls are read.

```bash
cosegvcf --vcf sib1.g.vcf.gz sib2.g.vcf.gz parents.vcf.gz --ped fam.ped --model AR_ANY \
         --shared-bed --callable-bed -o candidates.vcf.gz
# then, for example, keep candidates inside shared segments:
bedtools intersect -header -a candidates.vcf.gz -b candidates.shared_ibs2.bed
```

## Good practice

- Normalize and left-align all inputs against the same reference
  (`bcftools norm -f ref.fa -m -any`) so that the same variant is written the same way.
- Use the same chromosome naming (`chr1` or `1`) in every file; the program warns if it
  finds both.
- Prefer gVCFs (e.g. GATK HaplotypeCaller `-ERC GVCF`) or a jointly genotyped VCF: with
  variant-only VCFs a missing variant cannot be told apart from missing coverage.
- Candidate lists depend on data quality: check the depth and quality of every reported
  genotype, and confirm the final candidates by an independent method.

CoSegVCF is a research and laboratory aid. It does not classify variants and does not
replace clinical interpretation.

## Development

```bash
pip install -e ".[test]"
pytest
```

The simulation benchmark used in the paper is in [`benchmark/`](benchmark/).

## Citation

If you use CoSegVCF, please cite it as described in [`CITATION.cff`](CITATION.cff).

## License

MIT — see [`LICENSE`](LICENSE).
