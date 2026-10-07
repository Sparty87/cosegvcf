# Changelog

## 1.1.0 — 2026-10-05

- New optional BED exports (`--shared-bed`, `--callable-bed`, and the corresponding GUI
  options): identical-by-state segments shared by the affected individuals (IBS1, IBS2,
  shared homozygosity) and, with gVCF input, regions confidently genotyped in all affected.

## 1.0.0 — 2026-09-30

First public release.

- Graphical interface (Tkinter) and command-line interface sharing one engine.
- Nine inheritance models: AD, AR homozygous, compound heterozygous, combined recessive,
  XLR, XLD, de novo, mitochondrial, shared.
- gVCF reference blocks used to distinguish confident 0/0 calls from missing coverage.
- Compound heterozygous pairs phased through parental genotypes (`COSEG_CH_PHASE`).
- Multi-allelic splitting with Number=A/R INFO reduction.
- VCF, TSV and PED outputs with full parameter provenance in the VCF header.
- Test suite and simulation benchmark.
