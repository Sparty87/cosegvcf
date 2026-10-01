# Simulation benchmark

`simulate.py` generates families with a known causal event, writes one variant-only VCF
per individual (about 19,000 records each, exome-sized) and runs CoSegVCF on them.

Scenarios: dominant (7-member, 3-generation pedigree), homozygous and compound
heterozygous recessive (trio, quartet), de novo (trio, quartet), X-linked recessive
(5 members). Conditions: error-free genotypes; noisy (0.2% genotype errors, 2% low-depth
calls); noisy with missing genotypes ignored. Population AF threshold: 1%.

```bash
# full run as in the paper (50 families per scenario and condition, ~70 min on 2 cores)
seq 0 7 | xargs -P 2 -I{} python benchmark/simulate.py --reps 50 --only {} --out benchmark/results
python benchmark/simulate.py --merge --out benchmark/results
python benchmark/plot_figure.py          # needs matplotlib
```

`results/results.tsv` holds one row per simulated family, `results/summary.tsv` one row per
scenario and condition, and `results/figure2.{pdf,png}` is Figure 2 of the paper.
Seeds are fixed (`--seed`, default 20260930, offset by scenario index), so runs are
reproducible.
