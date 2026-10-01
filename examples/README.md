# Example: four-member family

Single-sample, variant-only VCFs of a father, a mother, an affected son (`PROBAND`) and an
unaffected daughter (`SISTER`), with one variant per inheritance pattern:

| Variant | Pattern | Found by |
| --- | --- | --- |
| chr1:100 C>T (GENEA) | proband 1/1, parents and sister 0/1 | `AR`, `AR_ANY` |
| chr1:200 G>A + chr1:300 G>C (GENEB) | paternal + maternal heterozygous in the proband | `CH`, `AR_ANY` |
| chr1:400 A>G (GENEC) | proband only | `DN`, `AD`, `SHARED` |
| chr1:500 A>G,T (GENED) | multi-allelic; T inherited from the unaffected father | none |
| chrX:5000000 C>T (GENEX) | hemizygous proband, carrier mother | `XLR` |
| chrM:150 C>T | proband, mother and sister | `MT` |

```bash
cosegvcf --vcf FATHER.vcf MOTHER.vcf PROBAND.vcf SISTER.vcf \
         --ped family.ped --model AR_ANY -o out_AR_ANY.vcf
```

In the GUI, add the four files, press **Import PED...** and select `family.ped`, then
choose a model and press **Run**.
