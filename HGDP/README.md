# Mapping data to reference population

## I. GIANT implementation

### GIANT reference (dated 8/8/2014)

* `refdata.*` (N=437, # of variants=79,399)
* `MDSplot.R`

and examples are in `EPIC-Omics.sh` and `HGI.sh` produces [`HGI.C1-C2.png`](HGI.C1-C2.png).

### 1000Genomes phase 3

* `keep.dat` contains a list of RSid's derived from HGI variants (# of variants=117,220).
* `keep.bed/bim/fam` contains genotypes for N=2,504 samples (too big to upload here).

The reference data were downloaded from files distributed with VEGAS2, 
[https://vegas2.qimrberghofer.edu.au/](https://vegas2.qimrberghofer.edu.au/), 
and extracted via the following script,

```bash
#!/usr/bin/bash

for pop in AFR AMR EAS EUR SAS
do
  echo ${pop}
  plink --bfile g1000p3_${pop} --extract keep.dat --make-bed --out ${pop}
  awk -v pop=${pop} '{$1=pop};1' ${pop}.fam > ${pop}-id.fam
done

echo AFR AMR EAS EUR SAS | tr ' ' '\n' | awk '{print $1".bed", $1".bim", $1"-id.fam"}' > keep.list
plink --merge-list keep.list --out keep
```

* `nonEUR.sh/sb/R` is the pipeline to produce [`nonEUR.C1-C2.png`](nonEUR.C1-C2.png).
* `keep.sh/sb/R` is the pipeline to produce [`INTERVAL.C1-C2.png`](INTERVAL.C1-C2.png) (not working properly).

## II. HGI implementation (dated 21/3/2021)

To project every GWAS participant into the same PC space, we used pre-computed PC loadings and
reference allele frequencies. For reference, we used unrelated samples from the 1000 Genomes Project and
the Human Genome Diversity Project (HGDP) and computed PC loadings and allele frequencies for the
117,221 SNPs that are i) available in every cohort, ii) MAF > 0.1% in the reference, and iii) LD pruned (r2
< 0.8; 500kb window). We then asked each cohort to project their samples using our automated script
provided here ([https://github.com/covid19-hg/pca_projection](https://github.com/covid19-hg/pca_projection)). It internally uses PLINK2 --score
function with variance-standardize option and reference allele frequencies (--read-freq); so that each
cohort specific genotype/dosage matrix is mean-centered and variance-standardized with regards to reference
allele frequencies, not cohort-specific allele frequencies. We further normalized the projected PC scores by
dividing by a square root of the number of variants used for projection to account for a subtle difference
due to missing variants.

Results: [hgdp_tgp_interval.PC1-2.pdf](hgdp_tgp_interval.PC1-2.pdf).
