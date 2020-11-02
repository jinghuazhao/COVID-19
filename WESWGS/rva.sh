#!/usr/bin/bash

# inputs
# - Filtered VCF file (see above section on Variant QC) OR GDS files if they already exist
# - Phenotype files containing two columns: sample ID and INT-transformed, covariate-adjusted and renormalised residuals without header
# - Genetic relatedness matrix (GRM) with sample IDs as row and column names

# outputs
# - A space-delimited file containing single variant scores
# - A binary file containing between-variant covariance matrices

# Meta-analysis

# Group file
#Column no.	Column name	Description
#1	CHR	Chromosome number (1-22)
#2	POS	Physical base-pair position on the chromosome (in b38 coordinates)
#3	REF	Reference allele
#4	ALT	Alternate allele
#5	AC	Allele count in genotypes for each ALT allele
#6	AN	Total number of alleles in called genotypes
