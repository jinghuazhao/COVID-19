#!/usr/bin/bash

# --- step1 ---

export TMPDIR=${HPC_WORK}/work
export WESWGS=${HOME}/COVID-19/WESWGS
export COHORT=INTERVAL
export PREFIX=
export PREFIX="singularity exec ${WESWGS}/burden_testing_latest.sif"

module load singularity
cd work
for panel in cvd2 cvd3 inf neu
do
  export panel=${panel}
  for weswgs in wes wgs
  do
    export weswgs=${weswgs}
    echo chr{1..22} chrX chrY | \
    tr ' ' '\n' | \
    parallel --env COHORT --env PREFIX --env panel --env weswgs -C' ' '${PREFIX} step1 ${COHORT}-${panel}-${weswgs} ${panel}-${weswgs}-{}.vcf.gz'
    for chr in chr{1..22} chrX chrY; do zcat ${COHORT}-${panel}-${weswgs}.${chr}.variantlist.gz; done | \
    sed 's/^chr//' | gzip -f > ${COHORT}-${panel}-${weswgs}.variantlist.gz
  done
done

${PREFIX} prepare-regions -o $(pwd)/geneset_data

# --- make (annotation) group files ---
#1. exon severe
#   variants with a "high" predicted consequence according to Ensembl (roughly equivalent to more severe than missense)
#2. exon CADD
#   all exonic variants (+50 bp outside of exons) weighted by CADD scores
#3. exon regulatory 
#   same as above, but weighted by Eigen scores (phred-scaled). Variants in regulatory regions that overlap with eQTL for that gene are also included.
#4. regulatory only
#   same as above but excluding exonic variants

export OPTS1="-g exon"
export OPTS2="-g exon -x 50 -s CADD"
export OPTS3="-g exon -x 50 -e promoter,enhancer,TF_bind -l promoter,enhancer,TF_bind -s EigenPhred"
export OPTS4="-e promoter,enhancer,TF_bind -l promoter,enhancer,TF_bind -s EigenPhred"

export chunks=2
for panel in cvd2 cvd3 inf neu
do
  for group in OPTS{1..4}
  do
    for i in $(seq $chunks)
    do
       ${PREFIX} make-group-file -C geneset_data/config.txt -i ${COHORT}-${panel}-wes.variantlist.gz ${!group} -o -w $(pwd) -d ${chunks} -c $i &
       ${PREFIX} make-group-file -C geneset_data/config.txt -i ${COHORT}-${panel}-wgs.variantlist.gz ${!group} -o -w $(pwd) -d ${chunks} -c $i &
    done
  done
done

find -name "group_file*.txt" -exec cat \{} \+ > group_file.txt
cut -f1 concat.group.file.txt | sort | uniq -c | awk '$1==1{print $2}'> singlesnp.genes.txt
fgrep -wvf singlesnp.genes.txt concat.group.file.txt > concat.group.file.filtered.txt


# --- step2 ---

# phenotype file

id         height
SAMPLE001  0.593
SAMPLE002  -0.135

# relatedness matrix files
# --- GRM generation ---

# GDS file
for i in chr{1..22} chrX chrY; do
  ${PREFIX} VCF2GDS ${COHORT}.vcf.gz ${COHORT}.${chr}.gds 10
done

# group file containing testing groups

${PREFIX} step2 --out test --threads 3

# <olink_protein>_<cohort>_<date_of_analysis>_<analyst_initials>.txt.bgz
# ACE2_INTERVAL_02112020_JHZ.txt.bgz

#Column no	Column name	Description
#1 	SNP	rsID, or NA if unavailable
#2 	CHR	Chromosome number
#3 	POS	Physical base pair position on the chromosome (in b38 coordinates)
#4 	N 	Number of non-missing observations
#5 	EFF_ALLELE 	Allele whose effect is reported (beta estimates)
#6	OTHER_ALLELE 	The other allele at the SNP
#7 	EFF_ALLELE_FREQ	Allele frequency of EFF_ALLELE
#8 	BETA	Effect size estimate, to at least 5 d.p.
#9	SE	Standard error of the beta estimates, to at least 5 d.p.
#10	P_LRT	P-value LRT
#11	P_SCORE	P-value score

# --- deprecated ---

function gihub()
{
cd burden_testing
for panel in cvd2 cvd3 inf neu
do
  export panel=${panel}
  for weswgs in wes wgs
  do
    export weswgs=${weswgs}
    echo chr{1..22} chrX chrY | \
    tr ' ' '\n' | \
    parallel --env COHORT --env panel --env weswgs -C' ' 'burden.1.6.5 step1 ${COHORT}-${panel}-${weswgs} ../work/${panel}-${weswgs}-{}.vcf.gz'
    for chr in chr{1..22} chrX chrY; do zcat ${COHORT}-${panel}-${weswgs}.${chr}.variantlist.gz; done | \
    sed 's/^chr//' | gzip -f > ${COHORT}-${panel}-${weswgs}.variantlist.gz
    tabix -f ${COHORT}-${panel}-${weswgs}.variantlist.gz
  # single_cohort_munge_variantlist ${COHORT}-${panel}-${weswgs}.variantlist.gz 1 1
  done
done

# --- geneset (annotation) data ---
# install axel, moreutils
# http://www.tucows.com/preview/231886/Axel
# git://git.joeyh.name/moreutils
# http://downloads.sourceforge.net/docbook2x/docbook2X-0.8.8.tar.gz
# ftp://xmlsoft.org/libxml2/
# ftp://xmlsoft.org/libxml2/libxslt-1.1.34.tar.gz
# ftp://xmlsoft.org/libxml2/libxml2-2.9.10.tar.gz
# ln -sf geneset_data/ensembl-vep/INSTALL.pl
}
