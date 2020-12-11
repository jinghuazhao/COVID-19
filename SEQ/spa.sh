#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export SEQ=${HOME}/COVID-19/SEQ
export COHORT=INTERVAL
export CONFIG=${SEQ}/geneset_data/config.txt
export STEP1="singularity exec ${SEQ}/burden_testing_latest.sif"
export STEP2="singularity exec --containall ${SEQ}/burden_testing_latest.sif"

cd work

# --- step1 ---

# obtain variant lists

for weswgs in wes wgs
do
  export weswgs=${weswgs}
  echo chr{1..22} chrX chrY | \
  tr ' ' '\n' | \
  parallel --env COHORT --env STEP1 --env weswgs -C' ' '${STEP1} step1 ${COHORT}-${weswgs} ${weswgs}-{}.vcf.gz'
  (
    echo -e "#CHROM\tPOS\tREF\tALT\tAC\tAN"
    for chr in chr{1..22} chrX chrY; do zcat ${COHORT}-${weswgs}.${chr}.variantlist.gz; done | \
    sed 's/^chr//' 
  ) | gzip -f > ${COHORT}-${weswgs}.variantlist.gz
done

join -v1 -t$'\t' <(gunzip -c ${COHORT}-wes.variantlist.gz | awk -vOFS="\t" 'NR>1 {snpid=$1":"$2"_"$3"_"$4;print snpid,$0}' | sort -k1,1) \
                 <(gunzip -c ${COHORT}-wgs.variantlist.gz | awk -vOFS="\t" 'NR>1 {snpid=$1":"$2"_"$3"_"$4;print snpid,$0}' | sort -k1,1) | \
sort -k2,2n -k3,3n | \
cut -f1 --complement | \
gzip -f > ${COHORT}-wes-wgs.variantlist.gz

(
  gunzip -c ${COHORT}-wgs.variantlist.gz | head -1
  gunzip -c ${COHORT}-wgs.variantlist.gz ${COHORT}-wes-wgs.variantlist.gz | sed '1d;s/X/23/;s/Y/24/' | sort -k1,1n -k2,2n | sed 's/23/X/;s/24/Y/'
) | \
gzip -f > ${COHORT}-wes+wgs.variantlist.gz

# prepare regions

${STEP1} prepare-regions -o $(pwd)/geneset_data

# make (annotation) group files
#
#1. exon severe
#   variants with a "high" predicted consequence according to Ensembl (roughly equivalent to more severe than missense)
#2. exon CADD
#   all exonic variants (+50 bp outside of exons) weighted by CADD scores
#3. exon regulatory
#   weighted by Eigen scores (phred-scaled). Variants in regulatory regions that overlap with eQTL for that gene are also included.
#4. regulatory only
#   excluding exonic variants

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
       ${STEP1} make-group-file -C ${CONFIG} -i ${COHORT}-${panel}-wes.variantlist.gz ${!group} -o -w $(pwd) -d ${chunks} -c $i &
       ${STEP1} make-group-file -C ${CONFIG} -i ${COHORT}-${panel}-wgs.variantlist.gz ${!group} -o -w $(pwd) -d ${chunks} -c $i &
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
  ${STEP2} VCF2GDS ${COHORT}.vcf.gz ${COHORT}.${chr}.gds 10
done

# group file containing testing groups

${STEP2} step2 -h

cd -

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

# http://www.tucows.com/preview/231886/Axel
# git://git.joeyh.name/moreutils
# sudo apt install docbook2x
# sudo apt install libxml2-utils
# git clone git://git.joeyh.name/moreutils
# cd moreutils
# make
# make install
# partial alternatives
# http://downloads.sourceforge.net/docbook2x/docbook2X-0.8.8.tar.gz
# ftp://xmlsoft.org/libxml2/
# ftp://xmlsoft.org/libxml2/libxslt-1.1.34.tar.gz
# ftp://xmlsoft.org/libxml2/libxml2-2.9.10.tar.gz
# ln -sf geneset_data/ensembl-vep/INSTALL.pl
}