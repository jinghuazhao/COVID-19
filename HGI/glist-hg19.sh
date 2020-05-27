#!/usr/bin/bash

function VEP()
{
# pbwt is incomplete!
  export TMPDIR=/rds/user/jhz22/hpc-work/work
  echo $(seq 22 -1 1) X | \
  tr ' ' '\n' | \
  parallel -j2 --env autosomes -C' ' '
    cd work
    (
      awk "BEGIN{print \"##fileformat=VCFv4.0\"}"
      awk -vOFS="\t" "BEGIN{print \"#CHROM\",\"POS\",\"ID\",\"REF\",\"ALT\",\"QUAL\",\"FILTER\",\"INFO\"}"
      if [ "{}" != "X" ]; then
         bgenix -g ${autosomes}/imputed/impute_{}_interval.bgen -list 2>&1 | \
         awk -vOFS="\t" "NR>9&&NF==7{print \$3+0,\$4,\$2,\$6,\$7,\".\",\".\",\".\"}"
      else
         export vcfgz=${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz;
         bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/INFO\n" ${vcfgz}  | \
         awk -v OFS="\t" "NR>1{print \$1,\$2,\$1 \":\" \$2 \"_\" \$3 \"/\" \$4, \$3, \$4, \$5, \$6, \$7}"
      fi
    ) | \
    gzip -f > INTERVAL-{}.vepinput.gz
  # Split large chromosomes into two chunks (at most comparable to chromosome 7)
    if [ {} -le 6 ] && [ "{}" != "X" ]; then
      gunzip -c INTERVAL-{}.vepinput.gz | \
      split -l 5000000 --numeric-suffixes=1 --additional-suffix=.vepinput - INTERVAL-{}.
      gzip -f INTERVAL-{}.01.vepinput
      (
        gunzip -c INTERVAL-{}.vepinput.gz | \
        awk "NR<3{print}"
        cat INTERVAL-{}.02.vepinput
        rm INTERVAL-{}.02.vepinput
      ) | \
      gzip -f > INTERVAL-{}.02.vepinput.gz
    fi
  '
  export VEP=${HPC_WORK}/ensembl-vep
  echo $(seq 22 -1 7) $(seq 6 -1 1 | parallel 'echo {}.01 {}.02') X | \
  tr ' ' '\n' | \
  parallel -j2 --env autosomes --env VEP -C' ' '
    cd work
    vep --input_file INTERVAL-{}.vepinput.gz --output_file INTERVAL-{}.tsv --cache --dir_cache ${VEP}/.vep --offline \
        --pick --force_overwrite --species homo_sapiens --assembly GRCh37 --tab
    cd -
  '
}

function Ensembl_GRCh37()
{
  echo $(seq 22 -1 1) X | \
  tr ' ' '\n' | \
  parallel -C' ' '
    (
      if [ "{}" != "X" ] && [ {} -lt 7 ]; then
        gunzip -c output/INTERVAL_annotation/INTERVAL-{}-1.txt.gz | \
        awk "NR>1"
        gunzip -c output/INTERVAL_annotation/INTERVAL-{}-2.txt.gz | \
        awk "NR>1"
      else
        gunzip -c output/INTERVAL_annotation/INTERVAL-{}.txt.gz | \
        awk "NR>1"
      fi
    ) | \
    cut -f1,2,4,6,7 | \
    awk -v OFS="\t" "\$5!=\"-\" {
      sub(/frameshift_variant|stop_gained|splice_acceptor_variant|splice_donor_variant/,\"pLoF\",\$3);
      split(\$2,a,\":\")
      split(a[2],b,\"-\")
      print a[1],b[1]-1,b[2],\$5 \":\" \$4 \":\" \$3
    }" > work/INTERVAL-{}.bed4
  '
}

function glist_annotate()
{
  export autosomes=/home/jhz22/rds/post_qc_data/interval/imputed/uk10k_1000g_b37
  seq 22 | \
  parallel -j4 -C' ' '
    qctool -g work/INTERVAL-{}.bgen -annotate-bed4 work/INTERVAL-{}.bed4 -osnp work/INTERVAL-{}.annotate
  '
  export X=/rds/project/jmmh2/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data
  qctool -g ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz -filetype vcf -annotate-bed4 work/INTERVAL-X.bed4 -osnp work/INTERVAL-X.annotate
  echo $(seq 22) X | \
  tr ' ' '\n' | \
  parallel -j4 -C' ' '
     awk "NR>9 && \$8!=\"NA\" && \$2!=\".\" && \$1!=\"#\"{print \$1}" work/INTERVAL-{}.annotate > work/INTERVAL-{}.incl
     export list=($(awk "NR>8 && \$8!=\"NA\"" work/INTERVAL-{}.annotate | cut -f8 | sort | uniq))
     (
       for g in ${list[@]}
       do
          awk -v g=${g} "\$8==g" work/INTERVAL-{}.annotate | \
          awk -vOFS="\t" "\$2!=\".\" {printf OFS \$2}" | \
          awk -v g=${g} -v OFS="\t" "{print g,\$0}"
       done
     ) > work/INTERVAL-{}.gene
  '
}

Ensembl_GRCh37
glist_annotate
