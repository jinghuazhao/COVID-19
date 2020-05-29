#!/usr/bin/bash

export ref=/home/jhz22/rds/post_qc_data/interval/reference_files/genetic/interval
export TMPDIR=/rds/user/jhz22/hpc-work/work

function do_vep()
{
  export VEP=${HPC_WORK}/ensembl-vep
  export chunk_size=100000
  seq 22 | \
  parallel -j2 --env ref --env VEP -C' ' '
    cd work
    (
      export n=$(wc -l $ref/impute_{}_interval.snpstats | cut -d" " -f1)
      export g=$(expr ${n} / ${chunk_size})
      for i in $(seq ${g})
      do
        (
          awk "BEGIN{print \"##fileformat=VCFv4.0\"}"
          awk -vOFS="\t" "BEGIN{print \"#CHROM\",\"POS\",\"ID\",\"REF\",\"ALT\",\"QUAL\",\"FILTER\",\"INFO\"}"
          if [ ${i} -eq ${g}]; then
             awk -v i=${i} -v n=${n} -v chunk_size=${chunk_size} "
                 NR==(i-1)*chunk_size,NR==n {print \$3,\$4,\$1,\$5,\$6,\".\",\".\",\$19}
                 " $ref/impute_{}_interval.snpstats
          else
             awk -v i=${i} -v chunk_size=${chunk_size} "
                 NR==(i-1)*chunk_size,NR==i*chunk_size {print \$3,\$4,\$1,\$5,\$6,\".\",\".\",\$19}
                 " $ref/impute_{}_interval.snpstats
          fi
        ) | \
        vep  --cache --dir_cache ${VEP}/.vep --offline --fork 4 --format vcf -o - --tab --pick --no_stats  \
             --species homo_sapiens --assembly GRCh37 --port 3337 | \
         grep -v '#'
      done
    ) | \
    gzip -f > INTERVAL-{}.vep.gz
    cd -
  '
}

function to_bed4()
{
  echo $(seq 22 -1 1) X | \
  tr ' ' '\n' | \
  parallel -C' ' '
    gunzip -c output/INTERVAL_annotation/INTERVAL-{}.vep.gz | \
    awk "NR>1" | \
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
  parallel -j4 -C' ' 'qctool -g work/INTERVAL-{}.bgen -annotate-bed4 work/INTERVAL-{}.bed4 -osnp work/INTERVAL-{}.annotate'
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
