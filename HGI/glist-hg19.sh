#!/usr/bin/bash

export ref=/home/jhz22/rds/post_qc_data/interval/reference_files/genetic/interval
export X=/rds/project/jmmh2/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data
export TMPDIR=/rds/user/jhz22/hpc-work/work
export chunk_size=10000

function local_vep_X()
# local annotation to guarantee success
{
  gunzip -c ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz | \
  bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/INFO\n" | \
  awk -v OFS="\t" "NR>1{print \$1,\$2,\$1 \":\" \$2 \"_\" \$3 \"/\" \$4, \$3, \$4, \$5, \$6, \$7}" > work/INTERVAL-X.query
  export n=$(wc -l work/INTERVAL-X.query | cut -d" " -f1)
  export g=$(expr ${n} / ${chunk_size})
  export s=$(expr $n - 1 - $(($g * $chunk_size)))
  (
    for i in $(seq ${g}); do
    (
      awk "BEGIN{print \"##fileformat=VCFv4.0\"}"
      awk -vOFS="\t" "BEGIN{print \"#CHROM\",\"POS\",\"ID\",\"REF\",\"ALT\",\"QUAL\",\"FILTER\",\"INFO\"}"
      awk -v i=${i} -v chunk_size=${chunk_size} -v OFS="\t" "NR==(i-1)*chunk_size+1,NR==i*chunk_size" work/INTERVAL-X.query
      if [ ${s} -gt 0 ] && [ ${i} -eq ${g} ]; then
         awk -v i=${i} -v chunk_size=${chunk_size} -v OFS="\t" -v n=${n} "NR==i*chunk_size+1,NR==n" work/INTERVAL-X.query
      fi
    ) | \
    vep  --cache --offline --format vcf -o - --tab --pick --no_stats --symbol \
         --species homo_sapiens --assembly GRCh37 --port 3337 | \
    (
      if [ ${i} -eq 1 ]; then cat; else grep -v "#"; fi
    ) 
    done
  ) | \
  gzip -f > work/INTERVAL-X.vep.gz
}

function local_vep()
{
  seq 22 | \
  parallel -j1 --env ref -C' ' '
    export n=$(awk "END{print NR-1}")
    export g=$(expr ${n} / ${chunk_size})
    export s=$(expr $n - $(($g * $chunk_size)))
    (
      for i in $(seq ${g}); do
        (
          awk "BEGIN{print \"##fileformat=VCFv4.0\"}"
          awk -vOFS="\t" "BEGIN{print \"#CHROM\",\"POS\",\"ID\",\"REF\",\"ALT\",\"QUAL\",\"FILTER\",\"INFO\"}"
          (
            sed "1d" ${ref}/impute_{}_interval.snpstats | \ 
            awk -v i=${i} -v chunk_size=${chunk_size} -v OFS="\t" "NR==(i-1)*chunk_size+1,NR==i*chunk_size {
                if(\$1==\".\") \$1=\$3+0 \":\" \$4 \"_\" \$5 \"/\" \$6; print \$3+0,\$4,\$1,\$5,\$6,\".\",\".\",\$19}"
            if [ ${s} -gt 0 ] && [ ${i} -eq ${g} ]; then
             sed "1d" ${ref}/impute_{}_interval.snpstats | \ 
             awk -v i=${i} -v chunk_size=${chunk_size} -v OFS="\t" -v n=${n} "NR==i*chunk_size+1,NR==n {
                   if(\$1==\".\") \$1=\$3+0 \":\" \$4 \"_\" \$5 \"/\" \$6; print \$3+0,\$4,\$1,\$5,\$6,\".\",\".\",\$19}"
            fi
          )
        ) | \
        vep  --cache --offline --format vcf -o - --tab --pick --no_stats --symbol \
             --species homo_sapiens --assembly GRCh37 --port 3337 | \
        (
          if [ ${i} -eq 1 ]; then cat; else grep -v "#"; fi
        )
      done
    ) | \
    gzip -f > work/INTERVAL-{}.vep.gz
  '
}

function to_bed4()
{
  echo $(seq 22 -1 1) X | \
  tr ' ' '\n' | \
  parallel -C' ' '
    gunzip -c work/INTERVAL_annotation/INTERVAL-{}.vep.gz | \
    awk "NR>43" | \
    cut -f1,2,4,7,18 | \
    awk -v OFS="\t" "\$3!=\"-\" {
      sub(/frameshift_variant|stop_gained|splice_acceptor_variant|splice_donor_variant/,\"pLoF\",\$4);
      split(\$2,a,\":\")
      split(a[2],b,\"-\")
      if (a[1]<10) printf 0
      print a[1],b[1]-1,b[1],\$3 \":\" \$5 \":\" \$4
    }" > work/INTERVAL-{}.bed4
  '
}

function glist_annotate()
{
  seq 22 | \
  parallel -j1 --env src -C' ' 'qctool -g work/INTERVAL-{}.bgen -annotate-bed4 work/INTERVAL-{}.bed4 -osnp work/INTERVAL-{}.annotate'
  sbatch glist-hg19.sb
  qctool -g work/INTERVAL-X-ploidy.vcf.gz -filetype vcf -annotate-bed4 work/INTERVAL-X.bed4 -osnp work/INTERVAL-X.annotate
  echo X | \
  parallel -j3 -C' ' '
     awk "NR>10 && \$8!=\"NA\" && \$2!=\".\" && \$1!=\"#\"{print \$2}" work/INTERVAL-{}.annotate > work/INTERVAL-{}.incl
     export list=($(awk "NR>10 && \$8!=\"NA\"" work/INTERVAL-{}.annotate | cut -f8 | sort | uniq))
     (
       for g in ${list[@]}
       do
          awk -v g=${g} "\$8==g" work/INTERVAL-{}.annotate | \
          awk -vOFS="\t" "\$2!=\".\" {printf OFS \$2}" | \
          awk -v g=${g} -v OFS="\t" "{print g \$0}"
       done
     ) > work/INTERVAL-{}.gene
     (
       for g in ${list[@]}
       do
          awk -v g=${g} "\$8==g" work/INTERVAL-{}.annotate | \
          awk -vOFS="\t" "\$2!=\".\" {split(\$2,a,\"_\");printf OFS a[1] \"_\" a[2] \"/\" a[3]}" | \
          awk -v g=${g} -v OFS="\t" "{print g \$0}"
       done
     ) > work/INTERVAL-{}.gene-snpid
     sed -i 's/\t/;/g;s/;/\t/' ../work/INTERVAL-{}.gene-snpid
  '
}
