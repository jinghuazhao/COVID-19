#!/usr/bin/bash

function glist_hg19()
{
  export glist=${INF}/csd3/glist-hg19
  (
    awk '$1!="X" && $1!="Y" && $1!="XY"' ${glist} | sort -k1,1n -k2,2n
    awk '$1=="X"' ${glist} | sort -k2,2n
    awk '$1=="XY"' ${glist} | sort -k2,2n
    awk '$1=="Y"' ${glist} | sort -k2,2n
  ) | \
  tr ' ' '\t' > work/glist-hg19.bed
}

function biomart()
{
  R --no-save -q <<\ \ END
    library(biomaRt)
    ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
    attr <- listAttributes(ensembl)
    gene <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'description', 'hgnc_symbol'), mart = ensembl)
    names(gene) <- c('ensGene','chrom','chromStart','chromEnd','desc','hgnc')
    select <- with(gene, chrom%in%c(paste(1:22),"X","XY","Y"))
    cols <- c(2:4,1,6,5)
    write.table(gene[select,cols],file="work/glist-hg19.biomart",quote=FALSE,row.names=FALSE,sep="\t")
  END
  export glist=work/glist-hg19.biomart
  (
    awk '$1!="X" && $1!="Y" && $1!="XY"' ${glist} | sort -k1,1n -k2,2n
    awk '$1=="X"' ${glist} | sort -k2,2n
    awk '$1=="XY"' ${glist} | sort -k2,2n
    awk '$1=="Y"' ${glist} | sort -k2,2n
  ) | \
  awk -vOFS='\t' '{print $1,$2,$3,$4 "_" $5 "_pLoF"}' > work/glist-hg19.bed
}

function gencode_v19()
{
  R --no-save -q <<\ \ END
    url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz"
    gtf <- rtracklayer::import("work//gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz")
    gencode <- as.data.frame(gtf)
    select <- with(gencode, seqnames%in%paste0("chr",c(paste(1:22),"X","XY","Y")))
    cols <- c(1:3,10,14)
    write.table(gencode[select,cols], file="work/glist-hg19.gencode",quote=FALSE,row.names=FALSE,sep="\t")
  END
  export glist=work/glist-hg19.gencode
  (
    awk '{sub(/chr/,"");if($1!="X" && $1!="Y" && $1!="XY") print}' ${glist} | sort -k1,1n -k2,2n
    awk '{sub(/chr/,"");if($1=="X") print}' ${glist} | sort -k2,2n
  ) | \
  awk -vOFS='\t' 'NR>1{print $1,$2,$3,$4 "_" $5}' > work/glist-hg19.bed
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
    }" > work/INTERVAL-{}.bed
  '
}

function glist_enshgnc()
{
  seq 22 | \
  parallel -j1 -C' ' '
    qctool -g work/INTERVAL-{}.bgen -annotate-bed4 work/INTERVAL-{}.bed -osnp work/INTERVAL-{}.annotate
  '
  export X=/rds/project/jmmh2/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data
  qctool -g ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz -filetype vcf  -annotate-bed4 work/INTERVAL-X.bed -osnp work/INTERVAL-X.annotate
  echo $(seq 22) X | \
  tr ' ' '\n' | \
  parallel -j1 -C' ' '
     awk "NR>9 && \$8!=\"NA\" && \$1!=\".\" && \$1!=\"#\"{print \$1}" work/INTERVAL-{}.annotate > work/INTERVAL-{}.incl
     export list=($(awk "NR>8 && \$8!=\"NA\"" work/INTERVAL-{}.annotate | cut -f8 | sort | uniq))
     (
       for g in ${list[@]}
       do
          awk -v g=${g} "\$8==g" work/INTERVAL-{}.annotate | \
          awk -vOFS="\t" "\$1!=\".\" {printf OFS \$1}" | \
          awk -v g=${g} -v OFS="\t" "{print g,\$0}"
       done
     ) > work/INTERVAL-{}.gene
  '
}

function snpstats()
{
  (
    cut -f1,7,8,15,18,19 $ref/impute_*_interval.snpstats | \
    head -1
    seq 22 | \
    parallel -j1 --env ref -C' ' '
      sed "1d" $ref/impute_{}_interval.snpstats | \
      cut -f1,3,4,7,8,15,18,19
    '
  ) | gzip -f > work/INTERVAL.snpstats.gz
}
