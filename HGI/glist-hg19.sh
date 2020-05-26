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

function annotate()
{
  export TMPDIR=/rds/user/jhz22/hpc-work/work
  echo $(seq 22 -1 1) X | \
  tr ' ' '\n' | \
  parallel -j2 --env autosomes -C' ' '
    cd work
    (
      awk "BEGIN{print \"##fileformat=VCFv4.0\"}"
      awk -vOFS="\t" "BEGIN{print \"#CHROM\",\"POS\",\"ID\",\"REF\",\"ALT\",\"QUAL\",\"FILTER\",\"INFO\"}"
      if [ "{}" == "X" ]; then
         export vcfgz=${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz; else export vcfgz=${autosomes}/{}.pbwt_reference_impute.vcf.gz
      fi
      bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/INFO\n" ${vcfgz}  | \
      awk -v OFS="\t" "NR>1{print \$1,\$2,\$1 \":\" \$2 \"_\" \$3 \"/\" \$4, \$3, \$4, \$5, \$6, \$7}"
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
  export ANNOVAR=${HPC_WORK}/annovar
  export LOFTEE=${HPC_WORK}/loftee
  export POLYPHEN=$HPC_WORK/polyphen-2.2.2
  export VEP=${HPC_WORK}/ensembl-vep
  export APV=VEP
  echo $(seq 22 -1 7) $(seq 6 -1 1 | parallel 'echo {}.01 {}.02') X | \
  tr ' ' '\n' | \
  parallel -j2 --env autosomes --env ANNOVAR --env POLYPHEN --env LOFTEE --env VEP --env APV -C' ' '
    cd work
    if [ ${APV} == "VEP" ]; then
   # VEP
      vep --input_file INTERVAL-{}.vepinput.gz --output_file INTERVAL-{}.tsv --cache --dir_cache ${VEP}/.vep --dir_plugins ${LOFTEE} --offline \
          --pick --force_overwrite --species homo_sapiens --assembly GRCh37 --tab
      #   --plugin LoF,loftee_path:${LOFTEE},human_ancestor_fa:human_ancestor.fa.gz,conservation_file:phylocsf_gerp.sql.gz
    elif [ ${APV} == "ANNOVAR" ]; then
   # ANNOVAR
      gunzip -c INTERVAL-${s}.vepinput.gz | \
      awk -v OFS="\t" "NR>2{print \$1,\$2,\$2,\$4,\$5}" > INTERVAL-${s}.avinput
      annotate_variation.pl -buildver hg19 INTERVAL-${s}.avinput ${ANNOVAR}/humandb/ -dbtype ensGene --outfile INTERVAL-${s}
      table_annovar.pl INTERVAL-${s}.avinput $ANNOVAR/test -buildver hg19 -out ${s} \
           -protocol ensGene,refGene,ccdsGene,wgEncodeGencodeBasicV19,cytoBand,exac03,avsnp147,dbnsfp30a,gwasCatalog \
           -operation g,g,g,g,r,f,f,f,r \
           -remove -nastring . -csvout -polish -xref $ANNOVA/example/gene_xref.txt
    else 
   # Polyphen-2
      gunzip -c INTERVAL-${s}.vepintput | \
      awk "NR>2{gsub(/_/,\" \",\$3);print \"chr\" \$3}" | \
      sort -k1,1 | \
      uniq > INTERVAL-${s}.pph.list
      mapsnps.pl -g hg19 -m -U -y INTERVAL-${s}.pph.input INTERVAL-${s}.pph.list 1>INTERVAL-${s}.pph.features 2>INTERVAL-${s}.log
      run_pph.pl INTERVAL-${s}.pph.input 1>$INTERVAL-{s}.pph.output 2>INTERVAL-${s}.pph.log
      run_weka.pl INTERVAL-${s}.pph.output >INTERVAL-${s}.pph.humdiv.output
      run_weka.pl -l $POLYPHEN/models/HumVar.UniRef100.NBd.f11.model INTERVAL-${s}.pph.output >INTERVAL-${s}.pph.humvar.output
    fi
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
    qctool -g ${autosomes}/imputed/impute_{}_interval.bgen -annotate-bed4 work/INTERVAL-{}.bed4 -osnp work/INTERVAL-{}.annotate
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
