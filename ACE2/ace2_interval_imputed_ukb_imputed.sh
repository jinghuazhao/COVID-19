# 23-4-2020 JHZ

export HPC_WORK=/rds/user/${USER}/hpc-work
export TMPDIR=${HPC_WORK}/work

export ANNOVAR=${HPC_WORK}/annovar
export LEFTEE=${HPC_WORK}/loftee
export POLYPHEN=${HPC_WORK}/polyphen-2.2.2
export VEP=${HPC_WORK}/ensembl-vep

function init()
{
# https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000130234;r=X:15579156-15620271
# ln -sf /rds/user/jhz22/hpc-work/ensembl-vep/clinvar_GRCh37.vcf.gz
# ln -sf /rds/user/jhz22/hpc-work/ensembl-vep/clinvar_GRCh37.vcf.gz.tbl
  R --no-save -q <<\ \ END
  snpstats <- read.table("ace2_interval_imputed_ukb_imputed_filtered_vars.txt",as.is=TRUE,col.names=c("chr","pos","ref","alt"))
  ord <- with(snpstats,order(chr,pos))
  ct <- snpstats[ord,]
  all <- within(ct,{snp <- paste0("chr",chr,":",pos,"_",ref,"_",alt); qual <- "."; filter <- "."; info <- "."})
  writetable <- function(d,f,...) write.table(d,file=f,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",...)
  for (f in c("ace2_interval_imputed_ukb_imputed"))
  {
    avinput <- paste0(f,".avinput")
    vars <- c("chr","pos","pos","ref","alt")
    writetable(all[vars],avinput)
    vepinput <- paste0(f,".vepinput")
    cat("##fileformat=VCFv4.0\n", file=vepinput)
    cat("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO\n",file=vepinput,append=TRUE,sep="\t")
    vars <- c("chr","pos","snp","ref","alt","qual","filter","info")
    writetable(all[vars],vepinput,append=TRUE)
  }
  END
# CADD
# https://cadd.gs.washington.edu/score
  cat ace2_interval_imputed_ukb_imputed.vepinput | gzip -f > ace2_interval_imputed_ukb_imputed.vcf.gz
# PROVEAN
# http://provean.jcvi.org/index.php
  sed '1,2d' ace2_interval_imputed_ukb_imputed.vepinput | awk -vOFS="," '{print $1,$2,$4,$5}' | xsel -i
}
for s in ace2_interval_imputed_ukb_imputed
do
   export s=${s}
 # ANNOVAR
   annotate_variation.pl -buildver hg19 ${s}.avinput ${ANNOVAR}/humandb/ -dbtype ensGene --outfile ${s}
   table_annovar.pl ${s}.avinput $ANNOVAR/test -buildver hg19 -out ${s} \
        -protocol ensGene,refGene,ccdsGene,wgEncodeGencodeBasicV19,cytoBand,exac03,avsnp147,dbnsfp30a,gwasCatalog \
        -operation g,g,g,g,r,f,f,f,r \
        -remove -nastring . -csvout -polish -xref $ANNOVA/example/gene_xref.txt
 # Polyphen-2
   awk 'NR>1{print "chr" $1 ":" $2, $3 "/" $4}' ${s}_filtered_vars.txt | \
   uniq > ${s}.pph.list
   mapsnps.pl -g hg19 -m -U -y ${s}.pph.input ${s}.pph.list 1>${s}.pph.features 2>${s}.log
   run_pph.pl ${s}.pph.input 1>${s}.pph.output 2>${s}.pph.log
   run_weka.pl ${s}.pph.output >${s}.pph.humdiv.output
   run_weka.pl -l $POLYPHEN/models/HumVar.UniRef100.NBd.f11.model ${s}.pph.output >${s}.pph.humvar.output
 # VEP
   export dbNSFP_1=clinvar_id,clinvar_clnsig,clinvar_review,clinvar_trait,1000Gp3_EUR_AF,CADD_phred,Eigen-PC-phred_coding,ExAC_NFE_AF,LRT_pred,
   export dbNSFP_2=FATHMM_pred,GERP++_RS,GTEx_V7_tissue,MutPred_protID,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,SIFT4G_pred,fathmm-MKL_coding_pred,
   export dbNSFP_3=rs_dbSNP151,fathmm-MKL_coding_pred,gnomAD_exomes_NFE_AF,gnomAD_genomes_NFE_AF
   export dbNSFP_fields=${dbNSFP_1}${dbNSFP_2}${dbNSFP_3}
   vep -i ${s}.vepinput -o ${s}.dbNSFP --cache --distance 500000 --force --offline --pick --tab \
       --plugin dbNSFP,${VEP}/dbNSFP4.0a/dbNSFP4.0a.gz,${dbNSFP_fields}
## Where the selected ClinVar INFO fields (from the ClinVar VCF file) are:
# - CLNSIG:     Clinical significance for this single variant
# - CLNREVSTAT: ClinVar review status for the Variation ID
# - CLNDN:      ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB
# Or the INFO fields you want in the ClinVar VCF file
  vep --i ukb-ACE2.vepinput --species homo_sapiens -o ukb-ACE2.clinvar --cache --offline --force_overwrite \
      --custom clinvar_GRCh37.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN --tab \
      --fields Uploaded_variation,Gene,Consequence,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN
  vep -i ${s}.vepinput -o ${s}.loftee --cache --distance 500000 --force --offline --pick --tab \
      --plugin LoF,loftee_path:${LOFTEE},human_ancestor_fa:human_ancestor.fa.gz,conservation_file:phylocsf_gerp.sql.gz
  R --no-save -q <<\ \ END
    loftee <- read.delim("ace2_interval_imputed_ukb_imputed.loftee",as.is=TRUE,skip=29)
    CSQ <- table(with(loftee,Consequence))
    lbls <- paste(names(CSQ), "\n", CSQ, sep="")
    pie(CSQ, labels = lbls, main="Pie Chart of Codons\n (with sample sizes)") 
  END
done
