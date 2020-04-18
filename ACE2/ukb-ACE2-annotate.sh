# 18-4-2020 JHZ

export UKB=/rds/project/jmmh2/rds-jmmh2-post_qc_data/uk_biobank/imputed/uk10k_hrc/HRC_UK10K
export HPC_WORK=/rds/user/${USER}/hpc-work
export TMPDIR=${HPC_WORK}/work

export ANNOVAR=${HPC_WORK}/annovar
export LEFTEE=${HPC_WORK}/loftee
export POLYPHEN=${HPC_WORK}/polyphen-2.2.2
export VEP=${HPC_WORK}/ensembl-vep

awk 'NR==10||(NR>10 && $4>=15579156 && $4<=15620192)' ${UKB}/ukb_imp_chrX_v3_snpstats.txt > ukb-ACE2.snpstats

R --no-save -q <<END
  snpstats <- read.table("ukb-ACE2.snpstats",as.is=TRUE,header=TRUE)
  ord <- with(snpstats,order(chromosome,position))
  ct <- snpstats[ord,]
  all <- within(ct,{chr <- chromosome; pos <- position; snp <- rsid; a1 <- alleleA; a2 <- alleleB; qual <- "."; filter <- ".";})
  writetable <- function(d,f,...) write.table(d,file=f,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",...)
  for (f in c("ukb-ACE2"))
  {
    avinput <- paste0(f,".avinput")
    vars <- c("chr","pos","pos","a1","a2")
    writetable(all[vars],avinput)
    vepinput <- paste0(f,".vepinput")
    cat("##fileformat=VCFv4.0\n", file=vepinput)
    cat("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO\n",file=vepinput,append=TRUE,sep="\t")
    vars <- c("chr","pos","snp","a1","a2","qual","filter","info")
    writetable(all[vars],vepinput,append=TRUE)
  }
END

ln -sf /rds/user/jhz22/hpc-work/ensembl-vep/clinvar_GRCh37.vcf.gz
ln -sf /rds/user/jhz22/hpc-work/ensembl-vep/clinvar_GRCh37.vcf.gz.tbl

export wd=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/ACE2

for s in ukb-ACE2
do
   export s=${s}
   # ANNOVAR
   annotate_variation.pl -buildver hg19 ${s}.avinput ${ANNOVAR}/humandb/ -dbtype ensGene --outfile ${s}
   table_annovar.pl ${s}.avinput $ANNOVAR/test -buildver hg19 -out ${s} \
        -protocol ensGene,refGene,ccdsGene,wgEncodeGencodeBasicV19,cytoBand,exac03,avsnp147,dbnsfp30a,gwasCatalog \
        -operation g,g,g,g,r,f,f,f,r \
        -remove -nastring . -csvout -polish -xref $ANNOVA/example/gene_xref.txt
   # Polyphen-2
   awk 'NR>1{print "chr" $3 ":" $4, $5 "/" $6}' ${s}.snpstats | \
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
   vep -i ${wd}/${s}.vepinput -o ${wd}/${s}.dbNSFP --cache --distance 500000 --force --offline --pick --tab \
       --plugin dbNSFP,${VEP}/dbNSFP4.0a/dbNSFP4.0a.gz,${dbNSFP_fields}
   cd ${LEFTEE}
   vep -i ${wd}/${s}.vepinput -o ${wd}/${s}.loftee --cache --distance 500000 --force --offline --pick --tab \
       --plugin LoF,loftee_path:.,human_ancestor_fa:human_ancestor.fa.gz
   cd -
## Where the selected ClinVar INFO fields (from the ClinVar VCF file) are:
# - CLNSIG:     Clinical significance for this single variant
# - CLNREVSTAT: ClinVar review status for the Variation ID
# - CLNDN:      ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB
# Or the INFO fields you want in the ClinVar VCF file
  vep --i ukb-ACE2.vepinput --species homo_sapiens -o ukb-ACE2.clinvar --cache --offline --force_overwrite \
      --custom clinvar_GRCh37.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN --tab \
      --fields Uploaded_variation,Gene,Consequence,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN
done

function IPD()
{
  bgenix -g ${UKB}/ukb_imp_chrX_v3.bgen -list -incl-range X:15579156-15620192 > ukb-ACE2
  bgenix -g ${UKB}/ukb_imp_chrX_v3.bgen -vcf -incl-range X:15579156-15620192  > ukb-ACE2.vcf

}
