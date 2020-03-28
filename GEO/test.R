# 28-3-2020 JHZ

library(rhdf5)
m <- function(rt,type="_filtered",suffix="_gene_bc_matrices_h5.h5")
{
  f <- paste0(rt,type,suffix)
  h5ls(f)
  h5readAttributes(f,"GRCh38")
  d <- h5read(f,name="GRCh38", read.attributes = TRUE)
  h5closeAll()
  d
}

rt <-c(
"GSM3489182_Donor_01",
"GSM3489182_Donor_01",
"GSM3489183_IPF_01",
"GSM3489183_IPF_01",
"GSM3489184_IPF_02",
"GSM3489184_IPF_02",
"GSM3489185_Donor_02",
"GSM3489185_Donor_02",
"GSM3489186_Cryobiopsy_01",
"GSM3489186_Cryobiopsy_01",
"GSM3489187_Donor_03",
"GSM3489187_Donor_03",
"GSM3489188_IPF_03",
"GSM3489188_IPF_03",
"GSM3489189_Donor_04",
"GSM3489189_Donor_04",
"GSM3489190_IPF_04",
"GSM3489190_IPF_04",
"GSM3489191_Donor_05",
"GSM3489191_Donor_05",
"GSM3489192_HP_01",
"GSM3489192_HP_01",
"GSM3489193_Donor_06",
"GSM3489193_Donor_06",
"GSM3489194_SSc-ILD_01",
"GSM3489194_SSc-ILD_01",
"GSM3489195_Donor_07",
"GSM3489195_Donor_07",
"GSM3489196_Myositis-ILD_01",
"GSM3489196_Myositis-ILD_01",
"GSM3489197_Donor_08",
"GSM3489197_Donor_08",
"GSM3489198_SSc-ILD_02",
"GSM3489198_SSc-ILD_02")

donorID <- grep ("Donor",rt)
filtered <- rt[donorID][2*(1:8)]

m_barcodes <- m_data <- m_gene_names <- m_genes <- m_indices <- m_indptr <- list()
for(i in 1:8)
{
 d[i] <- m(filtered[i])
 attach(d[i])
 m_barcodes[[i]] <- barcodes
 m_data[[i]] <- data
 m_gene_names[[i]] <- gene_names
 m_genes[[i]] <- gene_names
 m_indices[[i]] <- indices
 m_indptr[[i]] <- indptr
 detach(d[i])
 lapply(d[i],length)
}
lapply(m_data,length)
