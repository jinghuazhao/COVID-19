# 29-3-2020 JHZ

ls GSM*h5 | \
awk '/Donor/ && /filtered/{gsub(/_filtered_gene_bc_matrices_h5.h5/,"",$1);print}' > filtered.list

# Seurat tutorial
# https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
# https://www.dropbox.com/s/y6kwho066vugjrg/pbmc33k.R
# GSE122960
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122960
