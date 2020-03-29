# 29-3-2020 JHZ

ls GSM*h5 | \
awk '/Donor/ && /filtered/{gsub(/_filtered_gene_bc_matrices_h5.h5/,"",$1);print}' > filtered.list
