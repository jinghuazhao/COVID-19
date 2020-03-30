# 30-3-2020 JHZ

ls GSM*h5 | \
awk '/Donor/ && /filtered/{gsub(/_filtered_gene_bc_matrices_h5.h5/,"",$1);print}' > Donor.list

ls GSM*h5 | \
awk '/HP/ && /filtered/{gsub(/_filtered_gene_bc_matrices_h5.h5/,"",$1);print}' > HP.list

ls GSM*h5 | \
awk '/IPF/ && /filtered/{gsub(/_filtered_gene_bc_matrices_h5.h5/,"",$1);print}' > IPF.list

ls GSM*h5 | \
awk '/SSc/ && /filtered/{gsub(/_filtered_gene_bc_matrices_h5.h5/,"",$1);print}' > SSc.list

ls GSM*h5 | \
awk '/Cryobiopsy/ && /filtered/{gsub(/_filtered_gene_bc_matrices_h5.h5/,"",$1);print}' > Cryobiopsy.list

ls GSM*h5 | \
awk '/Myositis/ && /filtered/{gsub(/_filtered_gene_bc_matrices_h5.h5/,"",$1);print}' > Myositis.list
