#!/usr/bin/bash

git add .gitignore
git commit -m ".gitignore"
git add README.md
git commit -m "README"
git add doc
git commit -m "Documents and auxiliary files"
git add -f cvd2/README.txt cvd2/*sh cvd2/ps.R
git commit -m "CVD2"
git add -f ACE2/2805-6_2.metal ACE2/select.sh ACE2/ukb-ACE2-annotate.sh ACE2/ace2_interval_imputed_ukb_imputed.sh
git commit -m "SomaLogic"
git add -f GEO/GEO.R GEO/GEO2R.R GEO/test.ini GEO/test.R GEO/test.sh GEO/README.md
git commit -m "GEO test"
git add -f HGDP/README.md HGDP/4MDS*analysis*v1.3.pdf
git add -f HGDP/EPIC-Omics.sh HGDP/HGI.* HGDP/refdata.* HGDP/MDSplot.R HGDP/keep.* 
git add -f HGDP/nonEUR.sh HGDP/nonEUR.R HGDP/nonEUR.sb HGDP/nonEUR.C1-C2.png HGDP/INTERVAL.C1-C2.png HGDP/hgdp_tgp_interval.PC1-2.pdf
git commit -m "1000G mappings"
git add -f HGI/README.md HGI/INTERVAL.do HGI/INTERVAL.rec HGI/INTERVAL.sh HGI/autosomes.sb HGI/X.sb HGI/X-female.sb
git add -f HGI/glist-hg19.sb HGI/glist-hg19.sh HGI/L10RB_IFNAR2_variants.sh
git add -f HGI/20200731.sh HGI/20201201.sh HGI/20210317.sh HGI/bgen.sb HGI/bgen-X.sb
git add -f HGI/pca_project.* HGI/update_bgi.py HGI/update_bgi.py HGI/ethnic.do
git commit -m "HGI analysis"
git push
