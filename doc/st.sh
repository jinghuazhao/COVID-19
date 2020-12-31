#!/usr/bin/bash

git add .gitignore
git commit -m ".gitignore"
git add README.md
git commit -m "README"
git add doc
git commit -m "Documents and auxiliary files"
git add cvd2/README.txt cvd2/*sh cvd2/ps.R
git commit -m "CVD2"
git add ACE2/2805-6_2.metal ACE2/select.sh ACE2/ukb-ACE2-annotate.sh ACE2/ace2_interval_imputed_ukb_imputed.sh
git commit -m "SomaLogic"
git add GEO/GEO.R GEO/GEO2R.R GEO/test.ini GEO/test.R GEO/test.sh GEO/README.md
git commit -m "GEO test"
git add HGI/README.md HGI/INTERVAL.do HGI/INTERVAL.rec HGI/INTERVAL.sh HGI/autosomes.sb HGI/X.sb HGI/X-female.sb
git add HGI/glist-hg19.sb HGI/glist-hg19.sh HGI/L10RB_IFNAR2_variants.sh
git add HGI/20200731.sh HGI/20201201.sh HGI/bgen.sb HGI/bgen-X.sb
git commit -m "HGI analysis"
git add NGS/README.* NGS/ngs.png NGS/ngs.do NGS/ngs.R NGS/ngs.sb NGS/ngs.sh NGS/lod.sh NGS/pQTLtools.R NGS/utils
git commit -m "Olink/NGS pilot"
git add SEQ/INTERVAL.md SEQ/README.md SEQ/doc SEQ/weswgs.do SEQ/weswgs.R SEQ/weswgs.sh SEQ/wgs.wrap
git add SEQ/bgen.sb SEQ/spa.* SEQ/grm.wrap SEQ/rva.*
git commit -m "SCALLOP-seq meta-analysis"
git push
