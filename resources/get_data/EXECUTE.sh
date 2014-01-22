#!/bin/bash

# MMU.74
species=mmu
ensembl_v=74
outdir=/homes/mar/public_html/tools/lorem/data_for_download
path=http://www.ensembl.org/biomart/martservice?
c="./get_data.sh $species $ensembl_v $outdir $path"
bsub $c

# HSA.74
species=hsa
ensembl_v=74
outdir=/homes/mar/public_html/tools/lorem/data_for_download
path=http://www.ensembl.org/biomart/martservice?
c="./get_data.sh $species $ensembl_v $outdir $path"
bsub $c

# RNO.70
species=rno
ensembl_v=70
outdir=/homes/mar/public_html/tools/lorem/data_for_download
path=http://jan2013.archive.ensembl.org/biomart/martservice?
c="./get_data.sh $species $ensembl_v $outdir $path"
bsub $c

# DRE.70
species=dre
ensembl_v=70
outdir=/homes/mar/public_html/tools/lorem/data_for_download
path=http://jan2013.archive.ensembl.org/biomart/martservice?
c="./get_data.sh $species $ensembl_v $outdir $path"
bsub $c
