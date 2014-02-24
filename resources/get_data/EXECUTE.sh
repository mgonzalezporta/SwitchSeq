#!/bin/bash

# HSA.70
species=hsa
ensembl_v=70
outdir=/homes/mar/public_html/tools/switchseq/data_for_download
path=http://jan2013.archive.ensembl.org/biomart/martservice?
c="./get_data.sh $species $ensembl_v $outdir $path"
bsub $c