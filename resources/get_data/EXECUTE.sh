#!/bin/bash

# HSA.70
species=hsa
ensembl_v=73
outdir=/homes/mar/public_html/tools/switchseq/data_for_download
biomart_path=http://sep2013.archive.ensembl.org/biomart/martservice?
ensembl_api_path=/homes/mar/system/ensembl
bioperl=/homes/mar/system/bioperl-live
c="./get_data.sh $species $ensembl_v $outdir $biomart_path $ensembl_api_path $bioperl"
bsub $c
