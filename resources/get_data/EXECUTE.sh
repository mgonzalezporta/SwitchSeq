#!/bin/bash

species=mmu
ensembl_v=72
outdir=/homes/mar/public_html/tools/lorem/data_for_download
path=http://jun2013.archive.ensembl.org/biomart/martservice?
./get_data.sh $species $ensembl_v $outdir $path
