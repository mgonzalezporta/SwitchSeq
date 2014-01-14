#!/bin/bash

species=mmu
ensembl_v=72
outdir=/homes/mar/home_microarray/data/babraham_2013/lorem
path=http://jun2013.archive.ensembl.org/biomart/martservice?
./get_data.sh $species $ensembl_v $outdir $path
