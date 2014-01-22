#!/bin/bash

## arguments
species=$1
ensembl_v=$2
base_dir=$3
path=$4

## dir structure
outdir=$base_dir/$species.$ensembl_v
if [ ! -e $outdir ]; then mkdir $outdir; fi

## ensembl1
echo "Creating ensembl1.txt..."
in_xml=ensembl1.xml
out=$outdir/ensembl1.txt

echo "gene_id gene_name number_of_transcripts transcript_id transcript_biotype" > $out
perl query_biomart.pl $in_xml $species $path >> $out

## ensembl2
echo "Creating ensembl2.txt..."
in_xml=ensembl2.xml
out=$outdir/ensembl2.txt

echo "gene_id uniprot_id" > $out
perl query_biomart.pl $in_xml $species $path | grep -E '^[[:alnum:]]+\s[[:alnum:]]+' >> $out

## appris
echo "Retrieving appris_data.principal.txt..."
case $species in
	hsa)
		species_full=homo_sapiens
		;;
        mmu)
                species_full=mus_musculus
		;;
        rno)
                species_full=rattus_norvegicus
		;;
        dre)
                species_full=danio_rerio
		;;
esac

echo "s: $species sf: $species_full"

appris=http://appris.bioinfo.cnio.es/download/data/$species_full/appris_data.principal.txt
wget $appris 
mv appris_data.principal.txt $outdir

## protein sequences
echo "Retrieving protein sequences..."

if [ ! -e /homes/mar/system/ensembl.$ensembl_v ]
then
	url=http://www.ensembl.org/cvsdownloads/ensembl-$ensembl_v.tar.gz
	wget $url
	tar xzf ensembl-$ensembl_v.tar.gz
	rm ensembl-$ensembl_v.tar.gz
	mv ensembl /homes/mar/system/ensembl.$ensembl_v
fi

perl get_prot_seq.pl $ensembl_v $species $outdir
cd $outdir && tar czf prot_seq.tar.gz prot_seq/
