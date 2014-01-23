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
tmp=$out.tmp

echo "gene_id gene_name number_of_transcripts transcript_id transcript_biotype" > $out.tmp
perl query_biomart.pl $in_xml $species $path >> $out.tmp
# fix missing values
cat $out.tmp | sed 's/\t\t/\tNA\t/g' > $out
rm $out.tmp

## ensembl2
echo "Creating ensembl2.txt..."
in_xml=ensembl2.xml
out=$outdir/ensembl2.txt

echo "gene_id uniprot_id" > $out
perl query_biomart.pl $in_xml $species $path | grep -E '^[[:alnum:]]+\s[[:alnum:]]+' >> $out

## appris
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

if [ ! -z "$species_full" ]
then
	echo "Retrieving appris_data.principal.txt..."
	appris=http://appris.bioinfo.cnio.es/download/data/$species_full/appris_data.principal.txt
	wget -P $outdir $appris 
fi

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
