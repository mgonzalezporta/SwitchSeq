args=commandArgs(TRUE)
# input_file=as.character(args[1])
# out_dir=as.character(args[2])
# data_dir=as.character(args[3])
# species=as.character(args[4])
# ensembl_v=as.character(args[5])
# cond1=as.character(args[6])
# cond2=as.character(args[7])
# gId=as.character(args[8])

input_file="/Users/mar/Desktop/app/2013.05.30_45pairs.dexseq2.txt"
outdir="/Users/mar/Desktop/app/html_12/"
data_dir="/Users/mar/Desktop/app/data_10/"
species="hsa"
ensembl_v=66
cond1="3-45"
cond2="49-89"
ensembl_annot=paste(data_dir, "/", species, "/_ensembl", ensembl_v, ".annot_coding.1.txt", sep="")
gId="ENSG00000180104"

source("./PlotRNASeq/PlotRNASeq.R")

## create TranscriptExpressionSet object
tes=readTranscriptExpressionSet(gId=gId, infile=input_file, cond1=cond1, cond2=cond2)
tes=annotateTranscriptExpressionSet(gId=gId, infile=ensembl_annot, tes=tes)

## generate plots
outfile=getOutfile(gId=gId, plot_type="starplots", outdir=outdir)
plotStars(tes=tes, outfile="~/Desktop/test.stars.pdf")

outfile=getOutfile(gId=gId, plot_type="distrplots", outdir=outdir)
plotDistr(tes=tes, outfile="~/Desktop/test.distr.pdf")
