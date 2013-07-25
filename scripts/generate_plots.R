args=commandArgs(TRUE)
gId=as.character(args[1])
exp_data=as.character(args[2])
annot=as.character(args[3])
cond1=as.character(args[4])
cond2=as.character(args[5])
outdir=as.character(args[6])

# gId="ENSG00000180104"
# expdata="/Users/mar/Desktop/app/2013.05.30_45pairs.dexseq2.txt"
# annot="/Users/mar/Desktop/app/data_10/hsa/_ensembl66.annot_coding.1.txt"
# cond1="3-45"
# cond2="49-89"
# outdir="/Users/mar/Desktop/app/html_12/"

source("./PlotRNASeq/PlotRNASeq.R")

## load data
rpkms=readExpressionData(gId=gId, infile=expdata, cond1=cond1, cond2=cond2)
biotypes=readBiotypeData(gId=gId, infile=annot)

## create TranscriptExpressionSet object
tes=newTranscriptExpressionSet(
	gId=gId, 
	rpkms=rpkms, 
	biotypes=biotypes,
	cond1=cond1, 
	cond2=cond2
)

## generate plots
outfile=getOutfile(gId=gId, plot_type="starplots", outdir=outdir)
plotStars(tes=tes, outfile="~/Desktop/test.stars.pdf")

outfile=getOutfile(gId=gId, plot_type="distrplots", outdir=outdir)
plotDistr(tes=tes, outfile="~/Desktop/test.distr.pdf")
