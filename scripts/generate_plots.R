args=commandArgs(TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

# gId="ENSG00000124193"
# expdata="/Users/mar/github/lorem/test/test_data.txt"
# annot="/Users/mar/Desktop/app/data_10/hsa/_ensembl66.annot_coding.1.txt"
# cond1="3-6"
# cond2="7-10"
# outdir="./test"

source(paste(bin, "/scripts/PlotRNASeq/PlotRNASeq.R", sep=""))

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
plotStars(tes=tes, outfile=outfile)

outfile=getOutfile(gId=gId, plot_type="distrplots", outdir=outdir)
plotDistr(tes=tes, outfile=outfile)
