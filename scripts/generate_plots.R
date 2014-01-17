args=commandArgs(TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

source(paste(bin, "/scripts/PlotRNASeq/PlotRNASeq.R", sep=""), chdir=T)

## load data
rpkms=readExpressionData(gId=gId, infile=expdata, cond1=cond1, cond2=cond2)
biotypes=readBiotypeData(gId=gId, infile=annot)
if (filt != "NA") {
	significant_events=readSignificantEvents(gId=gId, infile=filt)
} else {
	significant_events="NA"
}

## create TranscriptExpressionSet object
tes=newTranscriptExpressionSet(
        gId=gId,
        rpkms=rpkms,
        biotypes=biotypes,
        cond1=cond1,
        cond2=cond2,
        significant_events=significant_events
)

## generate plots
outfile=getOutfile(gId=gId, plot_type="starplots", outdir=outdir)
plotStars(tes=tes, outfile=outfile)

outfile=getOutfile(gId=gId, plot_type="distrplots", outdir=outdir)
plotDistr(tes=tes, outfile=outfile)
