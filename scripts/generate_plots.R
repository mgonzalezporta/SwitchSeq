args=commandArgs(TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

library("tviz")

## functions
getOutfile=function(outdir=outdir, plot_type=plot_type, gId=gId) {
  # create plots directory
  outdir_plots=paste(outdir, "/data/plots", sep="")
  dir.create(outdir_plots, showWarnings = FALSE)

  # create plot type directory
  outdir_plots_type=paste(outdir_plots, plot_type, sep="/")
  dir.create(outdir_plots_type, showWarnings = FALSE)

  # create gId subdirectory
  s=unlist(strsplit(gId, ''))
  id=paste(s[c(1:(length(s)-3))], collapse='')
  outdir_plots_type_id=paste(outdir_plots_type, id, sep="/")
  dir.create(outdir_plots_type_id, showWarnings = FALSE)

  # outfile
  outfile=paste(outdir_plots_type_id, "/", gId, ".pdf", sep="")

  return(outfile)
}

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
