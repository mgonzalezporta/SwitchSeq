readExpressionData=function(gId=gId, infile=infile, cond1=cond1, cond2=cond2) {
  # header
  command=paste("head -n1", infile, sep=" ")
  header=read.table(pipe(command), as.is=T)

  # conditions
  cond1=.format_cond(cond1)
  cond2=.format_cond(cond2)

  # transcript RPKMs
  command=paste("grep", gId, infile, sep=" ")
  data=read.table(pipe(command), col.names=header,
    colClasses=c("NULL", "character", rep("numeric", length(header)-2)))
  rpkms=as.matrix(data[,c(cond1, cond2)])
  rownames(rpkms)=data[,1]

  return(rpkms)
}

readBiotypeData=function(gId=gId, infile=infile) {
  # biotype info
  command=paste("grep", gId, infile, sep=" ")
  tBiotype=read.table(pipe(command))[,4:5]
  colnames(tBiotype)=c("tId", "tBiotype")

  return(tBiotype)
}

readSignificantEvents=function(gId=gId, infile=significant_events) {
  # significant events, if any
  command=paste("grep", gId, infile, sep=" ")
  significant_events=as.character(read.table(pipe(command))[,2])

  return(significant_events)
}

newTranscriptExpressionSet=function(gId=gId, rpkms=rpkms, biotypes=biotypes, cond1=cond1, cond2=cond2, significant_events=significant_events) {
  # conditions
  cond1=.format_cond(cond1)
  cond2=.format_cond(cond2)

  cond1=cond1-1
  cond2=seq(from=length(cond1)+1, to=dim(rpkms)[2], by=1)
  conditions=list(cond1=cond1, cond2=cond2)

  # create TranscriptExpressionSet object
  tes=new("TranscriptExpressionSet",
    id=gId,
    conditions=conditions,
    rpkms=rpkms,
    biotypes=biotypes,
    significant_events=significant_events
  )

  return(tes)
}

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

## internal functions
.format_cond=function(cond) {
  x=as.numeric(strsplit(cond, "-")[[1]])
  cond=(x[1]-1):(x[2]-1)
  return(cond)
}
