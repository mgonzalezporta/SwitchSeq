## initialize
setMethod("initialize", "TranscriptExpressionSet",
  function(.Object, id=id, rpkms=rpkms, biotypes=biotypes, conditions=conditions, significant_events=significant_events){
    .Object@id=id
    .Object@rpkms=rpkms
    .Object@conditions=conditions
    .Object@gexp=colSums(.Object@rpkms)
    .Object@dominance=.calculate_dominance(.Object)
    .Object@.relexp=.calculate_relexp(.Object)
    .Object@.scaledexp=.calculate_scaledexp(.Object)
    .Object@biotypes=biotypes
    .Object@.cols=.get_transcript_colors(.Object)
    .Object@significant_events=significant_events
    return(.Object)
  })

## show + accessors
setMethod("show",
  "TranscriptExpressionSet",
  function(object) {
    cat("Object of class",class(object),"\n")
    cat("   Id:",id(object),"\n")
    cat("   RPKMs: data frame with", dim(object@rpkms)[1], "transcripts and", dim(object@rpkms)[2], "samples\n")
    cat("     transcripts:", head(rownames(object@rpkms), n=4), "(...)\n")
    cat("     samples:", head(colnames(object@rpkms), n=4), "(...)\n")
    cat("   Conditions:\n")
    cat("     condition 1:", head(object@conditions$cond1, n=4), "(...)\n")
    cat("     condition 2:", head(object@conditions$cond2, n=4), "(...)\n")
    cat("   Gene expression:", head(object@gexp, n=4), "(...)\n")
    cat("   Dominance:", head(object@dominance, n=4), "(...)\n")
    cat("   Unique transcript biotypes:", names(table(tes@biotypes$tBiotype)), "\n")
    cat("   Significant events:", object@significant_events, "\n")
  })

setMethod("id", "TranscriptExpressionSet", function(object) object@id)
setReplaceMethod("id",
                 signature(object="TranscriptExpressionSet",
                           value="character"),
                 function(object, value) {
                   object@id=value
                   return(object)
                 })

setMethod("conditions", "TranscriptExpressionSet", 
  function(object) object@conditions)

setMethod("rpkms", "TranscriptExpressionSet", 
  function(object) object@rpkms)

setMethod("gexp", "TranscriptExpressionSet", 
  function(object) object@gexp)

setMethod("dominance", "TranscriptExpressionSet", 
  function(object) object@dominance)

setMethod("biotypes", "TranscriptExpressionSet", 
  function(object) object@biotypes)

setMethod("significant_events", "TranscriptExpressionSet",
  function(object) object@significant_events)

## plots
setMethod("plotStars",
  signature(tes="TranscriptExpressionSet", 
            outfile="character"),
  function(tes, outfile) {
    # include tBiotype info + reorder transcripts
    data=.annotate_data(data=tes@.scaledexp, significant_events=tes@significant_events, biotypes=tes@biotypes)

    # colors
    col=unlist(tes@.cols)[seq(1, length(unlist(tes@.cols)), by=2)]

    # legend title
    legend_title=.get_legend_title_starplot(tes)

    # plot height
    nrows_layout=max(sapply(tes@conditions, length))
    if (nrows_layout<6) { 
      h=3 
    } else {
      h=nrows_layout*.5
    }

    # subplot location
    loc=matrix(ncol=2)
    for (i in seq(nrows_layout-1, by=-1, length.out=length(tes@conditions$cond1))) {
      loc=rbind(loc, c(0, i))
    }
    for (i in seq(nrows_layout-1, by=-1, length.out=length(tes@conditions$cond2))) {
      loc=rbind(loc, c(1, i))
    }
    loc=loc[complete.cases(loc),]

    # plot
    pdf(outfile, width=4, height=h)
      layout(t(1:2), widths=c(1,1.3))
      par(fg="white", col="gray35",
        mar=c(0, 0, 0, 0), oma=c(0, 0, 0, 0))
        
      stars(t(data), 
        draw.segments=T,
        col.segments=col,
        cex.main=0.6,
        cex=0.3,
        ncol=2,
        len=0.4,
        scale=F,
        locations=loc
      )

      # legend
      plot.new()
      par(oma=c(0,0,5,0), mar=c(0,0,0,0), xpd=NA)
      legend("center",
        title=legend_title, 
        title.adj=0,
        adj=0,
        legend=rownames(data), 
        fill=col, 
        border=col,
        cex=0.4
      )
    garbage <- dev.off()  # supress output
})

setMethod("plotDistr",
  signature(tes="TranscriptExpressionSet", 
            outfile="character"),
  function(tes, outfile) {
    if(length(tes@conditions$cond1)>=10 & length(tes@conditions$cond2)>=10) {
      .plot_boxplots(tes=tes, outfile=outfile)
    } else {
      .plot_segments(tes=tes, outfile=outfile)
    }
})

## internal methods
setMethod(".annotate_data", 
  signature(data="matrix",
	    significant_events="character",
            biotypes="data.frame"),
  function(data=data, significant_events=significant_events, biotypes=biotypes) {
    data=data[order(rownames(data)),]
  
    biotypes=biotypes[order(biotypes$tId),]

    biotypes$tId=as.character(biotypes$tId)
    se=biotypes$tId %in% significant_events
    if (sum(se) > 0) {
	biotypes[se,]$tId=paste(biotypes[se,]$tId, "**", sep=" ")
    }

    rn=paste(biotypes$tBiotype, " - ", biotypes$tId, sep="")
    rownames(data)=rn

    data=data[order(rownames(data)),]

    return(data)
}) 

setMethod(".plot_boxplots",
  signature(tes="TranscriptExpressionSet", 
            outfile="character"),
  function(tes, outfile) {

    # number of transcripts
    nOfT=dim(tes@.relexp)[1]

    # colors
    col=c(
      rep(list(c("white", "gray90")), 2),
      tes@.cols
    )

    # margins
    mar=c(
      list(c(0, 4, 0, 3)), 
      list(c(0, 2, 0, 5)), 
      rep(list(c(0, 0, 0, 0)), nOfT))

    # width
    if (nOfT>3) {
      #width=(4+nOfT-1)
      width=4+nOfT
    } else {
      width=7
    }

    # plot
    plot_count=0
    pdf(file=outfile, width=width, height=6)
      layout(t(1:(nOfT+3)), widths=c(0.5, rep(1, 2), rep(0.5, nOfT)))
      par(oma=c(13,0,2,0))

      # legend
      plot.new()
      legend="legend:\nlight colors (left) - condition 1\ndark colors (right) - condition 2"
      mtext(legend, side=1, las=2, line=5, cex=0.75, adj=0, padj=0.15)

      # gexp
      plot_count=plot_count+1
      .subplot_boxplots(
        ldata=list(tes@gexp[tes@conditions$cond1], tes@gexp[tes@conditions$cond2]),
        xlab=tes@id, 
        ylab="gene expression (RPKMs)", 
        type="gexp", 
        mar=mar, 
        col=col, 
        plot_count=plot_count
        )

      # dominance
      plot_count=plot_count+1
      .subplot_boxplots(
        ldata=list(tes@dominance[tes@conditions$cond1], tes@dominance[tes@conditions$cond2]),
        xlab=tes@id, 
        ylab="major transcript dominance", 
        type="dom", 
        mar=mar, 
        col=col, 
        plot_count=plot_count
        )
      
      # transcript relative abundances
      # annotate + reorder data first
      data=.annotate_data(data=tes@.relexp, significant_events=tes@significant_events, biotypes=tes@biotypes)

      for(i in 1:nOfT){
        plot_count=plot_count+1

	xlab=gsub(" - ", "\n", rownames(data)[i])
	if (sum("NA" %in% significant_events) == 0) {
		xlab=.expand_xlab(xlab)
	}

        .subplot_boxplots(
          ldata=list(data[i, tes@conditions$cond1], data[i, tes@conditions$cond2]),
          xlab=xlab,
          ylab="transcript relative abundance", 
          type="trel", 
          mar=mar, 
          col=col, 
          plot_count=plot_count
          )
      }

    garbage<-dev.off()  # supress output
})

setMethod(".subplot_boxplots", 
  signature(
    ldata="list",
    xlab="character",
    ylab="character",
    type="character",
    mar="list",
    col="list",
    plot_count="numeric"
    ), 
  function(ldata=ldata, xlab=xlab, ylab=ylab, type=type, 
    mar=mar, col=col, plot_count=plot_count) {

    # plot args
    par(bty="n", mar=mar[[plot_count]])

    # set ylim + yat
    ylim=c(0,1)
    yat=seq(0,1, length.out=6)
    if (type == "gexp") {
      ylim=round(c( min(sapply(ldata, min)), max(sapply(ldata, max)) ), -1)
      yat=round(seq(from=ylim[1], to=ylim[2], length.out=5), -1)
    }

    # boxplot
    boxplot(ldata, ylim=ylim, outline=T, pch=20, horizontal=F, xaxt="n", yaxt="n", col=col[[plot_count]])

    # axes and axes labels
    if (plot_count %in% 1:3){
        axis(2, las=0, at=yat)
        mtext(ylab, side=2, line=2.5, cex=0.75)
      }
    axis(1, labels=F, lwd.ticks=0, lwd=1, at=c(1,2))
    mtext(xlab, side=1, line=0.75, las=2, cex=0.75)
})

setMethod(".plot_segments",
  signature(tes="TranscriptExpressionSet", 
            outfile="character"),
  function(tes, outfile) {

    # number of transcripts
    nOfT=dim(tes@.relexp)[1]

    # colors
    col=c(
      rep(list(c("gray60", "gray30")), 2),
      tes@.cols
      )

    # margins
    mar=c(
      list(c(0, 4, 0, 3)), 
      list(c(0, 2, 0, 5)), 
      rep(list(c(0, 0, 0, 0)), nOfT))

    # width
    if (nOfT>3) {
      #width=(4+nOfT-1)
      width=4+nOfT
    } else {
      width=7
    }

    # plot
    plot_count=0
    pdf(file=outfile, width=width, height=6)
      layout(t(1:(nOfT+3)), widths=c(0.5, rep(1, 2), rep(0.5, nOfT)))
      par(oma=c(13,0,2,0))

      # legend
      plot.new()
      legend="legend:\nlight colors (left) - condition 1\ndark colors (right) - condition 2"
      mtext(legend, side=1, las=2, line=5, cex=0.75, adj=0, padj=0.15)

      # gexp
      plot_count=plot_count+1
      .subplot_segments(
        data1=tes@gexp[tes@conditions$cond1], 
        data2=tes@gexp[tes@conditions$cond2], 
        xlab=tes@id, 
        ylab="gene expression (RPKMs)", 
        type="gexp", 
        mar=mar, 
        col=col, 
        plot_count=plot_count
        )

      # dominance
      plot_count=plot_count+1
      .subplot_segments(
        data1=tes@dominance[tes@conditions$cond1], 
        data2=tes@dominance[tes@conditions$cond2], 
        xlab=tes@id, 
        ylab="major transcript dominance", 
        type="dom", 
        mar=mar, 
        col=col, 
        plot_count=plot_count
        )
      
      # transcript relative abundances
      # annotate + reorder data first
      data=.annotate_data(data=tes@.relexp, significant_events=tes@significant_events, biotypes=tes@biotypes)
      
      for (i in 1:nOfT) {
        plot_count=plot_count+1
        
        xlab=gsub(" - ", "\n", rownames(data)[i])
        if (sum("NA" %in% significant_events) == 0) {
                xlab=.expand_xlab(xlab)
        }

        .subplot_segments(
          data1=data[i, tes@conditions$cond1], 
          data2=data[i, tes@conditions$cond2], 
          xlab=xlab,
          ylab="transcript relative abundance",
          type="trel",
          mar=mar, 
          col=col, 
          plot_count=plot_count
          )
      }
    garbage<-dev.off()  # supress output
})

setMethod(".subplot_segments", 
  signature(
    data1="numeric",
    data2="numeric",
    xlab="character",
    ylab="character",
    type="character",
    mar="list",
    col="list",
    plot_count="numeric"
    ), 
  function(data1=data1, data2=data2, xlab=xlab, ylab=ylab, type=type, 
    mar=mar, col=col, plot_count=plot_count) {
    
    # plot args
    par(bty="n", mar=mar[[plot_count]]) 

    # data points
    x=c(0.1,0.2)
    y=c(mean(data1), mean(data2))
    sd=c(sd(data1), sd(data2))
    m=cbind(x,y,sd)

    # plot parameters
    col=col[[plot_count]]
    lwd=2
    epsilon=0.02
    xlim=c(0,0.3)

    # set ylim
    ylim=c(0,1)
    if (type == "gexp") {
      min_y=min(m[, "y"])
      max_y=max(m[, "y"])
      row_index_min=which(m[, "y"]==min_y)
      if (row_index_min==1) {
        row_index_max=2
      } else {
        row_index_max=1
      }
      ylim=c( min_y-m[row_index_min, "sd"]-1, max_y+m[row_index_max, "sd"]+1 )
    }
    
    # plot
    plot(x, y, xlim=xlim, ylim=ylim, col=col, xaxt="n", yaxt="n", ylab="", xlab="", lwd=lwd)
    segments(x, y-sd, x, y+sd, col=col, lwd=lwd)
    segments(x-epsilon,y-sd,x+epsilon,y-sd, col=col, lwd=lwd)
    segments(x-epsilon,y+sd,x+epsilon,y+sd, col=col, lwd=lwd)

    # axes and axes labels
    if (plot_count %in% 1:3){
        axis(2, las=0, pos=0)
        mtext(ylab, side=2, line=2.5, cex=0.75)
      }
    axis(1, labels=F, lwd.ticks=0, lwd=1, at=c(0.1,0.2))
    mtext(xlab, side=1, line=0.75, las=2, cex=0.75) 
})

setMethod(".distr_summary", "numeric", 
  function(object) {
    x=object
    x_m=formatC(round(mean(x), 2), 2, format="f")
    x_sd=formatC(round(sd(x), 2), 2, format="f")

    r=paste("mean=", x_m, ", sd=", x_sd, sep="")
    return(r)   
  })

setMethod(".get_legend_title_starplot", "TranscriptExpressionSet",
  function(object) {
    legend_title=paste(
      object@id, "\n\n",
      "gene expression summary (RPKMs):\n", 
      "   condition 1 (left) - ", .distr_summary(object@gexp[object@conditions$cond1]), "\n",
      "   condition 2 (right) - ", .distr_summary(object@gexp[object@conditions$cond2]), "\n\n",
      "dominance summary:\n",
      "   condition 1 (left) - ", .distr_summary(object@dominance[object@conditions$cond1]), "\n",
      "   condition 2 (right) - ", .distr_summary(object@dominance[object@conditions$cond2]), "\n\n",
      "legend:", sep="")

    return(legend_title)
  })

setMethod(".calculate_dominance", "TranscriptExpressionSet",
    function(object) {
          dominance=apply(object@rpkms, 2, .ratio_second)
          return(dominance)
    })

setMethod(".calculate_relexp", "TranscriptExpressionSet",
    function(object) {
      relexp=t(t(object@rpkms)/object@gexp)
      return(relexp)
    })

setMethod(".ratio_second", "numeric", 
  function(object) {
    x=object
    x=x[order(x, decreasing=T)]
    r=x[2]/x[1]

    return(r)   
  })

setMethod(".calculate_scaledexp", "TranscriptExpressionSet",
    function(object) {
      x=t( t(object@.relexp) * abs( sqrt( object@gexp/max(object@gexp ))))
      x[is.na(x)]=0
      r=x/max(x)

      return(r)
    })


setMethod(".get_transcript_colors", "TranscriptExpressionSet",
    function(object) {
      # define color intervals
      col_intervals=data.frame(
         tBiotype=c("nonsense_mediated_decay", "processed_transcript",
            "protein_coding", "retained_intron"),
         start=c(0.05, 0.30, 0.60, 0.85),
         end=c(0.15, 0.40, 0.70, 0.95)
      )

      # biotype frequencies
      freq=data.frame(table(object@biotypes$tBiotype))
      colnames(freq)=c("tBiotype", "n")

      # combine
      m=merge(freq, col_intervals, by="tBiotype")

      # obtain colors
      cols=apply(m, 1, function(x) {
        n=as.numeric(x["n"])
        start=as.numeric(x["start"])
        end=as.numeric(x["end"])

        col=c(
          rainbow( n, s=0.7, v=0.9, alpha=0.7, start=start, end=end),
          rainbow( n, s=0.7, v=0.5, alpha=0.7, start=start, end=end)
        )
        
        names(col)=rep(paste(x["tBiotype"], 1:n, sep="_"), 2)
        col=split(col, names(col))
        
        return(col)
        }
      )
      cols=unlist(cols, recursive=F)

      return(cols)
    })

setMethod(".expand_xlab", "character",
    function(object) {
	x=object
        xlab=gsub(" **", "\n**", x, fixed=T)

	return(xlab)
    })
