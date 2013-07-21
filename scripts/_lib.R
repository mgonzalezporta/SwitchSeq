new_switch_event=function(gId, input_file, ensembl_annot) {
	# load data
	command=paste("head -n1", input_file, sep=" ")
	header=read.table(pipe(command), as.is=T)

	command=paste("grep", gId, input_file, sep=" ")
	data=read.table(pipe(command), col.names=header,
		colClasses=c("NULL", "character", rep("numeric", length(header)-2)))

	# load transcript biotype info
	command=paste("grep", gId, ensembl_annot, sep=" ")
	tBiotype=read.table(pipe(command))[,4:5]
	colnames(tBiotype)=c("tId", "tBiotype")

	# format data
	data=merge(tBiotype, data, by="tId")
	header=paste(data$tBiotype, data$tId, sep=" - ")
	data=t(data[,-c(1:2)])
	colnames(data)=header
	data=data[,order(colnames(data))]

	return(data)
}

# scale_se=function(x, gexp) {
# 	for (i in 1:dim(x)[2]) {
# 		x[,i]=x[,i]*abs(sqrt(gexp/max(gexp)))
# 	}
# 	x[is.na(x)]=0
# 	x=x/max(x)

# 	return(x)
# }

plot_stars=function(gId, data, gexp, cond1, cond2, outfile) {

	# gexp info
	gexp_sum_cond1=summarise_gexp(gexp, cond1, "   condition 1 (left)")
	gexp_sum_cond2=summarise_gexp(gexp, cond2, "   condition 2 (right)")

	# dominance info
	dom_sum_cond1=summarise_dominance(data, cond1, "   condition 1 (left)")
	dom_sum_cond2=summarise_dominance(data, cond2, "   condition 2 (right)")

	# legend title
	legend_title=paste(gId, "\n\n",
		"gene expression summary (RPKMs):\n", 
		gexp_sum_cond1, "\n", 
		gexp_sum_cond2, "\n\n",
		"dominance summary:\n",
		dom_sum_cond1, "\n",
		dom_sum_cond2, "\n\n", 
		"legend:", sep="")

	# colors
	col=get_transcript_colors(data)
	col=unlist(col)[seq(1, length(unlist(col)), by=2)]

	# data normalisation
	data=data[c(cond1, cond2),]	## tmp
	gexp=gexp[c(cond1, cond2)]	## tmp
	data_scaled=scale_se(data, gexp)

	# plot height
	nrows_layout=max(length(cond1), length(cond2))
	if (nrows_layout<6) { 
		h=3 
	} else {
		h=nrows_layout*.5
	}

	# subplot location
	loc=matrix(ncol=2)
	for (i in seq(nrows_layout-1, by=-1, length.out=length(cond1))) {
		loc=rbind(loc, c(0, i))
	}
	for (i in seq(nrows_layout-1, by=-1, length.out=length(cond2))) {
		loc=rbind(loc, c(1, i))
	}
	loc=loc[complete.cases(loc),]

	# plot
	pdf(outfile, width=4, height=h)
		layout(t(1:2), widths=c(1,1.3))
		par(fg="white", col="gray35",
			mar=c(0, 0, 0, 0), oma=c(0, 0, 0, 0))
			
		stars(data_scaled, 
		    draw.segments = T,
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
			legend=colnames(data_scaled), 
			fill=col, 
			border=col,
			cex=0.4
		)
	garbage <- dev.off()	# supress output
}

plot_distrplot=function(gId, x, gexp, cond1, cond2, outfile) {
	if(length(cond1)>=10 & length(cond2)>=10) {
		plot_boxplots(gId, x, gexp, cond1, cond2, outfile)
	} else {
		plot_segments(gId, x, gexp, cond1, cond2, outfile)
	}
}

plot_segments=function(gId, data, gexp, cond1, cond2, outfile) {

	nOfT=dim(data)[2]
	plot_count=0

	col=c(
		rep(list(c("gray60", "gray30")), 2)
		)
	col=c(col, get_transcript_colors(data))

	mar=c(
		list(c(0, 4, 0, 3)), 
		list(c(0, 2, 0, 5)), 
		rep(list(c(0, 0, 0, 0)), nOfT))

	# plot
	pdf(file=outfile, width=(4+nOfT-1), height=6)
		layout(t(1:(nOfT+3)), widths=c(0.5, rep(1, 2), rep(0.5, nOfT)))
		par(oma=c(13,0,2,0))

		# legend
		plot.new()
		legend="legend:\nlight colors (left) - condition 1\ndark colors (right) - condition 2"
		#mtext(legend, side=1, las=2, line=-0.75, cex=0.75, adj=0)
		mtext(legend, side=1, las=2, line=5, cex=0.75, adj=0, padj=0.15)

		# gexp
		data_cond1=gexp[cond1]
		data_cond2=gexp[cond2]
		xlab=gId
		ylab="gene expression (RPKMs)"
		type="gexp"
		plot_count=plot_count+1
		subplot_segments(data_cond1, data_cond2, xlab, ylab, type, mar, col, plot_count)

		# dominance
		data_cond1=calculate_dominance(data[cond1,])
		data_cond2=calculate_dominance(data[cond2,])
		xlab=gId
		ylab="major transcript dominance"
		type="dom"
		plot_count=plot_count+1
		subplot_segments(data_cond1, data_cond2, xlab, ylab, type, mar, col, plot_count)
		
		# transcript relative abundances
		for(i in 1:nOfT){
			data_cond1=data[cond1,i]
			data_cond2=data[cond2,i]
			xlab=gsub(" - ", "\n", colnames(data)[i])
			ylab="transcript relative abundance"
			type="trel"
			plot_count=plot_count+1
			subplot_segments(data_cond1, data_cond2, xlab, ylab, type, mar, col, plot_count)
		}
	garbage<-dev.off()	# supress output
}

subplot_segments=function(data_cond1, data_cond2, xlab, ylab, type, mar, col, plot_count) {
	# plot args
	par(bty="n", mar=mar[[plot_count]])	

	# data points
	x=c(0.1,0.2)
	y=c(mean(data_cond1), mean(data_cond2))
	sd=c(sd(data_cond1), sd(data_cond2))
	m=cbind(x,y,sd)

	# plot parameters
	col=col[[plot_count]]
	lwd=2
	epsilon = 0.02
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

}

plot_boxplots=function(gId, data, gexp, cond1, cond2, outfile) {	
	# general agrs
	nOfT=dim(data)[2]
	plot_count=0

	col=c(
		rep(list(c("white", "gray90")), 2)
		)
	col=c(col, get_transcript_colors(data))

	mar=c(
		list(c(0, 4, 0, 3)), 
		list(c(0, 2, 0, 5)), 
		rep(list(c(0, 0, 0, 0)), nOfT))

	# plot
	pdf(file=outfile, width=(4+nOfT-1), height=6)
		layout(t(1:(nOfT+3)), widths=c(0.5, rep(1, 2), rep(0.5, nOfT)))
		par(oma=c(13,0,2,0))

		# legend
		plot.new()
		legend="legend:\nlight colors (left) - condition 1\ndark colors (right) - condition 2"
		#mtext(legend, side=1, las=2, line=-0.75, cex=0.75, adj=0)
		mtext(legend, side=1, las=2, line=5, cex=0.75, adj=0, padj=0.15)

		# gexp
		data_cond1=gexp[cond1]
		data_cond2=gexp[cond2]
		ldata=list(data_cond1, data_cond2)
		xlab=gId
		ylab="gene expression (RPKMs)"
		type="gexp"
		plot_count=plot_count+1
		subplot_boxplots(ldata, xlab, ylab, type, mar, col, plot_count)

		# dominance
		data_cond1=calculate_dominance(data[cond1,])
		data_cond2=calculate_dominance(data[cond2,])
		ldata=list(data_cond1, data_cond2)
		xlab=gId
		ylab="major transcript dominance"
		type="dom"
		plot_count=plot_count+1
		subplot_boxplots(ldata, xlab, ylab, type, mar, col, plot_count)

		# transcript relative abundances
		for(i in 1:nOfT){
			data_cond1=data[cond1,i]
			data_cond2=data[cond2,i]
			ldata=list(data_cond1, data_cond2)
			xlab=gsub(" - ", "\n", colnames(data)[i])
			ylab="transcript relative abundance"
			type="trel"
			plot_count=plot_count+1
			subplot_boxplots(ldata, xlab, ylab, type, mar, col, plot_count)
		}

	garbage<-dev.off()	# supress output
}

subplot_boxplots=function(ldata, xlab, ylab, type, mar, col, plot_count) {

	# plot args
	par(bty="n", mar=mar[[plot_count]])

	# set ylim + yat
	ylim=c(0,1)
	yat=seq(0,1, length.out=6)
	if (type == "gexp") {
		ylim=round(c( min(sapply(ldata, min)), max(sapply(ldata, max)) ), -1)
		yat=round(seq(from=ylim[1], to=ylim[2], length.out=5), -1)
		#if (yat[1] == 0) { yat[1]=1 }
	}

	# boxplot
	boxplot(ldata, ylim=ylim, outline=T, pch=20, horizontal=F, xaxt="n", yaxt="n", col=col[[plot_count]])



	# axes and axes labels
	if (plot_count %in% 1:3){
    	axis(2, las=0, at=yat)
    	#axis(2, at=yat)
        mtext(ylab, side=2, line=2.5, cex=0.75)
    }
	axis(1, labels=F, lwd.ticks=0, lwd=1, at=c(1,2))
	mtext(xlab, side=1, line=0.75, las=2, cex=0.75)
}

get_transcript_colors=function(x) {

	nOfT=dim(x)[2]

	biotypes=unlist(strsplit(colnames(x), " - "))
	biotypes=data.frame(table(biotypes[seq(1, length(biotypes), by=2)]))
	colnames(biotypes)=c("tBiotype", "n")
	
	col_intervals=data.frame(t(matrix(c(0.05, 0.15, 0.3, 0.4, 0.6, 0.7, 0.85, 0.95), nrow=2)))
	colnames(col_intervals)=c("start", "end")
	rownames(col_intervals)=NULL
	col_intervals$tBiotype=c("nonsense_mediated_decay", "processed_transcript",
		"protein_coding", "retained_intron")
	
	m=merge(biotypes, col_intervals, by="tBiotype")
	
	#col=t(data.frame(NA, NA))
	#rownames(col)=c("N", "T")
	result=list()
	for (i in 1:dim(m)[1]) {
		biotype=as.character(m[i,]$tBiotype)
        n=m[i,]$n
        start=m[i,]$start
        end=m[i,]$end

        col=c(
			rainbow( n, s = 0.7, v = 0.9, alpha=0.7, start=start, end=end),
        	rainbow( n, s = 0.7, v = 0.5, alpha=0.7, start=start, end=end)
        )
        names(col)=rep(paste(biotype, 1:n, sep="_"), 2)
        result=c(result, split(col, names(col)))
	}
	return(result)

}

# get_ratio_second=function(x) {
# 	x=x[order(x, decreasing=T)]
# 	r=x[2]/x[1]
# 	return(r)
# }

# calculate_dominance=function(x) {
# 	r=apply(x, 1, get_ratio_second)
# 	return(r)
# }

# format_cond=function(cond) {
# 	x=as.numeric(strsplit(cond, "-")[[1]])
# 	cond=(x[1]-2):(x[2]-2)
# 	return(cond)
# }

# normalise_se=function(x, gexp) {
# 	for (i in 1:dim(x)[2]) {
# 	        x[,i]=x[,i]/gexp
# 	}

# 	return(x)
# }

get_outfile=function(out_dir, plot, gId) {
	dir.create(paste(out_dir, "/plots", sep=""), showWarnings = FALSE)

	out_dir_plots=paste(out_dir, "plots", plot, sep="/")
	dir.create(out_dir_plots, showWarnings = FALSE)

	subdir=paste(unlist(strsplit(gId, ''))[c(1:12)], collapse='')
	out_subdir=paste(out_dir_plots, subdir, sep="/")
	dir.create(out_subdir, showWarnings = FALSE)
	
	outfile=paste(out_subdir, "/", gId, ".pdf", sep="")

	return(outfile)
}

summarise_gexp=function (gexp, cond, text) {
	x=gexp[cond]
	x_m=formatC(round(mean(x), 2), 2, format="f")
	x_sd=formatC(round(sd(x), 2), 2, format="f")

	result=paste(text, " - mean=", x_m, ", sd=", x_sd, sep="")
	return(result)
}

summarise_dominance=function(data, cond, text) {
	x=calculate_dominance(data[cond,])
	x_m=formatC(round(mean(x), 2), 2, format="f")
	x_sd=formatC(round(sd(x), 2), 2, format="f")

	result=paste(text, " - mean=", x_m, ", sd=", x_sd, sep="")
	return(result)
}