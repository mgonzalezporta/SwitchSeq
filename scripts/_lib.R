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



plot_stars=function(gId, x, gexp, cond1, cond2, outfile) {
	# data normalisation
	for (i in 1:dim(x)[2]) {
		x[,i]=x[,i]*abs(sqrt(gexp/max(gexp)))
	}
	x[is.na(x)]=0
	x=x/max(x)

	# gexp info
	gexp_sum_cond1=summarise_gexp(gexp, cond1, "   condition 1 (left)")
	gexp_sum_cond2=summarise_gexp(gexp, cond2, "   condition 2 (right)")

	# legend title
	legend_title=paste(gId, "\ngene expression summary (RPKMs):\n", 
		gexp_sum_cond1, "\n", gexp_sum_cond2, "\n\n", 
		"legend:", sep="")

	# colors
	nOfT=dim(x)[2]
	biotypes=unlist(strsplit(colnames(x), " - "))
	biotypes=data.frame(table(biotypes[seq(1, length(biotypes), by=2)]))
	colnames(biotypes)=c("tBiotype", "n")
	
	col_intervals=t(data.frame(
	        c("nonsense_mediated_decay", 0.05, 0.15),
	        c("processed_transcript", 0.3, 0.4),
	        c("protein_coding", 0.6, 0.7),
	        c("retained_intron", 0.85, 0.95)
	))
	rownames(col_intervals)=NULL
	colnames(col_intervals)=c("tBiotype", "start", "end")
	
	m=merge(biotypes, col_intervals, by="tBiotype")
	
	col=c()
	for (i in 1:dim(m)[1]) {
	        n=as.numeric(as.character(m[i,]$n))
	        start=as.numeric(as.character(m[i,]$start))
	        end=as.numeric(as.character(m[i,]$end))
	        col=c(col, rainbow( n, s = 0.7, v = 0.8, alpha=0.7, start=start ,end=end))
	}
	
	# plot layout
	nrows_layout=max(length(cond1), length(cond2))

	h=nrows_layout*.5
	if (h<4) { h=3 }

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
		palette(col)
		
		par(fg="white", col="gray35")
		par(mar=c(0, 0, 0, 0))
		par(oma=c(0, 0, 0, 0))
		par(mfrow=c(1,2))	
			
		stars(x, 
		    draw.segments = T,
			cex.main=0.6,
			cex=0.3,
			ncol=2,
			len=0.4,
			scale=F,
			locations=loc
		)

		plot.new()
		par(oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))

		legend("topleft",
			inset=c(0,0.18),
			title=legend_title, 
			title.adj=0,
			legend=colnames(x), 
			fill=col, 
			border=col,
			#title.col="gray5",
			cex=0.4
		)
	garbage <- dev.off()	# supress output
}

############################################
# Internal functions
summarise_gexp=function (gexp, cond, text) {
	x=gexp[cond]
	x_m=formatC(round(mean(x), 2), 2, format="f")
	x_sd=formatC(round(sd(x), 2), 2, format="f")

	result=paste(text, " - mean=", x_m, ", sd=", x_sd, sep="")
	return(result)
}