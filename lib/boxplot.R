#library(beeswarm)

# adapt with other gIds
# play with height

plot_boxplots=function(gId, x, gexp) {

	# gene exp
	a=gexp[seq(1, length(gexp), by=2)]
	b=gexp[seq(2, length(gexp), by=2)]
	
	lgexp=list(a, b)
	names(lgexp)=c("N", "T")
	
	# texp
	ltexp=list()
	for (i in 1:dim(x)[2]) {
		a=x[seq(1, dim(x)[1], by=2),i]
		b=x[seq(2, dim(x)[1], by=2),i]
		ltexp[[i]]=as.matrix(data.frame(a,b))
	
		# names
		colnames(ltexp[[i]])=c("N", "T")
		names(ltexp)[i]=gsub(" - ", " \n ", colnames(x)[i])
	}
	
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
	
	col=t(data.frame(NA, NA))
	rownames(col)=c("N", "T")
	for (i in 1:dim(m)[1]) {
	        n=as.numeric(as.character(m[i,]$n))
	        start=as.numeric(as.character(m[i,]$start))
	        end=as.numeric(as.character(m[i,]$end))

		tcol=t(data.frame(
			rainbow( n, s = 0.7, v = 0.9, alpha=0.7, start=start ,end=end),
                        rainbow( n, s = 0.7, v = 0.5, alpha=0.7, start=start ,end=end)
		))
		rownames(tcol)=c("N", "T")
		col=cbind(col, tcol)
	}
	col=col[,-1]

	#####
	
	# plot
	setwd("/nfs/gns/homes/mar/public_html/results_cagekid/switch_events/plots/boxplots")
        subDir=paste(unlist(strsplit(gId, ''))[c(1:11)], collapse='')
        dir.create(file.path('./', subDir), showWarnings = FALSE)
	setwd(subDir)
        outfile=paste(gId, ".pdf", sep="")

	pdf(file=outfile, width=6, height=(4+nOfT-1))
	layout(1:(nOfT+2), height=rep(2, nOfT+2))
	
	# par
	par(oma=c(4.5, 12, 0.5, 0.5), las=1, cex=1)

	# gexp
	par(mar=c(0.2, 0.2, 0.2, 0.2), bty="n")
	boxplot(lgexp, outline=T, pch=20, horizontal=T, xlab="", yaxt="n", col=c("white", "gray90"))
	#beeswarm(lgexp, add=TRUE, horizontal=T, col=col, pch=20)
	
	mtext("gene expression (FPKMs)", 1, 2.5)
	axis(2, labels=F, lwd.ticks=0, lwd=1, at=c(1,2))
	mtext(gId, 2, 0.75)
	
	frame()
	# texp
	par(mar=c(0.2, 0.2, 0.2, 0.2))
		# `mar` controls the space around each boxplot group
	
	# Plot all the boxes
	for(i in 1:length(ltexp)){
		boxplot(ltexp[[i]], ylim=c(0,1), axes=F, outline=T, pch=20, horizontal=TRUE, col=col[,i])
	#	beeswarm(ltexp[[i]][,1], add=TRUE, horizontal=T, col=col[1], pch=20, at=1)
	#	beeswarm(ltexp[[i]][,2], add=TRUE, horizontal=T, col=col[2], pch=20, at=2)
	
		# tIds
	        axis(2, labels=F, lwd.ticks=0, lwd=1, at=c(1,2))
		mtext(names(ltexp)[i], 2, 0.75)
	
	        # x axis
		if(i == length(ltexp)){
	        	axis(1, las=1)
		        mtext("transcript relative abundance", 1, 2.5)
	    	}
	}
	garbage <- dev.off()	# supress output

	# rotate plot
	c=paste("pdf90", outfile, sep=" ")
	system(c)
	c=paste("rm", outfile, sep=" ")
	system(c)
}

