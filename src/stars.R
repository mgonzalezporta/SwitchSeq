plot_stars=function(gId, x) {
	# data normalisation
	for (i in 1:dim(x)[2]) {
		x[,i]=x[,i]*abs(sqrt(gexp/max(gexp)))
	}
	x[is.na(x)]=0
	x=x/max(x)

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
	
	# optimise output dir structure
	setwd("./html_23/plots/stars")
	subDir=paste(unlist(strsplit(gId, ''))[c(1:12)], collapse='')
	dir.create(file.path('./', subDir), showWarnings = FALSE)
	outfile=paste(subDir, "/", gId, ".pdf", sep="")

	# plot layout
	cond1=1:3	##
	cond2=4:7	##
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
	# loc=t(matrix(
	# 	c(0,3, 0,2, 0,1, 1,3, 1,2, 1,1, 1,0),
	# 	nrow=2))

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
		par(oma=c(0, 0, 0.7, 0))
		legend("topleft",
			title=gId, legend=colnames(x), 
			fill=col, 
			border=col,
			title.col="gray5",
			cex=0.4
		)
	garbage <- dev.off()	# supress output
}
