args=commandArgs(TRUE)
plot=as.character(args[1])
input_file=as.character(args[2])
out_dir=as.character(args[3])
data_dir=as.character(args[4])
species=as.character(args[5])
ensembl_v=as.character(args[6])
cond1=as.character(args[7])
cond2=as.character(args[8])
gId=as.character(args[9])

# load data
command=paste("grep", gId, input_file, sep=" ")
f=pipe(command)
data=read.table(f, as.is=T)
command=paste("head -n1", input_file, sep=" ")
f=pipe(command)
header=read.table(f, as.is=T)
colnames(data)=header


ensembl_annot=paste(data_dir, "/", species, "/_ensembl", ensembl_v, ".annot_coding.txt", sep="")
command=paste("grep", gId, ensembl_annot, sep=" ")
f=pipe(command)
tBiotype=read.table(f)[,4:5]
colnames(tBiotype)=c("tId", "tBiotype")

# format data
x=data.frame(t(data))
colnames(x)=as.character(as.matrix(x[2,]))
x=x[-c(1:2),]
x=x[order(rownames(x)),]
# update colnames
tmp=data.frame(colnames(x))
colnames(tmp)="tId"
m=merge(tmp, tBiotype, by="tId", sort=F)
m$id=paste(m$tBiotype, m$tId, sep=" - ")
colnames(x)=m$id
x=x[,sort(colnames(x))]

# normalisation
for (i in 1:dim(x)[2]) {
        x[,i]=as.numeric(as.character(x[,i]))
}
gexp=apply(x, 1, sum)
for (i in 1:dim(x)[2]) {
        x[,i]=x[,i]/gexp
}

# format condition
tmp=as.numeric(strsplit(cond1, "-")[[1]])
cond1=(tmp[1]-2):(tmp[2]-2)
tmp=as.numeric(strsplit(cond2, "-")[[1]])
cond2=(tmp[1]-2):(tmp[2]-2)


#########
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

##

if (plot == "starplots") {
	source("./src/stars.R")
	outfile=get_outfile(out_dir, "starplots", gId)
	plot_stars(gId, x, cond1, cond2, outfile)
}

if (plot == "boxplots") {
	source("./src/boxplots.R")
	outfile=get_outfile(out_dir, "boxplots", gId)
	plot_boxplots(gId, x, gexp, out_dir_plots)
}

if (plot == "both") {
	source("./src/stars.R")
	outfile=get_outfile(out_dir, "starplots", gId)
	plot_stars(gId, x, outfile)

	source("./src/boxplots.R")
	outfile=get_outfile(out_dir, "boxplots", gId)
	plot_boxplots(gId, x, gexp, out_dir_plots)
}