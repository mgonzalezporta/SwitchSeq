args=commandArgs(TRUE)
# input_file=as.character(args[1])
# out_dir=as.character(args[2])
# data_dir=as.character(args[3])
# species=as.character(args[4])
# ensembl_v=as.character(args[5])
# cond1=as.character(args[6])
# cond2=as.character(args[7])
# gId=as.character(args[8])

# tmp
input_file="/Users/mar/Desktop/app/tExp.fpkms"
out_dir="/Users/mar/Desktop/app/html_12/"
data_dir="/Users/mar/Desktop/app/data_10/"
species="hsa"
ensembl_v=66
cond1="3-5"
cond2="6-9"
gId="ENSG00000106771"

# load data
command=paste("grep", gId, input_file, sep=" ")
f=pipe(command)
data=read.table(f, as.is=T)
command=paste("head -n1", input_file, sep=" ")
f=pipe(command)
header=read.table(f, as.is=T)
colnames(data)=header


ensembl_annot=paste(data_dir, "/", species, "/_ensembl", ensembl_v, ".annot_coding.1.txt", sep="")
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

source("./_lib.R")
outfile=get_outfile(out_dir, "starplots", gId)
plot_stars(gId, x, gexp, cond1, cond2, outfile)

outfile=get_outfile(out_dir, "distrplots", gId)
plot_distrplot(gId, x, gexp, out_dir_plots)

#########

