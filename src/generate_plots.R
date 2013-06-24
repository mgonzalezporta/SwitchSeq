args=commandArgs(TRUE)
gId=as.character(args[1])
input_file=as.character(args[2])

# load data
command=paste("grep", gId, input_file, sep=" ")
f=pipe(command)
data=read.table(f, as.is=T)
command=paste("head -n1", input_file, sep=" ")
f=pipe(command)
header=read.table(f, as.is=T)
colnames(data)=header

command=paste("grep", gId, "./test_23/hsa/_ensembl66.annot_coding.txt", sep=" ")
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

#########

source("./src/stars.R")
plot_stars(gId, x)

#source("/homes/mar/home_microarray/workspace/cagekid/scripts/analyses/switch_events/plots/boxplots.R")
#plot_boxplots(gId, x, gexp)
