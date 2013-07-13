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
#input_file="/Users/mar/Desktop/app/tExp.fpkms"
input_file="/Users/mar/Desktop/app/2013.05.30_45pairs.dexseq2.txt"
out_dir="/Users/mar/Desktop/app/html_12/"
data_dir="/Users/mar/Desktop/app/data_10/"
species="hsa"
ensembl_v=66
#cond1="3-5"
#cond2="6-9"
cond1="3-48"
cond2="49-92"
ensembl_annot=paste(data_dir, "/", species, "/_ensembl", ensembl_v, ".annot_coding.1.txt", sep="")
#gId="ENSG00000106771"
#gId="ENSG00000111615"
gId="ENSG00000180104"


source("./_lib.R")

# args
cond1=format_cond(cond1)
cond2=format_cond(cond2)

# data
se=new_switch_event(gId, input_file, ensembl_annot)
gexp=apply(se, 1, sum)
norm_se=normalise_se(se, gexp)

# plots
outfile=get_outfile(out_dir, "starplots", gId)
plot_stars(gId, x, gexp, cond1, cond2, outfile)

outfile=get_outfile(out_dir, "distrplots", gId)
plot_distrplot(gId, x, gexp, out_dir_plots)