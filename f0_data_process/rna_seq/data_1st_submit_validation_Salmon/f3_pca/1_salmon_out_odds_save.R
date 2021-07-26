#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("tximport")
#biocLite("tximportData")
library("tximport")
library("DESeq2")
library("tximportData")
library(gplots)
library(RColorBrewer)
source('plotPCAWithSampleNames.R')
#install.packages("ggrepel")

###################
#COPY REGION START 
################### 

file_dir = getwd()
setwd(file_dir)
tx2gene_dir="/nv/vol190/zanglab/zw5j/data/index/salmon_index/tx2gene_ensembl2symbol"
tx2gene <- read.csv(file.path(tx2gene_dir, "hg38_v92_tx2gene.csv"),header=TRUE,sep="\t")

# read the salmon-out of all samples
salmon_out_dir = "/nv/vol190/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/rna_seq/data_1st_submit_salmon/processed_results"
all_samples = list.files(salmon_out_dir)

#remove_samples <- c('R_7')
#all_samples = all_samples[!all_samples %in% remove_samples]
#kept_substring <- c('Jurkat')
#all_samples <- all_samples[grepl(kept_substring,all_samples)]


sample_files <- file.path(salmon_out_dir,all_samples, "quant.sf")
all(file.exists(sample_files))
names(sample_files) <- all_samples

# get the col match info of all samples
sampleTable_tmp = read.table(file.path(file_dir,"colNames.txt"),header=TRUE)
# I do not know why but after rewrite the sampleTable, the result could run correctly
sampleTable <- data.frame(condition=factor(sampleTable_tmp$condition),id=factor(sampleTable_tmp$id),name=factor(sampleTable_tmp$name))
rownames(sampleTable) <- sampleTable_tmp$id

###################
# COPY REGION END 
###################

txi=tximport(sample_files,type='salmon',tx2gene=tx2gene)
save(txi,file="f1_R_odds_saved/txi")
#load(file="f1_R_odds_saved/txi")
sampleTable <- sampleTable[colnames(txi$counts),]

ddsTxi=DESeqDataSetFromTximport(txi,colData=sampleTable,design=~condition)
dds=DESeq(ddsTxi)
save(dds,file = "f1_R_odds_saved/dds")
#load(file = "f1_R_odds_saved/dds")
#dir.create('f2_samlmon_pca')

log_dds<-rlog(dds)
save(log_dds,file = "f1_R_odds_saved/log_dds")
#load(file = "f1_R_odds_saved/log_dds")

