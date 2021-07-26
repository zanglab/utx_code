#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("tximport")
#biocLite("tximportData")
library("tximport")
library("DESeq2")
library("tximportData")
library(gplots)
library(RColorBrewer)
#source('plotPCAWithSampleNames.R')

###################
#COPY REGION START 
################### 

file_dir = getwd()
setwd(file_dir)
tx2gene_dir="/nv/vol190/zanglab/zw5j/data/index/salmon_index/tx2gene_ensembl2symbol"
tx2gene <- read.csv(file.path(tx2gene_dir, "hg38_v92_tx2gene.csv"),header=TRUE,sep="\t")

# read the salmon-out of all samples
# salmon_out_dir = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f1_TALL_CTCF_binding_vs_GeneExpr/f2_combine_patients/f0_process_rna/salmon_quant"
all_samples = list.files(salmon_out_dir)

#remove_samples <- c('R_7')
#all_samples = all_samples[!all_samples %in% remove_samples]

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

load(file="f1_R_odds_saved/txi")
sampleTable <- sampleTable[colnames(txi$counts),]
load(file = "f1_R_odds_saved/dds")
load(file = "f1_R_odds_saved/log_dds")
#dir.create('f2_samlmon_pca')

 

#####################
# identify DEG
#####################

func_deseq2<-function(sampleTable,txi,dds,treat,ctrl,csv_name)
{
counts_table=txi$counts
#counts_table <-counts_table[!rowSums(counts_table==0)>=10,]
# keep only those data in a/b group
filter_id = sampleTable[(sampleTable$condition==treat)|(sampleTable$condition==ctrl),]
counts_table = counts_table[,filter_id$id]

filtered_norm_counts<-counts_table[!rowSums(counts_table==0)>=1,]
filtered_norm_counts<-as.data.frame(filtered_norm_counts)
GeneID<-rownames(filtered_norm_counts)
filtered_norm_counts<-cbind(filtered_norm_counts,GeneID)
dim(filtered_norm_counts)
#head(filtered_norm_counts)

res=results(dds,contrast=c("condition",treat,ctrl)) # contrast=c("condition","treated","untreated") 
res_ordered<-res[order(res$padj),]
GeneID<-rownames(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_genes<-cbind(res_ordered,GeneID)
dim(res_genes)
res_genes_merged <- merge(filtered_norm_counts,res_genes,by=unique("GeneID"))
res_ordered<-res_genes_merged[order(res_genes_merged$padj),]
dim(res_ordered)
#head(res_ordered)
write.csv(res_ordered, file = file.path("f3_deseq_out/",csv_name),row.names=FALSE)
}


CtrlList = c('CD4','CD4','CD4','CD4','CD4','CD4')
TreatList = c('PD30','PD40','PD9','PDBG','CUTLL1','Jurkat')

for (i in 1:6){
treat = TreatList[i]
ctrl = CtrlList[i]
names_sep = c("treated_",treat,"_vs_ctrl_",ctrl,".csv")
func_deseq2(sampleTable,txi,dds,treat,ctrl,paste(names_sep,collapse=""))
}

#####record the gene level abundance/TMP
abundance_table=txi$abundance
filtered_norm_abundance<-abundance_table[!rowSums(abundance_table==0)>=120,]
filtered_norm_abundance<-as.data.frame(filtered_norm_abundance)
dim(filtered_norm_abundance)
head(filtered_norm_abundance)
write.csv(filtered_norm_abundance, file = file.path("f3_deseq_out/","gene_level_abundance_TPM_all.csv"),row.names=TRUE)

