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
salmon_out_dir = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f1_TALL_CTCF_binding_vs_GeneExpr/f2_combine_patients/f0_process_rna/salmon_quant"
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


####################
# PACA analysis
####################
library(RColorBrewer)
colors = brewer.pal(9, "Set1")
#colorCodes <- c('1'="red", 'a'="green", C="blue", D="yellow")
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label")
    code <- sampleTable[label,]$condition
    #code <- substr(label, 1, 1);#sprintf(code);sprintf('a');sprintf('a')
    ## use the following line to reset the label to one letter code
    # attr(x, "label") <- code
    attr(x, "nodePar") <- list(lab.col=colors[code])
  }
  return(x)
}


source('plotPCAWithSampleNames.R')
pdf("f2_samlmon_pca/all_samples_pca_log_allgenes.pdf", width=5, height=5)
plotPCAWithSampleNames(log_dds, intgroup="condition", ntop=40000,label_col=4)
dev.off()

dists <- dist(t(assay(log_dds)))
hc <- hclust(dists)
d <- dendrapply(as.dendrogram(hc), labelCol)
pdf("f2_samlmon_pca/all_samples_hclust_log_allgenes.pdf", width=6, height=4)
par(mar=c(15,2,0,0)+1)
plot(d)
dev.off()

###select several groups
#plotPCAWithSampleNames(log_dds[,!log_dds$condition==4], intgroup="condition", ntop=40000,label_col=3)
#plotPCAWithSampleNames(log_dds[,(log_dds$condition==1)|(log_dds$condition==2)|(log_dds$condition==3)], intgroup="condition", ntop=40000,label_col=3)
