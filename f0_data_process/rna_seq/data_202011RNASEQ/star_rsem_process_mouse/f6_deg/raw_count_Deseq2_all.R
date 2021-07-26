library("DESeq2")

raw_count_matrix = '../f3_expr/expected_count.csv'
rcm <- as.matrix(read.table(raw_count_matrix,sep=",",header=TRUE,row.names=1,check.names=FALSE))
# filter out genes with ≤ 10 reads on average across all samples
rcm<-rcm[,!grepl("ES",colnames(rcm))]
rcm = rcm[rowMeans(rcm)>0,]
dim(rcm)
head(rcm,2)

# prepare the matrix
coldata <- read.table("colname_all.txt",check.names=FALSE,header=TRUE)
# coldata <- data.frame(condition=factor(coldata$condition),id=factor(coldata$id))
rownames(coldata) <- coldata$id
coldata <- coldata[colnames(rcm),]

all(rownames(coldata) == colnames(rcm))
mode(rcm)<-"integer"

# Deseq2 
GeneID<-rownames(rcm)
dds <- DESeqDataSetFromMatrix(countData = rcm,
                            colData = coldata,
                            design = ~ condition)
dds <- DESeq(dds)

# # Deseq2 for each treat and control
CtrlList = c('UTX-WT','UTX-WT','UTX-WT')
TreatList = c('UTX-KO','UTX-DEL3-3','UTX-DEL3-72')

for (i in 1:3){
treat = TreatList[i]
ctrl = CtrlList[i]

res=lfcShrink(dds,contrast=c("condition",treat,ctrl))
# res_ordered<-res[order(res$padj),]
res_genes<-cbind(res,GeneID)

filter_id = coldata[(coldata$condition==treat)|(coldata$condition==ctrl),]
counts_table = rcm[,rownames(filter_id)]
head(counts_table)
counts_table<-cbind(counts_table,GeneID)
res_genes_merged <- merge(counts_table,res_genes,by=unique("GeneID"))
res_genes_merged <- res_genes_merged[order(res_genes_merged$padj),]
write.csv(res_genes_merged, file = file.path("f1_deseq2_out/",paste("treated_",treat,"_vs_ctrl_",ctrl,'.deseq2.csv',sep='')),row.names=FALSE)
}

