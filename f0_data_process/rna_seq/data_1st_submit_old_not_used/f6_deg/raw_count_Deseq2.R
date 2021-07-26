library("DESeq2")


raw_count_matrix = '../f3_expr/expected_count.csv'
rcm <- as.matrix(read.table(raw_count_matrix,sep=",",header=TRUE,row.names=1,check.names=FALSE))
# filter out genes with ≤ 10 reads on average across all samples
rcm = rcm[rowMeans(rcm)>0,]
dim(rcm)
head(rcm,2)


colname_txts <- c("colname.txt")
# Deseq2 for each patient sample
for (colname_txt in colname_txts){
coldata <- read.csv(colname_txt,row.names=1,check.names=FALSE,sep="\t")
rcm_tmp <- rcm[,rownames(coldata)]
rcm_tmp = rcm_tmp[rowMeans(rcm_tmp)>0,]
all(rownames(coldata) == colnames(rcm_tmp))

mode(rcm_tmp)<-"integer"
GeneID<-rownames(rcm_tmp)
dds <- DESeqDataSetFromMatrix(countData = rcm_tmp,
                            colData = coldata,
                            design = ~ condition)
dds <- DESeq(dds)
res=results(dds,contrast=c("condition","treat","control"))
# res_ordered<-res[order(res$padj),]
res_genes<-cbind(res,GeneID)
rcm_tmp<-cbind(rcm_tmp,GeneID)
res_genes_merged <- merge(rcm_tmp,res_genes,by=unique("GeneID"))
write.csv(res_genes_merged, file = file.path("f1_deseq2_out/",paste(colname_txt,'.deseq2.csv',sep='')),row.names=FALSE)
}