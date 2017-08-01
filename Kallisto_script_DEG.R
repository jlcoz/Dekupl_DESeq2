suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(pheatmap)))
suppressWarnings(suppressMessages(library(ggplot2)))

#INPUT
conditionA                              = snakemake@params$conditionA
conditionB                              = snakemake@params$conditionB
gene_counts                             = snakemake@input$gene_counts
pvalue_threshold                        = snakemake@params$pvalue_threshold
sample_conditions                       = snakemake@input$sample_conditions

#OUTPUT
differentially_expressed_genes          = snakemake@output$differentially_expressed_genes
differentially_expressed_genes_filtered = snakemake@output$differentially_expressed_genes_filtered
dist_matrix                             = snakemake@output$dist_matrix
norm_counts	                            = snakemake@output$norm_counts
pca_design                              = snakemake@output$pca_design

#LOG
output_log                              = snakemake@log[[1]]

#FUNCTION
logging <- function(str) {
  sink(file=paste(output_log), append=TRUE, split=TRUE)
  print(paste(Sys.time(),str))
  sink()
}

printing <- function(str) {
  sink(file=paste(output_log), append=TRUE, split=TRUE)
  print(paste(str))
  sink()
}

logging("Start Kallisto DE analysis")

# Load counts data
countsData = read.table(gene_counts,header=T,row.names=1)

# Load col data with sample specifications
colData = read.table(sample_conditions,header=T,row.names=1)
colData = colData[colnames(countsData),,drop=FALSE]


# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=countsData,
                              colData=colData,
                              design = ~ condition)

                              
dds <- DESeq(dds)

#WITHOUT FILTERING

#normalized counts 
NormCounts<- as.data.frame(counts(dds,normalized=TRUE))
  
#writing in a file normalized counts
normalized_counts<-data.frame(id=row.names(NormCounts),NormCounts,row.names=NULL)
write.table(normalized_counts,
            file=norm_counts, 
            sep="\t",
            row.names=F, 
            col.names=T, quote=F)

printing(resultsNames(dds))

results_Kallisto <- results(dds, contrast = c("condition", "A", "B"))

#write DEGs

write.table(data.frame(ID=rownames(results_Kallisto),
                       baseMean=results_Kallisto$baseMean,
                       log2FoldChange=results_Kallisto$log2FoldChange,
                       pvalue=results_Kallisto$pvalue,
                       padj=results_Kallisto$padj),
            file=differentially_expressed_genes,
            sep="\t",
            row.names=FALSE,
            quote=FALSE)

#WITH FILTER ON ADJUSTED PVALUE AND REMOVE NA 

results_Kallisto = na.omit(results_Kallisto)
results_Kallisto = results_Kallisto[results_Kallisto$padj<0.05,]

#Get DE gene_ID
gene_ID_DE = rownames(results_Kallisto)
normalized_counts_filtered = NormCounts[rownames(NormCounts)%in%gene_ID_DE,]

#writing in a file filtered norm counts
normalized_counts_filtered<-data.frame(id=rownames(NormCounts), NormCounts, row.names=NULL)

write.table(normalized_counts_filtered,
            file=norm_counts, 
            sep="\t",
            row.names=F, 
            col.names=T, quote=F)

#write DEGs filtered
write.table(data.frame(ID=rownames(results_Kallisto),
                       baseMean=results_Kallisto$baseMean,
                       log2FoldChange=results_Kallisto$log2FoldChange,
                       pvalue=results_Kallisto$pvalue,
                       padj=results_Kallisto$padj),
            file=differentially_expressed_genes_filtered,
            sep="\t",
            row.names=FALSE,
            quote=FALSE)

#PCA
rld<-rlog(dds)
sampleDists<-dist(t(assay(rld)))
sampleDistMatrix<-as.matrix(sampleDists)
rownames(sampleDistMatrix)<-colnames(rld)
colnames(sampleDistMatrix)<-colnames(rld)
colours=colorRampPalette(rev(brewer.pal(9,"Blues")) )(255)

pdf(dist_matrix,width=15,height=10)

pheatmap(sampleDistMatrix,
         main="Clustering of samples",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colours,
         fontsize = 14)
 
data<-plotPCA(rld,ntop=nrow(rld),returnData=TRUE)

write.table(data,pca_design,row.names=F, col.names=T, quote=F,sep="\t")

print(ggplot(data,aes(PC1,PC2,color=condition))+geom_point()+geom_text(aes(label=name),hjust=0,vjust=0))

dev.off()
