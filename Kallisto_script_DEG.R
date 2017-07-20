suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(pheatmap)))
suppressWarnings(suppressMessages(library(ggplot2)))

#INPUT
gene_counts                             = snakemake@input$gene_counts
sample_conditions                       = snakemake@input$sample_conditions
conditionA                              = snakemake@input$conditionA
conditionB                              = snakemake@input$conditionB
pvalue_threshold                        = snakemake@params$pvalue_threshold


#OUTPUT
differentially_expressed_genes          = snakemake@output$differentially_expressed_genes
differentially_expressed_genes_filtered = snakemake@output$differentially_expressed_genes_filtered
dist_matrix                             = snakemake@output$dist_matrix
norm_counts	                            = snakemake@output$norm_counts
pca_design                              = snakemake@output$pca_design
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

printing(colnames(countsData))
printing(rownames(colData))

colData = colData[colnames(countsData),,drop=FALSE]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=countsData,
                              colData=colData,
                              design = ~ condition)
                                       
dds <- DESeq(dds)

#normalized counts 
NormCount<- as.data.frame(counts(dds,normalized=TRUE))
  
#writing in a file normalized counts
normalized_counts<-data.frame(id=row.names(NormCount),NormCount,row.names=NULL)
write.table(normalized_counts,
            file=norm_counts, 
            sep="\t",
            row.names=F, 
            col.names=T, quote=F)

printing(resultsNames(dds))

results_Kallisto <- results(dds, contrast = c("condition", conditionA, conditionB))

# Write DEGs

#WITHOUT FILTERING
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

results_Kallisto = results_Kallisto[results_Kallisto$padj<pvalue_threshold,]

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
