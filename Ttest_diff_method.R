Ttest=snakemake@input$Ttest_filter
no_GENCODE=snakemake@input$counts
normalization_factor_path = snakemake@input$sample_conditions
pvalue_threshold=snakemake@config$Ttest$pvalue_threshold
log2fc_threshold=snakemake@config$Ttest$log2fc_threshold

output_diff_counts=snakemake@output$diff_counts
output_pvalue_all=snakemake@output$pvalue_all

system(paste(Ttest," -p ",pvalue_threshold, " -f ",log2fc_threshold," ",no_GENCODE," ",normalization_factor_path," A B | gzip -c > ",output_diff_counts,sep=""))
system(paste("mv raw_pvals.txt ",output_pvalue_all,sep=""))
