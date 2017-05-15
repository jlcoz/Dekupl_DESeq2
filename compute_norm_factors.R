if (!require("data.table")) {
  install.packages("data.table", repos="http://cran.rstudio.com/",dependencies=TRUE) 
  library("data.table")
}

  ## SET UP VARS
data=snakemake@input$gene_counts
output_tmp=snakemake@config$tmp_dir
output_norm_factors=snakemake@output

sampling_size = snakemake@input$sampling_size
seed = 30

dir.create(output_tmp, showWarnings = FALSE)

setwd(output_tmp)

  ## SAMPLING
system(paste("zcat ",data,
             " | awk -v sampling_size=",sampling_size," -v seed=",seed,
             " 'BEGIN{nb_kmers=0;}",
             "{if(nb_kmers==sampling_size)exit 1;if(NR % seed ==0 || NR ==1){print $0;nb_kmers++;}}' > selected_kmers",sep=""))

selected_kmers_counts <- data.frame(fread(paste("selected_kmers")))
 
loggeomeans <- rowMeans(log(selected_kmers_counts[,2:ncol(selected_kmers_counts)]))

normFactors <- apply(selected_kmers_counts[,2:ncol(selected_kmers_counts)], 2, function(cnts) { exp(median((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))})
 
 ## WRITING THE NORMALIZATION FACTORS
 ## SAMPLE NAMES \t NORMALIZATION FACTOR
write.table(data.frame(sample=names(normFactors),normalization_factor=as.vector(normFactors)),
             file = paste(output_norm_factors),
             sep="\t",
             quote=FALSE,
             row.names = FALSE)
