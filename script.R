start_of_analysis = Sys.time()

if (!require("data.table")) {
  install.packages("data.table", repos="http://cran.rstudio.com/",dependencies=TRUE) 
  library("data.table")
}

if (!require("foreach")) {
  install.packages("foreach", repos="http://cran.rstudio.com/",dependencies=TRUE) 
  library("foreach")
}

if (!require("doParallel")) {
  install.packages("doParallel", repos="http://cran.rstudio.com/",dependencies=TRUE) 
  library("doParallel")
}

if (!require("DESeq2")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("Biobase")
  biocLite("DESeq2")
  library("DESeq2")
}

data_path=snakemake@input$counts
normalization_factor_path = snakemake@input$sample_conditions
nb_conditionA=length(unlist(snakemake@config$samples)[unlist(snakemake@config$samples)=="A"])
nb_conditionB=length(unlist(snakemake@config$samples)[unlist(snakemake@config$samples)=="B"])
pvalue_threshold=snakemake@config$Ttest$pvalue_threshold
nb_core=snakemake@config$nb_threads

output_tmp=snakemake@config$tmp_dir
output_diff_counts=snakemake@output$diff_counts
output_pvalue_all=snakemake@output$pvalue_all

output_log=snakemake@log[[1]]

dir.create(output_tmp, showWarnings = FALSE)

split_lines = 1000000

registerDoParallel(cores=nb_core)

setwd(output_tmp)

# SPLIT THE MAIN FILE INTO CHUNKS WITH AUTOINCREMENTED NAMES AND SAVE THE HEADER INTO A FILE

system(paste("zcat ",data_path,
             " | awk -v split_lines=",split_lines,
             " -v output_tmp=",output_tmp,
             " 'NR%split_lines==1{OFS=\"\\t\";x=++i\"_subfile.txt\"}{OFS=\"\";print > output_tmp x}'",
             "; head -1 ",output_tmp,"1_subfile.txt >",output_tmp,"header_large_file.txt ; sed -i '1d' ",output_tmp,"1_subfile.txt; ",sep=""))

nb_line_last_file=system(paste("cd ",output_tmp," ; cat $(ls | sort -n | grep subfile |tail -1)|wc -l", sep=""), intern=TRUE)

#CONCATENATE THE LAST FILE CREATED WITH THE SECOND LAST ONE AND THEN REMOVE IT
#IN ORDER TO AVOID TOO SHORT FILES
#IF THE MERGED FILE IS TOO LARGE, IT WILL BE DIVIDED IN TWO

if(nb_line_last_file < (split_lines/2)){

  print(paste("The last file has",nb_line_last_file,"line(s) it will be concatenated to the second last one"))
  system(paste("cd ",output_tmp," ; file_number=$(ls | grep subfile|wc -l) ",
               "; file_number=$(echo $((file_number-1)))",
               "; last_2_files=$(ls | sort -n | grep subfile | tail -2)",
               "; cat $last_2_files > tmp_concat ",
               "&& mv tmp_concat ${file_number}_subfile.txt",
               "&& rm $((file_number+1))_subfile.txt ",
               sep=""))
               
  nb_line_last_file=system(paste("cd ",output_dir," ; cat $(ls | sort -n | grep subfile |tail -1)|wc -l", sep=""), intern=TRUE)
  
  print(paste("The last file has",nb_line_last_file,"line(s) it will be splitted in two"))
  
  system(paste("cd ",output_dir," ; file_number=$(ls | grep subfile | wc -l) ",
               "; last_file=$(ls | sort -n | grep subfile | tail -1)",
               "; split -n l/2 $last_file",
               "; mv xaa ${file_number}_subfile.txt",
               "; file_number=$(echo $((file_number+1)))",
               "; mv xab ${file_number}_subfile.txt",
               sep=""))
}

sink(output_log, append=TRUE, split=TRUE)
print(paste(Sys.time(),"Split done"))
sink()

## LOADING DATA
lst_files = system(paste("find",output_tmp,"-iname \"*_subfile.txt\" | sort -n"), intern = TRUE)
header = as.character(unlist(read.table(file = paste(output_tmp,"header_large_file.txt",sep=""), sep = "\t", header = FALSE)))


sink(output_log, append=TRUE, split=TRUE)
print(paste(Sys.time(),"Foreach on the", length(lst_files),"files"))
sink()

before_foreach<-Sys.time()

#DESeq2 ANALYSIS ON EACH CHUNKS
invisible(foreach(i=1:length(lst_files)) %dopar%{

  bigTab=data.frame(fread(paste(lst_files[i]),header=FALSE))
  #SET TAGS AS ROWNAMES
  rownames(bigTab)=bigTab[,1]
  #REMOVE THE TAG AS A COLUMN
  bigTab=bigTab[,2:ncol(bigTab)]
  names(bigTab)=header[2:length(header)]
  
  #FORMAT COLS DATA

  colDat <- data.frame(conds=factor(c(rep(0,nb_conditionA),rep(1,nb_conditionB))))

  countData = as.matrix(bigTab)

  dds <- DESeqDataSetFromMatrix(countData,
                              colData = colDat,
                              design = ~ conds)

  #LOADING PRIOR KNOWN NORMALISATION FACTORS

  size_factors = data.frame(fread(paste("cat ",normalization_factor_path," | awk '{print $1,$3}'")))

  NormCount_names = colnames(bigTab)
  rm(bigTab);gc()

  #RUN DESEQ2 AND COLLECT DESeq2 results
  
  #REPLACE SIZE FACTORS by SIZE FACTORS COMPUTED ON
  #THE ALL DATASET

  normFactors <- matrix(size_factors[,2],
                        ncol=ncol(dds),nrow=nrow(dds),
                        dimnames=list(1:nrow(dds),1:ncol(dds)),
                        byrow = TRUE)
  normalizationFactors(dds) <- normFactors
   
  #RUN DESeq2
  dds <- estimateDispersionsGeneEst(dds)
  dds <- estimateDispersionsFit(dds)
  dds <- estimateDispersionsMAP(dds)
  dds <- nbinomWaldTest(dds)
  resDESeq2 <- results(dds, pAdjustMethod = "none")
  
  #COLLECT COUNTS
  NormCount<- as.data.frame(counts(dds, normalized=TRUE))
  names(NormCount) <- NormCount_names
   
  # WRITE A TSV WITH THIS FORMAT FOR THE CURRENT CHUNK
  # Kmer_ID
  # meanA
  # meanB
  # log2FC
  # NormCount
  
  write.table(data.frame(ID=rownames(resDESeq2),
                         meanA=rowMeans(NormCount[,1:nb_conditionA]),
                         meanB=rowMeans(NormCount[,(as.numeric(nb_conditionA)+1):(as.numeric(nb_conditionA)+as.numeric(nb_conditionB))]),
                         log2FC=resDESeq2$log2FoldChange,
                         NormCount),
              file=paste(output_tmp,i,"_dataDESeq2_part_tmp", sep=""),
              sep="\t",quote=FALSE,
              row.names = FALSE,
              col.names = TRUE)
  
  # WRITE PVALUES FOR THE CURRENT CHUNK
  write.table(data.frame(ID=rownames(resDESeq2),pvalue=resDESeq2$pvalue),
                file=paste(output_tmp,i,"_pvalue_part_tmp",sep=""),
                sep="\t",quote=FALSE,
                row.names = FALSE,
                col.names = FALSE)

}) #END FOREACH

after_foreach<-Sys.time()

sink(output_log, append=TRUE, split=TRUE)
print(paste(Sys.time(),"Foreach done in : ",after_foreach-before_foreach))
sink()
  
  #MERGE ALL CHUNKS PVALUE INTO A FILE
system(paste("find ",output_tmp," -name '*_pvalue_part_tmp' | xargs cat > ",output_pvalue_all,sep=""))

sink(output_log, append=TRUE, split=TRUE)
print(paste(Sys.time(),"Pvalues merged into pvalueAll file"))
sink() 

  #MERGE ALL CHUNKS DESeq2 INTO A FILE
system(paste("find ",output_tmp," -name '*_dataDESeq2_part_tmp' | xargs cat | awk '{OFS=\"\\t\"}{if(NR==1){print}else{if($1 !~ \"ID\"){print}}}' > dataDESeq2All", sep=""))

sink(output_log, append=TRUE, split=TRUE)
print(paste(Sys.time(),"DESeq 2 results merged into dataDESeq2All file"))
sink()

  #SAVE THE COLNAMES INTO HEADER.TXT AND REMOVE THE COLNAME ID
system(paste("head -1 dataDESeq2All | cut -f2- > header_dataDESeq2All.txt ; sed -i 1d dataDESeq2All",sep=""))
  
  #CREATE AND WRITE THE ADJUSTED PVALUE UNDER THRESHOLD WITH THEIR ID
pvalueAll <- data.frame(fread(paste(output_pvalue_all),header=FALSE))
names(pvalueAll)=c("ID","pvalue")
adjPvalue <- p.adjust(as.numeric(as.character(pvalueAll[,"pvalue"])),"BH")
  
adjPvalue_dataframe = data.frame(ID=pvalueAll$ID,
                           pvalue=adjPvalue)

adjPvalue_dataframe = adjPvalue_dataframe[adjPvalue_dataframe$pvalue<pvalue_threshold,]

adjPvalue_dataframe = na.omit(adjPvalue_dataframe)

write.table(adjPvalue_dataframe,
            file="adj_pvalue",
            sep="\t",quote=FALSE,
            col.names = TRUE,
            row.names = FALSE)
  
print(paste(Sys.time(),"Pvalue are adjusted"))
  
  #SAVE THE HEADER
system(paste("head -1 adj_pvalue > header_adj_pvalue.txt ; sed -i 1d adj_pvalue",sep=""))
  
  #LEFT JOIN INTO dataDESeq2All
  #GET ALL THE INFORMATION (ID,MEAN_A,MEAN_B,LOG2FC,COUNTS) FOR DE KMERS
system(paste("sort -k1,1 adj_pvalue > sorted_adj_pvalue_tmp ; sort -k1,1 dataDESeq2All > sorted_dataDESeq2All_tmp ; join sorted_adj_pvalue_tmp sorted_dataDESeq2All_tmp | tr ' ' '\t' > dataDESeq2Filtered && rm sorted_dataDESeq2All_tmp sorted_adj_pvalue_tmp",sep = ""))

sink(output_log, append=TRUE, split=TRUE) 
print(paste(Sys.time(),"Pvalue and rest of data merged"))
sink()

  #CREATE THE FINAL HEADER USING ADJ_PVALUE AND DATADESeq2ALL ONES
system(paste("paste header_adj_pvalue.txt header_dataDESeq2All.txt > final_header.tmp ; cat final_header.tmp dataDESeq2Filtered > ",output_diff_counts," && rm final_header.tmp", sep = ""))

end_of_analysis=Sys.time()
sink(output_log, append=TRUE, split=TRUE)
print(paste(Sys.time()," Analysis done in :", difftime(end_of_analysis-start_of_analysis),sep="")
sink()
  #REMOVE THE CHUNKS FILE
system(paste("rm -rf",output_tmp))
