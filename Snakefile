#!/bin/python3
import os
import gzip
from snakemake.utils import R

__author__ = "Jérôme Audoux (jerome.audoux@inserm.fr)"


def getRAM():
    mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
    mem_gbytes = mem_bytes/(1024.**3)
    ramMax = mem_gbytes/2
    ramMax = int(ramMax)
    if ramMax > 20:
        ramMax = 20
    return ramMax

configfile: "/data/config.json"

# COMMON VARIABLES
SAMPLE_NAMES    = [i['name'] for i in config["samples"]]
CONDITION_COL   = "condition"
CONDITION_A     = config['Ttest']['condition']['A']
CONDITION_B     = config['Ttest']['condition']['B']
LIB_TYPE        = "rf"
MAX_CPU         = 1000
MAX_RAM         = getRAM()
R1_SUFFIX       = config['r1_suffix']
R2_SUFFIX       = config['r2_suffix']

if 'lib_type' in config:
  LIB_TYPE = config['lib_type']

# DIRECTORIES
BIN_DIR         = "/bin/dekupl/bin"
TMP_DIR         = "/data/" + config['tmp_dir']+"/dekupl_tmp"
FASTQ_DIR       = "/data/" + config['fastq_dir']
GENE_EXP_DIR    = "/data/gene_expression"
KALLISTO_DIR    = GENE_EXP_DIR + "/kallisto"
COUNTS_DIR      = "/data/kmer_counts"
KMER_DE_DIR     = "/data/" + CONDITION_A + "_vs_" + CONDITION_B + "_kmer_counts"
METADATA_DIR    = "/data/metadata"
REFERENCE_DIR   = "/data/references"

# FILES
RAW_COUNTS                  = COUNTS_DIR    + "/raw-counts.tsv.gz"
NO_GENCODE_COUNTS           = COUNTS_DIR    + "/noGENCODE-counts.tsv.gz"
DIFF_COUNTS                 = KMER_DE_DIR   + "/diff-counts.tsv.gz"
PVALUE_ALL                  = KMER_DE_DIR   + "/pvalue_all.tsv
MERGED_DIFF_COUNTS          = KMER_DE_DIR   + "/merged-diff-counts.tsv.gz"
ASSEMBLIES_FASTA            = KMER_DE_DIR   + "/merged-diff-counts.fa.gz"
ASSEMBLIES_BAM              = KMER_DE_DIR   + "/merged-diff-counts.bam"
SAMPLE_CONDITIONS           = METADATA_DIR  + "/sample_conditions.tsv"
SAMPLE_CONDITIONS_FULL      = METADATA_DIR  + "/sample_conditions_full.tsv"
GENOME_FASTA                = REFERENCE_DIR + "/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GSNAP_INDEX_DIR             = REFERENCE_DIR + "/GSNAP"
GSNAP_INDEX_NAME            = "Homo_sapiens.GRCh38.dna.primary_assembly"
GENCODE_FASTA               = REFERENCE_DIR + "/gencode.v24.transcripts.fa.gz"
GENCODE_COUNTS              = REFERENCE_DIR + "/gencode.v24.transcripts.tsv.gz"
TRANSCRIPT_TO_GENE_MAPPING  = REFERENCE_DIR + "/transcript_to_gene_mapping.tsv"
KALLISTO_INDEX              = REFERENCE_DIR + "/gencode.v24.transcripts-kallisto.idx"
TRANSCRIPT_COUNTS           = KALLISTO_DIR  + "/transcript_counts.tsv.gz"
GENE_COUNTS                 = KALLISTO_DIR  + "/gene_counts.tsv.gz"
DEGS                        = GENE_EXP_DIR  + "/" + CONDITION_A + "vs" + CONDITION_B + "-DEGs.tsv"
NORMALIZATION_FACTORS       = GENE_EXP_DIR  + "/normalization_factors.tsv"
NORMALIZED_COUNTS	        = GENE_EXP_DIR  + "/normalized_counts.tsv"
DIST_MATRIX                 = GENE_EXP_DIR  + "/clustering_of_samples.pdf"
PCA_DESIGN                  = GENE_EXP_DIR  + "/pca_design.tsv"
# binaries
REVCOMP         = "/bin/revCompFastq.pl"
DEKUPL_COUNTER  = BIN_DIR + "/dekupl-counter"
DIFF_FILTER     = BIN_DIR + "/diffFilter.pl"
TTEST_FILTER    = BIN_DIR + "/TtestFilter"
JOIN_COUNTS     = BIN_DIR + "/joinCounts"
MERGE_COUNTS    = BIN_DIR + "/mergeCounts.pl"
MERGE_TAGS      = BIN_DIR + "/mergeTags"
PIGZ            = "pigz"

# docker
KALLISTO        = "docker run --rm -v /data:/data -w /data ebio/kallisto kallisto"
JELLYFISH       = "docker run --rm -v /data:/data -w /data ebio/jellyfish"
JELLYFISH_COUNT = JELLYFISH + " /bin/bash -c \'jellyfish count"
JELLYFISH_DUMP  = JELLYFISH + " /bin/bash -c \' jellyfish dump"
GSNAP           = "docker run --rm -v /data:/data -w /data -i ebio/gsnap"
SAMTOOLS        = "docker run --rm -v /data:/data -w /data -i itsjeffreyy/samtools"

# plotter
PLOT_CHECKINGS  = BIN_DIR + "/PlotCheckings.R"
CHECKING_PLOTS  = KMER_DE_DIR   + "/fragment_number_per_library_raw.pdf"
LOGS            = "/data/Logs"

rule run:
  input: ASSEMBLIES_BAM

###############################################################################
#
# SOFTWARE INSTALLATION
#
rule compile_joinCounts:
  output: JOIN_COUNTS
  run:
    shell("cd share/joinCounts && make")
    shell("ln -s -f ../share/joinCounts/joinCounts bin/")

rule compile_mergeTags:
  output: MERGE_TAGS
  input: "share/mergeTags/mergeTags.c"
  run:
    shell("cd share/mergeTags && make")
    shell("ln -s -f ../share/mergeTags/mergeTags bin/")

rule compile_TtestFilter:
  input: "share/TtestFilter/TtestFilter.c"
  output: TTEST_FILTER
  run:
    shell("cd share/TtestFilter && make")
    shell("ln -s -f ../share/TtestFilter/TtestFilter bin/")

rule download_kallisto:
  output:
    kallisto_symlink = KALLISTO,
    kallisto_tarball = temp("share/kallisto.tar.gz")
  run:
    shell("wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz -O {output.kallisto_tarball}")
    shell("tar -xzf {output.kallisto_tarball} -C share")
    shell("ln -s ../share/kallisto_linux-v0.43.0/kallisto bin/kallisto")

###############################################################################
#
# DOWNLOAD REFERENCE FILES
#
# Download the gencode transcripts in fasta format
rule gencode_download:
  output: GENCODE_FASTA
  shell: "wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.transcripts.fa.gz -O {output}"

# Download the genome from Ensembl
rule genome_download:
  output: GENOME_FASTA
  shell: "wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O {output}"
  #shell: "wget ftp://ftp.ensembl.org/pub/release-86/fasta/ciona_intestinalis/dna/Ciona_intestinalis.KH.dna.chromosome.10.fa.gz -O {output}"

###############################################################################
#
# BUILD INDEXES FROM REFERENCE FILES
#
# Create a Kallisto index of the reference transrciptome
rule kallisto_index:
  input:
    transcripts   = GENCODE_FASTA
  output:
    KALLISTO_INDEX
  
  log: 
    exec_time= LOGS + "/kallisto_index_exec_time.log"  
      
  shell: """
        echo -e \"******\" >{log.exec_time}
        echo -e \"start of rule kallisto_index\" >>{log.exec_time}
        
        {KALLISTO} index -i {output} {input.transcripts}
        
        echo -e \"end of rule gencode_counts\" >>{log.exec_time} 
        echo -e \"******\" >>{log.exec_time}
        """
rule gsnap_index:
  input: GENOME_FASTA
  output: GSNAP_INDEX_DIR + "/" + GSNAP_INDEX_NAME
  shell: "{GSNAP} gmap_build -D {GSNAP_INDEX_DIR} -d {GSNAP_INDEX_NAME} --gunzip {input}"

###############################################################################
#
# UTILS
#
# Create a tabulated file with the sample name and conditions
rule sample_conditions:
  output: SAMPLE_CONDITIONS
  run:
    with open(output[0], "w") as f:
      f.write("\t".join(["sample",CONDITION_COL]) + "\n")
      for sample in config["samples"]:
        f.write("\t".join([sample["name"],sample[CONDITION_COL]]) + "\n")

rule sample_conditions_full:
  input:
    sample_conditions      = SAMPLE_CONDITIONS,
    normalization_factors  = NORMALIZATION_FACTORS
  output:
    SAMPLE_CONDITIONS_FULL
  shell: "join --header {input.sample_conditions} {input.normalization_factors} > {output}"


###############################################################################
#
# STEP 1: DIFFERENTIAL GENE EXPRESSION
#         Download kallisto, and quantify gene expression for all
#         the samples

# 1.3 Generic rule to quantify a sample with kallisto
rule kallisto_quantif:
  input:
    r1 = FASTQ_DIR + "/{sample}" + R1_SUFFIX,
    r2 = FASTQ_DIR + "/{sample}" + R2_SUFFIX,
    index = KALLISTO_INDEX
  output:
    KALLISTO_DIR + "/{sample}"
    
  log:
    exec_time = LOGS + "/{sample}_kallisto.log"
    
  threads: 1
  shell: """
         echo -e \"******\" >{log}
         echo -e \"start of rule kallito_quantif\" >>{log.exec_time}

         {KALLISTO} quant -i {input.index} -o {output} {input.r1} {input.r2} 2>>{log}
         
         echo -e \"end of rule kallito_quantif\" >>{log.exec_time}
         echo -e \"******\" >>{log.exec_time}
         """

# 1.4 Merge all transcripts counts from kallisto abundance files
rule transcript_counts:
  input:
    kallisto_outputs  = expand("{kallisto_dir}/{sample}", sample = SAMPLE_NAMES, kallisto_dir = KALLISTO_DIR)
  output:
    TRANSCRIPT_COUNTS
  run:
    extracted_counts  = expand("<(echo -e \'feature\t{sample}\' && tail -n+2 {kallisto_dir}/{sample}/abundance.tsv | cut -f1,4)", sample = SAMPLE_NAMES, kallisto_dir = KALLISTO_DIR)
    shell("bash -c \"{MERGE_COUNTS} {extracted_counts} | gzip -c > {output}\"")

# 1.5 Create a conversion table from transcript id to gene ids
rule transcript_to_gene_mapping:
  input: GENCODE_FASTA
  output: TRANSCRIPT_TO_GENE_MAPPING
  run:
    mapping = open(output[0], 'w')
    with gzip.open(input[0], 'rt') as f:
      for line in f:
        if line[0] == ">":
          fields = line[1:].split("|",2)
          mapping.write("\t".join([fields[0],fields[1]]) + "\n")

# 1.6 Convert transcript counts to gene counts
rule gene_counts:
  input:
    transcript_counts = TRANSCRIPT_COUNTS,
    transcript_to_gene_mapping = TRANSCRIPT_TO_GENE_MAPPING
  output:
    GENE_COUNTS
  run:
    # Load the conversion hash
    conversion_hash = {}
    with open(input['transcript_to_gene_mapping'], 'r') as f:
      for line in f:
        transcript_id, gene_id = line.split()
        conversion_hash[transcript_id] = gene_id
    # Summarize transcript into gene counts
    gene_counts = {}
    header = ""
    with gzip.open(input['transcript_counts'], 'rt') as f:
      header = f.readline().rstrip()
      for line in f:
        counts = line.split()
        transcript_id, trail = counts[0].split("|",1)
        gene_id = conversion_hash[transcript_id]
        counts[1:] = [ float(i) for i in counts[1:] ]
        if gene_id in gene_counts:
          gene_counts[gene_id] = [ sum(i) for i in zip(gene_counts[gene_id], counts[1:]) ]
        else:
          gene_counts[gene_id] = counts[1:]
    # print Gene counts
    with gzip.open(output[0], 'wb') as f:
      f.write(bytes(header + "\n",'UTF-8'))
      for gene_id in gene_counts:
        f.write(bytes(gene_id + "\t" + "\t".join([str(int(x)) for x in gene_counts[gene_id]]) + "\n",'UTF-8'))

# 1.7 Differential expression with DESEQ2
rule differential_gene_expression:
  input:
    gene_counts = GENE_COUNTS,
    sample_conditions = SAMPLE_CONDITIONS
  output:
    differentially_expressed_genes = DEGS,
    normalization_factors          = NORMALIZATION_FACTORS,
    dist_matrix			   = DIST_MATRIX,
    norm_counts		           = NORMALIZED_COUNTS,
    pca_design = PCA_DESIGN
    
  log : LOGS + "/DESeq2_diff_gene_exp.log"
  shell:
    """
    printf "
    library(DESeq2)
    library(RColorBrewer)
    library(pheatmap)
    library(ggplot2)
    
    write(date(),file=\\\"{log}\\\")

    # Load counts data
    countsData = read.table(\\\"{input.gene_counts}\\\",header=T,row.names=1)
                                                                    
    # Load col data with sample specifications
    colData = read.table(\\\"{input.sample_conditions}\\\",header=T,row.names=1)
    
    write(colnames(countsData),stderr())
    write(rownames(colData),stderr())

    colData = colData[colnames(countsData),,drop=FALSE]

    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData=countsData,colData=colData,design = ~ {CONDITION_COL})
    dds <- DESeq(dds)
    
    #normalized counts 
    NormCount<- as.data.frame(counts(dds, normalized=TRUE ))
  
    #writing in a file normalized counts
    normalized_counts<-data.frame(id=row.names(NormCount),NormCount,row.names=NULL)
    write.table(normalized_counts,file=\\\"{output.norm_counts}\\\", sep=\\\"\t\\\",row.names=F, col.names=T, quote=F)

    # Write normalization factors
    size_factors = data.frame(sample = names(sizeFactors(dds)), 
                              normalization_factor = sizeFactors(dds),
                              row.names=NULL)
    write.table(size_factors,
                file=\\\"{output.normalization_factors}\\\",
                sep=\\\"\t\\\",quote=FALSE, row.names = FALSE)

    write(resultsNames(dds),stderr())

    # Write DEGs
    res <- results(dds, contrast = c(\\\"{CONDITION_COL}\\\",\\\"{CONDITION_A}\\\",\\\"{CONDITION_B}\\\"))
    write.table(res,file=\\\"{output.differentially_expressed_genes}\\\",sep=\\\"\t\\\",quote=FALSE)


    rld<-rlog(dds)
    sampleDists<-dist(t(assay(rld) ) )
    sampleDistMatrix<-as.matrix( sampleDists )
    rownames(sampleDistMatrix)<-colnames(rld)
    colnames(sampleDistMatrix)<-colnames(rld)
    colours=colorRampPalette(rev(brewer.pal(9,\\\"Blues\\\")) )(255)
    
    pdf(\\\"{output.dist_matrix}\\\",width=15,height=10)
    
    pheatmap(sampleDistMatrix,
    main=\\\"clustering of samples\\\",
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colours,
    fontsize = 14)
     
    data<-plotPCA(rld,ntop=nrow(rld),returnData=TRUE)
    
    write.table(data,\\\"{output.pca_design}\\\",row.names=F, col.names=T, quote=F,sep=\\\"\t\\\")
    
    print(ggplot(data,aes(PC1,PC2,color=condition))+geom_point()+geom_text(aes(label=name),hjust=0,vjust=0))
    
    dev.off()
    
    write(date(),file=\\\"{log}\\\",append=T) " > tmp.txt

    Rscript tmp.txt
    """

###############################################################################
#
# STEP 2: KMER COUNTS
#         Compile DEkupl counter and count k-mers on all the samples
#
rule jellyfish_count:
  input:
    r1 = FASTQ_DIR + "/{sample}" + R1_SUFFIX,
    r2 = FASTQ_DIR + "/{sample}" + R2_SUFFIX
  output: COUNTS_DIR + "/{sample}.jf"
  threads: MAX_CPU
  run:
    if LIB_TYPE == 'rf':
      shell("{JELLYFISH_COUNT} -L 2 -m {config[kmer_length]} -s 10000 -t {threads} -o {output} -F 2 <(zcat {input.r1} | {REVCOMP}) <(zcat {input.r2})\'")
    elif LIB_TYPE == 'fr':
      shell("{JELLYFISH_COUNT} -L 2 -m {config[kmer_length]} -s 10000 -t {threads} -o {output} -F 2 <(zcat {input.r1}) <(zcat {input.r2} | {REVCOMP})\'")
    else:
      sys.exit('Unknown library type')

rule jellyfish_dump:
  input: COUNTS_DIR + "/{sample}.jf"
  output: COUNTS_DIR + "/{sample}.txt.gz"
  threads: MAX_CPU
  resources: ram=MAX_RAM
  shell: "{JELLYFISH_DUMP} -c {input} | sort -k 1 -S {resources.ram}G --parallel {threads}| {PIGZ} -p {threads} -c > {output}\'"

rule join_counts:
  input:
    fastq_files = expand("{counts_dir}/{sample}.txt.gz",counts_dir=COUNTS_DIR,sample=SAMPLE_NAMES),
    binary = JOIN_COUNTS
  params:
    sample_names = "\t".join(SAMPLE_NAMES)
  output: RAW_COUNTS
  run:
    shell("echo 'tag\t{params.sample_names}' | gzip -c > {output}")
    shell("""{JOIN_COUNTS} -r {config[dekupl_counter][min_recurrence]} \
          -a {config[dekupl_counter][min_recurrence_abundance]} \
          {input.fastq_files} | gzip -c >> {output}""")

###############################################################################
#
# STEP 3: FILTER-OUT KNOWN K-MERS
#         Download gencode transcripts set and remove the k-mer occuring this
#         set from the one found in the experimental data
#

# 3.2 Counts k-mer of all gencode transcript (for further filtration)
rule gencode_count:
  input: GENCODE_FASTA
  output: temp(GENCODE_FASTA + ".jf")
  threads: MAX_CPU
  shell: """{JELLYFISH_COUNT} -m {config[kmer_length]} \
            -s 10000 -t {threads} -o {output} <(zcat {input})\'"""

rule gencode_dump:
  input: GENCODE_FASTA + ".jf"
  output: GENCODE_COUNTS
  threads: MAX_CPU
  resources: ram=4
  shell: "{JELLYFISH_DUMP} -c {input} | sort -k 1 -S {resources.ram}G --parallel {threads}| {PIGZ} -p {threads} -c > {output}\'"

# 3.3 Filter counter k-mer that are present in the gencode set
rule filter_gencode_counts:
  input:
    counts = RAW_COUNTS,
    gencode_counts = GENCODE_COUNTS
  output: NO_GENCODE_COUNTS
  shell: "{DIFF_FILTER} {input.gencode_counts} {input.counts} | gzip -c > {output}"

###############################################################################
#
# STEP 4: SELECT DIFFERENTIALLY EXPRESSED K-MERS
#         Apply a T-test on all new k-mers to select only those that are
#         differentially expressed.
#
rule test_diff_counts:
  input:
    counts = NO_GENCODE_COUNTS,
    sample_conditions = SAMPLE_CONDITIONS_FULL
  output: 
    diff_counts = DIFF_COUNTS
    pvalue_all = PVALUE_ALL
  log = LOGS
  threads:6
  script: "./script.R"

rule merge_tags:
  input:
    counts = DIFF_COUNTS,
    binary = MERGE_TAGS
  output:
    MERGED_DIFF_COUNTS
  shell: "{MERGE_TAGS} -k {config[kmer_length]} {input.counts} | gzip -c > {output}"

rule plot_checkings:
  input: 
    raw_counts = RAW_COUNTS,
    no_gencode_counts = NO_GENCODE_COUNTS,
    diff_counts = DIFF_COUNTS,
    merged_diff_counts = MERGED_DIFF_COUNTS
    
  output: CHECKING_PLOTS
  
  log: 
      exec_time = LOGS + "/plot_checkings_exec_time.log"
  
  shell: """
  
         echo -e \"******\" >{log.exec_time}
         echo -e \"start of rule gencode_counts\" >>{log.exec_time}  
         chmod 755 {KMER_DE_DIR}
     {PLOT_CHECKINGS} {input.raw_counts} {input.no_gencode_counts} {input.diff_counts} {input.merged_diff_counts} {KMER_DE_DIR}
    
         echo -e \"end of rule gencode_counts\" >>{log.exec_time} 
         echo -e \"******\" >>{log.exec_time}
    
    """
###############################################################################
#
# STEP 4: ANNOTATED DIFF K-MERS
#
rule assembly_fasta:
  input: MERGED_DIFF_COUNTS
  output: ASSEMBLIES_FASTA
  shell: """
    zcat {input} | tail -n+2 | awk '{{print ">"$3"\\n"$2}}' | gzip -c > {output}
  """

rule gsnap_align:
  input: 
    fasta = ASSEMBLIES_FASTA, 
    gsnap = GSNAP_INDEX_DIR + "/" + GSNAP_INDEX_NAME,
    deseq = DEGS
  output: ASSEMBLIES_BAM
  threads: 10
  shell: """{GSNAP} gsnap -D {GSNAP_INDEX_DIR} -d {GSNAP_INDEX_NAME} --gunzip -t {threads} -A sam -B 2 -N 1 {input.fasta} | {SAMTOOLS} bash -c "samtools view -bS - | samtools sort - " > {output} && \
  {SAMTOOLS} samtools index {output}"""
