# Dekupl_DESeq2
Dekupl alternative version with DESeq2 or Ttest method in order to get differentially expressed kmers not seen in GENCODE data.
## How to configure
Configure the different features within the **config.json** file
  - ### About normalization factors
    - data_compute_norm_factors : Choose between **raw-counts** or **noGENCODE**
    - seed : Seed in order to compute the normalization factors
    - sampling_size : Number of kmers used for the sample

  - ### About differential expression analysis
    - diff_method : Choose between **Ttest** or **DESeq2** methods
    - chunk_size : Capped to 1 000 000 kmers, it's the size of each parts the noGENCODE file will be splitted
