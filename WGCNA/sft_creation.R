# Attach the DESeq2 library
library(DESeq2)
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
library(WGCNA) 
library(data.table)
library(dplyr)
library(genefilter) 
library(WGCNA)
allowWGCNAThreads() 
library(tidyverse)

# Unzip and read the GCT file
file_path <- "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"

# Use fread to read the file, skipping metadata lines
gene_data <- fread(file_path, skip = 2)

gene_data <- as.data.frame(gene_data)
rownames(gene_data) <- gene_data$Name

gene_data <- gene_data[, -c(1, 2)]  # Remove non-numeric columns with gene IDs and description of them

print(dim(gene_data))

# Use fread to read the file, skipping metadata lines
# This variable is a copy to update gene_data in the future (in the loop)
data <- fread(file_path, skip = 2)

data <- as.data.frame(data)
rownames(data) <- data$Name

data <- data[, -c(1, 2)] 

# Read metadata
metadata <- readr::read_tsv('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
pheno <- readr::read_tsv('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')

#Make new age categories
pheno <- pheno %>%
  mutate(
    new_AGE = case_when(
      AGE %in% c('20-29', '30-39') ~ '20–39',
      AGE %in% c('40-49', '50-59') ~ '40–59',
      AGE %in% c('60-69', '70-79') ~ '60–79',
      TRUE ~ 'Other' # Optional: Handle unexpected values
    )
  )

pheno <- pheno %>%
  mutate(individual = sub("^.*-(.*?)-.*", "\\1", SUBJID))

pheno$individual <- sub("^GTEX-", "", pheno$individual)
pheno <- pheno %>% select(-AGE, -DTHHRDY)

# metadata_article <- metadata[metadata$SMTS %in% c('Blood', 'Heart', 'Kidney', 'Liver', 'Lung', 'Muscle') | metadata$SMTSD %in% c('Brain - Hippocampus', 'Brain - Cerebellum', 'Brain - Caudate (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Cortex'),]

metadata_article <- metadata[metadata$SMTS %in% c('Heart', 'Kidney', 'Liver', 'Lung', 'Muscle') | metadata$SMTSD %in% c('Whole Blood', 'Brain - Hippocampus', 'Brain - Cerebellum', 'Brain - Caudate (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Cortex'),]

metadata_article <- metadata_article[metadata_article['SMAFRZE']=='RNASEQ', ]

for(i in c('Heart', 'Liver', 'Lung')) { 
    meta <- metadata_article[metadata_article$SMTS %in% c(i),]
    
    # take only sample names
    df <- data.frame(samples = c(meta$SAMPID))
    df <- df %>%
      mutate(individual = sub("^.*-(.*?)-.*", "\\1", samples))
    
    # construct meta_dds
    meta_dds <- merge(x = df, y = pheno, by='individual')
    meta_dds <- merge(x = meta_dds, y = metadata[c('SAMPID', 'SMTS')], by.x='samples', by.y='SAMPID')
    rownames(meta_dds) <- meta_dds$samples
    meta_dds <- meta_dds %>% select(-SUBJID, -samples)
    
    gene_data <- data
    # technical thing
    gene_data <- as.matrix(gene_data)
    
    #Filter into gene_data only samples from meta_dds
    # Identify selected columns based on matching names
    selected_columns <- which(colnames(gene_data) %in% rownames(meta_dds))
    
    # Subset gene_data using the selected column indices
    gene_data <- gene_data[, selected_columns]
    
    # Sort to make the order the same
    gene_data <- gene_data[, sort(colnames(gene_data))]
    meta_dds <- meta_dds[sort(rownames(meta_dds)), ]   
    
    # Filter: each sample contains >70 counts of all genes in sum
    gene_data <- gene_data[rowSums(gene_data[])>70,] 
    
    # Create a DESeq obj
    dds <- DESeqDataSetFromMatrix(round(gene_data),
                                  meta_dds,
                                  design = ~1)
    # Perform DESeq normalization (variance not stabilized yet)
    dds <- DESeq(dds)
    
    # Apply variance stabilizing transformation
    vsd <- vst(dds) # blind = FALSE
    
    # Extract the transformed data
    vst_data <- assay(vsd)
    normalized_counts <- t(vst_data)
    
    sft <- pickSoftThreshold(normalized_counts,
      dataIsExpr = TRUE,
      corFnc = cor,
      corOptions = list(method = 'spearman'),
      networkType = "signed"
    )
    
    sft_df <- data.frame(sft$fitIndices) %>%
      dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)
    
    saveRDS(sft, file = paste0("sft_object_", i, ".rds"))

    plot <- ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

    ggsave(paste0("plot_", i, ".png"), plot = plot)
    
}