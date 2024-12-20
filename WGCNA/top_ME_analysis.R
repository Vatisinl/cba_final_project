library(DESeq2)
library(tidyverse)
library(magrittr)
library(WGCNA) 
library(data.table)
library(dplyr)
library(genefilter) 
allowWGCNAThreads() 
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db) # Annotation database for human genes
library(enrichplot)   # For visualization
library(pathview)     # For KEGG pathway visualization
library(biomaRt)
library(ggpubr)
library(ComplexHeatmap)

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

#Powers are based on sft plots constructed by sft_creation.R
powers <- data.frame(tissue = c('Blood', 'Lung', 'Brain - Hippocampus', 'Brain - Cerebellum', 'Brain - Caudate (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Cortex', 'Kidney', 'Liver', 'Heart', 'Muscle'), 
                     power = c(8, 5, 5, 5, 4, 3, 4, 8, 4, 7, 6))

# Function to make fancy heatmap
make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = meta_dds,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.

  # Set up the module eigengene with its refinebio_accession_code
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("accession_code")
    
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(accession_code, new_AGE, individual) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "accession_code") %>%
    # Arrange by patient and time point
    dplyr::arrange(new_AGE, individual) %>%
    # Store sample
    tibble::column_to_rownames("accession_code")
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    time_point = col_annot_df$new_AGE,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in time_point
    col = list(new_AGE = c("20-39" = "#f1a340", "40-59" = "#998ec3"))
  )

  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)

  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()

  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()

  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )

  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
    name = module_name,
    # Supply color function
    col = color_func,
    # Supply column annotation
    bottom_annotation = col_annot,
    # We don't want to cluster samples
    cluster_columns = FALSE,
    # We don't need to show sample or gene labels
    show_row_names = FALSE,
    show_column_names = FALSE
  )

  # Return heatmap
  return(heatmap)
}

# Read matrisome genes list
matrisome <- read.csv('matrisome_all.tsv', sep = '\t')

# We will use it to annotate
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


for(i in c('Heart', 'Liver', 'Lung')) {
    meta <- metadata_article[metadata_article$SMTS %in% c(i),]

    print(paste0("Computations for ", i, ' have started'))

    normalized_counts <- readRDS(file = paste0("normalized_counts_object_", i, ".rds"))
    
        # take only sample names
    df <- data.frame(samples = c(meta$SAMPID))
    df <- df %>%
      mutate(individual = sub("^.*-(.*?)-.*", "\\1", samples))

    # construct meta_dds
    meta_dds <- merge(x = df, y = pheno, by='individual')
    meta_dds <- merge(x = meta_dds, y = metadata[c('SAMPID', 'SMTSD')], by.x='samples', by.y='SAMPID')
    rownames(meta_dds) <- meta_dds$samples
    meta_dds <- meta_dds %>% dplyr::select(-SUBJID, -samples)
    meta_dds <- meta_dds[sort(rownames(meta_dds)), ]  
    
    netwk <- readRDS(file = paste0("netwk_object_", i, ".rds"))

        # Convert labels to colors for plotting
    mergedColors = labels2colors(netwk$colors)
    
    png(paste0("dendro_plot_", i, ".png"), width = 700, height = 700)
    # Plot the dendrogram and the module colors underneath
    plotDendroAndColors(netwk$dendrograms[[1]],
        mergedColors[netwk$blockGenes[[1]]],
        "Module colors",
        dendroLabels = FALSE,
        hang = 0.03,
        addGuide = TRUE,
        guideHang = 0.05)
    dev.off()

    module_df <- data.frame(
        gene_id = names(netwk$colors),
        colors = labels2colors(netwk$colors))

    # Get Module Eigengenes per cluster
    MEs0 <- moduleEigengenes(normalized_counts, mergedColors)$eigengenes
    
    # Reorder modules so similar modules are next to each other
    MEs0 <- orderMEs(MEs0)
    module_order = names(MEs0) %>% gsub("ME","", .)
    
    # Add names
    MEs0$sample = row.names(MEs0)
    
    # tidy & plot data
    mME = MEs0 %>%
      pivot_longer(-sample) %>% 
      mutate(
        name = gsub("ME", "", name),
        name = factor(name, levels = module_order)
      )
    
    modules <- mME %>% ggplot(., aes(x=sample, y=name, fill=value)) 
      geom_tile() +
      theme_bw() +
      scale_fill_gradient2(
        low = "blue",
        high = "red",
        mid = "white",
        midpoint = 0,
        limit = c(-1,1)) +
      theme(axis.text.x = element_text(angle=90)) +
      labs(title = "Module-sample Relationships", y = "Modules", fill="corr")
    ggsave(paste0("modules_plot_", i, ".png"), plot = modules)

    module_eigengenes <- netwk$MEs
    moduleLabels <- netwk$colors
    
    # Create the design matrix from the `new_AGE` variable
    des_mat <- model.matrix(~ meta_dds$new_AGE)
    
    # lmFit() needs a transposed version of the matrix
    fit <- limma::lmFit(t(module_eigengenes), design = des_mat)
    # Apply empirical Bayes to smooth standard errors
    fit <- limma::eBayes(fit)
    # Apply multiple testing correction and obtain stats
    stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
      tibble::rownames_to_column("module")

    module_eigengenes$new_AGE <- meta_dds$new_AGE
    module_eigengenes$new_AGE <- as.factor(module_eigengenes$new_AGE)
    
    # Plot to compare samples of different ages: are they significant?
    significance <- ggplot(
      module_eigengenes,
      aes(
        x = new_AGE,
        y = .data[[stats_df[1,1]]],
        color = new_AGE
      )
    ) +
      # Boxplot with outlier points hidden
      geom_boxplot(width = 0.2, outlier.shape = NA) +
      # Sina plot to show individual data points
      ggforce::geom_sina(maxwidth = 0.3) +
      # Add pairwise comparisons using the Wilcoxon rank-sum test
      stat_compare_means(
        comparisons = list(c("20–39", "60–79"), c("40–59", "60–79"), c("20–39", "40–59")),  # specify the groups to compare
        method = "wilcox.test",  # specify the Wilcoxon rank-sum test
        label = "p.signif",  # label with significance level
        p.adjust.method = "BH"  # Adjust p-values for multiple comparisons
      ) +
      theme_classic()
    ggsave(paste0("significance_plot_", i, ".png"), plot = significance)

    meta_dds$accession_code <- rownames(meta_dds)
    module <- module_eigengenes %>%
      tibble::rownames_to_column("accession_code") %>%
      # Here we are performing an inner join with a subset of metadata
      dplyr::inner_join(meta_dds %>%
        dplyr::select(accession_code, new_AGE),
      by = c("accession_code" = "accession_code")
      )
    gene_module_key <- tibble::enframe(netwk$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
      dplyr::mutate(module = paste0("ME", module))

    # Take only the most significant module
    top_me <- gene_module_key %>%
  dplyr::filter(module == stats_df[1,1])
    # To remove gene versions
    top_me <- top_me %>%
      mutate(gene_wo_version = sub("\\..*", "\\1", gene))

    top_mod_heatmap <- make_module_heatmap(module_name = stats_df[1,1])

    png(paste0("top_module_heatmap_plot_", i, ".png"), width = 700, height = 700)
    ComplexHeatmap::draw(top_mod_heatmap)
    dev.off()

    # Retrieve corresponding gene names to genes from the most significant module
    annotated_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description"), 
                         filters = "ensembl_gene_id", 
                         values = top_me$gene_wo_version, 
                         mart = ensembl)

    top_me <- merge(x=top_me, y=annotated_genes, by.x='gene_wo_version', by.y='ensembl_gene_id')
    write.csv(top_me, paste0("top_module_genes_", i, ".csv"), row.names = FALSE)

    write.csv(intersect(x=top_me$hgnc_symbol, y=matrisome$Gene), paste0("interseсection_with_matrisome", i, ".csv"), row.names = FALSE)

    # GO enrichment
    ego <- enrichGO(gene = top_me$gene_wo_version,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05)
 
            # Check if the enrichGO result is empty before proceeding
    if (is.null(ego) || nrow(ego@result) == 0 ) {
        print(paste0("No GO enrichment results found. Skipping plotting for ", i))
    } else {
        if (nrow(ego@result[ego@result$p.adjust < 0.05, ]) == 0) {
            print("No significant GO terms found.")
        } else {
            # Create dot plot if results are available
            go_plot <- dotplot(ego, showCategory = 5)
        
            # Save plot as PNG
            png(paste0("enrichGO_plot_", i, ".png"), width = 700, height = 500)
            print(go_plot)
            dev.off()
            }

    }
        
    entrez_ids <- mapIds(org.Hs.eg.db, keys = top_me$gene_wo_version, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

    # Perform KEGG Enrichment
    kegg_results <- enrichKEGG(entrez_ids, organism = 'hsa')

    # The same test: if empty
        if (is.null(kegg_results) || nrow(kegg_results@result) == 0 ) {
        print(paste0("No KEGG enrichment results found. Skipping plotting for ", i))
    } else {
        if (nrow(kegg_results@result[kegg_results@result$p.adjust < 0.05, ]) == 0) {
            print("No significant KEGG terms found.")
        } else {
            # Create dot plot if results are available
            kegg_plot <- dotplot(kegg_results, showCategory = 5)
        
            # Save plot as PNG
            png(paste0("kegg_plot_", i, ".png"), width = 700, height = 500)
            print(kegg_plot)
            dev.off()
            }

    }
    
}
