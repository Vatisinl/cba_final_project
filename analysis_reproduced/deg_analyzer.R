#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)

gene_metadata_cols = c('Description', 'Source')

for (num_grs in c("3grs")) {
  for (sex in c('male', 'fem')) {
    deg_upreg_names_all <- data.frame(matrix(nrow=0, ncol=length(gene_metadata_cols)))
    colnames(deg_upreg_names_all) <- gene_metadata_cols
    deg_downreg_names_all <- data.frame(matrix(nrow=0, ncol=length(gene_metadata_cols)))
    colnames(deg_downreg_names_all) <- gene_metadata_cols
    
    search_pattern <- paste("degs_", num_grs, "_(blood|heart|kidney|brain|muscle|liver|lung)_", sex, ".csv", sep="")
    print(search_pattern)
    deg_files_grs <- list.files(pattern = search_pattern)
    
    for (test_filename in deg_files_grs) {
      #print(test_filename)
      degs_test <- read.delim(file=test_filename, sep = ',')
      degs_test <- degs_test[degs_test$gene_type != "Other", ]
    #  degs_test <- degs_test[degs_test$logCPM > 0 & degs_test$FDR <= 0.05, ]
      degs_test <- degs_test[degs_test$logCPM > 2, ]
      if(nrow(degs_test) == 0) {
        next
      }
      
      degs_test$Source <- test_filename
      degs_test_upreg <- degs_test[degs_test$logFC > 0, ]
      degs_test_downreg <- degs_test[degs_test$logFC < 0, ]
      
      deg_upreg_names_all <- as.data.frame(rbind(deg_upreg_names_all, degs_test_upreg[degs_test_upreg$gene_type == "ECM",gene_metadata_cols]))
      deg_downreg_names_all <- as.data.frame(rbind(deg_downreg_names_all, degs_test_downreg[degs_test_downreg$gene_type == "ECM",gene_metadata_cols]))
      
      save_filename_upreg <- paste("upreg", test_filename, sep = "_")
      save_filename_downreg <- paste("downreg", test_filename, sep = "_")
      
      print(paste("./StringDB/", num_grs, "/", save_filename_upreg, sep = ""))
      print(paste("./StringDB/", num_grs, "/", save_filename_downreg, sep = ""))
     # write.csv(degs_test_upreg, paste("./StringDB/", num_grs, "/", save_filename_upreg, sep = ""), row.names = FALSE)
    #  write.csv(degs_test_downreg, paste("./StringDB/", num_grs, "/", save_filename_downreg, sep = ""), row.names = FALSE)
   }
    
    agg_tbl_upreg <- as.data.frame(deg_upreg_names_all %>% group_by(Description) %>% 
      summarise(num_tissues = n()))
    agg_tbl_upreg <- agg_tbl_upreg[agg_tbl_upreg$num_tissues >= 3, ]
    agg_tbl_downreg <- as.data.frame(deg_downreg_names_all %>% group_by(Description) %>% 
                                     summarise(num_tissues = n()))
    agg_tbl_downreg <- agg_tbl_downreg[agg_tbl_downreg$num_tissues >= 3, ]
    
    #deg_upreg_names_all <- deg_upreg_names_all[!duplicated(deg_upreg_names_all), ]
    #deg_downreg_names_all <- deg_downreg_names_all[!duplicated(deg_downreg_names_all), ]
    print(paste("./StringDB/", num_grs, "/upreg_degs_", num_grs, "_", sex, ".csv", sep = ""))
    print(paste("./StringDB/", num_grs, "/downreg_degs_", num_grs, "_", sex, ".csv", sep = ""))
    #write.csv(deg_upreg_names_all, paste("./StringDB/", num_grs, "/upreg_degs_", num_grs, "_", sex, ".csv", sep = ""), row.names = FALSE)
    #write.csv(deg_downreg_names_all, paste("./StringDB/", num_grs, "/downreg_degs_", num_grs, "_", sex, ".csv", sep = ""), row.names = FALSE)
    write.csv(agg_tbl_upreg, paste("./StringDB/", num_grs, "/common_upreg_degs_", num_grs, "_", sex, ".csv", sep = ""), row.names = FALSE)
    write.csv(agg_tbl_downreg, paste("./StringDB/", num_grs, "/common_downreg_degs_", num_grs, "_", sex, ".csv", sep = ""), row.names = FALSE)
  }
}

common_upreg_male <- read.delim(file="./StringDB/3grs/common_upreg_degs_3grs_male.csv", sep = ',')
common_downreg_male <- read.delim(file="./StringDB/3grs/common_downreg_degs_3grs_male.csv", sep = ',')
common_degs_male <- rbind(common_upreg_male, common_downreg_male)

heatmap_data <- c()

for (sex in c('male', 'fem')) {
  common_upreg_sex <- read.delim(file=paste("./StringDB/", num_grs, "/common_upreg_degs_", num_grs, "_", sex, ".csv", sep = ""), sep = ',')
  common_downreg_sex <- read.delim(file=paste("./StringDB/", num_grs, "/common_downreg_degs_", num_grs, "_", sex, ".csv", sep = ""), sep = ',')
  common_degs_sex <- rbind(common_upreg_sex, common_downreg_sex)
  
  search_pattern <- paste("degs_", num_grs, "_(blood|heart|kidney|brain|muscle|liver|lung)_", sex, ".csv", sep="")
  print(search_pattern)
  deg_files_grs <- list.files(pattern = search_pattern)
  combined_df_sex <- data.frame(matrix(nrow=0, ncol=0))
  tissues <- c()
  
  for (test_filename in deg_files_grs) {
    cur_tissue <- str_extract(test_filename, "blood|heart|kidney|brain|muscle|liver|lung")
    tissues[length(tissues) + 1] <- cur_tissue
    test_df <- read.delim(file=test_filename, sep = ',')
    test_df <- test_df[test_df$Description %in% common_degs_sex$Description,]
    
    if(nrow(combined_df_sex) == 0) {
      combined_df_sex <- test_df[, c('Description', 'logFC')]
    } else {
      combined_df_sex <- combined_df_sex %>%
        # select(Description, logFC) %>%
        #  rename(C1 = logFC) %>%
        full_join(test_df %>% select(Description, logFC), by = "Description")
    }
  }
  colnames(combined_df_sex) <- c('Description', tissues)
  combined_df_sex[is.na(combined_df_sex)] <- 0
  
  heatmap_data[[length(heatmap_data)+1]] <- combined_df_sex
}

plotHeatmap = function(dat, plot_title = "") {
   heatmap_data <- as.matrix(dat[ , -1])  # Remove the gene column
   rownames(heatmap_data) <- dat$Description
   melted_data <- melt(heatmap_data)
  
  ggplot(melted_data, aes(Var2, Var1, fill=value)) + 
    ggtitle(plot_title) +
    geom_tile() + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Conditions", y = "Genes", fill = "Log2 Fold Change") +
    theme_minimal()
}

combined_male <- as.data.frame(heatmap_data[1])
combined_fem <- as.data.frame(heatmap_data[2])

plotHeatmap(combined_male, "Male Tissues")
plotHeatmap(combined_fem, "Female Tissues")

common_degs_all <- intersect(combined_male$Description, combined_fem$Description)

######################################################################################
####################################################################################

#test_df_plot <- test_df[test_df$gene_subtype == "Collagen", ]
#EnhancedVolcano(test_df_plot,
 #               lab = test_df_plot$Description,
#                x = 'logFC',
#                y = 'PValue',
#                ylim = c(0, 6),
 #               pCutoff = 0.05,
#                FCcutoff = 0)

#EnhancedVolcano(test_df,
 #                         lab = test_df$Description,
                        #  title = 'nTreg: Healthy vs XLA',
  #                        subtitle = "",
   #                       x = 'logFC',
    #                      y = 'log_fdr',
     #                     ylim = c(0, 5),
                          # xlim = c(-7, 6),
      #                    drawConnectors = TRUE,
                         # selectLab = c('HSPA6', 'LCN2', 'TNFSF13B', 'S100A11', 'IL1B', 'ENSG00000283579', 'BLVRB', 'CDA', 'ASGR1'),
       #                   FCcutoff = 0,
        #                  pCutoff = 0.05,
         #                 caption='',
          #                max.overlaps = Inf,
           #               borderWidth =0.5) + theme(panel.grid.major = element_line(size=0.5), axis.ticks=element_line(size=0.5))

