library(edgeR)
library(ggplot2)

metadat <- read.delim(file="~/Comp_Biol_of_Aging/Project/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
unique(metadat$SMTS)           
subj_metadat <- read.delim(file="~/Comp_Biol_of_Aging/Project/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

#dat.gct <- read.delim(file="~/Comp_Biol_of_Aging/Project/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz", skip=2)
dat.gct <- read.delim(file="~/Comp_Biol_of_Aging/Project/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", skip=2)

#########################################################################
#########################################################################
stratifyAge = function(age) {
  strata <- "other"
  if (age == "20-29" | age == "30-39") {
    strata <- "20-39"
  } else if (age == "40-49" | age == "50-59") {
    strata <- "40-59"
  } else if (age == "60-69" | age == "70-79") {
    strata <- "60-79"
  }
  strata
}

########################################################
analyze = function(gct_df, my_metadata, ecm_genes, matr_assoc_genes) {
  y <- DGEList(counts=gct_df, samples = my_metadata)
  # sum_counts <- apply(y_blood_m$counts, 2, sum)
  
  #Only keep the genes with counts > min_num_counts in at least min_num_samp samples (i.e., exclude the lowly expressed genes)
  keep <- filterByExpr(y, group=y$samples$age_strata)
  y <- y[keep,]
  
  y$samples$lib.size <- colSums(y$counts)  #Reset the library sizes
  y <- calcNormFactors(y, method="TMM")
  
  groups <- factor(y$samples$age_strata, levels = c("20-39", "40-59", "60-79"))
  design <- model.matrix(~groups)
  colnames(design) <- levels(groups)
  
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  
  qlf <- glmQLFTest(fit)
  topTags(qlf)
  result <- as.data.frame(topTags(qlf, n=nrow(y)))
  degs <- result[result$PValue<0.05,]
  
  degs <- merge(degs, ecm_genes, by.x = 'Description', by.y = 'Gene', all.x=T)
  degs <- merge(degs, matr_assoc_genes, by.x = 'Description', by.y = 'Gene', all.x=T)
  degs$gene_subtype <- ifelse(!is.na(degs$gene_subtype.x), degs$gene_subtype.x, ifelse(!is.na(degs$gene_subtype.y), degs$gene_subtype.y, "Undefined"))
  degs$gene_type <- ifelse(!is.na(degs$gene_subtype.x), "ECM", ifelse(!is.na(degs$gene_subtype.y), "Matr_Assoc", "Other"))
  degs <- degs[, -which(names(degs) %in% c("gene_subtype.x", "gene_subtype.y"))]
  
  degs <- degs[degs$logCPM > 0 & degs$gene_type != "Other", ] 
  
  res_list <- list()
  res_list$dge_obj <- y
  res_list$degs <- degs
  res_list
}

####################################################################
myPCA <- function(dge, my_metadata, tissue, degs_only=FALSE) {
  # Here we will use the log2 counts per million (logCPM)
  logCPM <- cpm(dge, log = TRUE)
  
  # Perform PCA using the prcomp function
  pca_results <- prcomp(t(logCPM), center = TRUE, scale. = TRUE)
  # Create a data frame with PCA results
  pca_data <- data.frame(pca_results$x)
  
  pca_data <- merge(pca_data, my_metadata, by.x = 0, by.y = "SAMPID")
  plot_title <- ifelse(my_metadata$SEX[1] == "M", 'Male', 'Female')

  # Plotting the PCA
  ggplot(pca_data, aes(x = PC1, y = PC2, color = age_strata)) +
    geom_point(size = 3) +
    labs(title = paste("PCA of Gene Expression Data for", plot_title, tissue, ifelse(degs_only, "(DEGs only)", "")),
         x = "Principal Component 1",
         y = "Principal Component 2") +
    theme_minimal() +
    scale_color_manual(values = c("blue", "red", "green"))
}
############################################################################
############################################################################

subj_metadat$age_strata <- sapply(subj_metadat$AGE, stratifyAge)

metadat_article <- metadat[metadat$SMTS %in% c('Heart', 'Kidney', 'Liver', 'Lung', 'Muscle') | metadat$SMTSD %in% c('Whole Blood', 'Brain - Hippocampus', 'Brain - Cerebellum', 'Brain - Caudate (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Cortex'),]
metadat_article <- metadat_article[metadat_article['SMAFRZE'] == 'RNASEQ', ]
#  metadat_article['SMAFRZE'] == 'WES' | metadat_article['SMAFRZE'] == 'WGS'
metadat_article <- metadat_article[metadat_article$SMNABTCHT!='RNA isolation_PAXgene Tissue miRNA', ]
metadat_article$donor <- unlist(lapply(strsplit(metadat_article$SAMPID, '-'), FUN=function(x){paste(x[1],x[2],sep="-")}))


as.matrix(sort(table(metadat_article['SMTSD']), decreasing=TRUE))

###############################################
###############################################
subsetMetadata = function(category, metadata_article, subj_metadata) {
  metadat_categ <- metadata_article[metadata_article$SMTS == category,]
 # metadat_categ <- metadata_article
  metadat_categ$SAMPID = gsub("-",".", metadat_categ$SAMPID)
  metadat_categ <- merge(metadat_categ[,c('SAMPID', 'SMTS', 'donor')], subj_metadata[,c('SUBJID', 'SEX', 'age_strata')], by.x = 'donor', by.y = 'SUBJID')
  metadat_categ$SEX <- ifelse(metadat_categ$SEX=="1", "M", "F")
  metadat_categ
}

metadat_full <- subsetMetadata('abc', metadat_article, subj_metadat)

subsetData = function(dat_gct, metadat_categ) {
  dat_gct_categ <- dat_gct[, colnames(dat_gct) %in% c('Description', 'Name') | colnames(dat_gct) %in% metadat_categ$SAMPID]
  dat_gct_categ
}
#####################################################################
#####################################################################
col_genes <- read.delim(file="~/Comp_Biol_of_Aging/Project/Collagen_genes.tsv", sep = '\t', header = T)
ecm_glyc_genes <- read.delim(file="~/Comp_Biol_of_Aging/Project/ECM_glycoproteins_genes.tsv", sep = '\t', header = T)
ecm_proteoglyc_genes <- read.delim(file="~/Comp_Biol_of_Aging/Project/ECM_proteoglyc_genes.tsv", sep = '\t', header = T)
ecm_affil_genes <- read.delim(file="~/Comp_Biol_of_Aging/Project/ECM-affiliated_genes.tsv", sep = '\t', header = T)
secreted_factor_genes <- read.delim(file="~/Comp_Biol_of_Aging/Project/Secreted_factor_genes.tsv", sep = '\t', header = T)
ecm_regulator_genes <- read.delim(file="~/Comp_Biol_of_Aging/Project/ECM-regulator_genes.tsv", sep = '\t', header = T)

col_genes <- col_genes[,c('Gene', 'Tissue')]
col_genes <- col_genes[!duplicated(col_genes), ]
col_genes$gene_subtype <- "Collagen"

ecm_glyc_genes <- ecm_glyc_genes[,c('Gene', 'Tissue')]
ecm_glyc_genes <- ecm_glyc_genes[!duplicated(ecm_glyc_genes), ]
ecm_glyc_genes$gene_subtype <- "ECM_glycoprotein"

ecm_proteoglyc_genes <- ecm_proteoglyc_genes[,c('Gene', 'Tissue')]
ecm_proteoglyc_genes <- ecm_proteoglyc_genes[!duplicated(ecm_proteoglyc_genes), ]
ecm_proteoglyc_genes$gene_subtype <- "Proteoglycans"

ecm_affil_genes <- ecm_affil_genes[,c('Gene', 'Tissue')]
ecm_affil_genes <- ecm_affil_genes[!duplicated(ecm_affil_genes), ]
ecm_affil_genes$gene_subtype <- "ECM affiliated"

secreted_factor_genes <- secreted_factor_genes[,c('Gene', 'Tissue')]
secreted_factor_genes <- secreted_factor_genes[!duplicated(secreted_factor_genes), ]
secreted_factor_genes$gene_subtype <- "Secreted factors"

ecm_regulator_genes <- ecm_regulator_genes[,c('Gene', 'Tissue')]
ecm_regulator_genes <- ecm_regulator_genes[!duplicated(ecm_regulator_genes), ]
ecm_regulator_genes$gene_subtype <- "ECM regulators"

ecm_genes <- rbind(col_genes, ecm_glyc_genes, ecm_proteoglyc_genes)
ecm_genes$Gene <- toupper(ecm_genes$Gene)
ecm_genes <- ecm_genes[!duplicated(ecm_genes), ]
ecm_genes_dedup <- ecm_genes[,c('Gene', 'gene_subtype')]
ecm_genes_dedup <- ecm_genes_dedup[!duplicated(ecm_genes_dedup),]

matr_associated_genes <- rbind(ecm_affil_genes, secreted_factor_genes, ecm_regulator_genes)
matr_associated_genes$Gene <- toupper(matr_associated_genes$Gene)
matr_associated_genes <- matr_associated_genes[!duplicated(matr_associated_genes), ]
matr_associated_genes_dedup <- matr_associated_genes[,c('Gene', 'gene_subtype')]
matr_associated_genes_dedup <- matr_associated_genes_dedup[!duplicated(matr_associated_genes_dedup),]


######################################################################
######################################################################
metadat_blood <- subsetMetadata('Blood', metadat_article, subj_metadat)
#metadat_blood <- metadat_blood[metadat_blood$age_strata != "40-59",]
table(metadat_blood$SEX)
############################################################################################
metadat_blood_m <- metadat_blood[metadat_blood$SEX == 'M',]
dat.gct_blood_m <- subsetData(dat.gct, metadat_blood_m)

res_m <- analyze(dat.gct_blood_m, metadat_blood_m, ecm_genes_dedup, matr_associated_genes_dedup)
deg_blood_m <- res_m$degs
y_blood_m <- res_m$dge_obj
##################################################################
##################################################################
metadat_blood_f <- metadat_blood[metadat_blood$SEX == 'F',]
dat.gct_blood_f <- subsetData(dat.gct, metadat_blood_f)

res_f <- analyze(dat.gct_blood_f, metadat_blood_f, ecm_genes_dedup, matr_associated_genes_dedup)
deg_blood_f <- res_f$degs
y_blood_f <- res_f$dge_obj
################################################################
y_blood_m_red <- y_blood_m[y_blood_m$genes$Name %in% deg_blood_m$Name, ]
y_blood_f_red <- y_blood_f[y_blood_f$genes$Name %in% deg_blood_f$Name, ]
myPCA(y_blood_m, metadat_blood_m, tissue = "Blood")
myPCA(y_blood_m_red, metadat_blood_m, tissue = "Blood", degs_only = TRUE)
myPCA(y_blood_f, metadat_blood_f, tissue = "Blood")
myPCA(y_blood_f_red, metadat_blood_f, tissue = "Blood", degs_only = TRUE)
#######################################################################
#######################################################################
write.csv(deg_blood_m, "./Results/degs_3grs_blood_male.csv", row.names = FALSE)
write.csv(deg_blood_f, "./Results/degs_3grs_blood_fem.csv", row.names = FALSE)
########################################################################
########################################################################
metadat_heart <- subsetMetadata('Heart', metadat_article, subj_metadat)
#metadat_heart <- metadat_heart[metadat_heart$age_strata != "40-59",]
table(metadat_heart$SEX)
########################################################################
metadat_heart_m <- metadat_heart[metadat_heart$SEX == 'M',]
dat.gct_heart_m <- subsetData(dat.gct, metadat_heart_m)
res_m <- analyze(dat.gct_heart_m, metadat_heart_m, ecm_genes_dedup, matr_associated_genes_dedup)
deg_heart_m <- res_m$degs
deg_heart_m <- deg_heart_m[deg_heart_m$logCPM > 0 & deg_heart_m$gene_type != "Other", ]
y_heart_m <- res_m$dge_obj
########################################################################
metadat_heart_f <- metadat_heart[metadat_heart$SEX == 'F',]
dat.gct_heart_f <- subsetData(dat.gct, metadat_heart_f)
res_f <- analyze(dat.gct_heart_f, metadat_heart_f, ecm_genes_dedup, matr_associated_genes_dedup)
deg_heart_f <- res_f$degs
y_heart_f <- res_f$dge_obj
################################################################
y_heart_m_red <- y_heart_m[y_heart_m$genes$Name %in% deg_heart_m$Name, ]
y_heart_f_red <- y_heart_f[y_heart_f$genes$Name %in% deg_heart_f$Name, ]
myPCA(y_heart_m, metadat_heart_m, tissue = "Heart")
myPCA(y_heart_m_red, metadat_heart_m, tissue = "Heart", degs_only = TRUE)
myPCA(y_heart_f, metadat_heart_f, tissue = "Heart")
myPCA(y_heart_f_red, metadat_heart_f, tissue = "Heart", degs_only = TRUE)
#####################################################################
write.csv(deg_heart_m,"./Results/degs_3grs_heart_male.csv", row.names = FALSE)
write.csv(deg_heart_f,"./Results/degs_3grs_heart_fem.csv", row.names = FALSE)
########################################################################
metadat_kidney <- subsetMetadata('Kidney', metadat_article, subj_metadat)
#metadat_kidney <- metadat_kidney[metadat_kidney$age_strata != "40-59",]
table(metadat_kidney$SEX)
########################################################################
metadat_kidney_m <- metadat_kidney[metadat_kidney$SEX == 'M',]
dat.gct_kidney_m <- subsetData(dat.gct, metadat_kidney_m)
res_m <- analyze(dat.gct_kidney_m, metadat_kidney_m, ecm_genes_dedup, matr_associated_genes_dedup)
deg_kidney_m <- res_m$degs
y_kidney_m <- res_m$dge_obj
deg_kidney_m <- deg_kidney_m[deg_kidney_m$logCPM > 0 & deg_kidney_m$gene_type != "Other", ]
########################################################################
metadat_kidney_f <- metadat_kidney[metadat_kidney$SEX == 'F',]
dat.gct_kidney_f <- subsetData(dat.gct, metadat_kidney_f)
res_f <- analyze(dat.gct_kidney_f, metadat_kidney_f, ecm_genes_dedup, matr_associated_genes_dedup)
deg_kidney_f <- res_f$degs
y_kidney_f <- res_f$dge_obj
deg_kidney_f <- deg_kidney_f[deg_kidney_f$logCPM > 0 & deg_kidney_f$gene_type != "Other", ]
################################################################
y_kidney_m_red <- y_kidney_m[y_kidney_m$genes$Name %in% deg_kidney_m$Name, ]
y_kidney_f_red <- y_kidney_f[y_kidney_f$genes$Name %in% deg_kidney_f$Name, ]
myPCA(y_kidney_m, metadat_kidney_m, tissue = "Kidney")
myPCA(y_kidney_m_red, metadat_kidney_m, tissue = "Kidney", degs_only = TRUE)
myPCA(y_kidney_f, metadat_kidney_f, tissue = "Kidney")
myPCA(y_kidney_f_red, metadat_kidney_f, tissue = "Kidney", degs_only = TRUE)
################################################################
write.csv(deg_kidney_m,"./Results/degs_3grs_kidney_male.csv", row.names = FALSE)
write.csv(deg_kidney_f,"./Results/degs_3grs_kidney_fem.csv", row.names = FALSE)
#####################################################################
########################################################################
metadat_lung <- subsetMetadata('Lung', metadat_article, subj_metadat)
#metadat_lung <- metadat_lung[metadat_lung$age_strata != "40-59",]
table(metadat_lung$SEX)
########################################################################
metadat_lung_m <- metadat_lung[metadat_lung$SEX == 'M',]
dat.gct_lung_m <- subsetData(dat.gct, metadat_lung_m)
res_m <- analyze(dat.gct_lung_m, metadat_lung_m, ecm_genes_dedup, matr_associated_genes_dedup)
deg_lung_m <- res_m$degs
y_lung_m <- res_m$dge_obj
deg_lung_m <- deg_lung_m[deg_lung_m$logCPM > 0 & deg_lung_m$gene_type != "Other", ]
########################################################################
metadat_lung_f <- metadat_lung[metadat_lung$SEX == 'F',]
dat.gct_lung_f <- subsetData(dat.gct, metadat_lung_f)
res_f <- analyze(dat.gct_lung_f, metadat_lung_f, ecm_genes_dedup, matr_associated_genes_dedup)
deg_lung_f <- res_f$degs
y_lung_f <- res_f$dge_obj
deg_lung_f <- deg_lung_f[deg_lung_f$logCPM > 0 & deg_lung_f$gene_type != "Other", ]
################################################################
y_lung_m_red <- y_lung_m[y_lung_m$genes$Name %in% deg_lung_m$Name, ]
y_lung_f_red <- y_lung_f[y_lung_f$genes$Name %in% deg_lung_f$Name, ]
myPCA(y_lung_m, metadat_lung_m, tissue = "Lung")
myPCA(y_lung_m_red, metadat_lung_m, tissue = "Lung", degs_only = TRUE)
myPCA(y_lung_f, metadat_lung_f, tissue = "Lung")
myPCA(y_lung_f_red, metadat_lung_f, tissue = "Lung", degs_only = TRUE)
##################################################################
write.csv(deg_lung_m,"./Results/degs_3grs_lung_male.csv", row.names = FALSE)
write.csv(deg_lung_f,"./Results/degs_3grs_lung_fem.csv", row.names = FALSE)
##################################################################
##################################################################
metadat_liver <- subsetMetadata('Liver', metadat_article, subj_metadat)
#metadat_liver <- metadat_liver[metadat_liver$age_strata != "40-59",]
table(metadat_liver$SEX)
########################################################################
metadat_liver_m <- metadat_liver[metadat_liver$SEX == 'M',]
dat.gct_liver_m <- subsetData(dat.gct, metadat_liver_m)
res_m <- analyze(dat.gct_liver_m, metadat_liver_m, ecm_genes_dedup, matr_associated_genes_dedup)
deg_liver_m <- res_m$degs
y_liver_m <- res_m$dge_obj
deg_liver_m <- deg_liver_m[deg_liver_m$logCPM > 0 & deg_liver_m$gene_type != "Other", ]
########################################################################
metadat_liver_f <- metadat_liver[metadat_liver$SEX == 'F',]
dat.gct_liver_f <- subsetData(dat.gct, metadat_liver_f)
res_f <- analyze(dat.gct_liver_f, metadat_liver_f, ecm_genes_dedup, matr_associated_genes_dedup)
deg_liver_f <- res_f$degs
y_liver_f <- res_f$dge_obj
deg_liver_f <- deg_liver_f[deg_liver_f$logCPM > 0 & deg_liver_f$gene_type != "Other", ]
################################################################
y_liver_m_red <- y_liver_m[y_liver_m$genes$Name %in% deg_liver_m$Name, ]
y_liver_f_red <- y_liver_f[y_liver_f$genes$Name %in% deg_liver_f$Name, ]
myPCA(y_liver_m, metadat_liver_m, tissue = "Liver")
myPCA(y_liver_m_red, metadat_liver_m, tissue = "Liver", degs_only = TRUE)
myPCA(y_liver_f, metadat_liver_f, tissue = "Liver")
myPCA(y_liver_f_red, metadat_liver_f, tissue = "Liver", degs_only = TRUE)
################################################################
write.csv(deg_liver_m, "./Results/degs_3grs_liver_male.csv", row.names = FALSE)
write.csv(deg_liver_f, "./Results/degs_3grs_liver_fem.csv", row.names = FALSE)
########################################################################
########################################################################
metadat_brain <- subsetMetadata('Brain', metadat_article, subj_metadat)
#metadat_brain <- metadat_brain[metadat_brain$age_strata != "40-59",]
table(metadat_brain$SEX)
########################################################################
metadat_brain_m <- metadat_brain[metadat_brain$SEX == 'M',]
dat.gct_brain_m <- subsetData(dat.gct, metadat_brain_m)
res_m <- analyze(dat.gct_brain_m, metadat_brain_m, ecm_genes_dedup, matr_associated_genes_dedup)
deg_brain_m <- res_m$degs
y_brain_m <- res_m$dge_obj
deg_brain_m <- deg_brain_m[deg_brain_m$logCPM > 0 & deg_brain_m$gene_type != "Other", ]
########################################################################
metadat_brain_f <- metadat_brain[metadat_brain$SEX == 'F',]
dat.gct_brain_f <- subsetData(dat.gct, metadat_brain_f)
res_f <- analyze(dat.gct_brain_f, metadat_brain_f, ecm_genes_dedup, matr_associated_genes_dedup)
deg_brain_f <- res_f$degs
y_brain_f <- res_f$dge_obj
deg_brain_f <- deg_brain_f[deg_brain_f$logCPM > 0 & deg_brain_f$gene_type != "Other", ]
################################################################
y_brain_m_red <- y_brain_m[y_brain_m$genes$Name %in% deg_brain_m$Name, ]
y_brain_f_red <- y_brain_f[y_brain_f$genes$Name %in% deg_brain_f$Name, ]
myPCA(y_brain_m, metadat_brain_m, tissue = "Brain")
myPCA(y_brain_m_red, metadat_brain_m, tissue = "Brain", degs_only = TRUE)
myPCA(y_brain_f, metadat_brain_f, tissue = "Brain")
myPCA(y_brain_f_red, metadat_brain_f, tissue = "Brain", degs_only = TRUE)
################################################################
write.csv(deg_brain_m, "./Results/degs_3grs_brain_male.csv", row.names = FALSE)
write.csv(deg_brain_f, "./Results/degs_3grs_brain_fem.csv", row.names = FALSE)
########################################################################
########################################################################
metadat_muscle <- subsetMetadata('Muscle', metadat_article, subj_metadat)
#metadat_muscle <- metadat_muscle[metadat_muscle$age_strata != "40-59", ]
table(metadat_muscle$SEX)
########################################################################
metadat_muscle_m <- metadat_muscle[metadat_muscle$SEX == 'M',]
dat.gct_muscle_m <- subsetData(dat.gct, metadat_muscle_m)
res_m <- analyze(dat.gct_muscle_m, metadat_muscle_m, ecm_genes_dedup, matr_associated_genes_dedup)
deg_muscle_m <- res_m$degs
y_muscle_m <- res_m$dge_obj
deg_muscle_m <- deg_muscle_m[deg_muscle_m$logCPM > 0 & deg_muscle_m$gene_type != "Other", ]
########################################################################
metadat_muscle_f <- metadat_muscle[metadat_muscle$SEX == 'F',]
dat.gct_muscle_f <- subsetData(dat.gct, metadat_muscle_f)
res_f <- analyze(dat.gct_muscle_f, metadat_muscle_f, ecm_genes_dedup, matr_associated_genes_dedup)
deg_muscle_f <- res_f$degs
y_muscle_f <- res_f$dge_obj
deg_muscle_f <- deg_muscle_f[deg_muscle_f$logCPM > 0 & deg_muscle_f$gene_type != "Other", ]
################################################################
y_muscle_m_red <- y_muscle_m[y_muscle_m$genes$Name %in% deg_muscle_m$Name, ]
y_muscle_f_red <- y_muscle_f[y_muscle_f$genes$Name %in% deg_muscle_f$Name, ]
myPCA(y_muscle_m, metadat_muscle_m, tissue = "Muscle")
myPCA(y_muscle_m_red, metadat_muscle_m, tissue = "Muscle", degs_only = TRUE)
myPCA(y_muscle_f, metadat_muscle_f, tissue = "Muscle")
myPCA(y_muscle_f_red, metadat_muscle_f, tissue = "Muscle", degs_only = TRUE)
################################################################
write.csv(deg_muscle_m, "./Results/degs_3grs_muscle_male.csv", row.names = FALSE)
write.csv(deg_muscle_f, "./Results/degs_3grs_muscle_fem.csv", row.names = FALSE)
########################################################################
########################################################################
