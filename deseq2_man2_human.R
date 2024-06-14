# BiocManager::install("ComplexHeatmap")

library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(gprofiler2)
library(GO.db)
library(igraph)
library(cowplot)
library(RColorBrewer)
library(RSkittleBrewer)
library(stringr)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(gridExtra)
library(grid)
library(ggpubr)
library(WGCNA)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)


library(stats)
library(genefilter)
library(tibble)
library(ggrepel)
library(pathview)
library(AnnotationDbi)
library(biomaRt)
library(scales)
library(enrichR)
library(GOSim)
library(foreach)
library(ggridges)
library(forcats)



setwd("/mmfs1/projects/changhui.yan/DeewanB/VOC_SRR/deseq2_human")
getwd()

# Load the saved R session
#load(file = "deseq2_human.RData")
load(file = "deseq2_human_wo_wuhan.RData")

##############################################################################
# Defining a function to remove outliers using IQR
remove_outliers <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  x <- ifelse(x < lower_bound | x > upper_bound, NA, x)
  return(x)
}

##############################################################################
# Defining filtering function which turns outliers into NAs to be removed
filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}


##############################################################################
# Defining custom Function to get GO terms for a gene
# Function to get GO terms for a gene
get_GO_terms <- function(gene_name) {
  tryCatch({
    # Convert gene name to ENTREZID using org.Hs.eg.db
    entrez_id <- mapIds(org.Hs.eg.db, gene_name, 
                        fromType = "SYMBOL", toType = "ENTREZID", 
                        column = "ENTREZID", keytype = "SYMBOL")
    
    if (is.na(entrez_id)) {
      return(NA)
    } else {
      # Get the GO terms associated with the ENTREZID using org.Hs.eg.db
      gene_info <- select(org.Hs.eg.db, keys = entrez_id, columns = "GO")
      if (length(gene_info$GO) == 0) {
        return(NA)
      } else {
        return(paste(gene_info$GO, collapse = "; "))
      }
    }
  }, error = function(e) {
    cat("Error occurred for gene:", gene_name, "\n")
    return(NA) # Return NA when an error occurs
  })
}
# Example showing how to apply get_GO_terms function to each gene in human_stringDB_genes$interacting_geneNames
#all_genes_host_control_immune_trial$GO <- sapply(all_genes_host_control_immune_trial$geneNames, get_GO_terms)






###################################################
# Custom function to determine the max value for a given GO term
get_max_Ycoordinate <- function(dataframe, go_term) {
  values <- as.numeric(dataframe$value[dataframe$GO_term == go_term])
  
  # Function to filter outliers
  filter_lims <- function(x) {
    l <- boxplot.stats(x)$stats[1]
    u <- boxplot.stats(x)$stats[5]
    
    for (i in 1:length(x)){
      x[i] <- ifelse(x[i] > l & x[i] < u, x[i], NA)
    }
    return(x)
  }
  
  filtered_values <- filter_lims(values)
  
  max_value <- max(filtered_values, na.rm = TRUE)
  return(max_value + (max_value * 0.1))
}


###################################################
# Custom function to determine the max value for a given GeneName
get_max_Ycoordinate_geneName <- function(dataframe, geneName) {
  max_value <- max(filter_lims(
    as.numeric(
      dataframe$value[dataframe$variable == geneName]
    )
  ), na.rm = TRUE)
  return(max_value + (max_value * 0.1))
}


###################################################
# Custom function to split x or y axis tick label into two lines if longer than 30 characters
custom_labeller <- function(labels) {
  labels <- lapply(labels, function(x) {
    if (nchar(x) > 20) {
      paste0(strwrap(x, width = 55), collapse = "\n")
    } else {
      x
    }
  })
  return(labels)
}



###################################################
# Function to average columns with shared names while preserving row names
average_columns_by_variant <- function(input_matrix) {
  # Identify unique column names
  unique_cols <- unique(colnames(input_matrix))
  
  # Initialize an empty matrix to store averaged results
  averaged_matrix <- matrix(NA, nrow = nrow(input_matrix), ncol = length(unique_cols))
  colnames(averaged_matrix) <- unique_cols
  rownames(averaged_matrix) <- rownames(input_matrix)  # Preserve row names
  
  # Loop through unique column names
  for (i in seq_along(unique_cols)) {
    colname <- unique_cols[i]
    # Find indices of columns with the same name
    selected_cols <- which(colnames(input_matrix) == colname)
    # Compute mean of selected columns
    col_mean <- rowMeans(input_matrix[, selected_cols, drop = FALSE])
    # Assign the mean values to the corresponding column in the result matrix
    averaged_matrix[, i] <- col_mean
  }
  
  # Return the averaged matrix
  return(averaged_matrix)
}


###################################################
# function to get GO terms for clustered genes:
get_main_GO_categories <- function(gene_ids) {
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(gene = gene_ids, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  
  # Extract and return the top GO terms
  top_GO_terms <- go_enrichment[1:10, c("ID", "Description")]
  return(top_GO_terms)
}

get_main_KEGG_categories <- function(gene_ids) {
  # Perform KEGG pathway enrichment analysis
  kegg_enrichment <- enrichKEGG(gene = gene_ids, organism = "hsa")
  
  # Extract and return the top KEGG pathways
  top_KEGG_pathways <- kegg_enrichment[1:10, c("ID", "Description")]
  return(top_KEGG_pathways)
}






##############################################################################################################################
##############################################################################################################################
################################################## HUMAN WITH HEALTHY CONTROLS ###############################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################



metadata_ncbi <- read.csv("SRP427025_metadata_from_ncbi.csv")
pheno_data_host <- metadata_ncbi[,c("Run","virus","batch","Timepoint")]
pheno_data_host <- pheno_data_host %>%
  filter(virus != "SA")
colnames(pheno_data_host) <- c("Run","variant","batch","Timepoint")


# Load the read count matrix
count_matrix <- read.csv("gene_count_matrix_human.csv", header=TRUE, row.names=1)
# Converting colnames from {geneName}|{geneNames} into {geneNames}
rownames(count_matrix) <- sub("^[^|]+\\|([^|]+)$", "\\1", rownames(count_matrix))

# Subset pheno_data_host based on the matched column names
pheno_data_host <- pheno_data_host %>% filter(pheno_data_host$Run %in% 
                                                intersect(pheno_data_host$Run, colnames(count_matrix)))


# Subsetting count matrix to only contain the SRR runs that were selected in metadata file
count_matrix <- count_matrix[, intersect(colnames(count_matrix), pheno_data_host$Run)]


count_matrix[is.na(count_matrix)] <- 0

# Design Matrix 
design_matrix <- as.data.frame(cbind(pheno_data_host$Run,
                                     pheno_data_host$variant,
                                     pheno_data_host$batch,
                                     pheno_data_host$Timepoint))
design_matrix %<>% {rownames(.) <- .$V1; .[, -1]} %>% {colnames(.) <- c("variant","batch","Timepoint"); .}

design_matrix_out <- t(design_matrix)
write.csv(design_matrix_out,"design_matrix_human.csv")
write.csv(count_matrix,"corrected_count_matrix_human.csv")

#################################################

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix) %in% rownames(design_matrix))

# Reorder the columns of count_matrix to match the row names of design_matrix
count_matrix <- count_matrix[, rownames(design_matrix)]
# order check
all(colnames(count_matrix) == rownames(design_matrix))


# Converting all read counts into integers
count_matrix[] <- lapply(count_matrix, as.integer)

design_matrix$variant <- factor(design_matrix$variant, levels = c("mock","Wuhan","delta","omicron"))
design_matrix$batch <- factor(design_matrix$batch)
design_matrix$Timepoint <- factor(design_matrix$Timepoint)



# Converting all read counts into integers
count_matrix_all <- count_matrix
count_matrix_all[] <- lapply(count_matrix, as.integer)

design_matrix$variant <- factor(design_matrix$variant, levels = c("mock","Wuhan","delta","omicron"))
design_matrix$group <- paste0(design_matrix$variant,"_",design_matrix$Timepoint)

dds_variant_time_all <- DESeqDataSetFromMatrix(countData = count_matrix_all,
                                             colData = design_matrix,
                                             design = ~ variant)


count_matrix_output <- count_matrix
# Update column names of the filtered matrix with the 'variant' values
colnames(count_matrix_output) <- design_matrix[colnames(count_matrix_output), "group"]

write.csv(count_matrix_output,"corrected_count_matrix_human_wt_group.csv")
write.csv(design_matrix,"design_matrix.csv")


######################################################################################################
##################################### DEG at each hpi#################################################
######################################################################################################

######### Expression at each timepoint:######### 

########### 1 hours ############
design_matrix_1 <- filter(design_matrix, Timepoint %in% c(1))

count_matrix_1 <- count_matrix[colnames(count_matrix) %in% rownames(design_matrix_1)]

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix_1) %in% rownames(design_matrix_1))

# Reorder the columns of count_matrix_1 to match the row names of design_matrix_1
count_matrix_1 <- count_matrix_1[, rownames(design_matrix_1)]
# order check
all(colnames(count_matrix_1) == rownames(design_matrix_1))

# Converting all read counts into integers
count_matrix_1[] <- lapply(count_matrix_1, as.integer)

design_matrix_1$variant <- factor(design_matrix_1$variant, levels = c("mock","Wuhan","delta","omicron"))

dds_variant_time_1 <- DESeqDataSetFromMatrix(countData = count_matrix_1,
                                             colData = design_matrix_1,
                                             design = ~ variant)

# set the factor level
dds_variant_time_1$variant <- relevel(dds_variant_time_1$variant, ref = "mock")
dds_variant_time_1$variant


# DESeq on 1h
deseq_variant_time_1 <- DESeq(dds_variant_time_1)


# Get differential expression results
result_wuhan_1 <- results(deseq_variant_time_1, contrast = c("variant","Wuhan", "mock"))
summary(result_wuhan_1)
head(result_wuhan_1)
result_wuhan_1 <- as.data.frame(result_wuhan_1)
result_wuhan_1 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_delta_1 <- results(deseq_variant_time_1, contrast = c("variant","delta","mock"))
summary(result_delta_1)
head(result_delta_1)
result_delta_1 <- as.data.frame(result_delta_1)
result_delta_1 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_omicron_1 <- results(deseq_variant_time_1, contrast = c("variant","omicron","mock"))
summary(result_omicron_1)
head(result_omicron_1)
result_omicron_1 <- as.data.frame(result_omicron_1)
result_omicron_1 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_wuhan_1 <- result_wuhan_1[, c("geneName", "baseMean","log2FoldChange", "padj")] %>% mutate(group = "Wuhan_1")
df_delta_1 <- result_delta_1[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Delta_1")
df_omicron_1 <- result_omicron_1[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Omicron_1")

log2FoldChange_table_1 <- rbind(df_wuhan_1, df_delta_1, df_omicron_1)

log2FoldChange_table_1 <- log2FoldChange_table_1 %>% filter(!is.na(log2FoldChange) & !is.na(padj))



# Filter for differentially expressed genes
#log2FC (1.5) = 2

diff_genes_1 <- log2FoldChange_table_1%>% 
  filter(abs(log2FoldChange) > 1 & padj < 0.05)

# Fetch Reactome gene sets
hs_immune_genes <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") 
# Filter rows in diff_genes_human_control based on ENTREZID in hs_immune_genes$entrez_gene
#selected_entrez_ids <- hs_immune_genes$entrez_gene

diff_genes_1_immune <- diff_genes_1[diff_genes_1$geneName %in% hs_immune_genes$gene_symbol, ]






########### 2 hours ############
design_matrix_2 <- filter(design_matrix, Timepoint %in% c(2))

count_matrix_2 <- count_matrix[colnames(count_matrix) %in% rownames(design_matrix_2)]

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix_2) %in% rownames(design_matrix_2))

# Reorder the columns of count_matrix_2 to match the row names of design_matrix_2
count_matrix_2 <- count_matrix_2[, rownames(design_matrix_2)]
# order check
all(colnames(count_matrix_2) == rownames(design_matrix_2))

# Converting all read counts into integers
count_matrix_2[] <- lapply(count_matrix_2, as.integer)

design_matrix_2$variant <- factor(design_matrix_2$variant, levels = c("mock","Wuhan","delta","omicron"))

dds_variant_time_2 <- DESeqDataSetFromMatrix(countData = count_matrix_2,
                                             colData = design_matrix_2,
                                             design = ~ variant)

# set the factor level
dds_variant_time_2$variant <- relevel(dds_variant_time_2$variant, ref = "mock")
dds_variant_time_2$variant


# DESeq on 1h
deseq_variant_time_2 <- DESeq(dds_variant_time_2)


# Get differential expression results
result_wuhan_2 <- results(deseq_variant_time_2, contrast = c("variant","Wuhan", "mock"))
summary(result_wuhan_2)
head(result_wuhan_2)
result_wuhan_2 <- as.data.frame(result_wuhan_2)
result_wuhan_2 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_delta_2 <- results(deseq_variant_time_2, contrast = c("variant","delta","mock"))
summary(result_delta_2)
head(result_delta_2)
result_delta_2 <- as.data.frame(result_delta_2)
result_delta_2 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_omicron_2 <- results(deseq_variant_time_2, contrast = c("variant","omicron","mock"))
summary(result_omicron_2)
head(result_omicron_2)
result_omicron_2 <- as.data.frame(result_omicron_2)
result_omicron_2 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_wuhan_2 <- result_wuhan_2[, c("geneName", "baseMean","log2FoldChange", "padj")] %>% mutate(group = "Wuhan_2")
df_delta_2 <- result_delta_2[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Delta_2")
df_omicron_2 <- result_omicron_2[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Omicron_2")

log2FoldChange_table_2 <- rbind(df_wuhan_2, df_delta_2, df_omicron_2)

log2FoldChange_table_2 <- log2FoldChange_table_2 %>% filter(!is.na(log2FoldChange) & !is.na(padj))



# Filter for differentially expressed genes
#log2FC (1.5) = 2

diff_genes_2 <- log2FoldChange_table_2%>% 
  filter(abs(log2FoldChange) > 1 & padj < 0.05)

# Fetch Reactome gene sets
hs_immune_genes <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") 
# Filter rows in diff_genes_human_control based on ENTREZID in hs_immune_genes$entrez_gene
#selected_entrez_ids <- hs_immune_genes$entrez_gene

diff_genes_2_immune <- diff_genes_2[diff_genes_2$geneName %in% hs_immune_genes$gene_symbol, ]








########### 3 hours ############
design_matrix_3 <- filter(design_matrix, Timepoint %in% c(3))

count_matrix_3 <- count_matrix[colnames(count_matrix) %in% rownames(design_matrix_3)]

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix_3) %in% rownames(design_matrix_3))

# Reorder the columns of count_matrix_3 to match the row names of design_matrix_3
count_matrix_3 <- count_matrix_3[, rownames(design_matrix_3)]
# order check
all(colnames(count_matrix_3) == rownames(design_matrix_3))

# Converting all read counts into integers
count_matrix_3[] <- lapply(count_matrix_3, as.integer)

design_matrix_3$variant <- factor(design_matrix_3$variant, levels = c("mock","Wuhan","delta","omicron"))

dds_variant_time_3 <- DESeqDataSetFromMatrix(countData = count_matrix_3,
                                             colData = design_matrix_3,
                                             design = ~ variant)

# set the factor level
dds_variant_time_3$variant <- relevel(dds_variant_time_3$variant, ref = "mock")
dds_variant_time_3$variant


# DESeq on 1h
deseq_variant_time_3 <- DESeq(dds_variant_time_3)


# Get differential expression results
result_wuhan_3 <- results(deseq_variant_time_3, contrast = c("variant","Wuhan", "mock"))
summary(result_wuhan_3)
head(result_wuhan_3)
result_wuhan_3 <- as.data.frame(result_wuhan_3)
result_wuhan_3 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_delta_3 <- results(deseq_variant_time_3, contrast = c("variant","delta","mock"))
summary(result_delta_3)
head(result_delta_3)
result_delta_3 <- as.data.frame(result_delta_3)
result_delta_3 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_omicron_3 <- results(deseq_variant_time_3, contrast = c("variant","omicron","mock"))
summary(result_omicron_3)
head(result_omicron_3)
result_omicron_3 <- as.data.frame(result_omicron_3)
result_omicron_3 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_wuhan_3 <- result_wuhan_3[, c("geneName", "baseMean","log2FoldChange", "padj")] %>% mutate(group = "Wuhan_3")
df_delta_3 <- result_delta_3[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Delta_3")
df_omicron_3 <- result_omicron_3[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Omicron_3")

log2FoldChange_table_3 <- rbind(df_wuhan_3, df_delta_3, df_omicron_3)

log2FoldChange_table_3 <- log2FoldChange_table_3 %>% filter(!is.na(log2FoldChange) & !is.na(padj))



# Filter for differentially expressed genes
#log2FC (1.5) = 2

diff_genes_3 <- log2FoldChange_table_3%>% 
  filter(abs(log2FoldChange) > 1 & padj < 0.05)

# Fetch Reactome gene sets
hs_immune_genes <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") 
# Filter rows in diff_genes_human_control based on ENTREZID in hs_immune_genes$entrez_gene
#selected_entrez_ids <- hs_immune_genes$entrez_gene

diff_genes_3_immune <- diff_genes_3[diff_genes_3$geneName %in% hs_immune_genes$gene_symbol, ]




######################################################################################################################

########### 4 hours ############
design_matrix_4 <- filter(design_matrix, Timepoint %in% c(4))

count_matrix_4 <- count_matrix[colnames(count_matrix) %in% rownames(design_matrix_4)]

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix_4) %in% rownames(design_matrix_4))

# Reorder the columns of count_matrix_4 to match the row names of design_matrix_4
count_matrix_4 <- count_matrix_4[, rownames(design_matrix_4)]
# order check
all(colnames(count_matrix_4) == rownames(design_matrix_4))

# Converting all read counts into integers
count_matrix_4[] <- lapply(count_matrix_4, as.integer)

design_matrix_4$variant <- factor(design_matrix_4$variant, levels = c("mock","Wuhan","delta","omicron"))

dds_variant_time_4 <- DESeqDataSetFromMatrix(countData = count_matrix_4,
                                             colData = design_matrix_4,
                                             design = ~ variant +batch)

# set the factor level
dds_variant_time_4$variant <- relevel(dds_variant_time_4$variant, ref = "mock")
dds_variant_time_4$variant


# DESeq on 4D
deseq_variant_time_4 <- DESeq(dds_variant_time_4)


# Get differential expression results
result_wuhan_4 <- results(deseq_variant_time_4, contrast = c("variant","Wuhan", "mock"))
summary(result_wuhan_4)
head(result_wuhan_4)
result_wuhan_4 <- as.data.frame(result_wuhan_4)
result_wuhan_4 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_delta_4 <- results(deseq_variant_time_4, contrast = c("variant","delta","mock"))
summary(result_delta_4)
head(result_delta_4)
result_delta_4 <- as.data.frame(result_delta_4)
result_delta_4 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_omicron_4 <- results(deseq_variant_time_4, contrast = c("variant","omicron","mock"))
summary(result_omicron_4)
head(result_omicron_4)
result_omicron_4 <- as.data.frame(result_omicron_4)
result_omicron_4 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_wuhan_4 <- result_wuhan_4[, c("geneName", "baseMean","log2FoldChange", "padj")] %>% mutate(group = "Wuhan_4")
df_delta_4 <- result_delta_4[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Delta_4")
df_omicron_4 <- result_omicron_4[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Omicron_4")

log2FoldChange_table_4 <- rbind(df_wuhan_4, df_delta_4, df_omicron_4)

log2FoldChange_table_4 <- log2FoldChange_table_4 %>% filter(!is.na(log2FoldChange) & !is.na(padj))



# Filter for differentially expressed genes
#log2FC (1.5) = 2

diff_genes_4 <- log2FoldChange_table_4%>% 
  filter(abs(log2FoldChange) > 2 & padj < 0.05)

# Fetch Reactome gene sets
hs_immune_genes <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") 
# Filter rows in diff_genes_human_control based on ENTREZID in hs_immune_genes$entrez_gene
#selected_entrez_ids <- hs_immune_genes$entrez_gene

diff_genes_4_immune <- diff_genes_4[diff_genes_4$geneName %in% hs_immune_genes$gene_symbol, ]





########### 5 hours ############
design_matrix_5 <- filter(design_matrix, Timepoint %in% c(5))

count_matrix_5 <- count_matrix[colnames(count_matrix) %in% rownames(design_matrix_5)]

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix_5) %in% rownames(design_matrix_5))

# Reorder the columns of count_matrix_5 to match the row names of design_matrix_5
count_matrix_5 <- count_matrix_5[, rownames(design_matrix_5)]
# order check
all(colnames(count_matrix_5) == rownames(design_matrix_5))

# Converting all read counts into integers
count_matrix_5[] <- lapply(count_matrix_5, as.integer)

design_matrix_5$variant <- factor(design_matrix_5$variant, levels = c("mock","Wuhan","delta","omicron"))

dds_variant_time_5 <- DESeqDataSetFromMatrix(countData = count_matrix_5,
                                             colData = design_matrix_5,
                                             design = ~ variant +batch)

# set the factor level
dds_variant_time_5$variant <- relevel(dds_variant_time_5$variant, ref = "mock")
dds_variant_time_5$variant


# DESeq on 1h
deseq_variant_time_5 <- DESeq(dds_variant_time_5)


# Get differential expression results
result_wuhan_5 <- results(deseq_variant_time_5, contrast = c("variant","Wuhan", "mock"))
summary(result_wuhan_5)
head(result_wuhan_5)
result_wuhan_5 <- as.data.frame(result_wuhan_5)
result_wuhan_5 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_delta_5 <- results(deseq_variant_time_5, contrast = c("variant","delta","mock"))
summary(result_delta_5)
head(result_delta_5)
result_delta_5 <- as.data.frame(result_delta_5)
result_delta_5 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_omicron_5 <- results(deseq_variant_time_5, contrast = c("variant","omicron","mock"))
summary(result_omicron_5)
head(result_omicron_5)
result_omicron_5 <- as.data.frame(result_omicron_5)
result_omicron_5 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_wuhan_5 <- result_wuhan_5[, c("geneName", "baseMean","log2FoldChange", "padj")] %>% mutate(group = "Wuhan_5")
df_delta_5 <- result_delta_5[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Delta_5")
df_omicron_5 <- result_omicron_5[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Omicron_5")

log2FoldChange_table_5 <- rbind(df_wuhan_5, df_delta_5, df_omicron_5)

log2FoldChange_table_5 <- log2FoldChange_table_5 %>% filter(!is.na(log2FoldChange) & !is.na(padj))



# Filter for differentially expressed genes
#log2FC (1.5) = 2

diff_genes_5 <- log2FoldChange_table_5%>% 
  filter(abs(log2FoldChange) > 2 & padj < 0.05)

# Fetch Reactome gene sets
hs_immune_genes <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") 
# Filter rows in diff_genes_human_control based on ENTREZID in hs_immune_genes$entrez_gene
#selected_entrez_ids <- hs_immune_genes$entrez_gene

diff_genes_5_immune <- diff_genes_5[diff_genes_5$geneName %in% hs_immune_genes$gene_symbol, ]






########### 6 hours ############
design_matrix_6 <- filter(design_matrix, Timepoint %in% c(6))

count_matrix_6 <- count_matrix[colnames(count_matrix) %in% rownames(design_matrix_6)]

# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix_6) %in% rownames(design_matrix_6))

# Reorder the columns of count_matrix_6 to match the row names of design_matrix_6
count_matrix_6 <- count_matrix_6[, rownames(design_matrix_6)]
# order check
all(colnames(count_matrix_6) == rownames(design_matrix_6))

# Converting all read counts into integers
count_matrix_6[] <- lapply(count_matrix_6, as.integer)

design_matrix_6$variant <- factor(design_matrix_6$variant, levels = c("mock","Wuhan","delta","omicron"))

dds_variant_time_6 <- DESeqDataSetFromMatrix(countData = count_matrix_6,
                                             colData = design_matrix_6,
                                             design = ~ variant +batch)

# set the factor level
dds_variant_time_6$variant <- relevel(dds_variant_time_6$variant, ref = "mock")
dds_variant_time_6$variant


# DESeq on 1h
deseq_variant_time_6 <- DESeq(dds_variant_time_6)


# Get differential expression results
result_wuhan_6 <- results(deseq_variant_time_6, contrast = c("variant","Wuhan", "mock"))
summary(result_wuhan_6)
head(result_wuhan_6)
result_wuhan_6 <- as.data.frame(result_wuhan_6)
result_wuhan_6 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_delta_6 <- results(deseq_variant_time_6, contrast = c("variant","delta","mock"))
summary(result_delta_6)
head(result_delta_6)
result_delta_6 <- as.data.frame(result_delta_6)
result_delta_6 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}


result_omicron_6 <- results(deseq_variant_time_6, contrast = c("variant","omicron","mock"))
summary(result_omicron_6)
head(result_omicron_6)
result_omicron_6 <- as.data.frame(result_omicron_6)
result_omicron_6 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_wuhan_6 <- result_wuhan_6[, c("geneName", "baseMean","log2FoldChange", "padj")] %>% mutate(group = "Wuhan_6")
df_delta_6 <- result_delta_6[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Delta_6")
df_omicron_6 <- result_omicron_6[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "Omicron_6")

log2FoldChange_table_6 <- rbind(df_wuhan_6, df_delta_6, df_omicron_6)

log2FoldChange_table_6 <- log2FoldChange_table_6 %>% filter(!is.na(log2FoldChange) & !is.na(padj))



# Filter for differentially expressed genes
#log2FC (1.5) = 2

diff_genes_6 <- log2FoldChange_table_6%>% 
  filter(abs(log2FoldChange) > 2 & padj < 0.05)

# Fetch Reactome gene sets
hs_immune_genes <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") 
# Filter rows in diff_genes_human_control based on ENTREZID in hs_immune_genes$entrez_gene
#selected_entrez_ids <- hs_immune_genes$entrez_gene

diff_genes_6_immune <- diff_genes_6[diff_genes_6$geneName %in% hs_immune_genes$gene_symbol, ]











######################################################################################################################


### 1 dpi
log2FoldChange_table_1$group <- factor(log2FoldChange_table_1$group, levels = c("Wuhan_1","Delta_1","Omicron_1"))
sig_alpha <- 0.05

day1 <- ggplot2::ggplot(log2FoldChange_table_1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5, na.rm = TRUE,
             aes(color = ifelse(log2FoldChange > 1 & -log10(padj) > -log10(sig_alpha), "red",
                                ifelse(log2FoldChange < -1 & -log10(padj) > -log10(sig_alpha), "dodgerblue", "darkgrey")))) +
  scale_color_manual(values = c("red" = "red", "dodgerblue" = "dodgerblue", "darkgrey" = "darkgrey"),
                     breaks = c("red", "dodgerblue"),
                     labels = c("Upregulated", "Downregulated"),
                     name = "Padj < 0.05 & FC > |2|") +
  theme_minimal() +
  xlim(c(-15, 15)) +
  ylim(c(-1, 6)) +
  ggtitle("log2 Fold Change for host genes in ALI cultures infected with SARS-CoV-2 VOCs compared to healthy controls at:
          1 dpi") +
  xlab("Log2 Fold Change") +
  ylab("-Log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5, size =14, face ="bold"),
        axis.text = element_text(size = 12, face = "italic"), 
        axis.title = element_text(size = 12, face = "italic"),
        strip.text = element_text(size = 12),
        panel.border = element_blank(),  # Add panel border
        strip.background = element_blank(),  # Remove strip background
        panel.spacing = unit(3, "lines"),  # Reduce panel spacing
        legend.position = "bottom") +  
  facet_wrap(~ group, nrow = 1)


### 2 dpi
log2FoldChange_table_2$group <- factor(log2FoldChange_table_2$group, levels = c("Wuhan_2","Delta_2","Omicron_2"))
sig_alpha <- 0.05

day2 <- ggplot2::ggplot(log2FoldChange_table_2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5, na.rm = TRUE,
             aes(color = ifelse(log2FoldChange > 2 & -log10(padj) > -log10(sig_alpha), "red",
                                ifelse(log2FoldChange < -1 & -log10(padj) > -log10(sig_alpha), "dodgerblue", "darkgrey")))) +
  scale_color_manual(values = c("red" = "red", "dodgerblue" = "dodgerblue", "darkgrey" = "darkgrey"),
                     breaks = c("red", "dodgerblue"),
                     labels = c("Upregulated", "Downregulated"),
                     name = "Padj < 0.05 & FC > |2|") +
  theme_minimal() +
  xlim(c(-15, 15)) +
  ylim(c(-1, 6)) +
  ggtitle("2 dpi") +
  xlab("Log2 Fold Change") +
  ylab("-Log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5, size =14, face ="bold"),
        axis.text = element_text(size = 12, face = "italic"), 
        axis.title = element_text(size = 12, face = "italic"),
        strip.text = element_text(size = 12),
        panel.border = element_blank(),  # Add panel border
        strip.background = element_blank(),  # Remove strip background
        panel.spacing = unit(3, "lines"),  # Reduce panel spacing
        legend.position = "bottom") +  
  facet_wrap(~ group, nrow = 1)



#### 3 dpi
log2FoldChange_table_3$group <- factor(log2FoldChange_table_3$group, levels = c("Wuhan_3","Delta_3","Omicron_3"))
sig_alpha <- 0.05

day3 <- ggplot2::ggplot(log2FoldChange_table_3, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5, na.rm = TRUE,
             aes(color = ifelse(log2FoldChange > 2 & -log10(padj) > -log10(sig_alpha), "red",
                                ifelse(log2FoldChange < -1 & -log10(padj) > -log10(sig_alpha), "dodgerblue", "darkgrey")))) +
  scale_color_manual(values = c("red" = "red", "dodgerblue" = "dodgerblue", "darkgrey" = "darkgrey"),
                     breaks = c("red", "dodgerblue"),
                     labels = c("Upregulated", "Downregulated"),
                     name = "Padj < 0.05 & FC > |2|") +
  theme_minimal() +
  xlim(c(-15, 15)) +
  ylim(c(-1, 6)) +
  ggtitle("3 dpi") +
  xlab("Log2 Fold Change") +
  ylab("-Log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5, size =14, face ="bold"),
        axis.text = element_text(size = 12, face = "italic"), 
        axis.title = element_text(size = 12, face = "italic"),
        strip.text = element_text(size = 12),
        panel.border = element_blank(),  # Add panel border
        strip.background = element_blank(),  # Remove strip background
        panel.spacing = unit(3, "lines"),  # Reduce panel spacing
        legend.position = "bottom") +  
  facet_wrap(~ group, nrow = 1)



volcano_dpi <-  grid.arrange(day1 + theme(legend.position = "none"),
                             day2 + theme(legend.position = "none"),
                             day3 + theme(legend.position = "bottom"),
                             nrow = 3, heights = c(1,1,1.1))



######################################################################################################################
######################################################################################################################
#################### Enrichment plots for top up and downregulated genes for each VOC  ###############################
######################################################################################################################
######################################################################################################################


# Fetch Reactome gene sets
hs_gene2name_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::mutate(gs_name = gsub("_", " ", gs_name)) %>%
  dplyr::mutate(gs_name = str_remove(gs_name, "^REACTOME ")) %>%  
  dplyr::mutate(gs_name = tolower(gs_name))

# Fetch Reactome gene sets
hs_gene2name_GO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% 
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::mutate(gs_name = gsub("_", " ", gs_name)) %>%
  dplyr::mutate(gs_name = str_remove(gs_name, "^GOBP ")) %>%  
  dplyr::mutate(gs_name = tolower(gs_name))


######################################################################################################################
######################################################################################################################
    ####################### GSEA for top up and downregulated genes for each VOC  ###############################
######################################################################################################################
######################################################################################################################
diff_genes_1 <- diff_genes_1 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_1$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))

diff_genes_1_immune <- diff_genes_1_immune %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_1_immune$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))
#######################################
# Beginning of infection: 
#######################################
### Omicron ###
diff_genes_1_Omicron <- diff_genes_1 %>%
  filter(group == "Omicron_1") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_1_Delta <- diff_genes_1 %>%
  filter(group == "Delta_1") %>%
  filter(!is.na(ENTREZID)) 

### Wuhan ###
diff_genes_1_Wuhan <- diff_genes_1 %>%
  filter(group == "Wuhan_1") %>%
  filter(!is.na(ENTREZID)) 


# Immune only
### Omicron ###
diff_genes_1_Omicron_immune <- diff_genes_1_immune %>%
  filter(group == "Omicron_1") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_1_Delta_immune <- diff_genes_1_immune %>%
  filter(group == "Delta_1") %>%
  filter(!is.na(ENTREZID)) 

### Wuhan ###
diff_genes_1_Wuhan_immune <- diff_genes_1_immune %>%
  filter(group == "Wuhan_1") %>%
  filter(!is.na(ENTREZID)) 



#############
### Wuhan ###
# Ranking genes for GSEA
result_wuhan_1 <- result_wuhan_1 %>% filter(!is.na(log2FoldChange))
result_wuhan_1 <- result_wuhan_1 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_wuhan_1$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))
gene_list_w_1 <- result_wuhan_1$log2FoldChange
names(gene_list_w_1) <- result_wuhan_1$ENTREZID
gene_list_w_1 = sort(gene_list_w_1, decreasing = TRUE)


#############
### Delta ###
# Ranking genes for GSEA
result_delta_1 <- result_delta_1 %>% filter(!is.na(log2FoldChange))
result_delta_1 <- result_delta_1 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_delta_1$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))
gene_list_d_1<- result_delta_1$log2FoldChange
names(gene_list_d_1) <- result_delta_1$ENTREZID
gene_list_d_1= sort(gene_list_d_1, decreasing = TRUE)


#############
### Omicron ###
# Ranking genes for GSEA
result_omicron_1 <- result_omicron_1 %>% filter(!is.na(log2FoldChange))
result_omicron_1 <- result_omicron_1 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_omicron_1$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))
gene_list_o_1<- result_omicron_1$log2FoldChange
names(gene_list_o_1) <- result_omicron_1$ENTREZID
gene_list_o_1= sort(gene_list_o_1, decreasing = TRUE)



# ClusterprofilR GSEA REACTOME pathway
# enricher_host_Wuhan_1 <- clusterProfiler::GSEA(gene_list_w_1,
#                                                pvalueCutoff = 0.05,
#                                                #OrgDb = org.Hs.eg.db,
#                                                pAdjustMethod = "BH",
#                                                minGSSize = 5,
#                                                #ont = "ALL",
#                                                gson = NULL,
#                                                TERM2GENE = hs_gene2name_reactome)
# 
# reactomeW_1 <- ggplot(data = enricher_host_Wuhan_1@result%>%
#                       top_n(10, rank), 
#                     mapping = aes(x = NES, y = Description, 
#                                   color = p.adjust, size = setSize / sum(setSize))) +
#   geom_point() + 
#   scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
#   labs(x = "Normalized Enrichment Score (NES)", y = "Reactome Pathways", color = "Adj. p-value", size = "Gene Ratio") +
#   theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
#         axis.text.y = element_text(size = 9, hjust = 1),
#         axis.text.x = element_text(size = 12, hjust = 1),
#         axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
#         axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
#         plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
#   scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
#   ggtitle("WuhanHu1") +
#   scale_y_discrete(labels = custom_labeller)
# reactomeW_1
# 
# # ClusterprofilR GSEA REACTOME pathway
# gse_w_1 <- clusterProfiler::GSEA(gene_list_w_1,
#                                pvalueCutoff = 0.05,
#                                pAdjustMethod = "BH",
#                                minGSSize = 5,
#                                gson = NULL,
#                                TERM2GENE = hs_gene2name_GO)
# 
# 
# GOplot_w_1 <- ggplot(data = gse_w_1@result%>%
#                      top_n(15, NES),
#                    mapping = aes(x = NES, y = Description, 
#                                  color = p.adjust, size = setSize )) +
#   geom_point() + 
#   scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
#   labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Gene Ontology Pathways", color = "Adj. p-value", size = "Gene Count") +
#   theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
#         axis.text.y = element_text(size = 9, hjust = 1),
#         axis.text.x = element_text(size = 12, hjust = 1),
#         axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
#         axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
#         plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
#   scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
#   ggtitle("WuhanHu1") +
#   scale_y_discrete(labels = custom_labeller)
# 
# GOplot_w_1




# ClusterprofilR GSEA REACTOME pathway
enricher_host_Delta_1 <- clusterProfiler::GSEA(gene_list_d_1,
                                               pvalueCutoff = 0.05,
                                               #OrgDb = org.Hs.eg.db,
                                               pAdjustMethod = "BH",
                                               minGSSize = 5,
                                               #ont = "ALL",
                                               gson = NULL,
                                               TERM2GENE = hs_gene2name_reactome)


reactomeD_1 <- ggplot(data = enricher_host_Delta_1@result%>%
                      top_n(10, rank), 
                    mapping = aes(x = NES, y = Description, 
                                  color = p.adjust, size = setSize / sum(setSize))) +
  geom_point() + 
  scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "Reactome Pathways", color = "Adj. p-value", size = "Gene Ratio") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Delta") +
  scale_y_discrete(labels = custom_labeller)

reactomeD_1




# ClusterprofilR GSEA REACTOME pathway
gse_d_1 <- clusterProfiler::GSEA(gene_list_d_1,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               minGSSize = 5,
                               gson = NULL,
                               TERM2GENE = hs_gene2name_GO)


GOplot_d_1 <- ggplot(data = gse_d_1@result%>%
                     top_n(15, NES), 
                   mapping = aes(x = NES, y = Description, 
                                 color = p.adjust, size = setSize)) +
  geom_point() + 
  scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Gene Ontology Pathways", color = "Adj. p-value", size = "Gene Count") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Delta") +
  scale_y_discrete(labels = custom_labeller)
GOplot_d_1






# ClusterprofilR GSEA REACTOME pathway
enricher_host_Omicron_1 <- clusterProfiler::GSEA(gene_list_o_1,
                                                 pvalueCutoff = 0.05,
                                                 #OrgDb = org.Hs.eg.db,
                                                 pAdjustMethod = "BH",
                                                 minGSSize = 5,
                                                 #ont = "ALL",
                                                 gson = NULL,
                                                 TERM2GENE = hs_gene2name_reactome)


reactomeO_1 <- ggplot(data = enricher_host_Omicron_1@result %>%
                      top_n(10, rank), 
                    mapping = aes(x = NES, y = Description, 
                                  color = p.adjust, size = setSize / sum(setSize)), showCategory = 10) +
  geom_point() + 
  labs(x = "Normalized Enrichment Score (NES)", y = "Reactome Pathways", color = "Adj. p-value", size = "Gene Ratio") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Omicron") +
  scale_y_discrete(labels = custom_labeller)
reactomeO_1

# ClusterprofilR GSEA REACTOME pathway
gse_o_1 <- clusterProfiler::GSEA(gene_list_o_1,
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               minGSSize = 5,
                               gson = NULL,
                               TERM2GENE = hs_gene2name_GO)


GOplot_o_1 <- ggplot(data = gse_o_1@result%>%
                     top_n(15, NES), 
                   mapping = aes(x = NES, y = Description, 
                                 color = p.adjust, size = setSize)) +
  geom_point() + 
  scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Gene Ontology Pathways", color = "Adj. p-value", size = "Gene Count") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Omicron") +
  scale_y_discrete(labels = custom_labeller)

GOplot_o_1


# Combine plots and labels using grid.arrange

# # Combine plots and labels using grid.arrange
# reactome_plot_1 <- grid.arrange(reactomeW_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                               reactomeD_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                               reactomeO_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                               nrow = 3, ncol = 1, heights = c(1.1,1.1,1.1))
# 
# 
# # Combine plots and labels using grid.arrange
# GO_plot_1 <-  grid.arrange(GOplot_w_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                          GOplot_d_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                          GOplot_o_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                          nrow = 3, ncol = 1, heights = c(1.1,1.1,1.1))
# 
# grid.arrange(
#   reactome_plot_1 ,
#   get_legend(reactomeD_1),
#   GO_plot_1 ,  
#   get_legend(GOplot_o_1),
#   ncol = 4, widths = c(1.1,0.14,1.1,0.14)
# )

GO_plot_1 <-  grid.arrange(GOplot_d_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
                           GOplot_o_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "right"),
                           nrow = 1, ncol = 2, widths = c(1.1,1.15))

reactome_plot_1 <-  grid.arrange(reactomeD_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
                                 reactomeO_1 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "right"),
                                 nrow = 1, ncol = 2, widths = c(1.1,1.15))


gsea_results_w_1 <- as.data.frame(gse_w_1@result)
gsea_results_w_1 <- as.data.frame(cbind(ID=gsea_results_w_1$ID,
                                        NES=gsea_results_w_1$NES,
                                        p.adjust=gsea_results_w_1$p.adjust))

gsea_results_d_1 <- as.data.frame(gse_d_1@result)
gsea_results_d_1 <- as.data.frame(cbind(ID=gsea_results_d_1$ID,
                                        NES=gsea_results_d_1$NES,
                                        p.adjust=gsea_results_d_1$p.adjust))

gsea_results_o_1 <- as.data.frame(gse_o_1@result)
gsea_results_o_1 <- as.data.frame(cbind(ID=gsea_results_o_1$ID,
                                        NES=gsea_results_o_1$NES,
                                        p.adjust=gsea_results_o_1$p.adjust))

######################################
# 2nd day of infection:
######################################

diff_genes_2 <- diff_genes_2 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_2$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))

diff_genes_2_immune <- diff_genes_2_immune %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_2_immune$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))
#######################################
# Beginning of infection: 
#######################################
### Omicron ###
diff_genes_2_Omicron <- diff_genes_2 %>%
  filter(group == "Omicron_2") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_2_Delta <- diff_genes_2 %>%
  filter(group == "Delta_2") %>%
  filter(!is.na(ENTREZID)) 

### Wuhan ###
diff_genes_2_Wuhan <- diff_genes_2 %>%
  filter(group == "Wuhan_2") %>%
  filter(!is.na(ENTREZID)) 


# Immune only
### Omicron ###
diff_genes_2_Omicron_immune <- diff_genes_2_immune %>%
  filter(group == "Omicron_2") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_2_Delta_immune <- diff_genes_2_immune %>%
  filter(group == "Delta_2") %>%
  filter(!is.na(ENTREZID)) 

### Wuhan ###
diff_genes_2_Wuhan_immune <- diff_genes_2_immune %>%
  filter(group == "Wuhan_2") %>%
  filter(!is.na(ENTREZID)) 

#############
### Wuhan ###
# Ranking genes for GSEA
result_wuhan_2 <- result_wuhan_2 %>% filter(!is.na(log2FoldChange))
result_wuhan_2 <- result_wuhan_2 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_wuhan_2$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))
gene_list_w_2 <- result_wuhan_2$log2FoldChange
names(gene_list_w_2) <- result_wuhan_2$ENTREZID
gene_list_w_2 = sort(gene_list_w_2, decreasing = TRUE)


#############
### Delta ###
# Ranking genes for GSEA
result_delta_2 <- result_delta_2 %>% filter(!is.na(log2FoldChange))
result_delta_2 <- result_delta_2 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_delta_2$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))
gene_list_d_2<- result_delta_2$log2FoldChange
names(gene_list_d_2) <- result_delta_2$ENTREZID
gene_list_d_2= sort(gene_list_d_2, decreasing = TRUE)


#############
### Omicron ###
# Ranking genes for GSEA
result_omicron_2 <- result_omicron_2 %>% filter(!is.na(log2FoldChange))
result_omicron_2 <- result_omicron_2 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_omicron_2$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))
gene_list_o_2<- result_omicron_2$log2FoldChange
names(gene_list_o_2) <- result_omicron_2$ENTREZID
gene_list_o_2= sort(gene_list_o_2, decreasing = TRUE)



# ClusterprofilR GSEA REACTOME pathway
# enricher_host_Wuhan_2 <- clusterProfiler::GSEA(gene_list_w_2,
#                                                pvalueCutoff = 0.05,
#                                                #OrgDb = org.Hs.eg.db,
#                                                pAdjustMethod = "BH",
#                                                minGSSize = 5,
#                                                #ont = "ALL",
#                                                gson = NULL,
#                                                TERM2GENE = hs_gene2name_reactome)
# 
# reactomeW_2 <- ggplot(data = enricher_host_Wuhan_2@result%>%
#                         top_n(10, rank), 
#                       mapping = aes(x = NES, y = Description, 
#                                     color = p.adjust, size = setSize / sum(setSize))) +
#   geom_point() + 
#   scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
#   labs(x = "Normalized Enrichment Score (NES)", y = "Reactome Pathways", color = "Adj. p-value", size = "Gene Ratio") +
#   theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
#         axis.text.y = element_text(size = 9, hjust = 1),
#         axis.text.x = element_text(size = 12, hjust = 1),
#         axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
#         axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
#         plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
#   scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
#   ggtitle("WuhanHu1") +
#   scale_y_discrete(labels = custom_labeller)
# reactomeW_2
# 
# # ClusterprofilR GSEA REACTOME pathway
# gse_w_2 <- clusterProfiler::GSEA(gene_list_w_2,
#                                  pvalueCutoff = 0.05,
#                                  pAdjustMethod = "BH",
#                                  minGSSize = 5,
#                                  gson = NULL,
#                                  TERM2GENE = hs_gene2name_GO)
# 
# 
# GOplot_w_2 <- ggplot(data = gse_w_2@result%>%
#                        top_n(15, NES),
#                      mapping = aes(x = NES, y = Description, 
#                                    color = p.adjust, size = setSize )) +
#   geom_point() + 
#   scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
#   labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Gene Ontology Pathways", color = "Adj. p-value", size = "Gene Count") +
#   theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
#         axis.text.y = element_text(size = 9, hjust = 1),
#         axis.text.x = element_text(size = 12, hjust = 1),
#         axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
#         axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
#         plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
#   scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
#   ggtitle("WuhanHu1") +
#   scale_y_discrete(labels = custom_labeller)
# 
# GOplot_w_2




# ClusterprofilR GSEA REACTOME pathway
enricher_host_Delta_2 <- clusterProfiler::GSEA(gene_list_d_2,
                                               pvalueCutoff = 0.05,
                                               #OrgDb = org.Hs.eg.db,
                                               pAdjustMethod = "BH",
                                               minGSSize = 5,
                                               #ont = "ALL",
                                               gson = NULL,
                                               TERM2GENE = hs_gene2name_reactome)


reactomeD_2 <- ggplot(data = enricher_host_Delta_2@result%>%
                        top_n(10, rank), 
                      mapping = aes(x = NES, y = Description, 
                                    color = p.adjust, size = setSize / sum(setSize))) +
  geom_point() + 
  scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "Reactome Pathways", color = "Adj. p-value", size = "Gene Ratio") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Delta") +
  scale_y_discrete(labels = custom_labeller)

reactomeD_2




# ClusterprofilR GSEA REACTOME pathway
gse_d_2 <- clusterProfiler::GSEA(gene_list_d_2,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 minGSSize = 5,
                                 gson = NULL,
                                 TERM2GENE = hs_gene2name_GO)


GOplot_d_2 <- ggplot(data = gse_d_2@result%>%
                       top_n(15, NES), 
                     mapping = aes(x = NES, y = Description, 
                                   color = p.adjust, size = setSize)) +
  geom_point() + 
  scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Gene Ontology Pathways", color = "Adj. p-value", size = "Gene Count") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Delta") +
  scale_y_discrete(labels = custom_labeller)
GOplot_d_2




# ClusterprofilR GSEA REACTOME pathway
enricher_host_Omicron_2 <- clusterProfiler::GSEA(gene_list_o_2,
                                                 pvalueCutoff = 0.05,
                                                 #OrgDb = org.Hs.eg.db,
                                                 pAdjustMethod = "BH",
                                                 minGSSize = 5,
                                                 #ont = "ALL",
                                                 gson = NULL,
                                                 TERM2GENE = hs_gene2name_reactome)


reactomeO_2 <- ggplot(data = enricher_host_Omicron_2@result %>%
                        top_n(10, rank), 
                      mapping = aes(x = NES, y = Description, 
                                    color = p.adjust, size = setSize / sum(setSize)), showCategory = 10) +
  geom_point() + 
  labs(x = "Normalized Enrichment Score (NES)", y = "Reactome Pathways", color = "Adj. p-value", size = "Gene Ratio") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Omicron") +
  scale_y_discrete(labels = custom_labeller)
reactomeO_2

# ClusterprofilR GSEA REACTOME pathway
gse_o_2 <- clusterProfiler::GSEA(gene_list_o_2,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 minGSSize = 5,
                                 gson = NULL,
                                 TERM2GENE = hs_gene2name_GO)


GOplot_o_2 <- ggplot(data = gse_o_2@result%>%
                       top_n(15, NES), 
                     mapping = aes(x = NES, y = Description, 
                                   color = p.adjust, size = setSize)) +
  geom_point() + 
  scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Gene Ontology Pathways", color = "Adj. p-value", size = "Gene Count") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Omicron") +
  scale_y_discrete(labels = custom_labeller)

GOplot_o_2


# Combine plots and labels using grid.arrange

# # Combine plots and labels using grid.arrange
# reactome_plot_2 <- grid.arrange(reactomeW_2 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                               reactomeD_2 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                               reactomeO_2 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                               nrow = 3, ncol = 1, heights = c(1.1,1.1,1.1))
# 
# 
# # Combine plots and labels using grid.arrange
# GO_plot_2 <-  grid.arrange(GOplot_w_2 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                          GOplot_d_2 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                          GOplot_o_2 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                          nrow = 3, ncol = 1, heights = c(1.1,1.1,1.1))
# 
# grid.arrange(
#   reactome_plot_2 ,
#   get_legend(reactomeD_2),
#   GO_plot_2 ,  
#   get_legend(GOplot_o_2),
#   ncol = 4, widths = c(1.1,0.14,1.1,0.14)
# )

GO_plot_2 <-  grid.arrange(GOplot_d_2 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
                           GOplot_o_2 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "right"),
                           nrow = 1, ncol = 2, widths = c(1.1,1.15))


gsea_results_w_2 <- as.data.frame(gse_w_2@result)
gsea_results_w_2 <- as.data.frame(cbind(ID=gsea_results_w_2$ID,
                                        NES=gsea_results_w_2$NES,
                                        p.adjust=gsea_results_w_2$p.adjust))

gsea_results_d_2 <- as.data.frame(gse_d_2@result)
gsea_results_d_2 <- as.data.frame(cbind(ID=gsea_results_d_2$ID,
                                        NES=gsea_results_d_2$NES,
                                        p.adjust=gsea_results_d_2$p.adjust))

gsea_results_o_2 <- as.data.frame(gse_o_2@result)
gsea_results_o_2 <- as.data.frame(cbind(ID=gsea_results_o_2$ID,
                                        NES=gsea_results_o_2$NES,
                                        p.adjust=gsea_results_o_2$p.adjust))
######################################
# 3rd day of infection:
######################################

diff_genes_3 <- diff_genes_3 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_3$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))

diff_genes_3_immune <- diff_genes_3_immune %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_3_immune$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))
#######################################
# Beginning of infection: 
#######################################
### Omicron ###
diff_genes_3_Omicron <- diff_genes_3 %>%
  filter(group == "Omicron_3") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_3_Delta <- diff_genes_3 %>%
  filter(group == "Delta_3") %>%
  filter(!is.na(ENTREZID)) 

### Wuhan ###
diff_genes_3_Wuhan <- diff_genes_3 %>%
  filter(group == "Wuhan_3") %>%
  filter(!is.na(ENTREZID)) 


# Immune only
### Omicron ###
diff_genes_3_Omicron_immune <- diff_genes_3_immune %>%
  filter(group == "Omicron_3") %>%
  filter(!is.na(ENTREZID)) 

### Delta ###
diff_genes_3_Delta_immune <- diff_genes_3_immune %>%
  filter(group == "Delta_3") %>%
  filter(!is.na(ENTREZID)) 

### Wuhan ###
diff_genes_3_Wuhan_immune <- diff_genes_3_immune %>%
  filter(group == "Wuhan_3") %>%
  filter(!is.na(ENTREZID)) 

#############
### Wuhan ###
# Ranking genes for GSEA
# result_wuhan_3 <- result_wuhan_3 %>% filter(!is.na(log2FoldChange))
# result_wuhan_3 <- result_wuhan_3 %>% 
#   mutate(ENTREZID = mapIds(org.Hs.eg.db, result_wuhan_3$geneName, 
#                            fromType = "SYMBOL", toType = "ENTREZID", 
#                            column = "ENTREZID", keytype = "SYMBOL"))
# gene_list_w_3 <- result_wuhan_3$log2FoldChange
# names(gene_list_w_3) <- result_wuhan_3$ENTREZID
# gene_list_w_3 = sort(gene_list_w_3, decreasing = TRUE)


#############
### Delta ###
# Ranking genes for GSEA
result_delta_3 <- result_delta_3 %>% filter(!is.na(log2FoldChange))
result_delta_3 <- result_delta_3 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_delta_3$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))
gene_list_d_3<- result_delta_3$log2FoldChange
names(gene_list_d_3) <- result_delta_3$ENTREZID
gene_list_d_3= sort(gene_list_d_3, decreasing = TRUE)


#############
### Omicron ###
# Ranking genes for GSEA
result_omicron_3 <- result_omicron_3 %>% filter(!is.na(log2FoldChange))
result_omicron_3 <- result_omicron_3 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, result_omicron_3$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))
gene_list_o_3<- result_omicron_3$log2FoldChange
names(gene_list_o_3) <- result_omicron_3$ENTREZID
gene_list_o_3= sort(gene_list_o_3, decreasing = TRUE)



# ClusterprofilR GSEA REACTOME pathway
# enricher_host_Wuhan_3 <- clusterProfiler::GSEA(gene_list_w_3,
#                                                pvalueCutoff = 0.05,
#                                                #OrgDb = org.Hs.eg.db,
#                                                pAdjustMethod = "BH",
#                                                minGSSize = 5,
#                                                #ont = "ALL",
#                                                gson = NULL,
#                                                TERM2GENE = hs_gene2name_reactome)
# 
# reactomeW_3 <- ggplot(data = enricher_host_Wuhan_3@result%>%
#                         top_n(10, rank), 
#                       mapping = aes(x = NES, y = Description, 
#                                     color = p.adjust, size = setSize / sum(setSize))) +
#   geom_point() + 
#   scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
#   labs(x = "Normalized Enrichment Score (NES)", y = "Reactome Pathways", color = "Adj. p-value", size = "Gene Ratio") +
#   theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
#         axis.text.y = element_text(size = 9, hjust = 1),
#         axis.text.x = element_text(size = 12, hjust = 1),
#         axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
#         axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
#         plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
#   scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
#   ggtitle("WuhanHu1") +
#   scale_y_discrete(labels = custom_labeller)
# reactomeW_3
# 
# # ClusterprofilR GSEA REACTOME pathway
# gse_w_3 <- clusterProfiler::GSEA(gene_list_w_3,
#                                  pvalueCutoff = 0.05,
#                                  pAdjustMethod = "BH",
#                                  minGSSize = 5,
#                                  gson = NULL,
#                                  TERM2GENE = hs_gene2name_GO)
# 
# 
# GOplot_w_3 <- ggplot(data = gse_w_3@result%>%
#                        top_n(15, NES),
#                      mapping = aes(x = NES, y = Description, 
#                                    color = p.adjust, size = setSize)) +
#   geom_point() + 
#   scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
#   labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Gene Ontology Pathways", color = "Adj. p-value", size = "Gene Count") +
#   theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
#         axis.text.y = element_text(size = 9, hjust = 1),
#         axis.text.x = element_text(size = 12, hjust = 1),
#         axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
#         axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
#         plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
#   scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
#   ggtitle("WuhanHu1") +
#   scale_y_discrete(labels = custom_labeller)
# 
# GOplot_w_3




# ClusterprofilR GSEA REACTOME pathway
enricher_host_Delta_3 <- clusterProfiler::GSEA(gene_list_d_3,
                                               pvalueCutoff = 0.05,
                                               #OrgDb = org.Hs.eg.db,
                                               pAdjustMethod = "BH",
                                               minGSSize = 5,
                                               #ont = "ALL",
                                               gson = NULL,
                                               TERM2GENE = hs_gene2name_reactome)


reactomeD_3 <- ggplot(data = enricher_host_Delta_3@result%>%
                        top_n(10, rank), 
                      mapping = aes(x = NES, y = Description, 
                                    color = p.adjust, size = setSize / sum(setSize))) +
  geom_point() + 
  scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "Reactome Pathways", color = "Adj. p-value", size = "Gene Ratio") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Delta") +
  scale_y_discrete(labels = custom_labeller)

reactomeD_3




# ClusterprofilR GSEA REACTOME pathway
gse_d_3 <- clusterProfiler::GSEA(gene_list_d_3,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 minGSSize = 5,
                                 gson = NULL,
                                 TERM2GENE = hs_gene2name_GO)


GOplot_d_3 <- ggplot(data = gse_d_3@result%>%
                       top_n(15, NES), 
                     mapping = aes(x = NES, y = Description, 
                                   color = p.adjust, size = setSize)) +
  geom_point() + 
  scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Gene Ontology Pathways", color = "Adj. p-value", size = "Gene Count") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Delta") +
  scale_y_discrete(labels = custom_labeller)
GOplot_d_3




# ClusterprofilR GSEA REACTOME pathway
enricher_host_Omicron_3 <- clusterProfiler::GSEA(gene_list_o_3,
                                                 pvalueCutoff = 0.05,
                                                 #OrgDb = org.Hs.eg.db,
                                                 pAdjustMethod = "BH",
                                                 minGSSize = 5,
                                                 #ont = "ALL",
                                                 gson = NULL,
                                                 TERM2GENE = hs_gene2name_reactome)


reactomeO_3 <- ggplot(data = enricher_host_Omicron_3@result %>%
                        top_n(10, rank), 
                      mapping = aes(x = NES, y = Description, 
                                    color = p.adjust, size = setSize / sum(setSize)), showCategory = 10) +
  geom_point() + 
  labs(x = "Normalized Enrichment Score (NES)", y = "Reactome Pathways", color = "Adj. p-value", size = "Gene Ratio") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Omicron") +
  scale_y_discrete(labels = custom_labeller)
reactomeO_3

# ClusterprofilR GSEA REACTOME pathway
gse_o_3 <- clusterProfiler::GSEA(gene_list_o_3,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 minGSSize = 5,
                                 gson = NULL,
                                 TERM2GENE = hs_gene2name_GO)


GOplot_o_3 <- ggplot(data = gse_o_3@result%>%
                       top_n(15, NES), 
                     mapping = aes(x = NES, y = Description, 
                                   color = p.adjust, size = setSize)) +
  geom_point() + 
  scale_y_discrete(labels = function(x) strwrap(x, width = 100)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Gene Ontology Pathways", color = "Adj. p-value", size = "Gene Count") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.text.x = element_text(size = 12, hjust = 1),
        axis.title.x =  element_text(size = 12, hjust = 0.5, face ="bold"),
        axis.title.y = element_text(size = 12, hjust = 0.5, face ="bold"),
        plot.background = element_rect(colour = "black", fill=NA, linewidth = 0.5)) +
  scale_color_viridis_c(limits = c(0, 0.05),  option = "plasma") +
  ggtitle("Omicron") +
  scale_y_discrete(labels = custom_labeller)

GOplot_o_3


# Combine plots and labels using grid.arrange

# # Combine plots and labels using grid.arrange
# reactome_plot_3 <- grid.arrange(reactomeW_3 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                               reactomeD_3 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                               reactomeO_3 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                               nrow = 3, ncol = 1, heights = c(1.1,1.1,1.1))
# 
# 
# # Combine plots and labels using grid.arrange
# GO_plot_3 <-  grid.arrange(GOplot_w_3 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                          GOplot_d_3 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                          GOplot_o_3 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
#                          nrow = 3, ncol = 1, heights = c(1.1,1.1,1.1))
# 
# grid.arrange(
#   reactome_plot_3 ,
#   get_legend(reactomeD_3),
#   GO_plot_3 ,  
#   get_legend(GOplot_o_3),
#   ncol = 4, widths = c(1.1,0.14,1.1,0.14)
# )

GO_plot_3 <-  grid.arrange(GOplot_d_3 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "none"),
                           GOplot_o_3 + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),legend.position = "right"),
                           nrow = 1, ncol = 2, widths = c(1.1,1.15))


grid.arrange(
  GO_plot_1,
  GO_plot_2,
  GO_plot_3 ,  
  ncol = 1, nrow= 3)

# gsea_results_w_3 <- as.data.frame(gse_w_3@result)
# gsea_results_w_3 <- as.data.frame(cbind(ID=gsea_results_w_3$ID,
#                                         NES=gsea_results_w_3$NES,
#                                         p.adjust=gsea_results_w_3$p.adjust))

gsea_results_d_3 <- as.data.frame(gse_d_3@result)
gsea_results_d_3 <- as.data.frame(cbind(ID=gsea_results_d_3$ID,
                                        NES=gsea_results_d_3$NES,
                                        p.adjust=gsea_results_d_3$p.adjust))

gsea_results_o_3 <- as.data.frame(gse_o_3@result)
gsea_results_o_3 <- as.data.frame(cbind(ID=gsea_results_o_3$ID,
                                        NES=gsea_results_o_3$NES,
                                        p.adjust=gsea_results_o_3$p.adjust))
















############ all three days ##############

# Extract unique gene names from each dataframe
gene_names_1 <- unique(diff_genes_1_heatmap$geneName)
gene_names_2 <- unique(diff_genes_2_heatmap$geneName)
gene_names_3 <- unique(diff_genes_3_heatmap$geneName)

length(unique(c(gene_names_1, gene_names_2, gene_names_3)))

matrix_all_times <- counts(dds_variant_time_1, normalized =T) [unique(c(gene_names_1, gene_names_2, gene_names_3)),]















#################################################################################
################ Pathway enrichr for exclusive DEGs ##############################


library(enrichR)
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")


# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
wuhan_genes <- unique(diff_genes_1_Wuhan$geneName)
enriched_wuhan <- enrichr(wuhan_genes, dbs)
enriched_wuhan_BP <- as.data.frame(enriched_wuhan[["GO_Biological_Process_2021"]])

w1 <- plotEnrich(
  enriched_wuhan_BP,
  showTerms = 10,
  numChar = 50,
  orderBy = "Adjusted.P.value",
  title = "Enriched pathways for DE host genes with Wuhan at 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 

# Extract unique gene names from top_diff_Delta_up and top_diff_Delta_down
delta_genes <- unique(diff_genes_1_Delta$geneName)
enriched_delta <- enrichr(delta_genes, dbs)
enriched_delta_BP <- as.data.frame(enriched_delta[["GO_Biological_Process_2021"]])

d1 <- plotEnrich(
  enriched_delta_BP,
  showTerms = 10,
  numChar = 50,
  orderBy = "Combined.Scores",
  title = "Enriched pathways for DE host genes with Delta at 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 


# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
omicron_genes <- unique(diff_genes_1_Omicron$geneName)
enriched_omicron <- enrichr(omicron_genes, dbs)
enriched_omicron_BP <- as.data.frame(enriched_omicron[["GO_Biological_Process_2021"]])

o1 <- plotEnrich(
  enriched_omicron_BP,
  showTerms = 10,
  numChar = 50,
  orderBy = "Combined.Scores",
  title = "Enriched pathways for DE host genes with Omicron at 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 



# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
wuhan_genes_2 <- unique(diff_genes_2_Wuhan$geneName)
enriched_wuhan_2 <- enrichr(wuhan_genes_2, dbs)
enriched_wuhan_BP_2 <- as.data.frame(enriched_wuhan_2[["GO_Biological_Process_2021"]])

w2 <- plotEnrich(
  enriched_wuhan_BP_2,
  showTerms = 10,
  numChar = 50,
  orderBy = "Adjusted.P.value",
  title = "Enriched pathways for DE host genes with Wuhan at 2 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 


# Extract unique gene names from top_diff_Delta_up and top_diff_Delta_down
delta_genes_2 <- unique(diff_genes_2_Delta$geneName)
enriched_delta_2 <- enrichr(delta_genes_2, dbs)
enriched_delta_BP_2 <- as.data.frame(enriched_delta_2[["GO_Biological_Process_2021"]])

d2 <- plotEnrich(
  enriched_delta_BP_2,
  showTerms = 10,
  numChar = 50,
  orderBy = "Combined.Scores",
  title = "Enriched pathways for DE host genes with Delta at 2 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12))


# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
omicron_genes_2 <- unique(diff_genes_2_Omicron$geneName)
enriched_omicron_2 <- enrichr(omicron_genes_2, dbs)
enriched_omicron_BP_2 <- as.data.frame(enriched_omicron_2[["GO_Biological_Process_2021"]])

o2 <- plotEnrich(
  enriched_omicron_BP_2,
  showTerms = 10,
  numChar = 50,
  orderBy = "Combined.Scores",
  title = "Enriched pathways for DE host genes with Omicron at 2 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12))




# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
wuhan_genes_3 <- unique(diff_genes_3_Wuhan$geneName)
enriched_wuhan_3 <- enrichr(wuhan_genes_3, dbs)
enriched_wuhan_BP_3 <- as.data.frame(enriched_wuhan_3[["GO_Biological_Process_2021"]])

w3 <- plotEnrich(
  enriched_wuhan_BP_3,
  showTerms = 10,
  numChar = 50,
  orderBy = "Adjusted.P.value",
  title = "Enriched pathways for DE host genes with Wuhan at 3 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 



# Extract unique gene names from top_diff_Delta_up and top_diff_Delta_down
delta_genes_3 <- unique(diff_genes_3_Delta$geneName)
enriched_delta_3 <- enrichr(delta_genes_3, dbs)
enriched_delta_BP_3 <- as.data.frame(enriched_delta_3[["GO_Biological_Process_2021"]])

d3 <- plotEnrich(
  enriched_delta_BP_3,
  showTerms = 10,
  numChar = 50,
  orderBy = "Combined.Scores",
  title = "Enriched pathways for DE host genes with Delta at 3 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12))


# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
omicron_genes_3 <- unique(diff_genes_3_Omicron$geneName)
enriched_omicron_3 <- enrichr(omicron_genes_3, dbs)
enriched_omicron_BP_3 <- as.data.frame(enriched_omicron_3[["GO_Biological_Process_2021"]])

o3 <- plotEnrich(
  enriched_omicron_BP_3,
  showTerms = 10,
  numChar = 50,
  orderBy = "Combined.Scores",
  title = "Enriched pathways for DE host genes with Omicron at 3 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 


 

# combined_wuhan <- grid.arrange(w1, w3, ncol = 1)

combined_delta <- grid.arrange(d1, d2, d3, ncol = 1)
combined_omicron <- grid.arrange(o1, o2, o3, ncol = 1)
combined_enrich <- grid.arrange(combined_delta, combined_omicron, ncol = 2)
  
combined_1 <- grid.arrange(w1, d1, o1, ncol = 1)
combined_2 <- grid.arrange(w2, d2, o2, ncol = 1)
combined_3 <- grid.arrange(w3, d3, o3, ncol = 1)







################################################################################
######################## Omicron V Delta at 1 dpi ##############################


result_oVd_1 <- results(deseq_variant_time_1, contrast = c("variant","omicron","delta"))
summary(result_oVd_1)
head(result_oVd_1)
result_oVd_1 <- as.data.frame(result_oVd_1)
result_oVd_1 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_oVd_1 <- result_oVd_1[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "OmicronVSDelta_1")

df_oVd_1 <- df_oVd_1 %>% filter(!is.na(log2FoldChange) & !is.na(padj))

# Filter for differentially expressed genes
#log2FC (1.5) = 2

diff_genes_oVd_1 <- df_oVd_1%>% 
  filter(abs(log2FoldChange) > 2 & padj < 0.05)

# filtering out immune genes
diff_genes_oVd_1 <- diff_genes_oVd_1[diff_genes_oVd_1$geneName %in% 
                                       hs_immune_genes$gene_symbol, ]

diff_genes_oVd_1 <- diff_genes_oVd_1 %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) 

fpkm_omicron_omicronVdelta_1 <- fpkm_human_t_1_melted %>%
  filter(Gene %in% diff_genes_oVd_1$geneName)

fpkm_omicron_omicronVdelta_1$Gene <- factor(
  fpkm_omicron_omicronVdelta_1$Gene,
  levels = diff_genes_oVd_1$geneName
)



# Set color palette
my_colors_c <- c("dodgerblue","forestgreen","goldenrod","#FC8D62")

##### Top GENES for highest FC in DE genes in Omicron

plot_omicronVdelta_1 <- fpkm_human_t_1_melted %>%
  filter(Gene %in% fpkm_omicron_omicronVdelta_1$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = diff_genes_oVd_1$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 6, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for DE host genes
       in human cells infected with Omicron at 1 DPI compared to Delta") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_omicronVdelta_1



######################## Omicron V Delta at 3 dpi ##############################


result_oVd_3 <- results(deseq_variant_time_3, contrast = c("variant","omicron","delta"))
summary(result_oVd_3)
head(result_oVd_3)
result_oVd_3 <- as.data.frame(result_oVd_3)
result_oVd_3 %<>% {.$geneName <- rownames(.); rownames(.) <- NULL; .}

df_oVd_3 <- result_oVd_3[, c("geneName", "baseMean", "log2FoldChange", "padj")] %>% mutate(group = "OmicronVSDelta_3")

df_oVd_3 <- df_oVd_3 %>% filter(!is.na(log2FoldChange) & !is.na(padj))

# Filter for differentially expressed genes
#log2FC (1.5) = 2

diff_genes_oVd_3 <- df_oVd_3%>% 
  filter(abs(log2FoldChange) > 2 & padj < 0.05)

# filtering out immune genes
diff_genes_oVd_3 <- diff_genes_oVd_3[diff_genes_oVd_3$geneName %in% 
                                       hs_immune_genes$gene_symbol, ]

diff_genes_oVd_3 <- diff_genes_oVd_3 %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) 

fpkm_omicron_omicronVdelta_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% diff_genes_oVd_3$geneName)

fpkm_omicron_omicronVdelta_3$Gene <- factor(
  fpkm_omicron_omicronVdelta_3$Gene,
  levels = diff_genes_oVd_3$geneName
)



# Set color palette
my_colors_c <- c("dodgerblue","forestgreen","goldenrod","#FC8D62")

##### Top GENES for highest FC in DE genes in Omicron

plot_omicronVdelta_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% fpkm_omicron_omicronVdelta_3$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = diff_genes_oVd_3$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 6, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for DE host genes
       in human cells infected with Omicron at 3 DPI compared to Delta") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_omicronVdelta_3








################################## VENN #########################################
################################################################################






# Function to calculate number of DEG (up and down) for each variant and print summary
calculate_and_print_summary <- function(df, name) {
  cat("Summary for", name, ":\n")
  df <- df %>% 
    distinct(geneName,group,log2FoldChange,.keep_all = TRUE)
  # Total number of DEG 
  total_rows <- table(df$group)
  cat("Total number of DEG for each variant:\n")
  print(total_rows)
  
  # # Count of positive and negative values in "log2FoldChange" column for each unique group
  # count_positive <- table(df$group[df$log2FoldChange > 0])
  # count_negative <- table(df$group[df$log2FoldChange < 0])
  # cat("Upregulated genes for each variant:\n")
  # print(count_positive)
  # cat("Downregulated genes for each variant :\n")
  # print(count_negative)
  
  cat("\n")
}


# Function to count shared genes between the diff_ dataframes based on group values
shared_genes <- function(df1, df2, df3 = NULL, variant1, variant2, variant3 = NULL) {
  # Subset data frames based on group values
  subset_df1 <- df1[df1$group == variant1, ]
  subset_df1 <- subset_df1 %>% 
    distinct(geneName, .keep_all = TRUE)
  
  subset_df2 <- df2[df2$group == variant2, ]
  subset_df2 <- subset_df2 %>% 
    distinct(geneName, .keep_all = TRUE)
  
  if (!is.null(df3) && !is.null(variant3)) {
    subset_df3 <- df3[df3$group == variant3, ]
    shared_genes <- Reduce(intersect, list(subset_df1$geneName, subset_df2$geneName, subset_df3$geneName))
    cat("Shared genes between", variant1, ",", variant2, "and", variant3, ":\n")
  } else {
    shared_genes <- intersect(subset_df1$geneName, subset_df2$geneName)
    cat("Shared genes between", variant1, "and", variant2, ":\n")
  }
  
  num_shared_genes <- length(shared_genes)
  cat(num_shared_genes, "shared genes found.\n")
  
  if (num_shared_genes > 0) {
    cat("Gene names:\n")
  } else {
    cat("No shared genes found.\n")
  }
  
  # Return the list of shared genes
  return(shared_genes)
}

# all DEG
calculate_and_print_summary(diff_genes_1, "diff_genes_1")
calculate_and_print_summary(diff_genes_2, "diff_genes_2")
calculate_and_print_summary(diff_genes_3, "diff_genes_3")

shared_DEG_hc_d_o <- shared_genes(diff_genes_3, diff_genes_3, NULL,"Omicron_3", "Delta_3", NULL)
shared_DEG_hc_w_o <- shared_genes(diff_genes_3, diff_genes_3, NULL, "Omicron_3", "Wuhan_3",NULL)
shared_DEG_hc_w_d <- shared_genes(diff_genes_3, diff_genes_3, NULL, "Wuhan_3", "Delta_3", NULL)
shared_DEG_hc <- shared_genes(diff_genes_3, diff_genes_3, diff_genes_3 , "Delta_3", "Omicron_3", "Wuhan_3")



library(VennDiagram)

######## Venn Diagram
DEG_counts_venn_1 <- list(Wuhan = 2180,
                            Delta = 2478,
                            Omicron = 1355,
                            `Delta+Wuhan+Omicron` = 750,
                            `Delta+Omicron` = 863,
                            `Omicron+Wuhan` = 954,
                            `Wuhan+Delta` = 1408)
dev.off()
# Create the Venn diagram
venn.plot_1 <- draw.triple.venn(
  area1 = DEG_counts_venn_1$Wuhan,
  area2 = DEG_counts_venn_1$Delta,
  area3 = DEG_counts_venn_1$Omicron,
  n12 = DEG_counts_venn_1$`Wuhan+Delta`,
  n13 = DEG_counts_venn_1$`Omicron+Wuhan`,
  n23 = DEG_counts_venn_1$`Delta+Omicron`,
  n123 = DEG_counts_venn_1$`Delta+Wuhan+Omicron`,
  category = c("Wuhan", "Delta", "Omicron"),
  fill = c("dodgerblue","goldenrod","#FC8D62"),
  alpha = 0.5,
  main = "Venn Diagram of Wuhan, Delta, and Omicron",
  cat.cex = 2,
  cex =2
)
venn.plot_1 <- arrangeGrob(venn.plot_1, left = textGrob("1 dpi", gp = gpar(fontsize = 18)))



######## Venn Diagram
DEG_counts_venn_2 <- list(Wuhan = 982,
                          Delta = 1160,
                          Omicron = 607,
                          `Delta+Wuhan+Omicron` = 186,
                          `Delta+Omicron` = 272,
                          `Omicron+Wuhan` = 213,
                          `Wuhan+Delta` = 589)
dev.off()
# Create the Venn diagram
venn.plot_2 <- draw.triple.venn(
  area1 = DEG_counts_venn_2$Wuhan,
  area2 = DEG_counts_venn_2$Delta,
  area3 = DEG_counts_venn_2$Omicron,
  n12 = DEG_counts_venn_2$`Wuhan+Delta`,
  n13 = DEG_counts_venn_2$`Omicron+Wuhan`,
  n23 = DEG_counts_venn_2$`Delta+Omicron`,
  n123 = DEG_counts_venn_2$`Delta+Wuhan+Omicron`,
  category = c("Wuhan", "Delta", "Omicron"),
  fill = c("dodgerblue","goldenrod","#FC8D62"),
  alpha = 0.5,
  main = "Venn Diagram of Wuhan, Delta, and Omicron",
  cat.cex = 2,
  cex =2
)
venn.plot_2 <- arrangeGrob(venn.plot_2, left = textGrob("2 dpi", gp = gpar(fontsize = 18)))


######## Venn Diagram
DEG_counts_venn_3 <- list(Wuhan = 748,
                          Delta = 1664,
                          Omicron = 727,
                          `Delta+Wuhan+Omicron` = 310,
                          `Delta+Omicron` = 514,
                          `Omicron+Wuhan` = 325,
                          `Wuhan+Delta` = 640)
dev.off()
# Create the Venn diagram
venn.plot_3 <- draw.triple.venn(
  area1 = DEG_counts_venn_3$Wuhan,
  area2 = DEG_counts_venn_3$Delta,
  area3 = DEG_counts_venn_3$Omicron,
  n12 = DEG_counts_venn_3$`Wuhan+Delta`,
  n13 = DEG_counts_venn_3$`Omicron+Wuhan`,
  n23 = DEG_counts_venn_3$`Delta+Omicron`,
  n123 = DEG_counts_venn_3$`Delta+Wuhan+Omicron`,
  category = c("Wuhan", "Delta", "Omicron"),
  fill = c("dodgerblue","goldenrod","#FC8D62"),
  alpha = 0.5,
  main = "Venn Diagram of Wuhan, Delta, and Omicron",
  cat.cex = 2,
  cex =2
)

venn.plot_3 <- arrangeGrob(venn.plot_3, left = textGrob("3 dpi", gp = gpar(fontsize = 18)))


# Immune venn

calculate_and_print_summary(diff_genes_1_immune, "diff_genes_1")
calculate_and_print_summary(diff_genes_2_immune, "diff_genes_2")
calculate_and_print_summary(diff_genes_3_immune, "diff_genes_3")

shared_DEG_hc_d_o <- shared_genes(diff_genes_1_immune, diff_genes_1_immune, NULL,"Omicron_1", "Delta_1", NULL)
shared_DEG_hc_w_o <- shared_genes(diff_genes_1_immune, diff_genes_1_immune, NULL, "Omicron_1", "Wuhan_1",NULL)
shared_DEG_hc_w_d <- shared_genes(diff_genes_1_immune, diff_genes_1_immune, NULL, "Wuhan_1", "Delta_1", NULL)
shared_DEG_hc <- shared_genes(diff_genes_1_immune, diff_genes_1_immune, diff_genes_1_immune , "Delta_1", "Omicron_1", "Wuhan_1")

######## Venn Diagram
DEG_counts_venn_1_imm <- list(Wuhan = 993,
                          Delta = 1419,
                          Omicron = 583,
                          `Delta+Wuhan+Omicron` = 403,
                          `Delta+Omicron` = 453,
                          `Omicron+Wuhan` = 439,
                          `Wuhan+Delta` = 799)
dev.off()
# Create the Venn diagram
venn.plot_1_imm <- draw.triple.venn(
  area1 = DEG_counts_venn_1_imm$Wuhan,
  area2 = DEG_counts_venn_1_imm$Delta,
  area3 = DEG_counts_venn_1_imm$Omicron,
  n12 = DEG_counts_venn_1_imm$`Wuhan+Delta`,
  n13 = DEG_counts_venn_1_imm$`Omicron+Wuhan`,
  n23 = DEG_counts_venn_1_imm$`Delta+Omicron`,
  n123 = DEG_counts_venn_1_imm$`Delta+Wuhan+Omicron`,
  category = c("Wuhan", "Delta", "Omicron"),
  fill = c("dodgerblue","goldenrod","#FC8D62"),
  alpha = 0.5,
  main = "Venn Diagram of Wuhan, Delta, and Omicron",
  cat.cex = 2,
  cex =2
)



######## Venn Diagram
DEG_counts_venn_2_imm <- list(Wuhan = 557,
                          Delta = 725,
                          Omicron = 404,
                          `Delta+Wuhan+Omicron` = 122,
                          `Delta+Omicron` = 189,
                          `Omicron+Wuhan` = 139,
                          `Wuhan+Delta` = 393)
dev.off()
# Create the Venn diagram
venn.plot_2_imm <- draw.triple.venn(
  area1 = DEG_counts_venn_2_imm$Wuhan,
  area2 = DEG_counts_venn_2_imm$Delta,
  area3 = DEG_counts_venn_2_imm$Omicron,
  n12 = DEG_counts_venn_2_imm$`Wuhan+Delta`,
  n13 = DEG_counts_venn_2_imm$`Omicron+Wuhan`,
  n23 = DEG_counts_venn_2_imm$`Delta+Omicron`,
  n123 = DEG_counts_venn_2_imm$`Delta+Wuhan+Omicron`,
  category = c("Wuhan", "Delta", "Omicron"),
  fill = c("dodgerblue","goldenrod","#FC8D62"),
  alpha = 0.5,
  main = "Venn Diagram of Wuhan, Delta, and Omicron",
  cat.cex = 2,
  cex =2
)


######## Venn Diagram
DEG_counts_venn_3_imm <- list(Wuhan = 562,
                          Delta = 1094,
                          Omicron = 531,
                          `Delta+Wuhan+Omicron` = 243,
                          `Delta+Omicron` = 396,
                          `Omicron+Wuhan` = 253,
                          `Wuhan+Delta` = 499)
dev.off()
# Create the Venn diagram
venn.plot_3_imm <- draw.triple.venn(
  area1 = DEG_counts_venn_3_imm$Wuhan,
  area2 = DEG_counts_venn_3_imm$Delta,
  area3 = DEG_counts_venn_3_imm$Omicron,
  n12 = DEG_counts_venn_3_imm$`Wuhan+Delta`,
  n13 = DEG_counts_venn_3_imm$`Omicron+Wuhan`,
  n23 = DEG_counts_venn_3_imm$`Delta+Omicron`,
  n123 = DEG_counts_venn_3_imm$`Delta+Wuhan+Omicron`,
  category = c("Wuhan", "Delta", "Omicron"),
  fill = c("dodgerblue","goldenrod","#FC8D62"),
  alpha = 0.5,
  main = "Venn Diagram of Wuhan, Delta, and Omicron",
  cat.cex = 2,
  cex =2
)

# Create two empty plots
empty_plot <- ggplot() + theme_void()

venn_all <- grid.arrange(venn.plot_1,empty_plot, venn.plot_2,empty_plot, venn.plot_3, ncol = 1,
                         nrow=5, heights=c(1,0.1,1,0.1,1)) 



dev.off()

venn_all_1 <- arrangeGrob(empty_plot, venn_all, empty_plot, ncol = 3, widths = c(.2,.8,.2))
venn_all_1 <- arrangeGrob(empty_plot, venn_all_1, nrow = 2, heights = c(.05,.9))
venn_all_1 <- arrangeGrob(venn_all_1, top = textGrob("All DEGs", gp = gpar(fontsize = 14, fontface = "bold")))

# Display the venn_all_1 plot using grid.draw
grid.draw(venn_all_1)

venn_all_imm <- grid.arrange(venn.plot_1_imm,empty_plot, venn.plot_2_imm,empty_plot, venn.plot_3_imm, ncol = 1,
                         nrow=5, heights=c(1,0.1,1,0.1,1)) 
dev.off()

venn_all_imm_1 <- arrangeGrob(empty_plot, venn_all_imm, empty_plot, ncol = 3, widths = c(.2,.8,.2))
venn_all_imm_1 <- arrangeGrob(empty_plot, venn_all_imm_1, nrow = 2, heights = c(.05,.9))

# Display the venn_all_1 plot using grid.draw
grid.draw(venn_all_imm_1)
# Add a title to the combined plot
venn_all_imm_1 <- arrangeGrob(venn_all_imm_1, top = textGrob("Immune DEGs", gp = gpar(fontsize = 14, fontface = "bold")))

combined_venn <- grid.arrange(venn_all_1, venn_all_imm_1,ncol = 2, widths = c(1.2,1))









################################################################################
################################################################################




allowWGCNAThreads()          # allow multi-threading (optional)

# 1. Fetch Data ------------------------------------------------
#data <- count_matrix
data <- count_matrix %>% 
  filter(rownames(.) %in% hs_immune_genes$gene_symbol)

wuhan_rows <- rownames(design_matrix)[design_matrix$variant == "Wuhan"]

data <- data %>%
        dplyr::select(-any_of(wuhan_rows))



# 2. WGCNA's QC step for outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK # allOK = TRUE means all the SRR runs passed and there are no outliers

table(gsg$goodGenes) # check outliers in rows (genes)
table(gsg$goodSamples) # check outliers in columns. (SRR runs)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,] # Selecting only genes that are not outliers
gsg <- goodSamplesGenes(t(data))
table(gsg$goodGenes) # check outliers in rows (genes)

dev.off()
# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


# pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


### NOTE: If there are batch effects observed, correct for them before moving ahead


# exclude outlier samples
#samples.to.be.excluded <- c('SRR23529187', 'SRR23529213', 'SRR23529216', 'SRR23529211', 'SRR23529212', 'SRR23529214', 'SRR23529215')
#samples.to.be.excluded <- c('SRR23529216', 'SRR23529215')
samples.to.be.excluded <- c('SRR23529181')

data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]


# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
colData <- pheno_data_host %>% 
  filter(!Run %in% wuhan_rows) %>%
  filter(!Run %in% samples.to.be.excluded) %>%
  dplyr::select(-batch) %>%
  column_to_rownames("Run")

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

#data.subset[] <- sapply(data.subset, as.numeric)
#data.subset <- round(data.subset)  # Round to nearest integer

# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model


# #varianceStabilzingTransformtion
#dds_norm <- varianceStabilizingTransformation(dds, blind = TRUE,fitType = "parametric")

## remove all genes with counts < 15 in more than 75% of samples (77*0.75=57)
#dds75 <- dds[rowSums(counts(dds) >= 15) >= 57,]
dds75 <- dds[rowSums(counts(dds) >= 15) >= 43,]
#dds75 <- dds[rowSums(counts(dds) >= 15) >= 62,]

nrow(dds75) # 11625 genes with counts less than 15 in over 75% of samples
# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>%
                      t()

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

# visualization to pick power

# Determining Maximum R square
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.9, color = 'red') +
  labs(x = 'Power', y = 'signed R^2') +
  theme_classic()

# Determining Minimum mean connectivity 
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 1)

# Determined soft threshold 
soft_power <- 10 # at 10, Rsquare value was just above 80% line and mean connectivity was lowest consecutively

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

# Assigining WGCNA correlation to cor variable to avoid namespace error
temp_cor <- cor
cor <- WGCNA::cor

# Blockwise module functions's memory estimate with respect to blocksize (4GB:8k, 16GB:20k, 32GB:30k)
bwnet <- blockwiseModules(norm.counts,
                          TOMType = "signed",
                          power = soft_power,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          minModuleSize = 30,
                          mergeCutHeight = 0.25,
                          maxBlockSize = 30000,
                          deepSplit = 1,
                          verbose = 3)

cor <- temp_cor # reassign the correlation funtion back the to original cor function 


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs

# get number of genes for each module
table(bwnet$colors)

# Plot dendrogram and module colors before and after merging genes
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors), # visualize genes before and after merging
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# several colors in the unmerged modules and fewer colors in merged modules indicate that the genes belonging to 
# different unmerged modules were similar and were merged into single module after merging and merged module 
# is chosen for further analysis

# grey module = all genes that doesn't fall into other modules were assigned to the grey module





# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

# binarize categorical variables disease status
traits <- colData %>% 
  mutate(mock = ifelse(grepl('mock', variant), 1, 0)) %>% 
  dplyr::select(3)

# binarize categorical variables variants
#colData$variant <- factor(colData$variant, levels = c("mock","Wuhan","delta","omicron"))
colData$variant <- factor(colData$variant, levels = c("mock","delta","omicron"))

variant.out <- binarizeCategoricalColumns(colData$variant,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)

traits <- cbind(traits, variant.out)
#colnames(traits) <- c("mock","Wuhan", "Delta", "Omicron")
colnames(traits) <- c("mock", "Delta", "Omicron")


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)


# visualize module-trait association as a heatmap
# adding binary train classes to ME dataframe
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[15:17],
             y = names(heatmap.data)[1:14],
             col = c("dodgerblue2", "skyblue", "white", "pink", "red2"),
             signifSymbols = c("***", "**", "*", ""),
             signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
             corFUN = "pearson",
             main = "Module-trait plot",
             cexMain = 2,
             rotLabX = 45)

# get number of genes for each module
table(bwnet$colors)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

module.gene.mapping <- as.data.frame(bwnet$colors)

# COPY OF MAIN MODULE MAP
module.gene.mapping1 <- module.gene.mapping %>%
  tibble::rownames_to_column(var = "RowNames") %>%
  setNames(c("Gene", "Module_col"))

EG_genelist_brown<- module.gene.mapping %>% 
   filter(`bwnet$colors` == 'brown') %>% 
   rownames() %>%
   as.data.frame(.)
EG_genelist_yellow<- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames() %>%
  as.data.frame(.)


#write.csv(EG_genelist_yellow,"EG_genelist_yellow.csv")



# Assuming module.gene.mapping1 is your dataframe with gene-module mapping
gene_lists_by_modules <- module.gene.mapping1 %>%
  group_by(Module_col) %>%
  summarize(Genes = list(unique(Gene))) %>%
  ungroup()



# Function to perform gene set enrichment analysis and return the most significant GO term
find_most_significant_GO_term <- function(gene_list) {
  # Convert gene symbols to Entrez IDs
  gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Perform enrichment analysis using clusterProfiler
  enrich_result <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  # Get the most significant GO term (lowest adjusted p-value)
  if (nrow(enrich_result) > 0) {
    most_significant_term <- enrich_result[5, ]
    return(most_significant_term)
  } else {
    return(NULL)
  }
}

# Apply the function to each list of genes in gene_lists_by_modules$Genes column
gene_lists_by_modules <- gene_lists_by_modules %>%
  mutate(Most_Significant_GO_Term = map(Genes, find_most_significant_GO_term))


# 6B. Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and 
# the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

# An eigengene is calculated as the first principal component (PC1) of the gene expression values within a module.
# The PC1 captures the major variation or pattern of gene expression across samples within that module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


# Calculate the gene significance and associated p-values
gene.signf.corr.delta <- cor(norm.counts, traits$Delta, use = 'p')
gene.signf.corr.delta.pvals <- corPvalueStudent(gene.signf.corr.delta, nSamples)
# The top 25 genes based on pvalue that are associated with delta
hub_genes_delta <- gene.signf.corr.delta.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  #rownames(.) %>%
  head(25)

# Calculate the gene significance and associated p-values
gene.signf.corr.omicron <- cor(norm.counts, traits$Omicron, use = 'p')
gene.signf.corr.omicron.pvals <- corPvalueStudent(gene.signf.corr.omicron, nSamples)
# The top 25 genes based on pvalue that are associated with omicron
hub_genes_omicron <- gene.signf.corr.omicron.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  #rownames(.) %>%
  head(50)


GS.trait_delta = as.data.frame(cor(norm.counts, traits$Delta, use = "p"))
GS.trait_omicron = as.data.frame(cor(norm.counts, traits$Omicron, use = "p"))
moduleColors <- labels2colors(bwnet$colors)
MEs = moduleEigengenes(norm.counts, moduleColors)$eigengenes
MM.trait = signedKME(norm.counts, MEs)

#par(mfrow = c(1, 2))
# verboseScatterplot(MM.trait$kMEyellow, GS.trait_delta$V1, xlab = "Module Membership", 
#                      ylab = "Gene Significance for Delta", main = paste("kME.", module,"vs. GS"), col = "yellow")
# abline(h = 0.60, col = "black", lty = 2)  # Adjust color and line type as needed
# abline(v = 0.80, col = "black", lty = 2)  # Adjust color and line type as needed
# 
# verboseScatterplot(MM.trait$kMElightgreen, GS.trait_omicron$V1, xlab = "Module Membership", 
#                    ylab = "Gene Significance for Omicron", main = paste("kME.", module,"vs. GS"), col = "lightgreen")  
# abline(h = 0.60, col = "black", lty = 2)  # Adjust color and line type as needed
# abline(v = 0.80, col = "black", lty = 2)  # Adjust color and line type as needed


genes_of_interest <- cbind(GS.trait_delta,MM.trait$kMEbrown,GS.trait_omicron,MM.trait$kMEyellow)
colnames(genes_of_interest) <- c("GS_delta","MM_brown","GS_omicron","MM_yellow")

# Hub genes for yellow module with delta as trait
genes_of_interest_delta <- genes_of_interest %>% 
  dplyr::select("GS_delta","MM_brown")
genes_of_interest_delta <- genes_of_interest_delta %>%
  filter(abs(GS_delta) >= 0.6,      
         abs(MM_brown) >= 0.8) 

# Hub genes for turquoise module with omicron as trait
genes_of_interest_omicron <- genes_of_interest %>% 
  dplyr::select("GS_omicron","MM_yellow")
genes_of_interest_omicron <- genes_of_interest_omicron %>%
  filter(abs(GS_omicron) >= 0.45,      
         abs(MM_yellow) >= 0.8)    


delta_hubs <- genes_of_interest_delta

hubgenes_delta_1dpi <- intersect(diff_genes_1_Delta$geneName,rownames(delta_hubs))

hubgenes_delta_2dpi <- intersect(diff_genes_2_Delta$geneName,rownames(delta_hubs))

hubgenes_delta_3dpi <- intersect(diff_genes_3_Delta$geneName,rownames(delta_hubs))

#hubgenes_delta_all <-  unique(c(hubgenes_delta_1dpi, hubgenes_delta_2dpi, hubgenes_delta_3dpi))


omicron_hubs <- genes_of_interest_omicron

hubgenes_omicron_1dpi <- intersect(diff_genes_1_Omicron$geneName,rownames(omicron_hubs))

hubgenes_omicron_2dpi <- intersect(diff_genes_2_Omicron$geneName,rownames(omicron_hubs))

hubgenes_omicron_3dpi <- intersect(diff_genes_3_Omicron$geneName,rownames(omicron_hubs))

#hubgenes_omicron_all <- unique(c(hubgenes_omicron_1dpi,hubgenes_omicron_3dpi,hubgenes_omicron_3dpi))


# Enrichment of hubgenes
library(enrichR)
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")

enriched_hub_omicron_1 <- enrichr(hubgenes_omicron_1dpi, dbs)
enriched_hub_omicron_1_BP <- as.data.frame(enriched_hub_omicron_1[["GO_Biological_Process_2021"]])

e_hub_o1 <- plotEnrich(
  enriched_hub_omicron_1_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for hubgenes for omicron at 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_hub_o1

enriched_hub_omicron_2 <- enrichr(hubgenes_omicron_2dpi, dbs)
enriched_hub_omicron_2_BP <- as.data.frame(enriched_hub_omicron_2[["GO_Biological_Process_2021"]])

e_hub_o2 <- plotEnrich(
  enriched_hub_omicron_2_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for hubgenes for omicron at 2 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_hub_o2

enriched_hub_omicron_3 <- enrichr(hubgenes_omicron_3dpi, dbs)
enriched_hub_omicron_3_BP <- as.data.frame(enriched_hub_omicron_3[["GO_Biological_Process_2021"]])

e_hub_o3 <- plotEnrich(
  enriched_hub_omicron_3_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for hubgenes for omicron at 3 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_hub_o3


enriched_hub_delta_1 <- enrichr(hubgenes_delta_1dpi, dbs)
enriched_hub_delta_1_BP <- as.data.frame(enriched_hub_delta_1[["GO_Biological_Process_2021"]])

e_hub_d1 <- plotEnrich(
  enriched_hub_delta_1_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for hubgenes for delta at 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_hub_d1

enriched_hub_delta_2 <- enrichr(hubgenes_delta_2dpi, dbs)
enriched_hub_delta_2_BP <- as.data.frame(enriched_hub_delta_2[["GO_Biological_Process_2021"]])

e_hub_d2 <- plotEnrich(
  enriched_hub_delta_2_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for hubgenes for delta at 2 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_hub_d2

enriched_hub_delta_3 <- enrichr(hubgenes_delta_3dpi, dbs)
enriched_hub_delta_3_BP <- as.data.frame(enriched_hub_delta_3[["GO_Biological_Process_2021"]])

e_hub_d3 <- plotEnrich(
  enriched_hub_delta_3_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for hubgenes for delta at 3 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_hub_d3


enrich_delta <- grid.arrange(e_hub_d1,
                               e_hub_d2,
                               e_hub_d3,
                               nrow=3)


enrich_omicron <- grid.arrange(e_hub_o1,
                               e_hub_o2,
                               e_hub_o3,
                               nrow=3)


######################################################################################################################
# Save the current R session
#save.image(file = "deseq2_human.RData")

######################################################################################################################
#save.image(file = "deseq2_human_wo_wuhan.RData")












##########################################################################
#########  BOX PLOT at 1 DPI (TOP 10 up and Top 10 down) #########
################### with highest Fold Change in fpkm #####################
#################### avg between control and VOCs ########################
##########################################################################

# top upregulated genes in Omicron but not DE in other two variants
# top_diff_Omicron_up <- diff_genes_1_Omicron %>%
#   filter(geneName %in% diff_genes_1_Delta$geneName & 
#            geneName %in% diff_genes_1_Wuhan$geneName & 
#            log2FoldChange > 0) %>%
#   arrange(desc(log2FoldChange)) 


hubgenes_delta_all
hubgenes_omicron_all

# Define the order of GO terms based on their ranking for each variant
top_diff_Omicron_up <- diff_genes_1_Omicron %>% 
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

top_diff_Omicron_down <- diff_genes_1_Omicron %>%
  filter(log2FoldChange < 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

# Define the order of GO terms based on their ranking for each variant
top_diff_Delta_up <- diff_genes_1_Delta %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)
top_diff_Delta_down <- diff_genes_1_Delta %>%
  filter(log2FoldChange < 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

# Define the order of GO terms based on their ranking for each variant
top_diff_Wuhan_up <- diff_genes_1_Wuhan %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)
top_diff_Wuhan_down <- diff_genes_1_Wuhan %>%
  filter(log2FoldChange < 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)




######### Boxplot prep 1 DPI ###########
normalized_counts_1 <- counts(dds_variant_time_1, normalized =T) 
# Convert to dataframe
fpkm_human_1 <- as.data.frame(normalized_counts_1)
# Transpose dataframe to create boxplot input
fpkm_human_t_1<- fpkm_human_1 %>% 
  t() %>% 
  as.data.frame() 

# Adding SRR values to new column to reshape count matrix
fpkm_human_t_1$Run <- rownames(fpkm_human_t_1)

# Melt the data
fpkm_human_t_1_melted <- reshape2::melt(fpkm_human_t_1, id.vars = "Run")

# Rename columns
colnames(fpkm_human_t_1_melted) <- c("Run", "Gene", "value")

# Add the variant column based on matching values between Run in fpkm_human_t_1_melted and metadata_ncbi
fpkm_human_t_1_melted$variant <- metadata_ncbi$virus[match(fpkm_human_t_1_melted$Run, metadata_ncbi$Run)]

# Reorder the 'variant' factor levels
fpkm_human_t_1_melted$variant <- factor(
  fpkm_human_t_1_melted$variant,
  levels = c("mock","Wuhan", "delta", "omicron")
)




# Reorder the levels of the variable factor for each variant
fpkm_omicron_up <- fpkm_human_t_1_melted %>%
  filter(Gene %in% top_diff_Omicron_up$geneName)
fpkm_omicron_up$Gene <- factor(
  fpkm_omicron_up$Gene,
  levels = top_diff_Omicron_up$geneName
)

fpkm_omicron_down <- fpkm_human_t_1_melted %>%
  filter(Gene %in% top_diff_Omicron_down$geneName)
fpkm_omicron_down$Gene <- factor(
  fpkm_omicron_down$Gene,
  levels = top_diff_Omicron_down$geneName
)


# Reorder the levels of the variable factor for each variant
fpkm_delta_up <- fpkm_human_t_1_melted %>%
  filter(Gene %in% top_diff_Delta_up$geneName)
fpkm_delta_up$Gene <- factor(
  fpkm_delta_up$Gene,
  levels = top_diff_Delta_up$geneName
)

fpkm_delta_down <- fpkm_human_t_1_melted %>%
  filter(Gene %in% top_diff_Delta_down$geneName)
fpkm_delta_down$Gene <- factor(
  fpkm_delta_down$Gene,
  levels = top_diff_Delta_down$geneName
)


# Reorder the levels of the variable factor for each variant
fpkm_wuhan_up <- fpkm_human_t_1_melted %>%
  filter(Gene %in% top_diff_Wuhan_up$geneName)
fpkm_wuhan_up$Gene <- factor(
  fpkm_wuhan_up$Gene,
  levels = top_diff_Wuhan_up$geneName
)

fpkm_wuhan_down <- fpkm_human_t_1_melted %>%
  filter(Gene %in% top_diff_Wuhan_down$geneName)
fpkm_wuhan_down$Gene <- factor(
  fpkm_wuhan_down$Gene,
  levels = top_diff_Wuhan_down$geneName
)


######### Top 10 1 dpi ############

# Set color palette
my_colors_c <- c("dodgerblue","forestgreen","goldenrod","#FC8D62")

##### Top GENES for highest FC in DE genes in Omicron

plot_top10_Omicron_up_1 <- fpkm_human_t_1_melted %>%
  filter(Gene %in% fpkm_omicron_up$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Omicron_up$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most upregulated host genes
       in human cells infected with Omicron at 1 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Omicron_up_1



plot_top10_Omicron_down_1 <- fpkm_human_t_1_melted %>%
  filter(Gene %in% fpkm_omicron_down$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Omicron_down$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most downregulated host genes
       in human cells infected with Omicron at 1 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Omicron_down_1


Omicron_combined_1 <- grid.arrange(plot_top10_Omicron_up_1, plot_top10_Omicron_down_1, ncol = 1, heights = c(1,1.2))




##### Top GENES for highest FC in DE genes in Delta

plot_top10_Delta_up_1 <- fpkm_human_t_1_melted %>%
  filter(Gene %in% fpkm_delta_up$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Delta_up$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most upregulated host genes
       in human cells infected with Delta at 1 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Delta_up_1



plot_top10_Delta_down_1 <- fpkm_human_t_1_melted %>%
  filter(Gene %in% fpkm_delta_down$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Delta_down$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most downregulated host genes
       in human cells infected with Delta at 1 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Delta_down_1


Delta_combined_1 <- grid.arrange(plot_top10_Delta_up_1, plot_top10_Delta_down_1, ncol = 1, heights = c(1,1.2))





##### Top GENES for highest FC in DE genes in Wuhan

plot_top10_Wuhan_up_1 <- fpkm_human_t_1_melted %>%
  filter(Gene %in% fpkm_wuhan_up$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Wuhan_up$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most upregulated host genes
       in human cells infected with Wuhan at 1 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Wuhan_up_1



plot_top10_Wuhan_down_1 <- fpkm_human_t_1_melted %>%
  filter(Gene %in% fpkm_wuhan_down$Gene) %>%
  filter(!Gene %in% c("RTEL1","RNF224")) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Wuhan_down$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most downregulated host genes
       in human cells infected with Wuhan at 1 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Wuhan_down_1



# Create two empty plots
empty_plot <- ggplot() + theme_void()
#plot_top10_Wuhan_down_1 <- arrangeGrob(empty_plot, plot_top10_Wuhan_down_1, empty_plot, ncol = 3, widths = c(.1,.8,.1))

Wuhan_combined_1 <- grid.arrange(plot_top10_Wuhan_up_1, plot_top10_Wuhan_down_1, ncol = 1, heights = c(1,1.2))






##########################################################################
#########  BOX PLOT at 2 DPI (TOP 10 up and Top 10 down) #########
################### with highest Fold Change in fpkm #####################
#################### avg between control and VOCs ########################
##########################################################################

# top upregulated genes in Omicron but not DE in other two variants
# top_diff_Omicron_up_2 <- diff_genes_2_Omicron %>%
#   filter(geneName %in% diff_genes_2_Delta$geneName & 
#            geneName %in% diff_genes_2_Wuhan$geneName & 
#            log2FoldChange > 0) %>%
#   arrange(desc(log2FoldChange)) 


# Define the order of GO terms based on their ranking for each variant
top_diff_Omicron_up_2 <- diff_genes_2_Omicron %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)
top_diff_Omicron_down_2 <- diff_genes_2_Omicron %>%
  filter(log2FoldChange < 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

# Define the order of GO terms based on their ranking for each variant
top_diff_Delta_up_2 <- diff_genes_2_Delta %>%
  filter(log2FoldChange > 0 & geneName != "RHOJ" & geneName != "ISY1-RAB43") %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)
top_diff_Delta_down_2 <- diff_genes_2_Delta %>%
  filter(log2FoldChange < 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

# Define the order of GO terms based on their ranking for each variant
top_diff_Wuhan_up_2 <- diff_genes_2_Wuhan %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)
top_diff_Wuhan_down_2 <- diff_genes_2_Wuhan %>%
  filter(log2FoldChange < 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)




######### Boxplot prep 2 DPI ###########
normalized_counts_2 <- counts(dds_variant_time_2, normalized =T) 
# Convert to dataframe
fpkm_human_2 <- as.data.frame(normalized_counts_2)
# Transpose dataframe to create boxplot input
fpkm_human_t_2<- fpkm_human_2 %>% 
  t() %>% 
  as.data.frame() 

# Adding SRR values to new column to reshape count matrix
fpkm_human_t_2$Run <- rownames(fpkm_human_t_2)

# Melt the data
fpkm_human_t_2_melted <- reshape2::melt(fpkm_human_t_2, id.vars = "Run")

# Rename columns
colnames(fpkm_human_t_2_melted) <- c("Run", "Gene", "value")

# Add the variant column based on matching values between Run in fpkm_human_t_2_melted and metadata_ncbi
fpkm_human_t_2_melted$variant <- metadata_ncbi$virus[match(fpkm_human_t_2_melted$Run, metadata_ncbi$Run)]

# Reorder the 'variant' factor levels
fpkm_human_t_2_melted$variant <- factor(
  fpkm_human_t_2_melted$variant,
  levels = c("mock","Wuhan", "delta", "omicron")
)




# Reorder the levels of the variable factor for each variant
fpkm_omicron_up_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% top_diff_Omicron_up_2$geneName)
fpkm_omicron_up_2$Gene <- factor(
  fpkm_omicron_up_2$Gene,
  levels = top_diff_Omicron_up_2$geneName
)


fpkm_omicron_down_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% top_diff_Omicron_down_2$geneName)
fpkm_omicron_down_2$Gene <- factor(
  fpkm_omicron_down_2$Gene,
  levels = top_diff_Omicron_down_2$geneName
)


# Reorder the levels of the variable factor for each variant
fpkm_delta_up_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% top_diff_Delta_up_2$geneName)
fpkm_delta_up_2$Gene <- factor(
  fpkm_delta_up_2$Gene,
  levels = top_diff_Delta_up_2$geneName
)

fpkm_delta_down_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% top_diff_Delta_down_2$geneName)
fpkm_delta_down_2$Gene <- factor(
  fpkm_delta_down_2$Gene,
  levels = top_diff_Delta_down_2$geneName
)


# Reorder the levels of the variable factor for each variant
fpkm_wuhan_up_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% top_diff_Wuhan_up_2$geneName)
fpkm_wuhan_up_2$Gene <- factor(
  fpkm_wuhan_up_2$Gene,
  levels = top_diff_Wuhan_up_2$geneName
)

fpkm_wuhan_down_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% top_diff_Wuhan_down_2$geneName)
fpkm_wuhan_down_2$Gene <- factor(
  fpkm_wuhan_down_2$Gene,
  levels = top_diff_Wuhan_down_2$geneName
)

######### Top 10 2 dpi ############

# Set color palette
my_colors_c <- c("dodgerblue","forestgreen","goldenrod","#FC8D62")

##### Top GENES for highest FC in DE genes in Omicron

plot_top10_Omicron_up_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% fpkm_omicron_up_2$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Omicron_up_2$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most upregulated host genes
       in human cells infected with Omicron at 2 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Omicron_up_2



plot_top10_Omicron_down_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% fpkm_omicron_down_2$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Omicron_down_2$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most downregulated host genes
       in human cells infected with Omicron at 2 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Omicron_down_2


Omicron_combined_2 <- grid.arrange(plot_top10_Omicron_up_2, plot_top10_Omicron_down_2, ncol = 1, heights = c(1,1.2))




##### Top GENES for highest FC in DE genes in Delta

plot_top10_Delta_up_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% fpkm_delta_up_2$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Delta_up_2$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most upregulated host genes
       in human cells infected with Delta at 2 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Delta_up_2



plot_top10_Delta_down_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% fpkm_delta_down_2$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Delta_down_2$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most downregulated host genes
       in human cells infected with Delta at 2 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Delta_down_2


Delta_combined_2 <- grid.arrange(plot_top10_Delta_up_2, plot_top10_Delta_down_2, ncol = 1, heights = c(1,1.2))





##### Top GENES for highest FC in DE genes in Wuhan

plot_top10_Wuhan_up_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% fpkm_wuhan_up_2$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Wuhan_up_2$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most upregulated host genes
       in human cells infected with Wuhan at 2 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Wuhan_up_2



plot_top10_Wuhan_down_2 <- fpkm_human_t_2_melted %>%
  filter(Gene %in% fpkm_wuhan_down_2$Gene) %>%
  #filter(!Gene %in% c("RTEL1","RNF224")) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Wuhan_down_2$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most downregulated host genes
       in human cells infected with Wuhan at 2 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Wuhan_down_2



# Create two empty plots
empty_plot <- ggplot() + theme_void()
#plot_top10_Wuhan_up_2 <- arrangeGrob(empty_plot, plot_top10_Wuhan_up_2, empty_plot, ncol = 3, widths = c(.15,.7,.15))

Wuhan_combined_2 <- grid.arrange(plot_top10_Wuhan_up_2, plot_top10_Wuhan_down_2, ncol = 1, heights = c(1,1.2))









##########################################################################
#########  BOX PLOT at 3 DPI (TOP 10 up and Top 10 down) #########
################### with highest Fold Change in fpkm #####################
#################### avg between control and VOCs ########################
##########################################################################

# top upregulated genes in Omicron but not DE in other two variants
# top_diff_Omicron_up_3 <- diff_genes_3_Omicron %>%
#   filter(geneName %in% diff_genes_3_Delta$geneName & 
#            geneName %in% diff_genes_3_Wuhan$geneName & 
#            log2FoldChange > 0) %>%
#   arrange(desc(log2FoldChange)) 


# Define the order of GO terms based on their ranking for each variant
top_diff_Omicron_up_3 <- diff_genes_3_Omicron %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)
top_diff_Omicron_down_3 <- diff_genes_3_Omicron %>%
  filter(log2FoldChange < 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

# Define the order of GO terms based on their ranking for each variant
top_diff_Delta_up_3 <- diff_genes_3_Delta %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)
top_diff_Delta_down_3 <- diff_genes_3_Delta %>%
  filter(log2FoldChange < 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

# Define the order of GO terms based on their ranking for each variant
top_diff_Wuhan_up_3 <- diff_genes_3_Wuhan %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)
top_diff_Wuhan_down_3 <- diff_genes_3_Wuhan %>%
  filter(log2FoldChange < 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)




######### Boxplot prep 3 DPI ###########
normalized_counts_3 <- counts(dds_variant_time_3, normalized =T) 
# Convert to dataframe
fpkm_human_3 <- as.data.frame(normalized_counts_3)
# Transpose dataframe to create boxplot input
fpkm_human_t_3<- fpkm_human_3 %>% 
  t() %>% 
  as.data.frame() 

# Adding SRR values to new column to reshape count matrix
fpkm_human_t_3$Run <- rownames(fpkm_human_t_3)

# Melt the data
fpkm_human_t_3_melted <- reshape2::melt(fpkm_human_t_3, id.vars = "Run")

# Rename columns
colnames(fpkm_human_t_3_melted) <- c("Run", "Gene", "value")

# Add the variant column based on matching values between Run in fpkm_human_t_3_melted and metadata_ncbi
fpkm_human_t_3_melted$variant <- metadata_ncbi$virus[match(fpkm_human_t_3_melted$Run, metadata_ncbi$Run)]

# Reorder the 'variant' factor levels
fpkm_human_t_3_melted$variant <- factor(
  fpkm_human_t_3_melted$variant,
  levels = c("mock","Wuhan", "delta", "omicron")
)




# Reorder the levels of the variable factor for each variant
fpkm_omicron_up_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% top_diff_Omicron_up_3$geneName)
fpkm_omicron_up_3$Gene <- factor(
  fpkm_omicron_up_3$Gene,
  levels = top_diff_Omicron_up_3$geneName
)


fpkm_omicron_down_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% top_diff_Omicron_down_3$geneName)
fpkm_omicron_down_3$Gene <- factor(
  fpkm_omicron_down_3$Gene,
  levels = top_diff_Omicron_down_3$geneName
)


# Reorder the levels of the variable factor for each variant
fpkm_delta_up_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% top_diff_Delta_up_3$geneName)
fpkm_delta_up_3$Gene <- factor(
  fpkm_delta_up_3$Gene,
  levels = top_diff_Delta_up_3$geneName
)

fpkm_delta_down_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% top_diff_Delta_down_3$geneName)
fpkm_delta_down_3$Gene <- factor(
  fpkm_delta_down_3$Gene,
  levels = top_diff_Delta_down_3$geneName
)


# Reorder the levels of the variable factor for each variant
fpkm_wuhan_up_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% top_diff_Wuhan_up_3$geneName)
fpkm_wuhan_up_3$Gene <- factor(
  fpkm_wuhan_up_3$Gene,
  levels = top_diff_Wuhan_up_3$geneName
)

fpkm_wuhan_down_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% top_diff_Wuhan_down_3$geneName)
fpkm_wuhan_down_3$Gene <- factor(
  fpkm_wuhan_down_3$Gene,
  levels = top_diff_Wuhan_down_3$geneName
)


######### Top 10 3 dpi ############

# Set color palette
my_colors_c <- c("dodgerblue","forestgreen","goldenrod","#FC8D62")

##### Top GENES for highest FC in DE genes in Omicron

plot_top10_Omicron_up_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% fpkm_omicron_up_3$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Omicron_up_3$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most upregulated host genes
       in human cells infected with Omicron at 3 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Omicron_up_3



plot_top10_Omicron_down_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% fpkm_omicron_down_3$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Omicron_down_3$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most downregulated host genes
       in human cells infected with Omicron at 3 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Omicron_down_3


Omicron_combined_3 <- grid.arrange(plot_top10_Omicron_up_3, plot_top10_Omicron_down_3, ncol = 1, heights = c(1,1.2))




##### Top GENES for highest FC in DE genes in Delta

plot_top10_Delta_up_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% fpkm_delta_up_3$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Delta_up_3$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most upregulated host genes
       in human cells infected with Delta at 3 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Delta_up_3



plot_top10_Delta_down_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% fpkm_delta_down_3$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Delta_down_3$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most downregulated host genes
       in human cells infected with Delta at 3 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Delta_down_3


Delta_combined_3 <- grid.arrange(plot_top10_Delta_up_3, plot_top10_Delta_down_3, ncol = 1, heights = c(1,1.2))





##### Top GENES for highest FC in DE genes in Wuhan

plot_top10_Wuhan_up_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% fpkm_wuhan_up_3$Gene) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Wuhan_up_3$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most upregulated host genes
       in human cells infected with Wuhan at 3 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Wuhan_up_3



plot_top10_Wuhan_down_3 <- fpkm_human_t_3_melted %>%
  filter(Gene %in% fpkm_wuhan_down_3$Gene) %>%
  #filter(!Gene %in% c("RTEL1","RNF224")) %>%
  group_by(Gene) %>%
  mutate(Gene = factor(Gene, levels = top_diff_Wuhan_down_3$geneName)) %>%
  #mutate(value2 = filter_lims(as.numeric(value))) %>%
  ggplot(aes(x = variant, y = value, color = variant)) +
  geom_boxplot(na.rm = TRUE, coef = 5, outlier.shape = NA, alpha = 0.5) +  # Set fill = NA for boxes
  scale_color_manual(values = my_colors_c) +  # Set point colors
  facet_wrap(~ Gene, scales = "free_y", ncol = 10, labeller = label_wrap_gen(width = 50)) +
  labs(x = NULL, y = "Normalized read count", title = "Normalized expression profile across all variants for most downregulated host genes
       in human cells infected with Wuhan at 3 DPI") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10, face = "bold"))
plot_top10_Wuhan_down_3



# Create two empty plots
empty_plot <- ggplot() + theme_void()
#plot_top10_Wuhan_up_3 <- arrangeGrob(empty_plot, plot_top10_Wuhan_up_3, empty_plot, ncol = 3, widths = c(.15,.7,.15))

Wuhan_combined_3 <- grid.arrange(plot_top10_Wuhan_up_3, plot_top10_Wuhan_down_3, ncol = 1, heights = c(1,1.2))











######## gene list for two relevant modules ##########
genes_module_color_df <- as.data.frame(bwnet$colors) 
genes_module_color_df <- genes_module_color_df %>% 
                          rownames_to_column(.) 
colnames(genes_module_color_df) <- c("gene", "color")

genes_module_color_df_yellow <- genes_module_color_df %>% filter(color == "yellow")
genes_module_color_df_turquoise <- genes_module_color_df %>% filter(color == "turquoise")

genes_module_color_df_yellow_turquoise <- rbind(genes_module_color_df_yellow,genes_module_color_df_turquoise)










######################################################################################################################
################## HEATMAP 1 DPI #####################################################################################
diff_genes_1_heatmap <- log2FoldChange_table_1 %>% filter(abs(log2FoldChange) > 2  & baseMean >50 & padj < 0.05)
#diff_genes_1_heatmap <- diff_genes_1_heatmap[diff_genes_1_heatmap$geneName %in% hs_immune_genes$gene_symbol, ]

#diff_genes_1_heatmap <- diff_genes_1_heatmap[diff_genes_1_heatmap$geneName %in% hs_immune_genes$gene_symbol, ]
dds_variant_time_1 <- estimateSizeFactors(dds_variant_time_1)


# Adding ENTREZ ID for differentially expressed geneNames
diff_genes_1_heatmap <- diff_genes_1_heatmap %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_1_heatmap$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))

matrix_1 <- counts(dds_variant_time_1, normalized =T) [unique(diff_genes_1_heatmap$geneName),]

# Update column names of the filtered matrix with the 'variant' values
colnames(matrix_1) <- design_matrix_1[colnames(matrix_1), "variant"]


# Scale function performs a z-score normalization on each row of matrix_1. 
# This means for each row, the mean of the row is subtracted from each element, 
# and the result is divided by the standard deviation of the row.
matrix.z_1 <- t(apply(matrix_1,1,scale))
colnames(matrix.z_1) <- colnames(matrix_1)



matrix.z_1 <- matrix.z_1[, order(match(colnames(matrix.z_1), c("mock", "Wuhan", "delta", "omicron")))] 

heatmap_1 <- Heatmap(matrix.z_1,
                     cluster_rows = T,
                     cluster_columns = F,
                     #column_order = order(colnames(matrix.z_1)),
                     column_labels = colnames(matrix.z_1),
                     column_names_gp = gpar(fontsize = 0),
                     column_names_side = "top",
                     column_names_rot = 45,
                     #column_split = 4,
                     #row_split = 3,
                     row_km = 3, # cluster based on kmean
                     column_split = colnames(matrix.z_1),
                     row_labels = rownames(matrix.z_1), 
                     row_names_gp = gpar(fontsize = 0),
                     name="norm. z-score")
heatmap_1



# Obtain the ordered row indices (genes) from heatmap_1
clusterlist_1 <- row_order(heatmap_1)

# Create a list to store genes for each cluster
cluster_genes_list_1 <- list()

# Extract gene names based on the ordered row indices per cluster
for (cluster_id in names(clusterlist_1)) {
  # Get the ordered indices of genes for the current cluster
  ordered_indices <- clusterlist_1[[cluster_id]]
  
  # Get the corresponding gene names from matrix.z_1
  cluster_genes <- rownames(matrix.z_1)[ordered_indices]
  
  # Store the genes in the cluster_genes_list_1 with cluster_id as the key
  cluster_genes_list_1[[cluster_id]] <- cluster_genes
}


# Create an empty data frame to store the gene IDs and clusters
clu_df_1 <- data.frame(GeneID = character(), Cluster = character(), stringsAsFactors = FALSE)
# Iterate over each cluster and its corresponding gene list
for (cluster_id in names(cluster_genes_list_1)) {
  # Get the gene IDs for the current cluster
  gene_ids <- cluster_genes_list_1[[cluster_id]]
  
  # Create a data frame for the current cluster
  cluster_df <- data.frame(GeneID = gene_ids, Cluster = cluster_id)
  
  # Append the cluster_df to clu_df_1
  clu_df_1 <- rbind(clu_df_1, cluster_df)
}


clu_df_1 <- clu_df_1 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, clu_df_1$GeneID, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))

length(clu_df_1$GeneID[clu_df_1$Cluster == "1"])
length(clu_df_1$GeneID[clu_df_1$Cluster == "2"])
length(clu_df_1$GeneID[clu_df_1$Cluster == "3"])
#length(clu_df_1$GeneID[clu_df_1$Cluster == "4"])

library(enrichR)
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")


# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
clu_df_1_genes_c1 <- unique(clu_df_1$GeneID [clu_df_1$Cluster == 1])

enriched_1_c1 <- enrichr(clu_df_1_genes_c1, dbs)
enriched_1_c1_BP <- as.data.frame(enriched_1_c1[["GO_Biological_Process_2021"]])

e_1_c1 <- plotEnrich(
  enriched_1_c1_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for Cluster 1 genes 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_1_c1

# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
clu_df_1_genes_c2 <- unique(clu_df_1$GeneID [clu_df_1$Cluster == 2])

enriched_1_c2 <- enrichr(clu_df_1_genes_c2, dbs)
enriched_1_c2_BP <- as.data.frame(enriched_1_c2[["GO_Biological_Process_2021"]])

e_1_c2 <- plotEnrich(
  enriched_1_c2_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for Cluster 2 genes 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_1_c2

# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
clu_df_1_genes_c3 <- unique(clu_df_1$GeneID [clu_df_1$Cluster == 3])

enriched_1_c3 <- enrichr(clu_df_1_genes_c3, dbs)
enriched_1_c3_BP <- as.data.frame(enriched_1_c3[["GO_Biological_Process_2021"]])

e_1_c3 <- plotEnrich(
  enriched_1_c3_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for Cluster 3 genes 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_1_c3




# Define cluster functions
cluster_functions_1 = list(
  `3` = strwrap(e_1_c3$data$Term, width = 60),
  `2` = strwrap(e_1_c2$data$Term, width = 60),
  `1` = strwrap(e_1_c1$data$Term, width = 60)
)                


# note how we set the width of this empty annotation
ha = rowAnnotation(foo = anno_empty(border = FALSE, 
                                    width = max_text_width(unlist(cluster_functions_1)) + unit(3, "mm")))

Heatmap(matrix.z_1,
        cluster_rows = T,
        cluster_columns = F,
        #column_order = order(colnames(matrix.z_1)),
        column_labels = colnames(matrix.z_1),
        column_names_gp = gpar(fontsize = 0),
        column_names_side = "top",
        column_names_rot = 45,
        #column_split = 4,
        #row_split = 4,
        row_km = 3, # cluster based on kmean
        column_split = colnames(matrix.z_1),
        row_labels = rownames(matrix.z_1), 
        row_names_gp = gpar(fontsize = 0),
        name="norm. z-score", 
        right_annotation = ha)

for(i in 1:3) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(cluster_functions_1[[i]], collapse = "\n"), x = unit(3, "mm"), just = "left")
  })
}


######################################################################################################################

################## HEATMAP 2 DPI #####################################################################################
diff_genes_2_heatmap <- log2FoldChange_table_2 %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
diff_genes_2_heatmap <- diff_genes_2_heatmap[diff_genes_2_heatmap$geneName %in% hs_immune_genes$gene_symbol, ]
dds_variant_time_2 <- estimateSizeFactors(dds_variant_time_2)

# Adding ENTREZ ID for differentially expressed geneNames
diff_genes_2_heatmap <- diff_genes_2_heatmap %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_2_heatmap$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))


matrix_2 <- counts(dds_variant_time_2, normalized =T) [genes_module_color_df_yellow_turquoise$gene,]

# Update column names of the filtered matrix with the 'variant' values
colnames(matrix_2) <- design_matrix_2[colnames(matrix_2), "variant"]

matrix.z_2 <- t(apply(matrix_2,1,scale))
colnames(matrix.z_2) <- colnames(matrix_2)

matrix.z_2 <- matrix.z_2[, order(match(colnames(matrix.z_2), c("mock", "Wuhan", "delta", "omicron")))] 

heatmap_2 <- Heatmap(matrix.z_2,
                     cluster_rows = T,
                     cluster_columns = F,
                     #column_order = order(colnames(matrix.z_2)),
                     column_labels = colnames(matrix.z_2),
                     column_names_gp = gpar(fontsize = 0),
                     column_names_side = "top",
                     column_names_rot = 45,
                     #column_split = 4,
                     #row_split = 3,
                     row_km = 2, # cluster based on kmean
                     column_split = colnames(matrix.z_2),
                     row_labels = rownames(matrix.z_2), 
                     row_names_gp = gpar(fontsize = 0),
                     name="norm. z-score")
heatmap_2

# Obtain the ordered row indices (genes) from heatmap_2
clusterlist_2 <- row_order(heatmap_2)

# Create a list to store genes for each cluster
cluster_genes_list_2 <- list()

# Extract gene names based on the ordered row indices per cluster
for (cluster_id in names(clusterlist_2)) {
  # Get the ordered indices of genes for the current cluster
  ordered_indices <- clusterlist_2[[cluster_id]]
  
  # Get the corresponding gene names from matrix.z_2
  cluster_genes <- rownames(matrix.z_2)[ordered_indices]
  
  # Store the genes in the cluster_genes_list_2 with cluster_id as the key
  cluster_genes_list_2[[cluster_id]] <- cluster_genes
}


# Create an empty data frame to store the gene IDs and clusters
clu_df_2 <- data.frame(GeneID = character(), Cluster = character(), stringsAsFactors = FALSE)
# Iterate over each cluster and its corresponding gene list
for (cluster_id in names(cluster_genes_list_2)) {
  # Get the gene IDs for the current cluster
  gene_ids <- cluster_genes_list_2[[cluster_id]]
  
  # Create a data frame for the current cluster
  cluster_df <- data.frame(GeneID = gene_ids, Cluster = cluster_id)
  
  # Append the cluster_df to clu_df_2
  clu_df_2 <- rbind(clu_df_2, cluster_df)
}


clu_df_2 <- clu_df_2 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, clu_df_2$GeneID, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))

length(clu_df_2$GeneID[clu_df_2$Cluster == "1"])
length(clu_df_2$GeneID[clu_df_2$Cluster == "2"])
#length(clu_df_2$GeneID[clu_df_2$Cluster == "3"])



#library(enrichR)
#dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")


# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
clu_df_2_genes_c1 <- unique(clu_df_2$GeneID [clu_df_2$Cluster == 1])

enriched_2_c1 <- enrichr(clu_df_2_genes_c1, dbs)
enriched_2_c1_BP <- as.data.frame(enriched_2_c1[["GO_Biological_Process_2021"]])

e_2_c1 <- plotEnrich(
  enriched_2_c1_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for Cluster 1 genes 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_2_c1

# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
clu_df_2_genes_c2 <- unique(clu_df_2$GeneID [clu_df_2$Cluster == 2])

enriched_2_c2 <- enrichr(clu_df_2_genes_c2, dbs)
enriched_2_c2_BP <- as.data.frame(enriched_2_c2[["GO_Biological_Process_2021"]])

e_2_c2 <- plotEnrich(
  enriched_2_c2_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for Cluster 2 genes 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_2_c2

# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
clu_df_2_genes_c3 <- unique(clu_df_2$GeneID [clu_df_2$Cluster == 3])

enriched_2_c3 <- enrichr(clu_df_2_genes_c3, dbs)
enriched_2_c3_BP <- as.data.frame(enriched_2_c3[["GO_Biological_Process_2021"]])

e_2_c3 <- plotEnrich(
  enriched_2_c3_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for Cluster 3 genes 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_2_c3




# Define cluster functions
cluster_functions_2 = list(
  #`3` = strwrap(e_2_c3$data$Term, width = 60),
  `1` = strwrap(e_2_c1$data$Term, width = 60),
  `2` = strwrap(e_2_c2$data$Term, width = 60)
)                


# note how we set the width of this empty annotation
ha = rowAnnotation(foo = anno_empty(border = FALSE, 
                                    width = max_text_width(unlist(cluster_functions_2)) + unit(3, "mm")))

Heatmap(matrix.z_2,
        cluster_rows = T,
        cluster_columns = F,
        #column_order = order(colnames(matrix.z_2)),
        column_labels = colnames(matrix.z_2),
        column_names_gp = gpar(fontsize = 0),
        column_names_side = "top",
        column_names_rot = 45,
        #column_split = 4,
        #row_split = 4,
        row_km = 2, # cluster based on kmean
        column_split = colnames(matrix.z_2),
        row_labels = rownames(matrix.z_2), 
        row_names_gp = gpar(fontsize = 0),
        name="norm. z-score", 
        right_annotation = ha)

for(i in 1:3) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(cluster_functions_2[[i]], collapse = "\n"), x = unit(3, "mm"), just = "left")
  })
}








################## HEATMAP 3 DPI #####################################################################################
diff_genes_3_heatmap <- log2FoldChange_table_3 %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
diff_genes_3_heatmap <- diff_genes_3_heatmap[diff_genes_3_heatmap$geneName %in% hs_immune_genes$gene_symbol, ]
dds_variant_time_3 <- estimateSizeFactors(dds_variant_time_3)

# Adding ENTREZ ID for differentially expressed geneNames
diff_genes_3_heatmap <- diff_genes_3_heatmap %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, diff_genes_3_heatmap$geneName, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))


matrix_3 <- counts(dds_variant_time_3, normalized =T) [unique(diff_genes_3_heatmap$geneName),]

# Update column names of the filtered matrix with the 'variant' values
colnames(matrix_3) <- design_matrix_3[colnames(matrix_3), "variant"]

matrix.z_3 <- t(apply(matrix_3,1,scale))
colnames(matrix.z_3) <- colnames(matrix_3)

matrix.z_3 <- matrix.z_3[, order(match(colnames(matrix.z_3), c("mock", "Wuhan", "delta", "omicron")))] 

heatmap_3 <- Heatmap(matrix.z_3,
                     cluster_rows = T,
                     cluster_columns = F,
                     #column_order = order(colnames(matrix.z_3)),
                     column_labels = colnames(matrix.z_3),
                     column_names_gp = gpar(fontsize = 0),
                     column_names_side = "top",
                     column_names_rot = 45,
                     #column_split = 4,
                     #row_split = 3,
                     row_km = 3, # cluster based on kmean
                     column_split = colnames(matrix.z_3),
                     row_labels = rownames(matrix.z_3), 
                     row_names_gp = gpar(fontsize = 0),
                     name="log2FoldChanges")
heatmap_3

# Obtain the ordered row indices (genes) from heatmap_3
clusterlist_3 <- row_order(heatmap_3)

# Create a list to store genes for each cluster
cluster_genes_list_3 <- list()

# Extract gene names based on the ordered row indices per cluster
for (cluster_id in names(clusterlist_3)) {
  # Get the ordered indices of genes for the current cluster
  ordered_indices <- clusterlist_3[[cluster_id]]
  
  # Get the corresponding gene names from matrix.z_3
  cluster_genes <- rownames(matrix.z_3)[ordered_indices]
  
  # Store the genes in the cluster_genes_list_3 with cluster_id as the key
  cluster_genes_list_3[[cluster_id]] <- cluster_genes
}


# Create an empty data frame to store the gene IDs and clusters
clu_df_3 <- data.frame(GeneID = character(), Cluster = character(), stringsAsFactors = FALSE)
# Iterate over each cluster and its corresponding gene list
for (cluster_id in names(cluster_genes_list_3)) {
  # Get the gene IDs for the current cluster
  gene_ids <- cluster_genes_list_3[[cluster_id]]
  
  # Create a data frame for the current cluster
  cluster_df <- data.frame(GeneID = gene_ids, Cluster = cluster_id)
  
  # Append the cluster_df to clu_df_3
  clu_df_3 <- rbind(clu_df_3, cluster_df)
}


clu_df_3 <- clu_df_3 %>% 
  mutate(ENTREZID = mapIds(org.Hs.eg.db, clu_df_3$GeneID, 
                           fromType = "SYMBOL", toType = "ENTREZID", 
                           column = "ENTREZID", keytype = "SYMBOL"))

length(clu_df_3$GeneID[clu_df_3$Cluster == "1"])
length(clu_df_3$GeneID[clu_df_3$Cluster == "2"])
length(clu_df_3$GeneID[clu_df_3$Cluster == "3"])



# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
clu_df_3_genes_c1 <- unique(clu_df_3$GeneID [clu_df_3$Cluster == 1])

enriched_3_c1 <- enrichr(clu_df_3_genes_c1, dbs)
enriched_3_c1_BP <- as.data.frame(enriched_3_c1[["GO_Biological_Process_2021"]])

e_3_c1 <- plotEnrich(
  enriched_3_c1_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for Cluster 1 genes 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_3_c1

# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
clu_df_3_genes_c2 <- unique(clu_df_3$GeneID [clu_df_3$Cluster == 2])

enriched_3_c2 <- enrichr(clu_df_3_genes_c2, dbs)
enriched_3_c2_BP <- as.data.frame(enriched_3_c2[["GO_Biological_Process_2021"]])

e_3_c2 <- plotEnrich(
  enriched_3_c2_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for Cluster 2 genes 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_3_c2

# Extract unique gene names from top_diff_Wuhan_up and top_diff_Wuhan_down
clu_df_3_genes_c3 <- unique(clu_df_3$GeneID [clu_df_3$Cluster == 3])

enriched_3_c3 <- enrichr(clu_df_3_genes_c3, dbs)
enriched_3_c3_BP <- as.data.frame(enriched_3_c3[["GO_Biological_Process_2021"]])

e_3_c3 <- plotEnrich(
  enriched_3_c3_BP,
  showTerms = 10,
  numChar = 60,
  orderBy = "Overlap",
  title = "Enriched pathways for Cluster 3 genes 1 DPI") + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(size = 12)) 
e_3_c3




# Define cluster functions
cluster_functions_3 = list(
  `1` = strwrap(e_3_c1$data$Term, width = 60),
  `2` = strwrap(e_3_c2$data$Term, width = 60),
  `3` = strwrap(e_3_c3$data$Term, width = 60)
)                


# note how we set the width of this empty annotation
ha = rowAnnotation(foo = anno_empty(border = FALSE, 
                                    width = max_text_width(unlist(cluster_functions_3)) + unit(3, "mm")))

Heatmap(matrix.z_3,
        cluster_rows = T,
        cluster_columns = F,
        #column_order = order(colnames(matrix.z_3)),
        column_labels = colnames(matrix.z_3),
        column_names_gp = gpar(fontsize = 0),
        column_names_side = "top",
        column_names_rot = 45,
        #column_split = 4,
        #row_split = 4,
        row_km = 3, # cluster based on kmean
        column_split = colnames(matrix.z_3),
        row_labels = rownames(matrix.z_3), 
        row_names_gp = gpar(fontsize = 0),
        name="norm. z-score", 
        right_annotation = ha,
        use_raster = FALSE)

for(i in 1:3) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(cluster_functions_3[[i]], collapse = "\n"), x = unit(3, "mm"), just = "left")
  })
}





