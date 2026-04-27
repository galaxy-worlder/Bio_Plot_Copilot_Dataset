#############
### Setup ###
#############

# install.packages("librarian")
librarian::shelf(tidyverse, biomaRt, ggtext, patchwork, GGally,
                 ggrepel, apeglm, ComplexHeatmap, EnhancedVolcano,
                 RColorBrewer, DESeq2, limma, gprofiler2, rrvgo,
                 AnnotationDbi, ggh4x, org.Hs.eg.db, circlize, scales,
                 ggseqlogo, universalmotif, patchwork, rstatix, ggpubr)

in_dir <- file.path("/path/to/nf-res/rnaseq")

work_dir <- "/path/to/workdir/rnaseq"

# set fig dir
fig_dir <- file.path(work_dir, "figs")
dir.create(fig_dir) # only run once

# out dir
out_dir <- file.path(work_dir, "output")
dir.create(out_dir) # only run once

########################
### Custom functions ###
########################

httr::set_config(httr::config(ssl_verifypeer = 0L)) # run if mart is unresponsive
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

ens_to_gen <- function(table) {
  gene_names_gen <- rownames(table)
  gene_names_ens <- sapply(strsplit(gene_names_gen, ".", fixed=T), function(x) x[1])
  table <- cbind(gene_names_ens, table)
  gene_IDs <- mapIds(org.Hs.eg.db, keys=table$gene_names_ens,
                     column="SYMBOL",keytype="ENSEMBL",multiVals="list")
  gene_IDs <- data.frame("gene_names_ens" = names(unlist(gene_IDs)),
                         "hgnc_symbol" = unlist(gene_IDs))
  table <- left_join(as.data.frame(table), gene_IDs, by="gene_names_ens")
  return(table)
}

add_diff_col <- function(table) {
  # add a column of NAs
  table$diffexpressed <- "NO"
  # if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
  table$diffexpressed[table$log2FoldChange > 1 & table$padj < 0.05] <- "UP"
  # if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
  table$diffexpressed[table$log2FoldChange < -1 & table$padj < 0.05] <- "DOWN"
  return(table)
}

add_label_col <- function(table) {
  table$delabel <- NA
  table$delabel[table$diffexpressed != "NO"] <- table$hgnc_symbol[table$diffexpressed != "NO"]
  return(table)
}

adjust_scale <- function(table) {
  table$pval_neglog10 <- ifelse(-log10(table$padj) >= 300, 300, -log10(table$padj))
  return(table)
}

add_top_label <- function(table) {
  ten_up <- table %>% arrange(padj, log2FoldChange) %>% dplyr::slice(1:10) %>% dplyr::select(hgnc_symbol)
  ten_down <- table %>% arrange(padj, -log2FoldChange) %>% dplyr::slice(1:10) %>% dplyr::select(hgnc_symbol)
  table <- table %>%
    mutate(top_genes = case_when(hgnc_symbol %in% ten_up$hgnc_symbol | hgnc_symbol %in% ten_down$hgnc_symbol ~ hgnc_symbol,
                                 TRUE ~ ""))
  return(table)
}

#################
### Load data ###
#################

# generate a list of all files within the QoRTS output directory
sample_dir <- file.path(in_dir, "qorts")
samples <- list.files(path = sample_dir, full.names = T)
files <- file.path(samples, "QC.geneCounts.formatted.for.DESeq.txt.gz")
files <- files %>%
  str_replace(., paste0(sample_dir, "/"), "")
names(files) <- dirname(files)

sampletype <- factor(c(rep("undiff_basal",4),rep("diff_basal",4),rep("diff_aicar",3), rep("diff_palm", 4)))
metadata <- data.frame(sampleName = names(files), fileName = files,
                       group = sampletype)

metadata <- metadata %>% relocate(fileName, .after=sampleName)
metadata$diff_state <- factor(c(rep("undiff", 4),rep("diff", 11)))
metadata$condition <- factor(c(rep("basal", 8),rep("aicar", 3),rep("palm", 4)))

# all samples
dds <- DESeqDataSetFromHTSeqCount(sampleTable=metadata,
                                  directory=paste0(work_dir, "rnaseq/qorts"),
                                  design=~group)
dds$group<-relevel(dds$group,ref="undiff_basal")

# subset samples
metadata_basal <- metadata %>% filter(condition == "basal")
metadata_diff <- metadata %>% filter(diff_state == "diff")

# basal undiff vs. diff
dds_basal <- DESeqDataSetFromHTSeqCount(sampleTable=metadata_basal, directory=paste0(work_dir, "rnaseq/qorts"), design=~diff_state)
dds_basal$diff_state<-relevel(dds_basal$diff_state,ref="undiff")

# diff conditions
dds_diff <- DESeqDataSetFromHTSeqCount(sampleTable=metadata_diff, directory=paste0(work_dir, "rnaseq/qorts"), design=~condition)
dds_diff$condition<-relevel(dds_diff$condition,ref="basal")

###########
### PCA ###
###########

# pre-filter counts
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize

dds <- dds[keep,]

# perform count normalization
dds <- estimateSizeFactors(dds)

# visualize via PCA
rld <- rlog(dds, blind=TRUE)
pca_vst <- vst(dds, blind=TRUE)
pca_df <- assay(pca_vst)

vst_df <- as.data.frame(pca_df)

vst_df <- ens_to_gen(vst_df)

colnames(vst_df)[2:16] <- c("ATS0263", "ATS0264", "ATS0265", "ATS0266",
                            "ATS0267", "ATS0268", "ATS0269", "ATS0270",
                            "ATS0271", "ATS0273", "ATS0274",
                            "ATS0275", "ATS0276", "ATS0277", "ATS0278")

# vst_df <- vst_df %>%
#   rownames_to_column(var = "ENSG_ID")
# write.table(vst_df, file.path(out_dir, "tovar-nishino_rnaseq_vst_counts.txt"),
#             sep = "\t", row.names = F, quote = F)

raw_df <- as.data.frame(assay(dds)) %>%
  rownames_to_column(var = "ENSG_ID")

colnames(raw_df) <- colnames(vst_df)

write.table(raw_df, file.path(out_dir, "tovar-nishino_rnaseq_raw_counts.txt"),
            sep = "\t", row.names = F, quote = F)

pca_res <- prcomp(t(pca_df))  # Transpose to have samples in rows
# Create a data frame with PCA results
pca_res_df <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  condition = factor(c(rep("undiff",4),rep("diff",4),rep("aicar",3),rep("palm",4)))
)


# Calculate the proportion of variance explained
pct_var_pca <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), digits = 2)


pca_res_df$label <- c("Undifferentiated +\nBasal", rep(NA, 6),
                      "Differentiated +\nBasal", NA,
                      "Differentiated +\nAICAR", NA,
                      "Differentiated +\nPalmitate",rep(NA, 3))
pca_res_df$color <- c(rep("#ff595e",4),rep("#ffca3a",4),rep("#8ac926",3), rep("#1982c4", 4))
pca_res_df$status <- c(rep("undiff",4),rep("diff",11))
rna_pca <- ggplot(pca_res_df, aes(PC1, PC2, color = condition, label=label, shape=status)) +
  geom_point(size=3, fill = pca_res_df$color, color = "black") +
  scale_shape_manual(values = c(21, 24)) +
  xlab(paste0("PC1: ",pct_var_pca[1],"% variance")) +
  ylab(paste0("PC2: ",pct_var_pca[2],"% variance")) + 
  scale_color_manual(values = c("#ff595e","#ffca3a","#8ac926","#1982c4")) +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, size = 5, lineheight = 1,
                  data = pca_res_df[c(1:4),], color = "black") + 
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_x = -0.5, size = 5, lineheight = 1,
                  data = pca_res_df[c(5:8),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_x = -0.05, size = 5, lineheight = 1,
                  data = pca_res_df[c(9:11),], color = "black") +
  geom_text_repel(min.segment.length = 0, family = "Helvetica",
                  box.padding = 1.5, nudge_x = 1, size = 5, lineheight = 1,
                  data = pca_res_df[c(12:15),], color = "black") +
  theme_bw(base_family = "Helvetica", base_size = 16) + theme(legend.position="none")

ggsave(rna_pca, filename = file.path(fig_dir, "rnaseq-pca.png"),
       units = "in", height = 4, width = 4.67, dpi = 600, device = ragg::agg_png())

#####################################
### DE analysis - undiff vs. diff ###
#####################################

smallestGroupSize <- 4
keep <- rowSums(counts(dds_basal) >= 10) >= smallestGroupSize

dds_basal <- dds_basal[keep,]

# perform count normalization
dds_basal <- estimateSizeFactors(dds_basal)

# estimate count dispersion
dds_basal <- estimateDispersions(dds_basal)

# norm count data for visualizations
vsd_basal <- vst(dds_basal, blind=FALSE)
mat_basal <- assay(vsd_basal)

# run deseq2
dds_basal <- DESeq(dds_basal)
resultsNames(dds_basal) # [1] "Intercept"                 "diff_state_diff_vs_undiff"

# set up contrasts to return DE results tables
contrast_diff_state <- c("diff_state", "diff", "undiff") # corrected previous error here that flipped the contrast
res_table_diff_state <- results(dds_basal,
                                contrast=contrast_diff_state)
res_table_diff_state <- ens_to_gen(res_table_diff_state)

# save DE results
write.table(res_table_diff_state, file.path(out_dir, "de_res_diff_state.tsv"),
            sep = "\t", quote = F, row.names = F)

#####################################
### Volcano plot - undiff vs diff ###
#####################################

# shrink estimates for viz
res_table_diff_lfc <- lfcShrink(dds_basal,
                                coef = "diff_state_diff_vs_undiff",
                                type="apeglm")

# convert to gene symbols
# gencode -> ensembl
res_table_diff_lfc <- ens_to_gen(res_table_diff_lfc)
res_table_diff_lfc <- add_label_col(add_diff_col(res_table_diff_lfc))
res_table_diff_lfc <- adjust_scale(res_table_diff_lfc)

diff_top_genes <- res_table_diff_lfc %>%
  dplyr::filter(padj < 0.05 & !is.na(hgnc_symbol) & log2FoldChange > 1) %>%
  arrange(padj) %>%
  dplyr::slice(1:10) %>%
  bind_rows(
    res_table_diff_lfc %>%
      dplyr::filter(padj < 0.05 & !is.na(hgnc_symbol) & log2FoldChange < -1) %>%
      arrange(padj) %>%
      dplyr::slice(1:10)
  )

vol_diff_state <- res_table_diff_lfc %>%
  ggplot(aes(x=log2FoldChange, y=pval_neglog10, color=diffexpressed)) +
  geom_text_repel(data = diff_top_genes, aes(x = log2FoldChange,
                                             y = pval_neglog10,
                                             label = hgnc_symbol),
                  max.overlaps = Inf, box.padding = 0.5,
                  color = "black", family = "Helvetica",
                  segment.size = 0.2,
                  nudge_y = -50) +
  geom_point(alpha = 0.25, size = 2) + 
  labs(x = expression(paste("log"[2],"(fold-change)")),
       y = expression(paste("-log"[10],"(p-value)"))) +
  theme_bw(base_family = "Helvetica",
           base_size = 14) +
  scale_color_manual(values=c("#ff595e", "gray", "#ffca3a")) +
  geom_vline(xintercept=c(-1, 1), col="black",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype="dashed") +
  guides(color = "none")

ggsave(vol_diff_state, filename = file.path(fig_dir, "vol-diff.png"), units = "in",
       height = 5, width = 6, dpi = 600, device = ragg::agg_png())


################################
### Heatmap - undiff vs diff ###
################################

uvd_set <- res_table_diff_lfc$gene_names_ens[abs(res_table_diff_lfc$log2FoldChange) >= 1 & res_table_diff_lfc$padj <= 0.05]
uvd_df <- vst_df %>%
  dplyr::select(-c(ATS0271:ATS0278)) %>%
  dplyr::filter(gene_names_ens %in% uvd_set)

uvd_mat <- as.matrix(uvd_df[,c(2:9)])
row.names(uvd_mat) <- uvd_df$gene_names_ens

uvd_mat_norm <- t(scale(t(uvd_mat)))
uvd_groups <- c(rep("Group 1", 4), rep("Group 2", 4))
cond_annot <- HeatmapAnnotation(
  condition = anno_block(gp = gpar(fill = c("#ff595e","#ffca3a"),
                                   border = NA,
                                   lty = "blank"),
                         labels = c("Undifferentiated + Basal",
                                    "Differentiated + Basal"),
                         labels_gp = gpar(fontsize = 12, fontfamily = "Helvetica", col = "white"))
)

my_palette <- rev(brewer.pal(11, "RdBu"))

uvd_breaks <- seq(-2.23, 2.12, length.out = length(my_palette))

# Plot heatmap
uvd_heatmap <- Heatmap(uvd_mat_norm,
                       top_annotation = cond_annot,
                       show_row_names = FALSE,
                       show_column_names = FALSE,
                       cluster_rows = TRUE,
                       cluster_columns = FALSE,
                       show_row_dend = FALSE,
                       column_order = rev(1:8),
                       column_split = uvd_groups,
                       column_title = NULL,
                       col = colorRamp2(diff_breaks, my_palette),
                       heatmap_legend_param = list(
                         title = "Norm. expr."
                       )
)

png(file.path(fig_dir, "uvd_heatmap.png"), units = "in",
    width = 5.6, height = 4, res = 600)
uvd_heatmap
dev.off()

#####################################
### DE analysis - basal vs others ###
#####################################

smallestGroupSize <- 3
keep <- rowSums(counts(dds_diff) >= 10) >= smallestGroupSize

dds_diff <- dds_diff[keep,]

# perform count normalization
dds_diff <- estimateSizeFactors(dds_diff)

# estimate count dispersion
dds_diff <- estimateDispersions(dds_diff)

# norm count data for visualizations
vsd_diff <- vst(dds_diff, blind=FALSE)
mat_diff <- assay(dds_diff)

# run deseq2
dds_diff <- DESeq(dds_diff)
resultsNames(dds_diff) # [1] "Intercept"                "condition_aicar_vs_basal" "condition_palm_vs_basal"

# set up contrasts to return DE results tables
contrast_aicar_state <- c("condition", "aicar", "basal")  # corrected previous error here that flipped the contrast
contrast_palm_state <- c("condition", "palm", "basal")  # corrected previous error here that flipped the contrast

res_table_aicar_state <- results(dds_diff,
                                 contrast=contrast_aicar_state)
res_table_aicar_state <- ens_to_gen(res_table_aicar_state)

res_table_palm_state <- results(dds_diff,
                                contrast=contrast_palm_state)
res_table_palm_state <- ens_to_gen(res_table_palm_state)

# save DE results
write.table(res_table_aicar_state, file.path(out_dir, "de_res_aicar_state.tsv"),
            sep = "\t", quote = F, row.names = F)
write.table(res_table_palm_state, file.path(out_dir, "de_res_palm_state.tsv"),
            sep = "\t", quote = F, row.names = F)

####################################
### Volcano plot - diff vs aicar ###
#################################### 

# lfc shrink for viz
res_table_aicar_lfc <- lfcShrink(dds_diff, coef = "condition_aicar_vs_basal", type="apeglm")

# convert to gene symbols
## aicar
# gencode -> ensembl

res_table_aicar_lfc <- ens_to_gen(res_table_aicar_lfc)
res_table_aicar_lfc <- add_label_col(add_diff_col(res_table_aicar_lfc))
res_table_aicar_lfc <- adjust_scale(res_table_aicar_lfc)

aicar_top_genes <- res_table_aicar_lfc %>%
  dplyr::filter(padj < 0.05 & !is.na(hgnc_symbol) & log2FoldChange > 1) %>%
  arrange(padj) %>%
  dplyr::slice(1:10) %>%
  bind_rows(
    res_table_aicar_lfc %>%
      dplyr::filter(padj < 0.05 & !is.na(hgnc_symbol) & log2FoldChange < -1) %>%
      arrange(padj) %>%
      dplyr::slice(1:10)
  )

# volcano plot
vol_aicar_state <- res_table_aicar_lfc %>%
  ggplot(aes(x=log2FoldChange, y=pval_neglog10, color=diffexpressed)) +
  geom_text_repel(data = subset(aicar_top_genes, log2FoldChange > 1), aes(x = log2FoldChange,
                                                                          y = pval_neglog10,
                                                                          label = hgnc_symbol),
                  max.overlaps = Inf, box.padding = 0.5,
                  color = "black", family = "Helvetica",
                  segment.size = 0.2, nudge_y = 20, nudge_x = 1.5) +
  geom_text_repel(data = subset(aicar_top_genes, log2FoldChange < -1), aes(x = log2FoldChange,
                                                                           y = pval_neglog10,
                                                                           label = hgnc_symbol),
                  max.overlaps = Inf, box.padding = 0.5,
                  color = "black", family = "Helvetica",
                  segment.size = 0.2, nudge_y = 35, nudge_x = -1) +
  geom_point(alpha = 0.25, size = 2) + 
  labs(x = expression(paste("log"[2],"(fold-change)")),
       y = expression(paste("-log"[10],"(p-value)"))) +
  theme_bw(base_family = "Helvetica",
           base_size = 14) +
  scale_color_manual(values=c("#ffca3a", "gray", "#8ac926")) +
  geom_vline(xintercept=c(-1, 1), col="black",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype="dashed") +
  guides(color = "none")

ggsave(vol_aicar_state, filename = file.path(fig_dir, "vol-aicar.png"), units = "in",
       height = 5, width = 6, dpi = 600, device = ragg::agg_png())

###################################
### Volcano plot - diff vs palm ###
###################################

# lfc shrink for viz
res_table_palm_lfc <- lfcShrink(dds_diff, coef = "condition_palm_vs_basal", type = "apeglm")

# convert to gene symbols
## palm
# gencode -> ensembl

res_table_palm_lfc <- ens_to_gen(res_table_palm_lfc)
res_table_palm_lfc <- add_label_col(add_diff_col(res_table_palm_lfc))
res_table_palm_lfc <- adjust_scale(res_table_palm_lfc)

palm_top_genes <- res_table_palm_lfc %>%
  dplyr::filter(padj < 0.05 & !is.na(hgnc_symbol) & log2FoldChange > 1) %>%
  arrange(padj) %>%
  dplyr::slice(1:10) %>%
  bind_rows(
    res_table_palm_lfc %>%
      dplyr::filter(padj < 0.05 & !is.na(hgnc_symbol) & log2FoldChange < -1) %>%
      arrange(padj) %>%
      dplyr::slice(1:10)
  )

# volcano plot
vol_palm_state <- res_table_palm_lfc %>%
  ggplot(aes(x=log2FoldChange, y=pval_neglog10, color=diffexpressed)) +
  geom_text_repel(data = subset(palm_top_genes, log2FoldChange > 1), aes(x = log2FoldChange,
                                                                         y = pval_neglog10,
                                                                         label = hgnc_symbol),
                  max.overlaps = Inf, box.padding = 0.5,
                  color = "black", family = "Helvetica",
                  segment.size = 0.2, nudge_y = 10, nudge_x = 1.5) +
  geom_text_repel(data = subset(palm_top_genes, log2FoldChange < -1), aes(x = log2FoldChange,
                                                                          y = pval_neglog10,
                                                                          label = hgnc_symbol),
                  max.overlaps = Inf, box.padding = 0.5,
                  color = "black", family = "Helvetica",
                  segment.size = 0.2, nudge_y = 35, nudge_x = -1) +
  geom_point(alpha = 0.25, size = 2) + 
  labs(x = expression(paste("log"[2],"(fold-change)")),
       y = expression(paste("-log"[10],"(p-value)"))) +
  theme_bw(base_family = "Helvetica",
           base_size = 14) +
  scale_color_manual(values=c("#ffca3a", "gray", "#1982c4")) +
  geom_vline(xintercept=c(-1, 1), col="black",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype="dashed") +
  guides(color = "none")

ggsave(vol_palm_state, filename = file.path(fig_dir, "vol-palm.png"), units = "in",
       height = 5, width = 6, dpi = 600, device = ragg::agg_png())


#########################################
### Pathway analysis - undiff vs diff ###
#########################################

filtered_basal <- res_table_diff_lfc %>%
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::distinct(gene_names_ens, .keep_all = TRUE)

lfc_vector_basal <- data.frame(log2FoldChange = filtered_basal$log2FoldChange,
                               gene_names_ens = filtered_basal$gene_names_ens,
                               symbol = filtered_basal$hgnc_symbol)

lfc_vector_basal <- lfc_vector_basal %>%
  arrange(-log2FoldChange)

# split into up and down genes 
basal_up <- lfc_vector_basal %>% dplyr::filter(log2FoldChange > 0)
basal_down <- lfc_vector_basal %>% dplyr::filter(log2FoldChange < 0) %>% arrange(log2FoldChange)

# generate separate enrichment results
gprof_basal_up <- gost(basal_up$symbol, organism = "hsapiens", ordered_query = TRUE, exclude_iea = TRUE)
gprof_basal_down <- gost(basal_down$symbol, organism = "hsapiens", ordered_query = TRUE, exclude_iea = TRUE)
pathwayres_basal_up <- gprof_basal_up$result
pathwayres_basal_down <- gprof_basal_down$result

# use rrvgo to summarize enrichment results - just go terms
go_sets <- c("GO:BP", "GO:MF")
pathwayres_basal_up_filt <- pathwayres_basal_up %>% dplyr::filter(source %in% go_sets)
pathwayres_basal_down_filt <- pathwayres_basal_down %>% dplyr::filter(source %in% go_sets)

# make MF sim matrix - basal
sim_mat_basal_mf_up <- calculateSimMatrix(pathwayres_basal_up_filt$term_id,
                                          orgdb="org.Hs.eg.db",
                                          method="Rel", ont = "MF")
basal_mf_up_scores <- setNames(-log10(pathwayres_basal_up_filt$p_value), pathwayres_basal_up_filt$term_id)
red_mat_basal_mf_up <- reduceSimMatrix(sim_mat_basal_mf_up,
                                       basal_mf_up_scores,
                                       threshold=0.7,
                                       orgdb="org.Hs.eg.db")

# make BP sim matrix - basal
sim_mat_basal_bp_up <- calculateSimMatrix(pathwayres_basal_up_filt$term_id,
                                          orgdb="org.Hs.eg.db",
                                          method="Rel", ont = "BP")
basal_bp_up_scores <- setNames(-log10(pathwayres_basal_up_filt$p_value), pathwayres_basal_up_filt$term_id)
red_mat_basal_bp_up <- reduceSimMatrix(sim_mat_basal_bp_up,
                                       basal_bp_up_scores,
                                       threshold=0.7,
                                       orgdb="org.Hs.eg.db")

sim_mat_basal_mf_down <- calculateSimMatrix(pathwayres_basal_down_filt$term_id,
                                            orgdb="org.Hs.eg.db",
                                            method="Rel", ont = "MF")
basal_mf_down_scores <- setNames(-log10(pathwayres_basal_down_filt$p_value), pathwayres_basal_down_filt$term_id)
red_mat_basal_mf_down <- reduceSimMatrix(sim_mat_basal_mf_down,
                                         basal_mf_down_scores,
                                         threshold=0.7,
                                         orgdb="org.Hs.eg.db")

sim_mat_basal_bp_down <- calculateSimMatrix(pathwayres_basal_down_filt$term_id,
                                            orgdb="org.Hs.eg.db",
                                            method="Rel", ont = "BP")
basal_bp_down_scores <- setNames(-log10(pathwayres_basal_down_filt$p_value), pathwayres_basal_down_filt$term_id)
red_mat_basal_bp_down <- reduceSimMatrix(sim_mat_basal_bp_down,
                                         basal_bp_down_scores,
                                         threshold=0.7,
                                         orgdb="org.Hs.eg.db")

red_mat_basal_down_filt <- red_mat_basal_bp_down %>%
  dplyr::filter(size > 10 & size < 1000) %>%
  distinct(parentTerm) %>%
  dplyr::slice(1:10) %>%
  bind_rows(
    red_mat_basal_mf_down %>%
      dplyr::filter(size > 10 & size < 1000) %>%
      distinct(parentTerm) %>%
      dplyr::slice(1:10)
  )
red_mat_basal_up_filt <- red_mat_basal_bp_up %>%
  dplyr::filter(size > 10 & size < 1000) %>%
  distinct(parentTerm) %>%
  dplyr::slice(1:10) %>%
  bind_rows(
    red_mat_basal_mf_up %>%
      dplyr::filter(size > 10 & size < 1000) %>%
      distinct(parentTerm) %>%
      dplyr::slice(1:10)
  )

red_basal_path <- red_mat_basal_bp_down %>%
  distinct(parentTerm) %>%
  bind_rows(
    red_mat_basal_mf_down %>%
      distinct(parentTerm) %>%
      bind_rows(
        red_mat_basal_bp_up %>%
          distinct(parentTerm) %>%
          bind_rows(
            red_mat_basal_mf_up %>%
              distinct(parentTerm)
          )
      )
  )

basal_path_out <- pathwayres_basal_up %>%
  dplyr::filter(term_name %in% red_basal_path$parentTerm) %>%
  mutate(direction = "UP") %>%
  bind_rows(
    pathwayres_basal_down %>%
      dplyr::filter(term_name %in% red_basal_path$parentTerm) %>%
      mutate(direction = "DOWN")
  )

write.table(basal_path_out[,-c(1,14)], file.path(out_dir, "path_res_diff_state.tsv"),
            sep = "\t", quote = F, row.names = F)

# plot top ten up and down pathways using reduced terms (aka up = diff, down = undiff)
basal_pathways <- pathwayres_basal_up_filt %>%
  dplyr::filter(term_name %in% red_mat_basal_up_filt$parentTerm & term_size > 10 & term_size < 1000) %>%
  mutate(group = "diff") %>%
  group_by(source) %>%
  arrange(p_value) %>%
  dplyr::slice(1:5) %>%
  bind_rows(
    pathwayres_basal_down_filt %>%
      dplyr::filter(term_name %in% red_mat_basal_down_filt$parentTerm & term_size > 10 & term_size < 1000) %>%
      mutate(group = "undiff") %>%
      group_by(source) %>%
      arrange(p_value) %>%
      dplyr::slice(1:5)
  )

basal_labs <- c("Undifferentiated + Basal", "Differentiated + Basal")
names(basal_labs) <- c("undiff","diff")
strip <- strip_themed(background_x = elem_list_rect(fill = c("#ff595e","#ffca3a")))
basal_path_plot <- basal_pathways %>%
  group_by(group) %>%
  mutate(term_name = fct_reorder(term_name, -log10(p_value)),
         group = factor(group, levels = c("undiff", "diff"))) %>%
  ggplot(aes(x = -log10(p_value), y = term_name)) +
  geom_bar(stat = "identity") +
  facet_wrap2(vars(group), scales = "free_y", ncol = 1,
              strip = strip, labeller = labeller(group = basal_labs)) +
  theme_linedraw(base_family = "Helvetica", base_size = 12) + labs(y = NULL, x = expression("-log"[10]~"(adjusted p-value)")) +
  theme(strip.background = element_rect(color = NA),
        strip.text = element_text(color = "white"),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.1),
        panel.grid.minor.x = element_line(linewidth = 0.1),
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(lineheight = 0.6),
        aspect.ratio = 0.8) +
  scale_y_discrete(labels = label_wrap(35)) +
  scale_x_continuous(breaks = seq(0,50,by=10), labels = seq(0,50,by=10)) +
  geom_vline(xintercept = -log10(0.05), color = "black", linewidth = 0.75)

ggsave(basal_path_plot, filename = file.path(fig_dir, "basal_pathways_plot.png"), units = "in",
       width = 4.25, height = 5, dpi = 600)

#########################################
### Pathway analysis - diff vs. aicar ###
#########################################

filtered_aicar <- res_table_aicar_lfc %>%
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::distinct(gene_names_ens, .keep_all = TRUE)

lfc_vector_aicar <- data.frame(log2FoldChange = filtered_aicar$log2FoldChange,
                               gene_names_ens = filtered_aicar$gene_names_ens,
                               symbol = filtered_aicar$hgnc_symbol)

lfc_vector_aicar <- lfc_vector_aicar %>%
  arrange(-log2FoldChange)

# split into up and down genes 
aicar_up <- lfc_vector_aicar %>% dplyr::filter(log2FoldChange > 0)
aicar_down <- lfc_vector_aicar %>% dplyr::filter(log2FoldChange < 0) %>% arrange(log2FoldChange)

# generate separate enrichment results
gprof_aicar_up <- gost(aicar_up$symbol, organism = "hsapiens", ordered_query = TRUE, exclude_iea = TRUE)
gprof_aicar_down <- gost(aicar_down$symbol, organism = "hsapiens", ordered_query = TRUE, exclude_iea = TRUE)
pathwayres_aicar_up <- gprof_aicar_up$result
pathwayres_aicar_down <- gprof_aicar_down$result

# use rrvgo to summarize enrichment results - just go terms
go_sets <- c("GO:BP", "GO:CC", "GO:MF")
pathwayres_aicar_up_filt <- pathwayres_aicar_up %>% dplyr::filter(source %in% go_sets)
pathwayres_aicar_down_filt <- pathwayres_aicar_down %>% dplyr::filter(source %in% go_sets)

# make MF sim matrix - aicar
sim_mat_aicar_mf_up <- calculateSimMatrix(pathwayres_aicar_up_filt$term_id,
                                          orgdb="org.Hs.eg.db",
                                          method="Rel", ont = "MF")
aicar_mf_up_scores <- setNames(-log10(pathwayres_aicar_up_filt$p_value), pathwayres_aicar_up_filt$term_id)
red_mat_aicar_mf_up <- reduceSimMatrix(sim_mat_aicar_mf_up,
                                       aicar_mf_up_scores,
                                       threshold=0.7,
                                       orgdb="org.Hs.eg.db")

# make BP sim matrix - basal
sim_mat_aicar_bp_up <- calculateSimMatrix(pathwayres_aicar_up_filt$term_id,
                                          orgdb="org.Hs.eg.db",
                                          method="Rel", ont = "BP")
aicar_bp_up_scores <- setNames(-log10(pathwayres_aicar_up_filt$p_value), pathwayres_aicar_up_filt$term_id)
red_mat_aicar_bp_up <- reduceSimMatrix(sim_mat_aicar_bp_up,
                                       aicar_bp_up_scores,
                                       threshold=0.7,
                                       orgdb="org.Hs.eg.db")

sim_mat_aicar_mf_down <- calculateSimMatrix(pathwayres_aicar_down_filt$term_id,
                                            orgdb="org.Hs.eg.db",
                                            method="Rel", ont = "MF")
aicar_mf_down_scores <- setNames(-log10(pathwayres_aicar_down_filt$p_value), pathwayres_aicar_down_filt$term_id)
red_mat_aicar_mf_down <- reduceSimMatrix(sim_mat_aicar_mf_down,
                                         aicar_mf_down_scores,
                                         threshold=0.7,
                                         orgdb="org.Hs.eg.db")

sim_mat_aicar_bp_down <- calculateSimMatrix(pathwayres_aicar_down_filt$term_id,
                                            orgdb="org.Hs.eg.db",
                                            method="Rel", ont = "BP")
aicar_bp_down_scores <- setNames(-log10(pathwayres_aicar_down_filt$p_value), pathwayres_aicar_down_filt$term_id)
red_mat_aicar_bp_down <- reduceSimMatrix(sim_mat_aicar_bp_down,
                                         aicar_bp_down_scores,
                                         threshold=0.7,
                                         orgdb="org.Hs.eg.db")

red_mat_aicar_down_filt <- red_mat_aicar_bp_down %>%
  dplyr::filter(size > 10 & size < 1000) %>%
  distinct(parentTerm) %>%
  dplyr::slice(1:15) %>%
  bind_rows(
    red_mat_aicar_mf_down %>%
      dplyr::filter(size > 10 & size < 1000) %>%
      distinct(parentTerm) %>%
      dplyr::slice(1:10)
  )
red_mat_aicar_up_filt <- red_mat_aicar_bp_up %>%
  dplyr::filter(size > 10 & size < 1000) %>%
  distinct(parentTerm) %>%
  dplyr::slice(1:10) %>%
  bind_rows(
    red_mat_aicar_mf_up %>%
      dplyr::filter(size > 10 & size < 1000) %>%
      distinct(parentTerm) %>%
      dplyr::slice(1:10)
  )

red_aicar_path <- red_mat_aicar_bp_down %>%
  distinct(parentTerm) %>%
  bind_rows(
    red_mat_aicar_mf_down %>%
      distinct(parentTerm) %>%
      bind_rows(
        red_mat_aicar_bp_up %>%
          distinct(parentTerm) %>%
          bind_rows(
            red_mat_aicar_mf_up %>%
              distinct(parentTerm)
          )
      )
  )

aicar_path_out <- pathwayres_aicar_up %>%
  dplyr::filter(term_name %in% red_aicar_path$parentTerm) %>%
  mutate(direction = "UP") %>%
  bind_rows(
    pathwayres_aicar_down %>%
      dplyr::filter(term_name %in% red_aicar_path$parentTerm) %>%
      mutate(direction = "DOWN")
  )

write.table(aicar_path_out[,-c(1,14)], file.path(out_dir, "path_res_aicar_state.tsv"),
            sep = "\t", quote = F, row.names = F)


# plot top ten up and down pathways using reduced terms (aka up = diff, down = undiff)
aicar_pathways <- pathwayres_aicar_up_filt %>%
  dplyr::filter(term_name %in% red_mat_aicar_up_filt$parentTerm & term_size > 10 & term_size < 1000) %>%
  mutate(group = "aicar") %>%
  group_by(source) %>%
  arrange(p_value) %>%
  dplyr::slice(1:5) %>%
  bind_rows(
    pathwayres_aicar_down_filt %>%
      dplyr::filter(term_name %in% red_mat_aicar_down_filt$parentTerm & term_size > 10 & term_size < 1000) %>%
      mutate(group = "diff") %>%
      group_by(source) %>%
      arrange(p_value) %>%
      dplyr::slice(1:5)
  )

aicar_labs <- c("Differentiated + Basal", "Differentiated + AICAR")
names(aicar_labs) <- c("diff","aicar")
strip <- strip_themed(background_x = elem_list_rect(fill = c("#ffca3a","#8ac926")))
aicar_path_plot <- aicar_pathways %>%
  group_by(group) %>%
  mutate(term_name = fct_reorder(term_name, -log10(p_value)),
         group = factor(group, levels = c("diff", "aicar"))) %>%
  ggplot(aes(x = -log10(p_value), y = term_name)) +
  geom_bar(stat = "identity") +
  facet_wrap2(vars(group), scales = "free_y", ncol = 1,
              strip = strip, labeller = labeller(group = aicar_labs)) +
  theme_linedraw(base_family = "Helvetica", base_size = 12) + labs(y = NULL, x = expression("-log"[10]~"(adjusted p-value)")) +
  theme(strip.background = element_rect(color = NA),
        strip.text = element_text(color = "white"),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.1),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(lineheight = 0.6),
        aspect.ratio = 0.8) +
  scale_y_discrete(labels = label_wrap(35)) +
  scale_x_continuous(breaks = seq(0,50,by=5), labels = seq(0,50,by=5)) +
  geom_vline(xintercept = -log10(0.05), color = "black", linewidth = 0.75)

ggsave(aicar_path_plot, filename = file.path(fig_dir, "aicar_pathways_plot.png"), units = "in",
       width = 4.25, height = 5, dpi = 600)

#########################################
### Pathway analysis - undiff vs diff ###
#########################################

filtered_palm <- res_table_palm_lfc %>%
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::distinct(gene_names_ens, .keep_all = TRUE)

lfc_vector_palm <- data.frame(log2FoldChange = filtered_palm$log2FoldChange,
                              gene_names_ens = filtered_palm$gene_names_ens,
                              symbol = filtered_palm$hgnc_symbol)

lfc_vector_palm <- lfc_vector_palm %>%
  arrange(-log2FoldChange)

# split into up and down genes 
palm_up <- lfc_vector_palm %>% dplyr::filter(log2FoldChange > 0)
palm_down <- lfc_vector_palm %>% dplyr::filter(log2FoldChange < 0) %>% arrange(log2FoldChange)

# generate separate enrichment results
gprof_palm_up <- gost(palm_up$symbol, organism = "hsapiens", ordered_query = TRUE, exclude_iea = TRUE)
gprof_palm_down <- gost(palm_down$symbol, organism = "hsapiens", ordered_query = TRUE, exclude_iea = TRUE)
pathwayres_palm_up <- gprof_palm_up$result
pathwayres_palm_down <- gprof_palm_down$result

# use rrvgo to summarize enrichment results - just go terms
go_sets <- c("GO:BP", "GO:CC", "GO:MF")
pathwayres_palm_up_filt <- pathwayres_palm_up %>% dplyr::filter(source %in% go_sets)
pathwayres_palm_down_filt <- pathwayres_palm_down %>% dplyr::filter(source %in% go_sets)

# make MF sim matrix - aicar
sim_mat_palm_mf_up <- calculateSimMatrix(pathwayres_palm_up_filt$term_id,
                                         orgdb="org.Hs.eg.db",
                                         method="Rel", ont = "MF")
palm_mf_up_scores <- setNames(-log10(pathwayres_palm_up_filt$p_value), pathwayres_palm_up_filt$term_id)
red_mat_palm_mf_up <- reduceSimMatrix(sim_mat_palm_mf_up,
                                      palm_mf_up_scores,
                                      threshold=0.7,
                                      orgdb="org.Hs.eg.db")

# make BP sim matrix - basal
sim_mat_palm_bp_up <- calculateSimMatrix(pathwayres_palm_up_filt$term_id,
                                         orgdb="org.Hs.eg.db",
                                         method="Rel", ont = "BP")
palm_bp_up_scores <- setNames(-log10(pathwayres_palm_up_filt$p_value), pathwayres_palm_up_filt$term_id)
red_mat_palm_bp_up <- reduceSimMatrix(sim_mat_palm_bp_up,
                                      palm_bp_up_scores,
                                      threshold=0.7,
                                      orgdb="org.Hs.eg.db")

sim_mat_palm_mf_down <- calculateSimMatrix(pathwayres_palm_down_filt$term_id,
                                           orgdb="org.Hs.eg.db",
                                           method="Rel", ont = "MF")
palm_mf_down_scores <- setNames(-log10(pathwayres_palm_down_filt$p_value), pathwayres_palm_down_filt$term_id)
red_mat_palm_mf_down <- reduceSimMatrix(sim_mat_palm_mf_down,
                                        palm_mf_down_scores,
                                        threshold=0.7,
                                        orgdb="org.Hs.eg.db")

sim_mat_palm_bp_down <- calculateSimMatrix(pathwayres_palm_down_filt$term_id,
                                           orgdb="org.Hs.eg.db",
                                           method="Rel", ont = "BP")
palm_bp_down_scores <- setNames(-log10(pathwayres_palm_down_filt$p_value), pathwayres_palm_down_filt$term_id)
red_mat_palm_bp_down <- reduceSimMatrix(sim_mat_palm_bp_down,
                                        palm_bp_down_scores,
                                        threshold=0.7,
                                        orgdb="org.Hs.eg.db")

red_mat_palm_down_filt <- red_mat_palm_bp_down %>%
  dplyr::filter(size > 10 & size < 1000) %>%
  distinct(parentTerm) %>%
  dplyr::slice(1:15) %>%
  bind_rows(
    red_mat_palm_mf_down %>%
      dplyr::filter(size > 10 & size < 1000) %>%
      distinct(parentTerm) %>%
      dplyr::slice(1:15)
  )
red_mat_palm_up_filt <- red_mat_palm_bp_up %>%
  dplyr::filter(size > 0 & size < 1000) %>%
  distinct(parentTerm) %>%
  dplyr::slice(1:15) %>%
  bind_rows(
    red_mat_palm_mf_up %>%
      dplyr::filter(size > 0 & size < 1000) %>%
      distinct(parentTerm) %>%
      dplyr::slice(1:15)
  )

red_palm_path <- red_mat_palm_bp_down %>%
  distinct(parentTerm) %>%
  bind_rows(
    red_mat_palm_mf_down %>%
      distinct(parentTerm) %>%
      bind_rows(
        red_mat_palm_bp_up %>%
          distinct(parentTerm) %>%
          bind_rows(
            red_mat_palm_mf_up %>%
              distinct(parentTerm)
          )
      )
  )

palm_path_out <- pathwayres_palm_up %>%
  dplyr::filter(term_name %in% red_palm_path$parentTerm) %>%
  mutate(direction = "UP") %>%
  bind_rows(
    pathwayres_palm_down %>%
      dplyr::filter(term_name %in% red_palm_path$parentTerm) %>%
      mutate(direction = "DOWN")
  )

write.table(palm_path_out[,-c(1,14)], file.path(out_dir, "path_res_palm_state.tsv"),
            sep = "\t", quote = F, row.names = F)


# plot top ten up and down pathways using reduced terms (aka up = diff, down = undiff)
palm_pathways <- pathwayres_palm_up_filt %>%
  dplyr::filter(term_name %in% red_mat_palm_up_filt$parentTerm & term_size > 3 & term_size < 1000 &
                  !str_detect(term_name, "sinoatrial") & !str_detect(term_name, "cardiac")) %>%
  mutate(group = "palm") %>%
  group_by(source) %>%
  arrange(p_value) %>%
  dplyr::slice(1:5) %>%
  bind_rows(
    pathwayres_palm_down_filt %>%
      dplyr::filter(term_name %in% red_mat_palm_down_filt$parentTerm & term_size > 3 & term_size < 500) %>%
      mutate(group = "diff") %>%
      group_by(source) %>%
      arrange(p_value) %>%
      dplyr::slice(1:5)
  )

palm_labs <- c("Differentiated + Basal", "Differentiated + Palmitate")
names(palm_labs) <- c("diff","palm")
strip <- strip_themed(background_x = elem_list_rect(fill = c("#ffca3a","#1982c4")))
palm_path_plot <- palm_pathways %>%
  group_by(group) %>%
  mutate(term_name = fct_reorder(term_name, -log10(p_value)),
         group = factor(group, levels = c("diff", "palm"))) %>%
  ggplot(aes(x = -log10(p_value), y = term_name)) +
  geom_bar(stat = "identity") +
  facet_wrap2(vars(group), scales = "free_y", ncol = 1,
              strip = strip, labeller = labeller(group = palm_labs)) +
  theme_linedraw(base_family = "Helvetica", base_size = 12) + labs(y = NULL, x = expression("-log"[10]~"(adjusted p-value)")) +
  theme(strip.background = element_rect(color = NA),
        strip.text = element_text(color = "white"),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.1),
        panel.grid.minor.x = element_line(linewidth = 0.1),
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(lineheight = 0.6),
        aspect.ratio = 0.8) +
  scale_y_discrete(labels = label_wrap(35)) +
  scale_x_continuous(breaks = seq(0,40,by=10), labels = seq(0,40,by=10)) +
  geom_vline(xintercept = -log10(0.05), color = "black", linewidth = 0.75)


ggsave(palm_path_plot, filename = file.path(fig_dir, "palm_pathways_plot.png"), units = "in",
       width = 4.25, height = 5, dpi = 600)

#####################
### Joint heatmap ###
#####################

diff_set <- unique(c(res_table_aicar_lfc$gene_names_ens[abs(res_table_aicar_lfc$log2FoldChange) >= 1 & res_table_aicar_lfc$padj <= 0.05],
                     res_table_palm_lfc$gene_names_ens[abs(res_table_palm_lfc$log2FoldChange) >= 1 & res_table_palm_lfc$padj <= 0.05]))
diff_df <- vst_df %>%
  dplyr::select(-c(ATS0263:ATS0266)) %>%
  dplyr::filter(gene_names_ens %in% diff_set)

diff_mat <- as.matrix(diff_df[,c(2:12)])
row.names(diff_mat) <- diff_df$gene_names_ens

diff_mat_norm <- t(scale(t(diff_mat)))
diff_groups <- c(rep("Group 1", 4), rep("Group 2", 3),
                 rep("Group 3", 4))
diff_annot <- HeatmapAnnotation(
  condition = anno_block(gp = gpar(fill = c("#ffca3a", "#8ac926","#1982c4"),
                                   border = NA,
                                   lty = "blank"),
                         labels = c("Differentiated + Basal",
                                    "Diff. + AICAR",
                                    "Differentiated + Palmitate"),
                         labels_gp = gpar(fontsize = 12, fontfamily = "Helvetica", col = "white"))
)

diff_breaks <- seq(-2.78, 2.89, length.out = length(my_palette))

# Plot heatmap
diff_heatmap <- Heatmap(diff_mat_norm,
                        top_annotation = diff_annot,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = TRUE,
                        cluster_columns = FALSE,
                        show_row_dend = FALSE,
                        column_order = rev(1:11),
                        column_split = diff_groups,
                        col = colorRamp2(diff_breaks, my_palette),
                        column_title = NULL,
                        use_raster = FALSE,
                        heatmap_legend_param = list(
                          title = "Norm. expr."
                        )
)

png(file.path(fig_dir, "diff_heatmap.png"), units = "in",
    width = 6.5, height = 4, res = 600)
diff_heatmap
dev.off()

#####################
### RREB1 heatmap ###
#####################

rreb1_df <- vst_df %>%
  filter(hgnc_symbol == "RREB1")

rreb1_mat <- as.matrix(rreb1_df[,c(2:16)])
rreb1_mat <- as_tibble(t(scale(t(rreb1_mat)))) %>%
  pivot_longer(ATS0263:ATS0278) %>%
  mutate(gene = "RREB1")

rreb1_labels <- data.frame(
  x = c(2.5, 6.5, 10, 13.5),  # Midpoints of each group
  label = c("Undiff. +\nBasal", "Diff. +\nBasal",
            "Diff. +\nAICAR", "Diff. +\nPalmitate")
)

rreb1_heatmap <- ggplot(data = rreb1_mat, mapping = aes(x = name,
                                                        y = gene,
                                                        fill = value)) +
  geom_tile() +
  xlab(label = "Sample") +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFFF",
                       high = "#FF0000",
                       name = "Norm. Expr.") +
  labs(y = NULL, x = NULL) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 12),
        aspect.ratio = 0.075,
        line = element_blank(),
        legend.position = "top") +
  coord_cartesian(clip = "off") +
  geom_text(data = rreb1_labels,
            aes(x = x, y = -1, label = label),
            inherit.aes = FALSE, size = 4, vjust = 0.5,
            family = "Helvetica")

ggsave(rreb1_heatmap, filename = file.path(fig_dir, "rreb1_heatmap.png"),
       units = "in", width = 8, height = 2, res = 600, device = ragg::agg_png)


plag1_df <- vst_df %>%
  filter(hgnc_symbol %in% c("SP1","SP3","SP8"))

plag1_mat <- as.matrix(plag1_df[,c(2:16)])
plag1_mat <- as_tibble(t(scale(t(plag1_mat)))) %>%
  pivot_longer(ATS0263:ATS0278) %>%
  mutate(gene = "PLAG1")

plag1_labels <- data.frame(
  x = c(2.5, 6.5, 10, 13.5),  # Midpoints of each group
  label = c("Undiff. +\nBasal", "Diff. +\nBasal",
            "Diff. +\nAICAR", "Diff. +\nPalmitate")
)

plag1_heatmap <- ggplot(data = plag1_mat, mapping = aes(x = name,
                                                        y = gene,
                                                        fill = value)) +
  geom_tile() +
  xlab(label = "Sample") +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFFF",
                       high = "#FF0000",
                       name = "Norm. Expr.") +
  labs(y = NULL, x = NULL) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 12),
        aspect.ratio = 0.075,
        line = element_blank(),
        legend.position = "top") +
  coord_cartesian(clip = "off") +
  geom_text(data = plag1_labels,
            aes(x = x, y = -1, label = label),
            inherit.aes = FALSE, size = 4, vjust = 0.5,
            family = "Helvetica")

ggsave(rreb1_heatmap, filename = file.path(fig_dir, "rreb1_heatmap.png"),
       units = "in", width = 8, height = 2, res = 600, device = ragg::agg_png)

# download meme files for motifs from jaspar db
rreb1 <- read_meme("/path/to/MA0073.1.meme")
rreb1_logo <- view_motifs(jaspar[[11]]) + 
  annotate('rect', xmin = 12.5, xmax = 13.5, ymin = 0, ymax = 0.9, alpha = .1, col='black', fill='yellow') +
  annotate('segment', x = 12.5, xend = 13.5, y=1.2, yend=1.2, size=1) + 
  annotate('text', x=13, y=1.3, label='rs3810155 (C/G)',
           family = "Helvetica") +
  theme(axis.title.y = element_text(family = "Helvetica"),
        axis.text = element_text(family = "Helvetica")) +
  labs(y = "Bits")

ggsave(rreb1_logo, filename = file.path(fig_dir, "rreb1_logo.png"),
       units = "in", width = 7, height = 3, res = 600, device = ragg::agg_png)

sp_df <- vst_df %>%
  filter(hgnc_symbol %in% c("SP1", "SP2", "SP3"))

sp_mat <- as.matrix(sp_df[,c(2:16)])
sp_mat <- as_tibble(t(scale(t(sp_mat)))) %>%
  mutate(gene = sp_df$hgnc_symbol) %>%
  pivot_longer(ATS0263:ATS0278)

sp_labels <- data.frame(
  x = c(2.5, 6.5, 10, 13.5),  # Midpoints of each group
  label = c("Undiff. +\nBasal", "Diff. +\nBasal",
            "Diff. +\nAICAR", "Diff. +\nPalmitate")
)

sp_heatmap <- ggplot(data = sp_mat, mapping = aes(x = name,
                                                  y = gene,
                                                  fill = value)) +
  geom_tile() +
  xlab(label = "Sample") +
  scale_fill_gradientn(colors = c("#075AFF",
                                  "#FFFFFF",
                                  "#FF0000"),
                       name = "Norm. Expr.",
                       values = rescale(c(-2.5, 0, 0.8)),
                       guide = "colorbar",
                       limits = c(-2.5, 0.8),
                       na.value = "#FFFFFF") +
  labs(y = NULL, x = NULL) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 12),
        aspect.ratio = 0.075,
        line = element_blank(),
        legend.position = "top") +
  coord_cartesian(clip = "off") +
  geom_text(data = plag1_labels,
            aes(x = x, y = -1, label = label),
            inherit.aes = FALSE, size = 4, vjust = 0.5,
            family = "Helvetica")

ggsave(sp_heatmap, filename = file.path(fig_dir, "sp_heatmap.png"),
       units = "in", width = 8, height = 2, res = 600, device = ragg::agg_png)

sp1 <- read_meme(file.path("/path/to/MA0079.5.meme"))
sp1_logo <- view_motifs(sp1) + 
  annotate('rect', xmin = 2.5, xmax = 3.5, ymin = 0, ymax = 2, alpha = .1, col='black', fill='yellow') +
  theme(axis.title.y = element_text(family = "Helvetica"),
        axis.text = element_text(family = "Helvetica"),
        plot.title = element_text(family = "Helvetica", hjust = 0.5, size = 14)) +
  labs(y = "Bits", title = "SP1 motif (MA0079.5)")

ggsave(sp1_logo, filename = file.path(fig_dir, "sp1_logo.png"),
       units = "in", width = 7, height = 2, res = 600, device = ragg::agg_png)

sp4 <- read_meme(file.path("/path/to/MA0685.2.meme"))
sp4_logo <- view_motifs(sp4) + 
  annotate('rect', xmin = 2.5, xmax = 3.5, ymin = 0, ymax = 1.99, alpha = .1, col='black', fill='yellow') +
  theme(axis.title.y = element_text(family = "Helvetica"),
        axis.text = element_text(family = "Helvetica"),
        plot.title = element_text(family = "Helvetica", hjust = 0.5, size = 14)) +
  labs(y = "Bits", title = "SP4 motif (MA0685.2)")


ggsave(sp4_logo, filename = file.path(fig_dir, "sp4_logo.png"),
       units = "in", width = 7, height = 2, res = 600, device = ragg::agg_png)

logo_plot <- sp4_logo / sp1_logo
ggsave(logo_plot, filename = file.path(fig_dir, "sp_logo.png"),
       units = "in", width = 4, height = 4.28, res = 600, device = ragg::agg_png)

undiff_motif_df <- vst_df %>%
  filter(hgnc_symbol %in% c("KLF3", "RFX2", "KLF6",
                            "RREB1", "SP8", "RFX5",
                            "JUND", "SP4", "SP1",
                            "HSF2"))

undiff_mat <- as.matrix(undiff_motif_df[,c(2:16)])
undiff_mat <- as_tibble(t(scale(t(undiff_mat)))) %>%
  mutate(gene = undiff_motif_df$hgnc_symbol) %>%
  pivot_longer(ATS0263:ATS0278) %>%
  mutate(cond = case_when(name %in% c("ATS0263", "ATS0264", "ATS0265", "ATS0266") ~ "undiff",
                          name %in% c("ATS0267", "ATS0268", "ATS0269", "ATS0270") ~ "diff",
                          name %in% c("ATS0271", "ATS0273", "ATS0274") ~ "aicar",
                          name %in% c("ATS0275", "ATS0276", "ATS0277", "ATS0278") ~ "palm"),
         rep = case_when(name %in% c("ATS0263", "ATS0267", "ATS0271", "ATS0275") ~ "rep1",
                         name %in% c("ATS0264", "ATS0268", "ATS0273", "ATS0276") ~ "rep2",
                         name %in% c("ATS0265", "ATS0269", "ATS0274", "ATS0277") ~ "rep3",
                         name %in% c("ATS0266", "ATS0270", "ATS0278") ~ "rep4"),
         cond = factor(cond, levels = c("undiff","diff","aicar","palm")),
         gene = factor(gene, levels = rev(c("KLF3", "RFX2", "KLF6", "RREB1",
                                            "RFX5", "JUND", "SP4",
                                            "SP1", "HSF2"))))

undiff_labels <- c("Undifferentiated +\nBasal", "Differentiated +\nBasal",
                   "Differentiated +\nAICAR", "Differentiated +\nPalmitate")
names(undiff_labels) <- c("undiff", "diff", "aicar", "palm")

rect_annot_df <- undiff_mat %>%
  group_by(cond) %>%
  summarize(
    xmin = 0.5,
    xmax = n()/9,
    ymin = 1.5,
    ymax = 3.5
  )

undiff_heatmap <- ggplot(data = undiff_mat, mapping = aes(x = rep,
                                                          y = gene,
                                                          fill = value)) +
  geom_tile() +
  facet_grid(.~cond, space = "free_x",
             labeller = labeller(cond = undiff_labels),
             switch = "x") +
  xlab(label = "Sample") +
  scale_fill_gradientn(colors = c("#075AFF",
                                  "#FFFFFF",
                                  "#FF0000"),
                       name = "Norm. Expr.",
                       values = rescale(c(-3.15, 0, 1.9)),
                       guide = "colorbar",
                       limits = c(-3.15, 1.9),
                       na.value = "#FFFFFF") +
  labs(y = NULL, x = NULL) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 12, face = "italic"),
        line = element_blank(),
        legend.position = "top",
        strip.background = element_rect(color = NA),
        strip.text = element_text(color = "black"),
        panel.spacing = unit(0.1, "lines")) +
  geom_rect(data = rect_annot_df,
            aes(xmin = xmin,
                xmax = xmax + 0.5,
                ymin = ymin,
                ymax = ymax
                ),
            inherit.aes = FALSE, color = "black", fill = NA)

ggsave(undiff_heatmap, filename = file.path(fig_dir, "undiff_heatmap.png"),
       units = "in", width = 7, height = 4.2, res = 600, device = ragg::agg_png)

# pairwise corr between expression of SP1, SP4 and activity at rs490972 ref allele
rs490_df <- read.delim("/path/to/workdir/mpra/output/rs490972_test-df.tsv")
rs490_df <- rs490_df %>%
  mutate(rep = case_when(replicate %in% c("ratio1", "ratio5", "ratio9", "ratio13") ~ "rep1",
                         replicate %in% c("ratio2", "ratio6", "ratio10", "ratio14") ~ "rep2",
                         replicate %in% c("ratio3", "ratio7", "ratio11", "ratio15") ~ "rep3",
                         replicate %in% c("ratio4", "ratio8", "ratio12", "ratio16") ~ "rep4")) %>%
  dplyr::select(-c(.groups, value, refname)) %>%
  mutate(condition = as.factor(condition))

undiff_motif_df <- vst_df %>%
  filter(hgnc_symbol %in% c("KLF3", "RFX2", "KLF6",
                            "RREB1", "SP8", "RFX5",
                            "JUND", "SP4", "SP1",
                            "HSF2", "SLC44A2", "KRI1")) %>%
  dplyr::rename(gene = hgnc_symbol) %>%
  pivot_longer(ATS0263:ATS0278) %>%
  mutate(cond = case_when(name %in% c("ATS0263", "ATS0264", "ATS0265", "ATS0266") ~ "undiff",
                          name %in% c("ATS0267", "ATS0268", "ATS0269", "ATS0270") ~ "diff",
                          name %in% c("ATS0271", "ATS0273", "ATS0274") ~ "aicar",
                          name %in% c("ATS0275", "ATS0276", "ATS0277", "ATS0278") ~ "palm"),
         rep = case_when(name %in% c("ATS0263", "ATS0267", "ATS0271", "ATS0275") ~ "rep1",
                         name %in% c("ATS0264", "ATS0268", "ATS0273", "ATS0276") ~ "rep2",
                         name %in% c("ATS0265", "ATS0269", "ATS0274", "ATS0277") ~ "rep3",
                         name %in% c("ATS0266", "ATS0270", "ATS0278") ~ "rep4"),
         cond = factor(cond, levels = c("undiff","diff","aicar","palm")),
         gene = factor(gene, levels = rev(c("KLF3", "RFX2", "KLF6", "RREB1",
                                            "RFX5", "JUND", "SP4",
                                            "SP1", "HSF2", "SLC44A2", "KRI1"))))

rs490_df <- undiff_motif_df %>%
  dplyr::select(-c(name, gene_names_ens)) %>%
  inner_join(rs490_df, by = c("rep" = "rep", "cond" = "condition"))

rs490_cor <- rs490_df %>%
  group_by(allele, gene) %>%
  summarize(coef = cor.test(value, plot_value, method = "spearman")$estimate,
            pval = cor.test(value, plot_value, method = "spearman")$p.value)

rs490_cor_sp <- rs490_cor %>%
  filter(gene %in% c("SP1","SP4") & allele == "R") %>%
  mutate(value = c(9.1, 5.8),
         plot_value = c(0.625, 0.625),
         rlabel = paste0("R^2 == ", round(coef, 3)),
         plabel = paste0("p == ", round(pval, 3)))

sp_cor_plot <- rs490_df %>%
  filter(gene %in% c("SP1","SP4") & allele == "R") %>%
  mutate(cond = factor(cond, levels = c("undiff", "diff", "aicar", "palm")),
         gene = factor(gene, levels = c("SP4", "SP1"))) %>%
  ggplot(aes(x = value, y = plot_value, fill = cond)) +
  facet_wrap(~gene, ncol = 1, scales = "free") +
  geom_point(pch = 21, size = 3) +
  scale_fill_manual(values = c("#ff595e","#ffca3a","#8ac926","#1982c4"),
                    labels = c("Undiff", "Diff", "AICAR", "Palmitate"),
                    name = "Condition") +
  theme_bw(base_size = 12, base_family = "Helvetica") + labs(x = "VST-normalized gene expression",
                                                             y = parse(text = "'rs490972-G MPRA activity, ' * log[2]~(RNA/DNA)")) +
  theme(strip.text = element_text(face = "bold.italic", size = 12, family = "Helvetica"),
        strip.background = element_rect(linewidth=0.5),
        axis.title = element_text(family = "Helvetica")) +
  geom_text(rs490_cor_sp, mapping = aes(x = value, y = plot_value, label = rlabel),
            inherit.aes = FALSE, parse = TRUE, hjust = 0, family = "Helvetica") +
  geom_text(rs490_cor_sp, mapping = aes(x = value, y = plot_value-0.2, label = plabel),
            inherit.aes = FALSE, parse = TRUE, hjust = 0, family = "Helvetica")

ggsave(sp_cor_plot, filename = file.path(fig_dir, "sp_cor_plot.png"),
       units = "in", height = 4.54, width = 4.25, dpi = 600, device = ragg::agg_png())


# pairwise corr between expression of RREB1 and activity at rs3810155?

rs381_df <- read.delim("/path/to/workdir/mpra/output/rs3810155_test-df.tsv")
rs381_df <- rs381_df %>%
  mutate(rep = case_when(replicate %in% c("ratio1", "ratio5", "ratio9", "ratio13") ~ "rep1",
                         replicate %in% c("ratio2", "ratio6", "ratio10", "ratio14") ~ "rep2",
                         replicate %in% c("ratio3", "ratio7", "ratio11", "ratio15") ~ "rep3",
                         replicate %in% c("ratio4", "ratio8", "ratio12", "ratio16") ~ "rep4")) %>%
  dplyr::select(-c(.groups, value, refname)) %>%
  mutate(condition = as.factor(condition))

rs381_df <- undiff_motif_df %>%
  dplyr::select(-c(name, gene_names_ens)) %>%
  inner_join(rs381_df, by = c("rep" = "rep", "cond" = "condition"))

rs381_cor <- rs381_df %>%
  group_by(allele, gene) %>%
  summarize(coef = cor.test(value, plot_value, method = "spearman")$estimate,
            pval = cor.test(value, plot_value, method = "spearman")$p.value)

rs381_cor_rreb <- rs381_cor %>%
  filter(gene == "RREB1" & allele == "R") %>%
  mutate(value = 8.35,
         plot_value = 1.55,
         rlabel = paste0("R^2 == ", round(coef, 3)),
         plabel = paste0("p == ", round(pval, 3)))

rs381_df %>%
  filter(gene == "RREB1" & allele == "R") %>%
  mutate(cond = factor(cond, levels = c("undiff", "diff", "aicar", "palm"))) %>%
  ggplot(aes(x = value, y = plot_value, fill = cond)) +
  geom_point(pch = 21, size = 3) +
  scale_fill_manual(values = c("#ff595e","#ffca3a","#8ac926","#1982c4"),
                    labels = c("Undiff", "Diff", "AICAR", "Palmitate"),
                    name = "Condition") +
  theme_bw(base_size = 12, base_family = "Helvetica") + labs(x = "VST-normalized counts",
                                                             y = parse(text = "log[2]~(RNA/DNA) * ', rs3810155-G'")) +
  theme(axis.title = element_text(family = "Helvetica")) +
  geom_text(rs381_cor_rreb, mapping = aes(x = value, y = plot_value, label = rlabel),
            inherit.aes = FALSE, parse = TRUE, hjust = 0, family = "Helvetica") +
  geom_text(rs381_cor_rreb, mapping = aes(x = value, y = plot_value-0.15, label = plabel),
            inherit.aes = FALSE, parse = TRUE, hjust = 0, family = "Helvetica")

ggsave(sp_cor_plot, filename = file.path(fig_dir, "sp_cor_plot.png"),
       units = "in", height = 4.25, width = 4.25, dpi = 600, device = ragg::agg_png())

rs490_df %>%
  filter(gene == "SP1" & allele == "R") %>%
  ggplot(aes(x = value, y = plot_value, color = cond)) +
  geom_point() + geom_abline(slope = 1, linetype = "dashed")
