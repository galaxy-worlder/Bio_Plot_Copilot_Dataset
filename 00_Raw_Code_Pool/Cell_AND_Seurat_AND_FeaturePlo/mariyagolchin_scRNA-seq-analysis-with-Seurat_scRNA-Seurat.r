
setRepositories(ind = 1:4)  # This will select the first four repositories (CRAN and CRAN mirrors)

install.packages("tidyverse", repos = "https://cran.r-project.org")
library(tidyverse)
library(Seurat)
set.seed(42)
#The seed is given an integer value to ensure that the results of pseudo-random generation are reproducible
# 1. Import data ===================================================
# /home/golchinpour/projects/1-scRNA-seq-seurat/40k_NSCLC_DTC_3p_HT_nextgem_donor_1_count_sample_feature_bc_matrix.h5
# nsclc_sparcematrix 
nsclc_sm <- Read10X_h5(filename = "40k_NSCLC_DTC_3p_HT_nextgem_donor_1_count_sample_feature_bc_matrix.h5")
str(nsclc_sm) # Check the multiple modalities (list of matrixes) - we're interested in Gene expression
cts <- nsclc_sm$`Gene Expression`
# count matrix-sparce matrix
cts[345:350, 1380:1382]
# row:genes
# colomns:cells

# 6 x 3 sparse Matrix of class "dgCMatrix"
#            CTAGGTACAGGACAGT-1 CTAGGTAGTACTGCCG-1 CTAGGTAGTTCGGCCA-1
# TMEM51                      .                  .                  .
# C1orf195                    .                  .                  .
# FHAD1                       .                  .                  .
# AL031283.2                  .                  .                  .
# AL031283.3                  .                  .                  .
# AL031283.1                  .                  .                  .
# > 
> cts[c('EFHD2','CASP9'), 1380:1382]
# tow gene in row
# 3 sample in colomns
# ================================================
# 2 x 3 sparse Matrix of class "dgCMatrix"
#       CTAGGTACAGGACAGT-1 CTAGGTAGTACTGCCG-1 CTAGGTAGTTCGGCCA-1
# EFHD2                  .                  5                  .
# CASP9                  1                  .                  .
# > 

# 2. Create your Seurat object (raw counts) ===========================================
nsclc_seu <- CreateSeuratObject(counts = cts, project = 'NSCLC', min.cells = 3, min.features = 200)

# min.cells = 3 ==> cell that have at least non zero value for 3 genes
# سلولی که حداقل برای سه تا ژن مقدار غیر صفر داشته باشه

#min.features = 200 ==> featuears that are genes, are expressed in at least 200 cells
فیچر ها که ژن هامون هستند خداقل تو 
200 تا سلول بیان شدن
str(nsclc_seu)

# =========================================================
# 1. QC metrices
# =========================================================
# 1- cell with low genes
# 2- cell with high genes
# 3- mitochondrial gene

# mitochondrial gene transcribe in mitocondria, not in sitoplasem
# پس اگه درصد ژن های میتوکندری زیاد باشه، یعنی سلول باکیفیت پایین هست

# پس مطمعن میشیم که سلول های با درصد بالای ژن های میتوکوندریال نداریم.

# PercentageFeatureSet درصد ژن های با یه کاندیشن خاص رو با این تابع در میارن

# mitochondrial gene for human start with "MT"
# calculate mitochondrial gene that store in metadata
nsclc_seu[['percent_mt']] <- PercentageFeatureSet(nsclc_seu, pattern = '^MT-')
> head(nsclc_seu@meta.data)
#                    orig.ident nCount_RNA nFeature_RNA percent_mt
# AAACCCAGTTCTCACC-1      NSCLC       2975         1194  3.9663866
# AAACGAACAACATACC-1      NSCLC      11225         2825  0.9977728
# AAACGAATCATTACCT-1      NSCLC       1700          855  2.1176471
# AAACGCTAGTTAGTAG-1      NSCLC        957          627  8.0459770
# AAACGCTTCCATCTCG-1      NSCLC      30735         5762  3.5464454
# AAAGAACAGCACCGAA-1      NSCLC       2139         1103  2.9920524
# > 

# Create the violin plot
vln_plot <- VlnPlot(nsclc_seu, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)

# Save the plot
ggsave("1-vln_plot.png", plot = vln_plot, width = 12, height = 6, dpi = 300)

# =====================================================
# Create the feature scatter plot with a linear regression line
scatter_plot <- FeatureScatter(nsclc_seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm') +
  theme_bw()
# adds a linear regression line to the scatter plot.
# Save the scatter plot
ggsave("2-FeatureScatter_plot.png", plot = scatter_plot, width = 12, height = 6, dpi = 300)

========================================================
## Filtering 
=========================================================
#cell with higher percent of MT >5 filter
nsclc_seu <- subset(nsclc_seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mt < 5)
#
=============================================================
## Normalisation 
=============================================================
nsclc_seu <- NormalizeData(nsclc_seu, normalization.method = 'LogNormalize', scale.factor = 10000)
nsclc_seu <- NormalizeData(nsclc_seu) # دیفالت پارامترای بالا رو داره

=============================================================
## Identify highly-variable features ===========================================
=============================================================
nsclc_seu <- FindVariableFeatures(nsclc_seu, selection.method =  'vst', nfeatures = 2000)
# identify the most variable genes in a single-cell RNA sequencing dataset
#  to identify genes that show the most variation in expression across all the cells.
#  selection method: dispersion or vst
# 'vst' (Variance Stabilizing Transformation)
#  is a method that normalizes the data to account for the relationship between mean expression and variance

# بای دیفالت 2000 فیچر برمیگردونه ولی میتونیم تعدادش تغییر بدیم
# Identify the top 10 HVGs
# Get the top 10 variable features
# Get the top 10 variable features


> head(VariableFeatures(nsclc_seu))
[1] "IGHG4" "MT1G"  "IGHA1" "IGLC2" "IGHGP" "DCN"  
> 
top10 <- head(VariableFeatures(nsclc_seu), 10)
# > top10
#  [1] "IGHG4"   "MT1G"    "IGHA1"   "IGLC2"   "IGHGP"   "DCN"     "JCHAIN" 
#  [8] "LYZ"     "IGHG3"   "RARRES2"
# Create a VariableFeaturePlot
top10_plot <- VariableFeaturePlot(nsclc_seu)

# Label the top 10 variable features on the plot
top10_plot <- LabelPoints(plot = top10_plot, points = top10, repel = TRUE)

# Set the background to white
top10_plot <- top10_plot + theme_bw()

# Save the plot with a white background
ggsave("3-VariableFeaturePlot_top10_white.png", plot = top10_plot, width = 12, height = 6, dpi = 300)

# =================================================================
        #   scale data      
        #   to remove unwanted sorces of variations, such as bacheffect, biological variation,....
        #  a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
# =================================================================

all_genes <- rownames(nsclc_seu)
nsclc_seu <- ScaleData(nsclc_seu, features = all_genes)
# nsclc_seu@assays$RNA

============================================================================
## Dimetionality reduction:PCA
============================================================================

# Run PCA
nsclc_seu <- RunPCA(nsclc_seu, features = VariableFeatures(nsclc_seu))
# visualize PCA results a few different ways
print(nsclc_seu[["pca"]], dims = 1:5, nfeatures = 5)
#  first 5 principal components (PCs 1 to 5).
# The first principal component (PC1) explains the greatest amount of variance,
#  the second principal component (PC2) explains the next greatest amount of variance,

# How to Interpret Each PC:
# Each PC is a linear combination of genes:
#  The genes listed under each PC (positive or negative) 
#  are those that contribute most to defining that component. 
#  A positive contribution means that the gene is positively correlated with the PC,
#   while a negative contribution means the gene is negatively correlated.

# Principal components are orthogonal: Each PC is independent (uncorrelated) with the others,
#  meaning that PC1 captures the maximum variance along one axis, 
#  and PC2 captures the next largest amount of variance in a different, independent direction, and so on.

Viz_nsclc_seu <- VizDimLoadings(nsclc_seu, dims = 1:2, reduction = "pca")
ggsave("4-VizDimLoadings.png", plot = Viz_nsclc_seu, width = 12, height = 8, dpi = 300)

#####
dimplot_nsclc_seu<- DimPlot(nsclc_seu, reduction = "pca") + NoLegend()
ggsave("5-dimplot.png", plot = dimplot_nsclc_seu, width = 12, height = 8, dpi = 300)

# PCA Heatmap
# ????????????????????????????//
pca_heatmap <- DimHeatmap(nsclc_seu, dims = 1:15, cells = 50, balanced = TRUE)
ggsave("6-PCA_Heatmap.png", plot = pca_heatmap, width = 12, height = 8, dpi = 300)


===================================================================
#####  Determine the ‘dimensionality’ of the dataset
# Elbow Plot

#  Seurat clusters cells based on their PCA scores, with each PC 
#  essentially representing a ‘metafeature’ that combines information across a correlated feature set.

# Elbow plot’: a ranking of principle components based on the percentage
#  of variance explained by each one (ElbowPlot() function).
#   In this example, we can observe an ‘elbow’ around PC9-10, 
#   suggesting that the majority of true signal is captured in the first 10 PCs.



elbow_plot <- ElbowPlot(nsclc_seu)+ theme_bw()
ggsave("7-ElbowPlot.png", plot = elbow_plot, width = 12, height = 8, dpi = 300)

===================================================
##Cluster the cells
===================================================
nsclc_seu <- FindNeighbors(nsclc_seu, dims = 1:15)
nsclc_seu <- FindClusters(nsclc_seu, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
# View(nsclc_seu@meta.data)
# Create the DimPlot
# DimPlot is a function in Seurat used to create dimensionality reduction plots.
# By default, DimPlot uses the UMAP or t-SNE coordinates

# group.by = 'RNA_snn_res.0.1'
# This argument specifies the grouping variable for coloring the points in the plot.
#  The cells will be colored based on this grouping factor.

# It represents the clustering results based on the shared nearest neighbor (SNN) method for clustering. 
# The value 0.1 indicates the resolution used for clustering:
# A higher resolution (e.g., 0.5 or 1.0) typically results in more, smaller clusters.
# A lower resolution (e.g., 0.1) results in fewer, larger clusters.
# This parameter assigns each cell to a cluster, and the DimPlot function will color the cells according to their cluster assignment.

dimplot5 <- DimPlot(nsclc_seu, group.by = 'RNA_snn_res.0.5', label = TRUE)
ggsave("8-5-DimPlot_RNA_snn_res_0.5.png", plot = dimplot5, width = 12, height = 8, dpi = 300)

dimplot1 <- DimPlot(nsclc_seu, group.by = 'RNA_snn_res.0.1', label = TRUE)
ggsave("8-1-DimPlot_RNA_snn_res_0.1.png", plot = dimplot1, width = 12, height = 8, dpi = 300)
# Save the DimPlot as a PNG file


# The clusters can be found using the Idents() function.
head(Idents(nsclc_seu), 5)


saveRDS(nsclc_seu, file = 'nsclc_seu.RDS')

# =================================================================
# Finding differentially expressed features (cluster biomarkers)
# =================================================================
# Seurat can help you find markers that define clusters via differential expression (DE).
# By default, it identifies positive and negative markers of a single cluster (specified in ident.1),
#  compared to all other cells. FindAllMarkers() automates this process for all clusters

# find all markers of cluster 2
cluster2.markers <- FindMarkers(nsclc_seu, ident.1 = 2)
head(cluster2.markers, n = 5)
> head(cluster2.markers, n = 5)
#                p_val avg_log2FC pct.1 pct.2    p_val_adj
# BTG1    2.665721e-34   1.295503 0.987 0.883 6.261512e-30
# CLEC2B  1.949733e-30   1.446552 0.895 0.518 4.579729e-26
# TNFAIP3 3.445764e-30   1.640095 0.850 0.513 8.093755e-26
# YPEL5   1.161465e-28   1.446788 0.889 0.537 2.728165e-24
# TUBA4A  1.223130e-28   2.065607 0.739 0.386 2.873010e-24
# > 

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(nsclc_seu, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

> head(cluster5.markers, n = 5)
#                p_val avg_log2FC pct.1 pct.2    p_val_adj
# CREM    3.619695e-56   4.216724 0.911 0.152 8.502301e-52
# TNFRSF9 8.233142e-53   4.342475 0.741 0.043 1.933883e-48
# DUSP4   2.829214e-48   3.214177 0.893 0.189 6.645542e-44
# TNFAIP3 1.082674e-41   2.466625 0.946 0.265 2.543093e-37
# CXCL13  1.427097e-40   3.152509 0.723 0.082 3.352108e-36
# > 

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
nsclc_seu.markers <- FindAllMarkers(nsclc_seu, only.pos = TRUE)
nsclc_seu.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(nsclc_seu, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


# We include several tools for visualizing marker expression. 
# VlnPlot() (shows expression probability distributions across clusters),
#  and FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. 

VlnPlot_cb <- VlnPlot(nsclc_seu, features = c("CREM", "BTG1"))
ggsave("9-VlnPlot_cb.png", plot = VlnPlot_cb, width = 12, height = 8, dpi = 300)

FeaturePlot_cb <- FeaturePlot(nsclc_seu, features = c("CREM", "BTG1"))
ggsave("10-FeaturePlot_cb.png", plot = FeaturePlot_cb, width = 12, height = 8, dpi = 300)

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

nsclc_seu.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

DoHeatmap <- DoHeatmap(nsclc_seu, features = top10$gene) + NoLegend()
ggsave("11-DoHeatmap.png", plot = DoHeatmap, width = 12, height = 8, dpi = 300)

# ???? for our data???
# # Assigning cell type identity to clusters
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#     "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

reff: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
