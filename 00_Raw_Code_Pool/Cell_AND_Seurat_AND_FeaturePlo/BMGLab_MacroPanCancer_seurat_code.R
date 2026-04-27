library(Seurat)
library(dplyr)
library(hdf5r)
library(DoMultiBarHeatmap)

# Read data
cell_annotation.data <- read.delim("path/to/file", row.names = 1, header = TRUE)
umi_matrix.data <- read.delim("path/to/file",sep = "\t", header = TRUE)
#tpm_matrix.data <- read.delim("../Desktop/LAB/DATASETS/M7/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt", row.names = 1, sep = "\t", header = TRUE)

# Create seurat object
seurat_obj <- CreateSeuratObject(counts = umi_matrix.data, meta.data = cell_annotation.data , min.cells = 3, min.features = 200)


# If there are existing idents subset them before 
dim(seurat_obj)
Idents(object = seurat_obj) <- seurat_obj$Cell_subtype
seurat_obj <- subset(x = seurat_obj, idents = c("mo-Mac", "Alveolar Mac"))
dim(seurat_obj)

# Normalization (dont use if data is tpm)
seurat_obj <- NormalizeData(seurat_obj)

# If tpm used use this code
#seurat_obj <- SetAssayData(object = seurat_obj, slot="data", new.data = log1p(tpm_matrix.data))

# Find variable features 
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Özellik sıralamasını ve boyutlandırmayı yap
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.4))
seurat_obj <- RunUMAP(seurat_obj, dims = (1:10))
DimPlot(seurat_obj, reduction = "umap", label = T)

# macrophage control
VlnPlot(seurat_obj, features = c("LYZ", "CD68"))
VlnPlot(seurat_obj, features = c("CD163", "IL4I1"))
VlnPlot(seurat_obj, features = c("CD14", "CD64"))
#VlnPlot(seurat_obj, features = c("CD71", "CCR5"))
#VlnPlot(seurat_obj, features = c("CD16", "CD206"))

# dentritic cell control
VlnPlot(seurat_obj, features = c("IL3RA", "IRF4"))

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat_obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# find all markers of any cluster
cluster9.markers <- FindMarkers(seurat_obj, ident.1 = 9, min.pct = 0.25)
head(cluster9.markers, n = 5)
cluster10.markers <- FindMarkers(seurat_obj, ident.1 = 10, min.pct = 0.25)
head(cluster10.markers, n = 5)
cluster7.markers <- FindMarkers(seurat_obj, ident.1 = 7, min.pct = 0.25)
head(cluster7.markers, n = 5)
cluster11.markers <- FindMarkers(seurat_obj, ident.1 = 11, min.pct = 0.25)
head(cluster11.markers, n = 5)

# rename clusters according to markers
new.cluster.ids <- c("0","1", "2", "3","4","5","6","7","8", "Macrophages", "10")
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

#subset
Idents(object = seurat_obj) <- seurat_obj$seurat_clusters
seurat_obj <- RenameIdents(seurat_obj, '9' = "Head and Neck Macrophages")
seurat_obj <- subset(x = seurat_obj, subset = (HN2_seurat_M))
dim(seurat_obj)

# Identification of M1/M2 macrophages
## Markers for M1 macrophages
FeaturePlot(GSE130148_Lung, features = c("CD80", "CD86"), label = F)
## Markers for M2 macrophages
FeaturePlot(GSE130148_Lung, features = c("CD163", "CD206", "CD68"), label = F)

dim(seurat_obj)
saveRDS(seurat_obj, file = "HN2_seurat.RDS")
