install.packages('Seurat')
library(Seurat)
library(dplyr)

rna_path <-"C:/Users/zjohn/Desktop/single_nuc" 

b1.data <- Read10X(data.dir = paste(rna_path, "/filtered_feature_bc_matrix/", sep=""))
b2.data <- Read10X(data.dir = paste(rna_path, "/filtered_feature_bc_matrix/", sep=""))
c1.data <- Read10X(data.dir = paste(rna_path, "/filtered_feature_bc_matrix/", sep=""))

b1 <- CreateSeuratObject(counts = b1.data, project = "behav")
b2 <- CreateSeuratObject(counts = b2.data, project = "behav2")
c1 <- CreateSeuratObject(counts = c1.data, project = "control")

b1$cond <- "BHVE"
b2$cond <- "BHVE"
c1$cond <- "CTRL"

# Normalize Data without scale factor per https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
b1 <- NormalizeData(b1)
b2 <- NormalizeData(b2)
c1 <- NormalizeData(c1)

# Normalize Data with scale factor per https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
b1 <- NormalizeData(b1, normalization.method = "LogNormalize", scale.factor = 100000)
b2 <- NormalizeData(b2, normalization.method = "LogNormalize", scale.factor = 100000)
c1 <- NormalizeData(c1, normalization.method = "LogNormalize", scale.factor = 100000)

combined <- merge(x=c1, y=c(b1,b2), merge.data = TRUE, add.cell.ids = c("CTRL", "BHVE", "BHVE"))

combined <- subset(combined, subset = nFeature_RNA > 200)
combined <- subset(combined, subset = nFeature_RNA < 1500)
combined <- subset(combined, subset = nCount_RNA > 200)
combined <- subset(combined, subset = nCount_RNA < 1500)

# combined <- FindVariableFeatures(object = combined, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
combined <- FindVariableFeatures(object = combined, selection.method = "vst", loess.span = 0.3, dispersion.function = "default", clip.max = "auto", num.bin = 50, binning.method = "equal_width", verbose = TRUE)
plot1 <- VariableFeaturePlot(combined)
plot1

all.genes <- rownames(combined)
length(all.genes)
combined <- ScaleData(combined, features = all.genes)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))

#clustering
combined <- FindNeighbors(combined, dims = 1:25)
combined <- FindClusters(combined, resolution = 1.2)

combined <- RunUMAP(combined, dims = 1:20)
DimPlot(combined, reduction = "umap", pt.size = 0.5, label = TRUE)



#from here onward still playing around, some commands may not work or be really time intensive

my.data=fetch.data(nbt,c("ident","PC1","nGene","orig.ident","PAX6","DLX2","ACTB"))

#find the top 20 distinguishing genes for each cluster
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(head(combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC),400), "test.csv") #write them into a .csv file
top5 <- combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top3 <- combined.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
top2 <- combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


my.data=FetchData(combined,)
                   
doHeatMap(nbt,features = combined.markers,slim.col.label = TRUE,remove.key = TRUE,cexRow=0.1)

#another command for finding markers that distinguish all clusters
markers.all=FindAllMarkers(combined,test.use = "poisson", logfc.threshold=0.1, min.pct=.25, label = TRUE, verbose=TRUE)
?FindAllMarkers


#picking top marker genes by cluster
top5 <- subset(markers.all %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)) #this is how to get genes with biggest expression differences
top5 <- subset(markers.all %>% group_by(cluster) %>% top_n(n = -4, wt = p_val_adj)) #this is how to get most significant genes
top3 <- subset(markers.all %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC))
top3 <- subset(markers.all %>% group_by(cluster) %>% top_n(n = -1, wt = p_val_adj))
?top_n
top3
top2 <- markers.all %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#all markers found by FindAllMarkers
#DoHeatmap(combined, features = markers.all$gene, cells = NULL, group.by = "ident",
#          group.bar = TRUE, group.colors = NULL, disp.min = -2.5,
#          disp.max = NULL, slot = "scale.data", assay = NULL, label = TRUE,
#          size = 5.5, hjust = 0, angle = 45, raster = TRUE,
#          draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
#          combine = TRUE)

#plot top 5 for each cluster
DoHeatmap(combined, features = top5$gene, cells = NULL, group.by = "ident",
          group.bar = TRUE, group.colors = NULL, disp.min = -2.5,
          disp.max = NULL, slot = "scale.data", assay = NULL, label = TRUE,
          size = 5.5, hjust = 0, angle = 45, raster = TRUE,
          draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
          combine = TRUE)

#plot top 3 for each cluster
DoHeatmap(combined, features = top3$gene, cells = NULL, group.by = "ident",
          group.bar = TRUE, group.colors = NULL, disp.min = -2.5,
          disp.max = NULL, slot = "scale.data", assay = NULL, label = TRUE,
          size = 5.5, hjust = 0, angle = 45, raster = TRUE,
          draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
          combine = TRUE)

top5_list <- unique(top5$gene)
top2_list <- unique(top2$gene)
DotPlot(combined,features = top5_list)
DotPlot(combined,features = top2_list, angle = 45)
?DotPlot


print(head(combined.markers,40))
nbt = set.all.ident(nbt, "orig.ident")
dot.plot(nbt,genes.plot = rownames(combined.markers)[1:200],cex.use=4)

# Run the standard workflow for visualization and clustering
combined <- ScaleData(object = combined, vars.to.regress = NULL)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

unique(combined.markers.test$gene)


# t-SNE and Clustering
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "umap", dims = 1:2)
combined <- FindClusters(combined, resolution = 0.005)

DimPlot(combined, reduction = "umap", split.by = "cond", label = TRUE)
FeaturePlot(combined, features = c("slc17a6b"), split.by = "cond", reduction = "umap", pt.size = 1.5, label=TRUE, order = TRUE)

# Find clusters
folder <- "Dim25"
saveRDS(combined, file = paste(rna_path, "/results/clusters/", folder, "/combined.RDS", sep=""))
print("Finding DEG between clusters")
Idents(object = combined) <- "seurat_clusters"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  print(i)
  nk.markers <- FindMarkers(combined, ident.1 = i, verbose = FALSE)
  nk.markers$gene_name <- row.names(nk.markers)
  sig_nk.markers <- nk.markers[which(nk.markers$p_val_adj < 0.05 & abs(nk.markers$avg_logFC) > 2),]
  write.table(nk.markers, file = paste(rna_path, "/results/clusters/", folder, "/all_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE)
  write.table(sig_nk.markers, file = paste(rna_path, "/results/clusters/", folder, "/sig_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE)
  write.table(sig_nk.markers$gene_name, file = paste(rna_path, "/results/clusters/", folder, "/genes_sig_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE)
  sig_nk_pos.markers <- nk.markers[which(nk.markers$p_val_adj < 0.05 & nk.markers$avg_logFC > 2),]
  write.table(sig_nk_pos.markers, file = paste(rna_path, "/results/clusters/", folder, "/sig_pos_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE)
  write.table(sig_nk_pos.markers$gene_name, file = paste(rna_path, "/results/clusters/", folder,"/genes_sig_pos_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE)
}
print("Done finding clusters")