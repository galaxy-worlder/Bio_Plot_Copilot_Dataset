library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(scMC)
library(cowplot)
library(dplyr)


setwd("path to working directory")

dirs <- list.dirs(path = 'path/to folder/', recursive = F, full.names = F)
dirs

##loop through each folder to fetch the files
for(x in dirs){
  name <- gsub('_processed','', x)
  
  cts <- ReadMtx(mtx = paste0('/path/to folder/',x,'/matrix.mtx'),
                 features = paste0('/path/to folder/',x,'/genes.tsv'),
                 cells = paste0('/path/to folder/',x,'/barcodes.tsv'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts, project = "Vitiligo", min.cells = 3, min.features = 200))
}

ls()

sample_names <- c("Pt1", "Pt2", "Pt3", "Pt4", "Pt5", "Pt6")  # Update with your actual sample names
object.list <- lapply(sample_names, function(x) get(x))
names(object.list) <- sample_names


for (i in 1:length(object.list)) {
  object.list[[i]][["percent.mt"]] <- PercentageFeatureSet(object.list[[i]], pattern = "^MT-")
}

#Filter cells based on QC metrics
for (i in 1:length(object.list)) {
  object.list[[i]] <- subset(object.list[[i]], subset = nFeature_RNA < 7000 & nCount_RNA < 40000 & percent.mt < 10)
}

view(object.list[[1]]@meta.data)

for (i in 1:length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], verbose = T)
  object.list[[i]] <- FindVariableFeatures(object.list[[i]], selection.method = "mean.var.plot", 
                                           verbose = T, mean.cutoff = c(0.01, 5), 
                                           dispersion.cutoff = c(0.25, Inf))
  object.list[[i]] <- ScaleData(object.list[[i]], verbose = T)
}


object.list <- identifyNeighbors(object.list)

object.list <- identifyClusters(object.list, mode = "separate",resolution = 1)

features.integration = identifyIntegrationFeatures(object.list)

object.list <- identifyConfidentCells(object.list, features.integration,quantile.cutoff = 0.5)

object.list <- identifyMarkers(object.list, test.use = "bimod")


features <- c('KRT15','KRT5','KRT14','TYMS','TOP2A','KRT1','KRT2','FLG','S100A8','CXCL10','PMEL','PTPRC','CD207','CD3D')




for (i in 1:length(object.list)) {
  object1 <- RunUMAP(object.list[[i]], reduction = 'pca', dims = 1:40)
  
  # Overlay feature plots
  gg <- FeaturePlot(object1, features = features, combine = FALSE, cols = c("lightgrey", "blue"))
  gg <- patchwork::wrap_plots(plots = gg, ncol = 4)
  cowplot::save_plot(filename = paste0("overlayMarkers_", i, "_umap.pdf"), plot = gg, base_width = 8, base_height = 11)
  
  # DimPlot for clusters
  p <- DimPlot(object1, reduction = "umap", label = TRUE)
  cowplot::save_plot(filename = paste0("data_", i, "_umap.pdf"), plot = p, base_width = 4, base_height = 3)
}

#step 5. Learn technical variation between any two datasets
structured_mat <- learnTechnicalVariation(object.list, features.integration, similarity.cutoff = 0.65)


# Get barcodes from the first dataset in the list
head(Cells(object.list[[1]]))


# Iterate through each Seurat object in the list
for (i in 1:length(object.list)) {
  sample_name <- names(object.list)[i]  # Extract the sample name
  object.list[[i]] <- RenameCells(object.list[[i]], 
                                  new.names = paste0(sample_name, "_", colnames(object.list[[i]])))
}

combined <- merge(x = object.list[[1]], y = object.list[-1])
view(combined@meta.data)

#create a sample column
combined$sample <- rownames(combined@meta.data)

#split sample column
combined@meta.data <- separate(combined@meta.data, col = 'sample', into = c('Patient', 'Barcode'), 
                                    sep = '_')

combined$condition <- ifelse(grepl("-1$", combined@meta.data$Barcode),"Non-lesional", "lesional")

unique(combined@meta.data$Patient)
unique(combined@meta.data$condition)
view(combined@meta.data)

DefaultAssay(combined) <- "RNA"

#Join layers to make all data accessible
combined <- JoinLayers(combined, overwrite = TRUE)

VariableFeatures(combined) <- features.integration
combined <- integrateData(combined, structured_mat, lambda = 1)

########### Part III: Run the standard workflow for visualization and clustering ###########
nPC = 40
combined <- FindNeighbors(combined, reduction = "scMC", dims = 1:nPC)
combined <- FindClusters(combined, algorithm = 4, resolution = 1)
levels(Idents(combined))
combined <- RunUMAP(combined, reduction='scMC', dims = 1:nPC)

combined <- ScaleData(combined, feature = rownames(combined), verbose = FALSE)


## Visualization
DimPlot(combined, reduction = "umap", label = TRUE)
DimPlot(combined, reduction = "umap", split.by = "condition", combine = T)

pdf("combined_umap.pdf",height = 8,width= 12)
p1 <- DimPlot(combined, reduction = "umap", group.by = "Patient")
p2 <-DimPlot(combined, reduction = "umap", label = TRUE)
gg <- patchwork::wrap_plots(p1, p2, ncol = 2)
gg
dev.off()

table(combined@meta.data$seurat_clusters, combined@meta.data$Patient )


features = c('CXCL10','PMEL','MLANA','CD74','IL32','CD3D',"DCT", 
             "CSF1R",'CD3E', "CSF1R", "CD207", "CLEC10A")


gg <- FeaturePlot(combined, features = features)
gg <- patchwork::wrap_plots(plots = gg, ncol = 4)
gg

color.use <- scPalette(nlevels(Idents(combined)))
gg <- StackedVlnPlot(combined, features = features, colors.use = color.use)

pdf("Marker eexpression  per cluster.pdf", width=12, height = 9)
gg

dev.off()


gg <- StackedVlnPlot(combined, features = features, colors.use = color.use,split.by = "condition")
gg

####only T cell and MElanocyte cluster-16


combined.sub <- subset(combined, idents = c('16'))
combined.sub <- FindNeighbors(combined.sub, reduction = "scMC", dims = 1:nPC)
combined.sub <- FindClusters(combined.sub, algorithm = 4, resolution = 1)
levels(Idents(combined.sub))

DimPlot(combined.sub, reduction = "umap", label = TRUE)
m=combined.sub@meta.data

features = c('CXCL10','PMEL','MLANA','CD74','IL32',"CD207","CD3D","DCT", "CSF1R","CD3E", "CSF1R", "CD207", "CLEC10A")

gg <- StackedVlnPlot(combined.sub, features = features, colors.use = color.use)

pdf("Marker expression  per cluster_Cluster 16.pdf", width=12, height = 9)
gg

dev.off()


###Cluuster 2 and 4 from cluster 16

combined.sub2 <- subset(combined.sub, idents = c('2','4'))


levels(Idents(combined.sub2))

table(combined.sub2@meta.data$seurat_clusters, combined.sub2@meta.data$condition)

new.cluster.ids <- c("T.cells","T.cells")
names(new.cluster.ids) <- levels(combined.sub2)
combined.sub2 <- RenameIdents(combined.sub2, new.cluster.ids)

DimPlot(combined.sub2, reduction = "umap", label = TRUE)

FeaturePlot(combined.sub2, features = c("PCMTD2"), split.by = "condition")

DotPlot(combined.sub2, features = "PCMTD2",
        split.by = "condition") + RotatedAxis()

features = "PCMTD2"

VlnPlot(combined.sub2, features = features, split.by = "condition")


av=(AverageExpression(object = combined.sub2, group.by = c('identity', 'condition'))$RNA)

head(av)

i=which(rownames(av)=="PCMTD2")

av[i,]

saveRDS(combined, file = "vitiligo_final.rds")

