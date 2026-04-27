## Processing single cell data from PanglaoDB repository using Seurat
## Bulk download can be done from here:
## https://panglaodb.se/bulk.html

#######################################################################################
##                                                                                   ##
##                      Initial settings                                             ##
##                                                                                   ##
#######################################################################################

## Dependencies
library(Seurat)
library(miceadds)   ## for load.Rdata
library(Matrix)
library(tidyverse)
#library(future)

#plan("multiprocess", workers = 2)

## Setting data input
##  path to github repository directory
path2project <- '/media/ag-cherrmann/cramirez/sars-cov2/'     ## change in cluster
## directory to save seurat objects
panglaodb_seurat_dir <- '/media/ag-cherrmann/cramirez/sars-cov2/analysis/panglaodb_seurat'

#######################################################################################
##                                                                                   ##
##                     Auxiliary functions                                           ##
##                                                                                   ##
#######################################################################################

## Seurat standard workflow
st_workflow <- function(
        seurat_object,
        n_features = 3000,
        n_pca_dims = 15
){
        cat('Normalizing and finding variable features\n')
        seurat.p <- NormalizeData(seurat_object) %>%
                FindVariableFeatures(selection.method = 'vst',
                                     nfeatures = n_features)
        cat('Scaling and projection\n')
        seurat.p <- ScaleData(seurat.p, 
                              verbose = FALSE) %>% 
                RunPCA(npcs = n_pca_dims, 
                       verbose = FALSE) %>%
                RunUMAP(reduction = "pca", 
                        dims = 1:n_pca_dims) 
        return(seurat.p)
}

## define seurat object
define_seurat <- function(file.name, sample_name) {
        load.Rdata(file.name, objname = 'sparseMat')
        colnames(sparseMat) <- paste0(sample_name, '-', 
                                      colnames(sparseMat))
        sample.seu <- CreateSeuratObject(
                counts = sparseMat, 
                project = sample_name, 
                assay = 'RNA', 
                min.cells = 30, 
                min.features = 300 
        )
        return(sample.seu)
}

#######################################################################################
## loading metadata
## Metadata can be obtained by cloning the pangladb repository:
## https://github.com/oscar-franzen/PanglaoDB
## Run 1_download_data.sh in order to get the metadata.txt file in the correct dir
metadata <- read.table(paste0(path2project, 
                              'data/panglaodb/PanglaoDB/data/metadata.txt'), 
                       sep = ',', colClasses = 'character')
colnames(metadata) <- c('SRA_accession', 'SRS_accession',
                        'Tissue_origin', 'scRNA-seq_protocol',
                        'Species', 'Sequencing_instrument',
                        'Number_of_expressed_genes', 'Median_number_of_expressed_genes',
                        'Number_of_clusters_in_sample', 'is_tumor',
                        'is_primary', 'is_cell_line')


write_rds(metadata, 
         paste0(path2project, '/analysis/panglaoDB_metadata_full.rds'))

#####################################################################################
##                                                                                 ##
##                 Transform GE matrix to seurat object                            ##
##                                                                                 ##
#####################################################################################

## **1. Records with no SRS ID were removed (30)**
## **2. Only count matrices were kept (45 records given in RPKM were left)**
## matrices to Seurat 
#for (i in 1:nrow(metadata)){
#        file.name <- paste0('var/www/html/SRA/SRA.final/',
#                            metadata$SRA_accession[i],
#                            '_',
#                            metadata$SRS_accession[i], 
#                            '.sparse.RData')
#        file.name.full <- paste0(path2project,
#                                 '/data/',         ## change this in cluster
#                                 file.name)
#        cat(file.name.full, '\n')
#        if (file.exists(file.name.full)) {
#                sample.seu <- define_seurat(
#                        file.name = file.name.full,
#
#                        sample_name = paste0(metadata$SRA_accession[i],
#                                             '_',
#                                             metadata$SRS_accession[i])
#                )
#                sample.seu$'SRA_accession' <- metadata$SRA_accession[i]
#                sample.seu$'SRS_accession' <- metadata$SRS_accession[i]
#                sample.seu$'Tissue_origin' <- metadata$Tissue_origin[i]
#                sample.seu$'scRNA_seq_protocol' <- metadata$`scRNA-seq_protocol`[i]
#                sample.seu$'Sequencing_instrument' <- metadata$Sequencing_instrument[i]
#                sample.seu$'Number_of_expressed_genes' <- metadata$Number_of_expressed_genes[i]
#                
#                cat('Finishing ', file.name, 'sparse to seurat conversion\n')
#                seurat.file <- paste0(panglaodb_seurat_dir, '/',
#                                      metadata$SRA_accession[i], 
#                                      '_', metadata$SRS_accession[i], 
#                                      '_seurat.rds')
#                cat('Saving file to ', seurat.file, '\n')
#                saveRDS(sample.seu, seurat.file)
#        }
#}

cat('Subsetting human profiles\n')
species <- grepl('Homo sapiens', metadata$Species)
tissue <- grepl('lung|kidney|colon|liver|oma', tolower(metadata$Tissue_origin))
metadata.h <- metadata[species & tissue, ]
write_rds(metadata.h, 
         paste0(path2project, '/analysis/panglaoDB_metadata_selected_tissues.rds'))

#cat('Reading human datasets into list\n')
#seurat.list <- list()
#file.names <- c()
#for (i in 1:nrow(metadata.h)){
#        file.name <- paste0(metadata.h$SRA_accession[i],
#                            '_',
#                            metadata.h$SRS_accession[i], 
#                            '_seurat.rds')
#        file.names <- c(file.names, file.name)
#        file.name.full <- paste0(path2project,
#                                 'analysis/panglaodb_seurat/',         ## change in cluster
#                                 file.name)
#        if ( file.exists(file.name.full) ) {
#            cat('Reading ', file.name.full, 'seurat file\n')
#            sample.seu <- readRDS(file.name.full)
#            seurat.list <- c(seurat.list, sample.seu)
#        }
#}
#cat('Saving seurat object')
seurat.list.path <- paste0(path2project, 'analysis/seurat.list.rds')
#saveRDS(seurat.list, seurat.list.path)

#####################################################################################
##                                                                                 ##
##              Coembedding samples                                                ##
##                                                                                 ##
####################################################################################
#seurat.list <- read_rds(seurat.list.path)
#cat('Processing seurat objects\n')
#for (i in 1:length(seurat.list)) {
#        cat('Normalizing', i, 'file\n')
#        seurat.list[[i]] <- NormalizeData(seurat.list[[i]], 
#                                          verbose = T)
#        #seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]], 
#        #                                         selection.method = "vst", 
#        #                                         nfeatures = 8000, verbose = T)
#}
#cat('Saving seurat object')
seurat.list.path <- paste0(path2project, 'analysis/seurat.norm.list.rds')
#saveRDS(seurat.list, seurat.list.path)


#cat('Taking the intersection of variable features in all samples')
#var_feat <- lapply(seurat.list, VariableFeatures)
#var_fea_int <- Reduce(intersect, var_feat)

#cat('Merging seurat objects\n')
#merge.seu <- merge(x= seurat.list[[1]], y = seurat.list[2:length(seurat.list)])
merge.rds <- paste0(path2project, '/analysis/panglaodb.merge.seu.rds')
#saveRDS(merge.seu, merge.rds)


#cat('Processing merged seurat\n')
#merge.seu <- read_rds(merge.rds)
#merge.seu <- ScaleData(merge.seu, do.scale=TRUE)
#merge.seu <- FindVariableFeatures(merge.seu, 
#                                  selection.method = "vst", 
#                                  nfeatures = 1000, verbose = T)
#merge.seu <- RunPCA(merge.seu, verbose = FALSE)
#merge.seu <- RunUMAP(merge.seu, dims = 1:30, metric = 'cosine')
#merge.seu <- FindNeighbors(merge.seu, reduction = 'umap', dims = 1:2)
#merge.seu <- FindClusters(merge.seu, resolution = 0.05)
merge.p.rds <- paste0(path2project, '/analysis/panglaodb.merge.seu.processed.rds')
#saveRDS(merge.seu, merge.p.rds)
#cat('Merged processed data saved to ', merge.p.rds, '\n')

cat('Reading merged seurat file\n')
merge.seu <- read_rds(merge.p.rds)
cat('Adding abbreviations\n')
abbs <- read.csv(paste0(path2project, '/analysis/panglaodb_tissue_abbs.csv'),
                 stringsAsFactors = FALSE)
abbs <- plyr::mapvalues(x = merge.seu$Tissue_origin,
                        from = abbs$names, to = abbs$abbs)
merge.seu$"abbs" <- abbs

#cat('Subseting interesting genes\n')
interesting_genes <- readRDS(paste0(path2project, 
                                    '/analysis/interesting_genes.rds'))
interesting_genes <- gsub('_', '-', interesting_genes)
#merge.seu.sub.mtx <- merge.seu@assays$RNA@counts[interesting_genes, ]
#merge.sub.seu <- CreateSeuratObject(
#      counts=merge.seu.sub.mtx,
#      project = 'sars-cov2', 
#      assay = 'RNA', 
#      min.cells = 0, 
#      min.features = 0
#)
#merge.sub.seu$'cell_type' <- merge.seu$abbs
#write_rds(merge.sub.seu,
#          paste0(path2project, 'analysis/merge.sub.seu.rds'))

#cat('Plotting figures\n')
#jpeg(paste0(path2project, '/docs/panglao_seurat_cluster.jpg'))
#DimPlot(merge.seu, group.by = 'seurat_clusters', reduction = 'umap')
#dev.off()
#jpeg(paste0(path2project, '/docs/panglao_sraid.jpg'))
#DimPlot(merge.seu, group.by = 'SRA_accession', reduction = 'umap')
#dev.off()
#jpeg(paste0(path2project, '/docs/panglao_cell_type.jpg'))
#DimPlot(merge.seu, group.by = 'abbs', reduction = 'umap', label=TRUE) + NoLegend()
#dev.off()
#jpeg(paste0(path2project, '/docs/panglao_markers_ACE2.jpg'))
#FeaturePlot(merge.seu, 
#            features = c('ACE2-ENSG00000130234.10'),
#            reduction = 'umap') 
#dev.off()
#jpeg(paste0(path2project, '/docs/panglao_markers_TMPRSS2.jpg'))
#FeaturePlot(merge.seu, 
#            features = c('TMPRSS2-ENSG00000184012.11'),
#                         reduction = 'umap') 
#dev.off()
#jpeg(paste0(path2project, '/docs/panglao_markers_NEU1.jpg'))
#FeaturePlot(merge.seu, 
#            features = c('NEU1-ENSG00000204386.10'), 
#            reduction = 'umap') 
#dev.off()
#jpeg(paste0(path2project, '/docs/panglao_markers_NPTX1.jpg'))
#FeaturePlot(merge.seu, 
#            features = c('NPTX1-ENSG00000171246.5'), 
#            reduction = 'umap') 
#dev.off()
#jpeg(paste0(path2project, '/docs/panglao_markers_BRD2.jpg'))
#FeaturePlot(merge.seu, 
#            features = c('BRD2-ENSG00000204256.12'), 
#            reduction = 'umap') 
#dev.off()
jpeg(paste0(path2project, '/docs/panglao_markers_ATP1B1P1.jpg'))
FeaturePlot(merge.seu, 
            features = c('ATP1B1-ENSG00000143153.12'), 
            reduction = 'umap') 
dev.off()

#interesting_genes <- gsub('_', '-', interesting_genes)
for (i in interesting_genes){
        gene.n <- gsub('-.*', '', i)
        cat('Plotting to: ', paste0(path2project, '/docs/panglao_markers/', gene.n, '.jpg'), '\n')
        jpeg(paste0(path2project, '/docs/panglao_markers/', gene.n, '.jpg'))
         p <- FeaturePlot(merge.seu, features = i, reduction = 'umap')
         plot(p)
        dev.off()
}
#jpeg(paste0(path2project, '/docs/panglao_interactome_1.jpg'))
#FeaturePlot(merge.seu, features = interesting_genes[1:4], reduction = 'umap', ncol=2)
#jpeg(paste0(path2project, '/docs/panglao_interactome_2.jpg'))
#FeaturePlot(merge.seu, features = interesting_genes[5:8], reduction = 'umap', ncol=2)
#dev.off()
#jpeg(paste0(path2project, '/docs/panglao_interactome_3.jpg'))
#FeaturePlot(merge.seu, features = interesting_genes[9:12], reduction = 'umap', ncol=2)
#dev.off()
#jpeg(paste0(path2project, '/docs/panglao_interactome_4.jpg'))
#FeaturePlot(merge.seu, 
#            features = interesting_genes[13:16], 
#            reduction = 'umap', ncol=2)
#dev.off()
