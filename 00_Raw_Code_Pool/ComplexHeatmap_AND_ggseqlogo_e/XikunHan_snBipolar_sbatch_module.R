#!/usr/bin/env Rscript


f_log <- function(mes) {
  cat("\n", sprintf('%s: %s', date(), mes), "\n")
}

f_log("Start job")

library(data.table)
library(optparse)


option_list <- list(
  make_option(c("--cell_type"), type="character", default=NULL, help="cell_type"),
  make_option(c("--cohort"), type="character", default=NULL, help="cohort"),
  make_option(c("--dir_main"), type="character", default=NULL, help="dir_main"),
  make_option(c("--sub_outpath"), type="character", default=NULL, help="sub_outpath")
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)


cell_type <- opt$cell_type
cohort <- opt$cohort
sub_outpath <- opt$sub_outpath
dir_main <- opt$dir_main

setwd(sub_outpath)


# start ####

if(1) {
  
  library(data.table)
  library(optparse)
  
  # single-cell analysis package
  library(SingleCellExperiment)
  library(Seurat)
  library(SeuratObject)
  
  library(coin)
  library(UCell)
  library(qlcMatrix)
  library(rstatix)
  library(corrplot)
  
  
  # plotting and data science packages
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(igraph)
  library(ggpubr)
  
  
  library(readr)
  
  # co-expression network analysis packages:
  library(WGCNA)
  library(hdWGCNA) # devtools::install_github('smorabit/hdWGCNA', ref='dev')
  
  
  library(ggrepel)
  library(ggforestplot)
  
  
  library(enrichR)
  library(GeneOverlap)
  library(SuperExactTest)
  
  library(ggseqlogo)
  
  
  # devtools::install_github("XikunHan/geneticToolBox")
  library(geneticToolBox)
  
  
  # single-cell analysis package
  library(Seurat)
  
  # plotting and data science packages
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  
  # co-expression network analysis packages:
  library(WGCNA)
  library(hdWGCNA)
  
  # network analysis & visualization package:
  library(igraph)
  
  # packages for TF motif analysis
  library(JASPAR2020)
  library(motifmatchr)
  library(TFBSTools)
  library(EnsDb.Hsapiens.v86)
  library(GenomicRanges)
  
  library(ComplexHeatmap)
  
  # BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
  library(BSgenome.Hsapiens.UCSC.hg38)
  
}


theme_set(theme_cowplot())
set.seed(12345)

# run

if(1) {

if(1) {
  
  
  if(1) {
    sce <- readr::read_rds("./data/sce_snBD_McLean_MtSinai.rds")
    table(sce$Cohort)
    sce <- sce[, sce$Cohort == cohort]
    sce
  }


  set.seed(123)


  seurat_obj <- CreateSeuratObject(counts = counts(sce), min.cells = 3, min.genes = 500,
    meta.data = as.data.frame(sce@colData),
    project = "BD")

  seurat_obj


  seurat_obj$cell_group <- ifelse(grepl("^In", seurat_obj$Celltype_inferred), "In",
    ifelse(grepl("^Ex", seurat_obj$Celltype_inferred), "Ex", seurat_obj$Celltype_inferred))

  table(seurat_obj$Celltype_inferred, seurat_obj$cell_group)

  rm(sce)
  gc()

  # names(seurat_obj)
  table(grepl("^RPL|^PRS", rownames(seurat_obj)))
  seurat_obj <- seurat_obj[!grepl("^RPL|^PRS", rownames(seurat_obj)), ]
  seurat_obj


  f_log("Run normlization and UMAP")
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, features=VariableFeatures(seurat_obj))
  seurat_obj <- RunUMAP(seurat_obj, dims=1:20)
  seurat_obj <- FindNeighbors(seurat_obj)
  seurat_obj <- FindClusters(seurat_obj)


  write_rds(seurat_obj, file='hdWGCNA_process_1.rds')



  DimPlot(seurat_obj, group.by = "Celltype_inferred")
  DimPlot(seurat_obj, group.by = "cell_group")
  DimPlot(seurat_obj, group.by = "SubID")


  f_log("Run SetupForWGCNA")
  
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction", 
    fraction = 0.05,
    wgcna_name = "BD"
  )

  head(seurat_obj@meta.data)

  f_log("Run MetacellsByGroups")
  
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("SubID", "cell_group"), 
    k = 25, 
    ident.group = 'cell_group' 
  )

  # normalize metacell expression matrix:
  seurat_obj <- NormalizeMetacells(seurat_obj)

  metacell_obj <- GetMetacellObject(seurat_obj)

  write_rds(seurat_obj, file='hdWGCNA_before-RunUMAPMetacells.rds')

  
}



# run
if(1) {
  
  f_log("Run normalize Metacells")
  
  seurat_obj$SubID <- as.factor(seurat_obj$SubID)
  seurat_obj <- NormalizeMetacells(seurat_obj)  
  seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
  seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj)) 
  seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='SubID')
  seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:15)
  
  
  
  seurat_obj
  
  
  p1 <- DimPlotMetacells(seurat_obj, group.by='cell_group') + umap_theme() + ggtitle("Cell Type")
  
  
  p2 <- DimPlotMetacells(seurat_obj, group.by='SubID') + umap_theme() + ggtitle("SubID")
  
  p1 | p2
  
  p1
  p2
  
  write_rds(seurat_obj, file='hdWGCNA_after-RunUMAPMetacells.rds')
}




# run
if(1) {
  #Co-expression of specific celltype

  f_log("Run ExIn SetDatExpr")
  

    seurat_obj <- SetDatExpr(
      seurat_obj,
      group_name = cell_type, # the name of the group of interest in the group.by column
      group.by='cell_group' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    )
    



## Select soft-power threshold
  f_log("Run TestSoftPowers")
  
seurat_obj <- TestSoftPowers(
  seurat_obj,
  setDatExpr = FALSE, 
)


plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)

seurat_obj <- ConstructNetwork(
  seurat_obj, 
  setDatExpr=FALSE,
  overwrite_tom = TRUE
)

table(seurat_obj$Phenotype)


seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="SubID"
)


hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized=FALSE)


seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_group', group_name = cell_type)


seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = cell_type
)

write_rds(seurat_obj, file= sprintf('hdWGCNA_%s_post-eigengenes.rds', cell_type))


}




# run
if(1) {
  
  f_log("Run ModuleFeaturePlot")
  
  plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features='hMEs', 
    order=TRUE 
  )
  
  plot_list
  
  modules <- GetModules(seurat_obj)
  write.csv(modules, sprintf("%s_modules.subset-first.csv", cell_type))
  
  seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = 25,
    method='UCell'
  )
  write_rds(seurat_obj, file= sprintf('hdWGCNA_%s_post-eigengenes_module.rds', cell_type))
} 
  
}


if(1) {
  
  pdf(sprintf("p_module_Dendrogram_%s_%s.pdf", cohort, cell_type), height = 10, width = 8)
  PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
  dev.off()
  
  png(sprintf("p_module_Dendrogram_%s_%s.png", cohort, cell_type), height = 10, width = 8, units = "in", res = 500)
  PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
  dev.off()
  
  
  pdf("p1.pdf", height = 10, width = 8)
  ModuleCorrelogram(seurat_obj)
  
  
  MEs <- GetMEs(seurat_obj, harmonized=TRUE)
  colnames(MEs)
  mods <- colnames(MEs); mods <- mods[mods != 'grey']
  
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
  seurat_obj@meta.data
  
  p <- DotPlot(seurat_obj, features= mods, group.by = 'cell_group')
  p <- p +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')
  
  
  group1 <- seurat_obj@meta.data %>% subset(cell_group == cell_type & Phenotype == "BD") %>% rownames
  group2 <- seurat_obj@meta.data %>% subset(cell_group == cell_type & Phenotype == "CON") %>% rownames
  
  DMEs <- FindDMEs(
    seurat_obj,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use='wilcox',
    wgcna_name='BD'
  )
  
if(1) {
  
  p<- hdWGCNA::PlotDMEsLollipop(
    seurat_obj, 
    DMEs, 
    wgcna_name='BD', 
    pvalue = "p_val_adj"
  )
  
  p
  str(p)
  
  modules <- GetModules(seurat_obj) %>% subset(module != "grey") %>% mutate(module = droplevels(module))
  modules
  mod_colors <- dplyr::select(modules, c(module, color)) %>%  distinct
  mod_colors
  rownames(mod_colors) <- mod_colors$module
  mod_colors[rownames(DMEs[order(DMEs$avg_log2FC),]), ]$color
  
  
  p <-  p +
    theme(
      axis.text.y = element_text(color= mod_colors[rownames(DMEs[order(DMEs$avg_log2FC),]), ]$color, size=14, face="bold")
    )
  p
  
  geneticToolBox::ggsave_all(sprintf("hdWGCNA_DME_%s_%s", cohort, cell_type), p, width = 8, height = 6)
}
  
  
  DMEs$p_val_adj <- ifelse(DMEs$p_val_adj <= .Machine$double.xmin , 1e-300, DMEs$p_val_adj)
  DMEs
  
  
  PlotDMEsVolcano2 <- function (seurat_obj, DMEs, plot_labels = TRUE, mod_point_size = 4, 
            label_size = 4, show_cutoff = TRUE, wgcna_name = NULL) 
  {
    DMEs <- na.omit(DMEs)
    lowest <- DMEs %>% subset(p_val_adj != 0) %>% top_n(-1, 
                                                        wt = p_val_adj) %>% .$p_val_adj
    DMEs$p_val_adj <- ifelse(DMEs$p_val_adj == 0, lowest, DMEs$p_val_adj)
    max_fc <- max(abs(DMEs$avg_log2FC))
    max_fc <- DMEs %>% subset(abs(avg_log2FC) != Inf) %>% .$avg_log2FC %>% abs %>% max
    DMEs$avg_log2FC <- ifelse(DMEs$avg_log2FC == -Inf, -1 * 
                                max_fc, DMEs$avg_log2FC)
    DMEs$avg_log2FC <- ifelse(DMEs$avg_log2FC == Inf, max_fc, 
                              DMEs$avg_log2FC)
    print(DMEs)
    modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != 
                                                               "grey") %>% mutate(module = droplevels(module))
    module_colors <- modules %>% dplyr::select(c(module, color)) %>% 
      distinct
    mods <- levels(modules$module)
    mods <- mods[mods %in% DMEs$module]
    mod_colors <- module_colors$color
    names(mod_colors) <- as.character(module_colors$module)
    DMEs$anno <- ifelse(DMEs$p_val_adj < 0.05, DMEs$module, 
                        paste0(DMEs$module, " (ns)"))
    xmax <- max_fc
    print(DMEs)
    print(xmax)
    p <- DMEs %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), 
                             fill = module, color = module))
    if (show_cutoff) {
      p <- p + geom_vline(xintercept = 0, linetype = "dashed", 
                          color = "grey75", alpha = 0.8) + geom_rect(data = DMEs[1, 
                          ], aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05)), 
                          fill = "grey75", alpha = 0.8, color = NA)
    }
    p <- p + geom_point(size = mod_point_size, pch = 21, color = "black")
    if (plot_labels) {
      p <- p + geom_text_repel(aes(label = anno), color = "black", 
                               min.segment.length = 0, max.overlaps = Inf, size = label_size)
    }
    p <- p + scale_fill_manual(values = mod_colors) + scale_color_manual(values = mod_colors) + 
      xlim((-1 * xmax) - 0.1, xmax + 0.1) + xlab(bquote("Average log"[2] ~ 
                                                          "(Fold Change)")) + ylab(bquote("-log"[10] ~ "(Adj. P-value)")) + 
      theme(panel.border = element_rect(color = "black", fill = NA, 
                                        size = 1), panel.grid.major = element_blank(), axis.line = element_blank(), 
            plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
      NoLegend()
    p + NoLegend()
  }
  
  p <- PlotDMEsVolcano2(
    seurat_obj,
    DMEs,
    wgcna_name = 'BD'
  )
  
  p
  
  
  mods_sort <- paste0(cell_type, sort(as.integer(stringr::str_remove_all(mods, cell_type))))
  mods_sort
   
  p <- DotPlot(seurat_obj, features= mods_sort, group.by = 'Phenotype')
  print(p)
  
  p <- p +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')
  print(p)
  dev.off()
  
  ModuleNetworkPlot2 <- function(seurat_obj, mods = "all", outdir = "ModuleNetworks", 
            plot_size = c(6, 6), wgcna_name = NULL, label_center = FALSE, 
            edge.alpha = 0.25, vertex.label.cex = 1, vertex.size = 6, n_genes = 25, n_hubs = 25,
            n_conns = 500,
            ...) 
  {
    if (is.null(wgcna_name)) {
      wgcna_name <- seurat_obj@misc$active_wgcna
    }
    MEs <- GetMEs(seurat_obj, wgcna_name)
    modules <- GetModules(seurat_obj, wgcna_name)
    if (mods == "all") {
      mods <- levels(modules$module)
      mods <- mods[mods != "grey"]
    }
    if (!all(paste0("kME_", as.character(mods)) %in% colnames(modules))) {
      stop("Eigengene-based connectivity (kME) not found. Did you run ModuleEigengenes and ModuleConnectivity?")
    }
    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }
    cat(paste0("Writing output files to ", outdir))
    TOM <- GetTOM(seurat_obj, wgcna_name)
    hub_list <- lapply(mods, function(cur_mod) {
      cur <- subset(modules, module == cur_mod)
      cur <- cur[, c("gene_name", paste0("kME_", cur_mod))] %>% 
        top_n(n_hubs)
      colnames(cur)[2] <- "var"
      cur %>% arrange(desc(var)) %>% .$gene_name
    })
    names(hub_list) <- mods
    for (cur_mod in mods) {
      print(cur_mod)
      cur_color <- modules %>% subset(module == cur_mod) %>% 
        .$color %>% unique
      
      cur_kME <- paste0("kME_", cur_mod)
      cur_genes <- hub_list[[cur_mod]]
      matchind <- match(cur_genes, colnames(TOM))
      reducedTOM = TOM[matchind, matchind]
      orderind <- order(reducedTOM, decreasing = TRUE)
      connections2keep <- orderind[1:n_conns]
      reducedTOM <- matrix(0, nrow(reducedTOM), ncol(reducedTOM))
      reducedTOM[connections2keep] <- 1
      if (label_center) {
        cur_genes[11:25] <- ""
      }
      gA <- graph.adjacency(as.matrix(reducedTOM[1:10, 1:10]), 
                            mode = "undirected", weighted = TRUE, diag = FALSE)
      gB <- graph.adjacency(as.matrix(reducedTOM[11:n_genes, 
                                                 11:n_genes]), mode = "undirected", weighted = TRUE, 
                            diag = FALSE)
      layoutCircle <- rbind(layout.circle(gA)/2, layout.circle(gB))
      g1 <- graph.adjacency(as.matrix(reducedTOM), mode = "undirected", 
                            weighted = TRUE, diag = FALSE)
      pdf(paste0(outdir, "/", cur_mod, ".pdf"), width = plot_size[1], 
          height = plot_size[2], useDingbats = FALSE)
      plot(g1, edge.color = adjustcolor(cur_color, alpha.f = 0.25), 
           edge.alpha = edge.alpha, vertex.color = cur_color, 
           vertex.label = as.character(cur_genes), vertex.label.dist = 1.1, 
           vertex.label.degree = -pi/4, vertex.label.color = "black", 
           vertex.label.family = "Helvetica", vertex.label.font = 3, 
           vertex.label.cex = vertex.label.cex, vertex.frame.color = "black", 
           layout = jitter(layoutCircle), vertex.size = vertex.size, 
           main = paste(cur_mod))
      dev.off()
    }
  }
  
  
  ModuleNetworkPlot2(seurat_obj, n_conns = 50)
  
  graphics.off()
  
  modules <- GetModules(seurat_obj)
  
  fwrite(modules, sprintf("df_GetModule_BD_exp_%s_%s.csv", cohort, cell_type))
  
  if(0) {
    modules <- fread(sprintf("df_GetModule_BD_exp_%s_%s.csv", cohort, cell_type))
    modules <- as.data.frame(modules)
    str(modules)
    
    }
  
  df_m <- modules[modules$module != "grey", ]
  df_m$kME_grey <- NULL
  df_m
  
  skimr::skim(df_m)
  
  module_heatmap <- function(df_m, n_hubs = 4) {
    
    rownames(df_m) <- df_m$gene_name
    mods <- as.character(unique(df_m$module))
    mods

    hub_list <- lapply(mods, function(cur_mod){
      cur <- subset(df_m, module == cur_mod)
      cur <- cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
        top_n(n_hubs)
      colnames(cur)[2] <- 'var'
      cur %>% arrange(desc(var)) %>% .$gene_name
    })
    names(hub_list) <- mods
    
    
    
    df_s <- df_m[df_m$gene_name %in% unlist(hub_list), ]
    
    
    library(ComplexHeatmap)
    
    df_s
    df_s_value <- df_s[,grep("kME", names(df_s), value = T) ]
if(0) {
  ht1 = Heatmap(df_s_value ,  cluster_rows = T, cluster_columns = F,split = factor(df_s$module),
                column_names_gp = gpar(fontsize = 9))
  
  
}
    
    col <- df_s$color
    names(col) <- df_s$module
    col
    
    
    df_col <- df_s[!duplicated(df_s$module), ]
    str( df_col)
    df_col <- df_col[order(df_col$module), ]
    col <- df_col$color
    names(col) <- df_col$module
    col
    
    
    ha = rowAnnotation(Module = df_s$module, 
                       col = list(Module = col))
    
    
    ht1 = Heatmap(df_s_value , name = "Hub Genes", cluster_rows = T, cluster_columns = T,
                  right_annotation = ha, 
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 9, col = col))
    ht1
    }
 
  ht <- module_heatmap(df_m)
  
  ht
  
  
  
  pdf(sprintf("p_GetModule_BD_exp_%s_%s.pdf", cohort, cell_type), width = 8, height = 10)
  ht
  dev.off()
  
  pdf("p2.pdf", height = 10, width = 8)
  
  HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3, n_other=5,
    edge_prop = 0.75,
    mods = 'all'
  )
  
  
  g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)
  g
  
  seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs = 10, 
    n_neighbors=15,
    min_dist=0.3,
    spread=2
  )
  
  umap_df <- GetModuleUMAP(seurat_obj)
  
  ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
    geom_point(
      color=umap_df$color, 
      size=umap_df$kME*2 
    ) +
    umap_theme()
  
  graphics.off()
  
  
  
  hub_genes <- GetHubGenes(seurat_obj, 100)
  df_hubs <- subset(hub_genes, gene_name %in% df_gwas$gene_name)
  df_hubs
  
  
  label_genes <- df_hubs$gene_name
  label_genes 
  
  
  
  
  pdf(sprintf("hubgene_umap_igraph_BD_exp_%s_%s.pdf", cohort, cell_type), width=10, height=10)
  
  p <- ModuleUMAPPlot(
    seurat_obj,
    edge.alpha=0.2,
    sample_edges=TRUE,
    edge_prop=0.1, 
    label_hubs=0,
    keep_grey_edges=FALSE,
    vertex.label.cex = 0.5,
    label_genes = label_genes
  )
  p
  
  dev.off()
  
  
  setDT(df_hubs)[, list(module, gene_name)][, paste0(.SD,collapse=", "), by = "module"]
  
  
  library(GeneOverlap)
  
  
  modules <- GetModules(seurat_obj)
  table( modules$module)
  mods <- levels(modules$module)
  genome.size <- nrow(modules)
  genome.size 
  
  overlap_df <- do.call(rbind, lapply(mods, function(cur_mod){
    
    cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name
    
    cur_overlap <- testGeneOverlap(newGeneOverlap(
      cur_genes,
      df_gwas$gene_name,
      genome.size=genome.size
    ))
    
    cur_overlap <- data.frame(
      'odds.ratio' = cur_overlap@odds.ratio,
      'pval' = cur_overlap@pval,
      'Jaccard' = cur_overlap@Jaccard,
      'size_intersection' = length(cur_overlap@intersection),
      'module' = cur_mod
    )
    
    cur_overlap
    
  })) %>% as.data.frame()
  
  
  
  
  overlap_df <- overlap_df %>% mutate(fdr=p.adjust(pval, method='fdr'))
  overlap_df <- overlap_df %>% subset(module != 'grey')
  
  overlap_df$shape <- ifelse(overlap_df$fdr < 0.05, 21, 4)
  overlap_df <- overlap_df %>% arrange(odds.ratio, descending=TRUE)
  overlap_df$module <- factor(as.character(overlap_df$module), levels=as.character(overlap_df$module))
  
  mod_colors <- dplyr::select(modules, c(module, color)) %>%
    distinct
  cp <- mod_colors$color; names(cp) <- mod_colors$module
  
  p <- overlap_df %>%
    ggplot(aes(y=module, x=odds.ratio, size= size_intersection, color=module)) +
    geom_segment(aes(y=module, yend=module, x=0, xend=odds.ratio), size=0.5, color='grey') +
    geom_point() +
    geom_point(shape=overlap_df$shape, color='black', fill=NA) +
    scale_color_manual(values=cp, guide='none') +
    ylab('') + xlab("Odds ratio") +
    scale_x_continuous(breaks = c(0, 1, 2,3)) +
    labs(size='Size\nintersection') +
    ggtitle('Overlap with GWAS genes') +
    theme(
      panel.border = element_rect(size=1, color='black', fill=NA),
      axis.line.y = element_blank(),
      axis.line.x = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain')
    )
  p
  
  pdf(paste0('GWAS_overlap2.pdf'), width=8, height=8)
  p
  dev.off()
  
  
  
  
  mod_scores <-  GetModuleScores(seurat_obj)
  mod_scores
  
  seurat_obj@meta.data <- cbind(
    seurat_obj@meta.data,
    mod_scores
  )
  
  seurat_obj@meta.data
  
  # plot with Seurat's DotPlot function
  p <- DotPlot(
    seurat_obj,
    features = colnames(mod_scores),
    group.by = 'cell_group'
  )
  
  # flip the x/y axes, rotate the axis labels, and change color scheme:
  p <- p +
    RotatedAxis() +
    scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
    xlab('') + ylab('')
  
  
  p
  
  
  
  
  p <- DotPlot(
    seurat_obj,
    features = colnames(mod_scores),
    group.by = 'Phenotype'
  )
  
  
  # flip the x/y axes, rotate the axis labels, and change color scheme:
  p <- p +
    RotatedAxis() +
    scale_color_gradient2(high='darkorchid4', mid='grey95', low='seagreen') +
    xlab('') + ylab('')
  p
  
  
  # motif analysis  ####
  
  # single-cell analysis package
  library(Seurat)
  
  # plotting and data science packages
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  
  # co-expression network analysis packages:
  library(WGCNA)
  library(hdWGCNA)
  
  # network analysis & visualization package:
  library(igraph)
  
  # packages for TF motif analysis
  library(JASPAR2020)
  library(motifmatchr)
  library(TFBSTools)
  library(EnsDb.Hsapiens.v86)
  library(GenomicRanges)
  
  theme_set(theme_cowplot())
  set.seed(12345)
  
  
  
  pfm_core <- TFBSTools::getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )
  pfm_core
  
  
  library(BSgenome.Hsapiens.UCSC.hg38)
  
  MotifScan2 <- function (seurat_obj, species_genome, pfm, EnsDb, wgcna_name = NULL) 
  {
    if (is.null(wgcna_name)) {
      wgcna_name <- seurat_obj@misc$active_wgcna
    }
    motif_df <- data.frame(motif_name = purrr::map(1:length(pfm), 
                                                   function(i) {
                                                     pfm[[i]]@name
                                                   }) %>% unlist, motif_ID = purrr::map(1:length(pfm), 
                                                                                        function(i) {
                                                                                          pfm[[i]]@ID
                                                                                        }) %>% unlist)
    
    gene.promoters <- ensembldb::promoters(EnsDb, filter = ~gene_biotype == 
                                             "protein_coding") %>% subset(seqnames %in% c(1:100))
    gene.coords <- ensembldb::genes(EnsDb, filter = ~gene_biotype == 
                                      "protein_coding") %>% subset(seqnames %in% c(1:100))
    gene.promoters$symbol <- gene.coords$symbol[match(gene.promoters$gene_id, 
                                                      names(gene.coords))]
    gene.promoters <- keepSeqlevels(gene.promoters, value = levels(droplevels(seqnames(gene.promoters))))
    old_levels <- levels(seqnames(gene.promoters))
    new_levels <- ifelse(old_levels %in% c("X", "Y"), old_levels, 
                         paste0("chr", old_levels))
    gene.promoters <- renameSeqlevels(gene.promoters, new_levels)
    genome(seqinfo(gene.promoters)) <- species_genome
    my_promoters <- GRanges(seqnames = droplevels(seqnames(gene.promoters)), 
                            IRanges(start = start(gene.promoters), end = end(gene.promoters)), 
                            symbol = gene.promoters$symbol, genome = species_genome)
    print("Matching motifs...")
    motif_ix <- motifmatchr::matchMotifs(pfm, my_promoters, 
                                         genome = species_genome)
    tf_match <- motifmatchr::motifMatches(motif_ix)
    rownames(tf_match) <- my_promoters$symbol
    colnames(tf_match) <- motif_df$motif_name
    gene_list <- rownames(seurat_obj)
    gene_list <- gene_list[gene_list %in% rownames(tf_match)]
    tf_match <- tf_match[gene_list, ]
    print("Getting TF target genes...")
    tfs <- motif_df$motif_name
    tf_targets <- list()
    n_targets <- list()
    for (cur_tf in tfs) {
      tf_targets[[cur_tf]] <- names(tf_match[, cur_tf][tf_match[, 
                                                                cur_tf]])
      n_targets[[cur_tf]] <- length(tf_targets[[cur_tf]])
    }
    n_targets <- unlist(n_targets)
    motif_df$n_targets <- n_targets
    seurat_obj <- SetMotifMatrix(seurat_obj, tf_match)
    seurat_obj <- SetMotifs(seurat_obj, motif_df)
    seurat_obj <- SetMotifTargets(seurat_obj, tf_targets)
    seurat_obj <- SetPFMList(seurat_obj, pfm)
    seurat_obj
  }
  
  # run the motif scan with these settings for the mouse dataset
  seurat_obj <- MotifScan2(
    seurat_obj,
    species_genome = 'hg38',
    pfm = pfm_core,
    EnsDb = EnsDb.Hsapiens.v86
  )
  
  
  target_genes <- GetMotifTargets(seurat_obj)
  length(target_genes) 
  str(target_genes)
  target_genes
  
}



