## Analysis of HEI-193 snARC-seq combined
library(Signac)
library(Seurat)
library(plyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(ggseqlogo)
library(ggforce)
library(biovizBase)
library(qlcMatrix)
library(ggplot2)
library(ggthemes)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(future)
library(pheatmap)
library(rGREAT)
library(DESeq2)
library(MAST)
library(cicero)
library(SeuratWrappers)
library(monocle3)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(genomation)
library(reticulate)
library(SCpubr)

setwd('/raleighlab/data1/liuj/schwannoma/snarc-seq/analysis/snarc_rna_3')

# load prior snarc data from snarc_rna_2
load("/raleighlab/data1/liuj/schwannoma/snarc-seq/analysis/scnarc_rna_2.RData")

# Retain only rho screen hits
pertubationsToKeep <- c(screenTbl[screenTbl$thresh.rho > 4,"gene"],"non-targeting")
pertubationsToKeep <- pertubationsToKeep[grep("NF2",invert=TRUE,pertubationsToKeep)]

dat.keep <- subset(dat.subNo12Cicero, idents=as.character(unique(Idents(dat.subNo12Cicero)))[grep(paste(pertubationsToKeep,collapse="|",sep=""),as.character(unique(Idents(dat.subNo12Cicero))))])

# generate perturbation signal from RNA only to start 
##############################################
# Normalize to Average of 0Gy both reps only #
##############################################
DefaultAssay(dat.keep) <- "RNA"
dat.keep <- NormalizeData(dat.keep, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE) # CPM normalize

# generate pseudobulk NT vector for all
ntcells = subset(dat.keep, idents=('non-targeting_NoRT'))
df = GetAssayData(ntcells, slot = 'data')
ntVector<- apply(data.frame(df),1,mean)

df = GetAssayData(dat.keep, slot = 'data')
dfT0GyNorm = round((data.frame(df)+0.01)/(ntVector+0.01),5)

# Normalize to sgNT WITHIN conditions 


#######################################################################
# get scores for signatures from single cell HEI-193 RNA-seq clusters #
top50 <- read.table('/raleighlab/data1/liuj/schwannoma/scrna_cells/analysis/plots/top50markers.txt',header=TRUE,sep='\t') # Tim's scRNA data only
#top50 <- read.table('/raleighlab/data1/liuj/schwannoma/scrna_cells/analysis/combined/vs_perturb_scRNA_cells_top50markers_harmony_integration_res0.4.txt',header=TRUE,sep='\t') # integrated with perturb-seq
#load('/raleighlab/data1/liuj/schwannoma/scrna/analysis_scrna/top50markers_20221107.Rds') # human schwannoma cell states

# use sgRNA level or gene level with radiation condition
sigRes <- data.frame(matrix(data = 0,nrow = length(as.character(unique(dat.keep$sgRNACond))), ncol = length(unique(top50$cluster))))
row.names(sigRes) <-as.character(unique(dat.keep$sgRNACond))

#colnames(sigRes) <- c("Immunogenic schwannoma", "Myelin remodeling schwannoma", "Microglia", "Vascular remodeling schwannoma ", "T cells", "Neural regeneration schwannoma",
#                       "Axon injury schwannoma", "Vascular endothelia", "Macrophages", "Complement activation schwannoma", "Stromal cells", "Progenitor schwannoma", "Vascular smooth muscle")
colnames(sigRes) <- c("Fe Metabolism","S phase","Gluthathione","EMT","Cell Stress (ATF3, TK1)","G2/M","TGFb","Cell Stress (ATF3, SGK1)","RNA splicing","Mitochondria",
                      "Ribosome","Cell membrane","p53 pathway","Interferon","Ribosome PLP2+")

for (i in 1:nrow(sigRes)){
  for (j in 1:ncol(sigRes)){
    #goi = pull(top50[top50$cluster == j-1,"gene"][1:10,]) # top n marker genes for clinical samples
    goi = top50[top50$cluster == j-1,"gene"][1:10]
    goi = goi[grep("AC092069.1",invert = TRUE,goi)]
    cells = row.names(sigRes)[i]
    cells = colnames(subset(dat.keep,idents=cells))
    cells = gsub("-",".",cells)
    
   #pseudobulk = apply(dfT.gemnorm[,row.names(sgRNATblConcordant[sgRNATblConcordant$geneCond == cells,])],1,mean) # pseudobulk using mean across all sgRNA conds 
    pseudobulk = apply(dfT0GyNorm[,cells],1,mean) # pseudobulk using mean across all sgRNA conds, using 0 Gy normalized NTs 
    score = mean(pseudobulk[goi],na.rm = TRUE) # signature using mean 

    sigRes[i,j] <- score
  }
}

sigRes <- sigRes[complete.cases(sigRes),]
sigRes <- sigRes[order(row.names(sigRes)),]

#sigResClinical <- sigRes[,grep("schwannoma",colnames(sigRes))] # pull out schwannoma cell states only from clinical data 

sigRes <- cbind(sigRes,sigResClinical)

sgRNAClassRT <- data.frame(row.names = row.names(sigRes), genecond = row.names(sigRes), cond = str_split_fixed(row.names(sigRes), "_", 2)[,2])
statesClass <- data.frame(row.names = colnames(sigRes), cluster = colnames(sigRes),type = c(rep("states",15),rep("types",7)))


#####################
# cluster by groups #

row.names(metaDataChromVar) <- metaDataChromVar$id
metaDataChromVar2 <- metaDataChromVar
row.names(metaDataChromVar2) <- gsub("NoRT","5x1.8",row.names(metaDataChromVar2))
metaDataChromVar <- rbind(metaDataChromVar,metaDataChromVar2)
metaDataChromVar <- metaDataChromVar[row.names(sigRes),]

pdf("vs_snarc_RNA_modules_res03_and_clinical_0GyNorm_split.pdf",width = 16, height=7)

Heatmap(as.matrix(t(log(sigRes[complete.cases(sigRes),],2))),
        top_annotation = HeatmapAnnotation(cond = sgRNAClassRT$cond,
                                           GO = metaDataChromVar$GO, 
                                           col = list(cond = c("NoRT" = "#1B9E77", "5x1.8" = "#D95F02"),
                                           GO = c("Histone demethylase" = "#9E0142", "Histone acetyltransferase" = "#D53E4F",
                                                                                                "Histone deacetylase" = "#F46D43", "SAGA complex" = "#FDAE61",
                                                                                                "HBO1 complex" = "#FEE08B", "Histone methyltransferase" = "#FFFFBF",
                                                                                                "NuA4 complex" = "#E6F598", "NuRD complex" = "#ABDDA4","PRC2 complex" = "#66C2A5",
                                                                                                "SWI SNF" = "#3288BD", "Glycosylase" = "#5E4FA2", "non-targeting" = "gray50")),
                                          #gamma = anno_barplot(metaDataChromVar$gamma, bar_width = 0.9,ylim = c(-1,1)), 
                                          rho = anno_barplot(metaDataChromVar$rho, bar_width = 0.9)),
        column_split = factor(sgRNAClassRT$cond, levels = c("NoRT","5x1.8")),
        row_split = statesClass$type,
        cluster_column_slice = FALSE,
        width = ncol(t(sigRes))*unit(4, "mm"),
        height = nrow(t(sigRes))*unit(4, "mm"),
        col = colorRamp2(seq(-3, 3, length = 3), c("#2268AD", "#EEEEEE", "#B31B2C"))
)
dev.off()
 
#########################################################################
# Analysis of gene activity scores separated by perturbation/conditions #
# compute gene activities
DefaultAssay(dat.keep) <- "ATAC"
gene.activities <- GeneActivity(dat.keep)
dat.keep[['gene.activities']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(dat.keep) <- "gene.activities"
dat.keep <- NormalizeData(dat.keep, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE) # CPM normalize
gene.activities = GetAssayData(dat.keep, slot = 'data')

# generate pseudobulk NT vector for all gene activities 
ntVectorGA<- apply(gene.activities[,colnames(ntcells)],1,mean)

GAT0GyNorm = round((data.frame(gene.activities)+0.01)/(ntVectorGA+0.01),5)

#######################################################################
# get scores for signatures from single cell HEI-193 RNA-seq clusters #
top50 <- read.table('/raleighlab/data1/liuj/schwannoma/scrna_cells/analysis/plots/top50markers.txt',header=TRUE,sep='\t') # Tim's scRNA data only
#load('/raleighlab/data1/liuj/schwannoma/scrna/analysis_scrna/top50markers_20221107.Rds') # human schwannoma cell states

# use sgRNA level or gene level with radiation condition
sigResGA <- data.frame(matrix(data = 0,nrow = length(as.character(unique(dat.keep$sgRNACond))), ncol = length(unique(top50$cluster))))
row.names(sigResGA) <-as.character(unique(dat.keep$sgRNACond))

#colnames(sigResGA) <- c("Immunogenic schwannoma", "Myelin remodeling schwannoma", "Microglia", "Vascular remodeling schwannoma ", "T cells", "Neural regeneration schwannoma",
#                       "Axon injury schwannoma", "Vascular endothelia", "Macrophages", "Complement activation schwannoma", "Stromal cells", "Progenitor schwannoma", "Vascular smooth muscle")
colnames(sigResGA) <- c("Fe Metabolism","S phase","Gluthathione","EMT","Cell Stress (ATF3, TK1)","G2/M","TGFb","Cell Stress (ATF3, SGK1)","RNA splicing","Mitochondria",
                     "Ribosome","Cell membrane","p53 pathway","Interferon","Ribosome PLP2+")

for (i in 1:nrow(sigResGA)){
  for (j in 1:ncol(sigResGA)){
    #goi = pull(top50[top50$cluster == j-1,"gene"][1:10,]) # top n marker genes for clinical samples
    goi = top50[top50$cluster == j-1,"gene"][1:10]
    goi = goi[grep("AC092069.1",invert = TRUE,goi)]
    cells = row.names(sigResGA)[i]
    cells = colnames(subset(dat.keep,idents=cells))
    cells = gsub("-",".",cells)
    
    pseudobulk = apply(GAT0GyNorm[,cells],1,mean) # pseudobulk using mean across all sgRNA conds, using 0 Gy normalized NTs 
    score = mean(pseudobulk[goi],na.rm = TRUE) # signature using mean 
    
    sigResGA[i,j] <- score
  }
}

sigResGA <- sigResGA[complete.cases(sigResGA),]
sigResGA <- sigResGA[order(row.names(sigResGA)),]

#sigResGAClinical <- sigResGA[,grep("schwannoma",colnames(sigResGA))] # pull out schwannoma cell states only from clinical data 

sigResGA <- cbind(sigResGA,sigResGAClinical)

sgRNAClassRT <- data.frame(row.names = row.names(sigResGA), genecond = row.names(sigResGA), cond = str_split_fixed(row.names(sigResGA), "_", 2)[,2])
statesClass <- data.frame(row.names = colnames(sigResGA), cluster = colnames(sigResGA),type = c(rep("states",15),rep("types",7)))


#####################
# plot 

pdf("vs_snarc_GeneActivity_modules_res03_and_clinical_0GyNorm_split.pdf",width = 16, height=7)

Heatmap(as.matrix(t(log(sigResGA[complete.cases(sigResGA),],2))),
        top_annotation = HeatmapAnnotation(cond = sgRNAClassRT$cond,
                                           GO = metaDataChromVar$GO, 
                                           col = list(cond = c("NoRT" = "#1B9E77", "5x1.8" = "#D95F02"),
                                                      GO = c("Histone demethylase" = "#9E0142", "Histone acetyltransferase" = "#D53E4F",
                                                             "Histone deacetylase" = "#F46D43", "SAGA complex" = "#FDAE61",
                                                             "HBO1 complex" = "#FEE08B", "Histone methyltransferase" = "#FFFFBF",
                                                             "NuA4 complex" = "#E6F598", "NuRD complex" = "#ABDDA4","PRC2 complex" = "#66C2A5",
                                                             "SWI SNF" = "#3288BD", "Glycosylase" = "#5E4FA2", "non-targeting" = "gray50")),
                                           #gamma = anno_barplot(metaDataChromVar$gamma, bar_width = 0.9,ylim = c(-1,1)), 
                                           rho = anno_barplot(metaDataChromVar$rho, bar_width = 0.9)),
        column_split = factor(sgRNAClassRT$cond, levels = c("NoRT","5x1.8")),
        row_split = statesClass$type,
        cluster_column_slice = FALSE,
        width = ncol(t(sigResGA))*unit(4, "mm"),
        height = nrow(t(sigResGA))*unit(4, "mm"),
        col = colorRamp2(seq(-6, 6, length = 3), c("#2268AD", "#EEEEEE", "#B31B2C"))
)
dev.off()

#########################################################
# Combined snARC-RNA and gene activity scores from ATAC #
colnames(sigRes) <- paste(colnames(sigRes), "RNA",sep=' ') 
colnames(sigResGA) <- paste(colnames(sigResGA), "ATAC",sep=' ') 

sigResBoth <-cbind(sigResGA, sigRes)
statesClassBoth <- data.frame(row.names = colnames(sigResBoth), cluster = colnames(sigResBoth),type = c(rep("states - ATAC",15),rep("types - ATAC",7),rep("states - RNA",15),rep("types - RNA",7)))

pdf("vs_snarc_RNA_GeneActivity_modules_res03_and_clinical_0GyNorm_split.pdf",width = 16, height=12)

Heatmap(as.matrix(t(log(sigResBoth[complete.cases(sigResBoth),],2))),
        top_annotation = HeatmapAnnotation(cond = sgRNAClassRT$cond,
                                           GO = metaDataChromVar$GO, 
                                           col = list(cond = c("NoRT" = "#1B9E77", "5x1.8" = "#D95F02"),
                                                      GO = c("Histone demethylase" = "#9E0142", "Histone acetyltransferase" = "#D53E4F",
                                                             "Histone deacetylase" = "#F46D43", "SAGA complex" = "#FDAE61",
                                                             "HBO1 complex" = "#FEE08B", "Histone methyltransferase" = "#FFFFBF",
                                                             "NuA4 complex" = "#E6F598", "NuRD complex" = "#ABDDA4","PRC2 complex" = "#66C2A5",
                                                             "SWI SNF" = "#3288BD", "Glycosylase" = "#5E4FA2", "non-targeting" = "gray50")),
                                           #gamma = anno_barplot(metaDataChromVar$gamma, bar_width = 0.9,ylim = c(-1,1)), 
                                           rho = anno_barplot(metaDataChromVar$rho, bar_width = 0.9)),
        column_split = factor(sgRNAClassRT$cond, levels = c("NoRT","5x1.8")),
        row_split = factor(statesClassBoth$type, levels = c("states - ATAC", "types - ATAC", "states - RNA", "types - RNA")),
        cluster_column_slice = FALSE,
        width = ncol(t(sigResBoth))*unit(4, "mm"),
        height = nrow(t(sigResBoth))*unit(4, "mm"),
        col = colorRamp2(seq(-6, 6, length = 3), c("#2268AD", "#EEEEEE", "#B31B2C"))
)
dev.off()

#########################################################################################
# Normalize to sgNT WITHIN conditions - GEM normalization of RNA and ATAC gene activity #
# generate pseudobulk NT vector for all subset by radiation 
DefaultAssay(dat.keep) <- "RNA"

ntcells = subset(dat.keep, idents=('non-targeting_NoRT'))
df = GetAssayData(ntcells, slot = 'data')
ntVector<- apply(data.frame(df),1,mean)

ntcells9 = subset(dat.keep, idents=('non-targeting_5x1.8'))
df = GetAssayData(ntcells9, slot = 'data')
ntVector9 <- apply(data.frame(df),1,mean)

# subset by noRT and RT
df = subset(dat.keep, idents=as.character(unique(Idents(dat.keep)))[grep("NoRT",as.character(unique(Idents(dat.keep))))])
df = GetAssayData(df, slot = 'data')
dfTGEMNorm = round((data.frame(df)+0.01)/(ntVector+0.01),5)

df = subset(dat.keep, idents=as.character(unique(Idents(dat.keep)))[grep("5x1.8",as.character(unique(Idents(dat.keep))))])
df = GetAssayData(df, slot = 'data')
dfTGEMNorm = cbind(dfTGEMNorm,round((data.frame(df)+0.01)/(ntVector9+0.01),5))

# GEM normalize Gene Activities
# generate pseudobulk NT vectors
ntVectorGA9<- apply(gene.activities[,colnames(ntcells9)],1,mean)

df = gene.activities[,names(Idents(dat.keep)[grep("NoRT",Idents(dat.keep))])]
GATGEMNorm = round((data.frame(df)+0.01)/(ntVectorGA+0.01),5)

df = gene.activities[,names(Idents(dat.keep)[grep("5x1.8",Idents(dat.keep))])]
GATGEMNorm = cbind(GATGEMNorm,round((data.frame(df)+0.01)/(ntVectorGA9+0.01),5))

############################
# snARC-RNA GEM normalized #
#load('/raleighlab/data1/liuj/schwannoma/scrna/analysis_scrna/top50markers_20221107.Rds') # human schwannoma cell states
top50 <- read.table('/raleighlab/data1/liuj/schwannoma/scrna_cells/analysis/plots/top50markers.txt',header=TRUE,sep='\t') # Tim's scRNA data only

# use sgRNA level or gene level with radiation condition
sigResGEM <- data.frame(matrix(data = 0,nrow = length(as.character(unique(dat.keep$sgRNACond))), ncol = length(unique(top50$cluster))))
row.names(sigResGEM) <-as.character(unique(dat.keep$sgRNACond))

#colnames(sigResGEM) <- c("Immunogenic schwannoma", "Myelin remodeling schwannoma", "Microglia", "Vascular remodeling schwannoma ", "T cells", "Neural regeneration schwannoma",
#                       "Axon injury schwannoma", "Vascular endothelia", "Macrophages", "Complement activation schwannoma", "Stromal cells", "Progenitor schwannoma", "Vascular smooth muscle")
colnames(sigResGEM) <- c("Fe Metabolism","S phase","Gluthathione","EMT","Cell Stress (ATF3, TK1)","G2/M","TGFb","Cell Stress (ATF3, SGK1)","RNA splicing","Mitochondria",
                      "Ribosome","Cell membrane","p53 pathway","Interferon","Ribosome PLP2+")

for (i in 1:nrow(sigResGEM)){
  for (j in 1:ncol(sigResGEM)){
    #goi = pull(top50[top50$cluster == j-1,"gene"][1:10,]) # top n marker genes for clinical samples
    goi = top50[top50$cluster == j-1,"gene"][1:10]
    goi = goi[grep("AC092069.1",invert = TRUE,goi)]
    cells = row.names(sigResGEM)[i]
    cells = colnames(subset(dat.keep,idents=cells))
    cells = gsub("-",".",cells)
    
    pseudobulk = apply(dfTGEMNorm[,cells],1,mean) # pseudobulk using mean normalized by NTs WITHIN conditions 
    score = mean(pseudobulk[goi],na.rm = TRUE) # signature using mean 
    
    sigResGEM[i,j] <- score
  }
}

sigResGEM <- sigResGEM[complete.cases(sigResGEM),]
sigResGEM <- sigResGEM[order(row.names(sigResGEM)),]

#sigResGEMClinical <- sigResGEM[,grep("schwannoma",colnames(sigResGEM))] # pull out schwannoma cell states only from clinical data 

sigResGEM <- cbind(sigResGEM,sigResGEMClinical)

# cluster by groups #

pdf("vs_snarc_RNA_modules_res03_and_clinical_GEMNorm_split.pdf",width = 16, height=7)

Heatmap(as.matrix(t(log(sigResGEM[complete.cases(sigResGEM),],2))),
        top_annotation = HeatmapAnnotation(cond = sgRNAClassRT$cond,
                                           GO = metaDataChromVar$GO, 
                                           col = list(cond = c("NoRT" = "#1B9E77", "5x1.8" = "#D95F02"),
                                                      GO = c("Histone demethylase" = "#9E0142", "Histone acetyltransferase" = "#D53E4F",
                                                             "Histone deacetylase" = "#F46D43", "SAGA complex" = "#FDAE61",
                                                             "HBO1 complex" = "#FEE08B", "Histone methyltransferase" = "#FFFFBF",
                                                             "NuA4 complex" = "#E6F598", "NuRD complex" = "#ABDDA4","PRC2 complex" = "#66C2A5",
                                                             "SWI SNF" = "#3288BD", "Glycosylase" = "#5E4FA2", "non-targeting" = "gray50")),
                                           #gamma = anno_barplot(metaDataChromVar$gamma, bar_width = 0.9,ylim = c(-1,1)), 
                                           rho = anno_barplot(metaDataChromVar$rho, bar_width = 0.9)),
        column_split = factor(sgRNAClassRT$cond, levels = c("NoRT","5x1.8")),
        row_split = statesClass$type,
        cluster_column_slice = FALSE,
        width = ncol(t(sigResGEM))*unit(4, "mm"),
        height = nrow(t(sigResGEM))*unit(4, "mm"),
        col = colorRamp2(seq(-2, 2, length = 3), c("#2268AD", "#EEEEEE", "#B31B2C"))
)
dev.off()

################################
# Gene activity GEM normalized #
top50 <- read.table('/raleighlab/data1/liuj/schwannoma/scrna_cells/analysis/plots/top50markers.txt',header=TRUE,sep='\t') # Tim's scRNA data only
#load('/raleighlab/data1/liuj/schwannoma/scrna/analysis_scrna/top50markers_20221107.Rds') # human schwannoma cell states

# use sgRNA level or gene level with radiation condition
sigResGEMGA <- data.frame(matrix(data = 0,nrow = length(as.character(unique(dat.keep$sgRNACond))), ncol = length(unique(top50$cluster))))
row.names(sigResGEMGA) <-as.character(unique(dat.keep$sgRNACond))

#colnames(sigResGEMGA) <- c("Immunogenic schwannoma", "Myelin remodeling schwannoma", "Microglia", "Vascular remodeling schwannoma ", "T cells", "Neural regeneration schwannoma",
#                       "Axon injury schwannoma", "Vascular endothelia", "Macrophages", "Complement activation schwannoma", "Stromal cells", "Progenitor schwannoma", "Vascular smooth muscle")
colnames(sigResGEMGA) <- c("Fe Metabolism","S phase","Gluthathione","EMT","Cell Stress (ATF3, TK1)","G2/M","TGFb","Cell Stress (ATF3, SGK1)","RNA splicing","Mitochondria",
                       "Ribosome","Cell membrane","p53 pathway","Interferon","Ribosome PLP2+")

for (i in 1:nrow(sigResGEMGA)){
  for (j in 1:ncol(sigResGEMGA)){
    #goi = pull(top50[top50$cluster == j-1,"gene"][1:10,]) # top n marker genes for clinical samples
    goi = top50[top50$cluster == j-1,"gene"][1:10]
    goi = goi[grep("AC092069.1",invert = TRUE,goi)]
    cells = row.names(sigResGEMGA)[i]
    cells = colnames(subset(dat.keep,idents=cells))
    cells = gsub("-",".",cells)
    
    pseudobulk = apply(GATGEMNorm[,cells],1,mean) # pseudobulk using mean across normalized WITHIN conditions by sgNTC
    score = mean(pseudobulk[goi],na.rm = TRUE) # signature using mean 
    
    sigResGEMGA[i,j] <- score
  }
}

sigResGEMGA <- sigResGEMGA[complete.cases(sigResGEMGA),]
sigResGEMGA <- sigResGEMGA[order(row.names(sigResGEMGA)),]

#sigResGEMGAClinical <- sigResGEMGA[,grep("schwannoma",colnames(sigResGEMGA))] # pull out schwannoma cell states only from clinical data 

sigResGEMGA <- cbind(sigResGEMGA,sigResGEMGAClinical)

sgRNAClassRT <- data.frame(row.names = row.names(sigResGEMGA), genecond = row.names(sigResGEMGA), cond = str_split_fixed(row.names(sigResGEMGA), "_", 2)[,2])
statesClass <- data.frame(row.names = colnames(sigResGEMGA), cluster = colnames(sigResGEMGA),type = c(rep("states",15),rep("types",7)))

pdf("vs_snarc_GeneActivity_modules_res03_and_clinical_GEMNorm_split.pdf",width = 16, height=7)

Heatmap(as.matrix(t(log(sigResGEMGA[complete.cases(sigResGEMGA),],2))),
        top_annotation = HeatmapAnnotation(cond = sgRNAClassRT$cond,
                                           GO = metaDataChromVar$GO, 
                                           col = list(cond = c("NoRT" = "#1B9E77", "5x1.8" = "#D95F02"),
                                                      GO = c("Histone demethylase" = "#9E0142", "Histone acetyltransferase" = "#D53E4F",
                                                             "Histone deacetylase" = "#F46D43", "SAGA complex" = "#FDAE61",
                                                             "HBO1 complex" = "#FEE08B", "Histone methyltransferase" = "#FFFFBF",
                                                             "NuA4 complex" = "#E6F598", "NuRD complex" = "#ABDDA4","PRC2 complex" = "#66C2A5",
                                                             "SWI SNF" = "#3288BD", "Glycosylase" = "#5E4FA2", "non-targeting" = "gray50")),
                                           #gamma = anno_barplot(metaDataChromVar$gamma, bar_width = 0.9,ylim = c(-1,1)), 
                                           rho = anno_barplot(metaDataChromVar$rho, bar_width = 0.9)),
        column_split = factor(sgRNAClassRT$cond, levels = c("NoRT","5x1.8")),
        row_split = statesClass$type,
        cluster_column_slice = FALSE,
        width = ncol(t(sigResGEMGA))*unit(4, "mm"),
        height = nrow(t(sigResGEMGA))*unit(4, "mm"),
        col = colorRamp2(seq(-4, 4, length = 3), c("#2268AD", "#EEEEEE", "#B31B2C"))
)
dev.off()

###########################################################################
# Combined snARC-RNA and gene activity scores from ATAC  - GEM normalized #
colnames(sigResGEM) <- paste(colnames(sigResGEM), "RNA",sep=' ') 
colnames(sigResGEMGA) <- paste(colnames(sigResGEMGA), "ATAC",sep=' ') 

sigResGEMBoth <-cbind(sigResGEMGA, sigResGEM)
statesClassBoth <- data.frame(row.names = colnames(sigResGEMBoth), cluster = colnames(sigResGEMBoth),type = c(rep("states - ATAC",15),rep("types - ATAC",7),rep("states - RNA",15),rep("types - RNA",7)))

pdf("vs_snarc_RNA_GeneActivity_modules_res03_and_clinical_GEMNorm_split.pdf",width = 16, height=12)

Heatmap(as.matrix(t(log(sigResGEMBoth[complete.cases(sigResGEMBoth),],2))),
        top_annotation = HeatmapAnnotation(cond = sgRNAClassRT$cond,
                                           GO = metaDataChromVar$GO, 
                                           col = list(cond = c("NoRT" = "#1B9E77", "5x1.8" = "#D95F02"),
                                                      GO = c("Histone demethylase" = "#9E0142", "Histone acetyltransferase" = "#D53E4F",
                                                             "Histone deacetylase" = "#F46D43", "SAGA complex" = "#FDAE61",
                                                             "HBO1 complex" = "#FEE08B", "Histone methyltransferase" = "#FFFFBF",
                                                             "NuA4 complex" = "#E6F598", "NuRD complex" = "#ABDDA4","PRC2 complex" = "#66C2A5",
                                                             "SWI SNF" = "#3288BD", "Glycosylase" = "#5E4FA2", "non-targeting" = "gray50")),
                                           #gamma = anno_barplot(metaDataChromVar$gamma, bar_width = 0.9,ylim = c(-1,1)), 
                                           rho = anno_barplot(metaDataChromVar$rho, bar_width = 0.9)),
        column_split = factor(sgRNAClassRT$cond, levels = c("NoRT","5x1.8")),
        row_split = factor(statesClassBoth$type, levels = c("states - ATAC", "types - ATAC", "states - RNA", "types - RNA")),
        cluster_column_slice = FALSE,
        width = ncol(t(sigResGEMBoth))*unit(4, "mm"),
        height = nrow(t(sigResGEMBoth))*unit(4, "mm"),
        col = colorRamp2(seq(-4, 4, length = 3), c("#2268AD", "#EEEEEE", "#B31B2C"))
)
dev.off()

#################################################
# Normalize RT to noRT of the same perturbation #

df = GetAssayData(dat.keep, slot = 'data', assay = 'RNA')

# load genes and calculate the means 
top50 <- read.table('/raleighlab/data1/liuj/schwannoma/scrna_cells/analysis/plots/top50markers.txt',header=TRUE,sep='\t') # Tim's scRNA data only
#load('/raleighlab/data1/liuj/schwannoma/scrna/analysis_scrna/top50markers_20221107.Rds') # human schwannoma cell states

# use sgRNA level or gene level with radiation condition
sigResPerturbNormRNA <- data.frame(matrix(data = 0,nrow = length(as.character(unique(dat.keep$sgRNACond))), ncol = length(unique(top50$cluster))))
row.names(sigResPerturbNormRNA) <-as.character(unique(dat.keep$sgRNACond))
sigResPerturbNormRNA <- sigResPerturbNormRNA[grep("5x1.8",row.names(sigResPerturbNormRNA)),]
row.names(sigResPerturbNormRNA) <- gsub("_5x1.8","",row.names(sigResPerturbNormRNA))

#colnames(sigResPerturbNormRNA) <- c("Immunogenic schwannoma", "Myelin remodeling schwannoma", "Microglia", "Vascular remodeling schwannoma ", "T cells", "Neural regeneration schwannoma",
#                       "Axon injury schwannoma", "Vascular endothelia", "Macrophages", "Complement activation schwannoma", "Stromal cells", "Progenitor schwannoma", "Vascular smooth muscle")
colnames(sigResPerturbNormRNA) <- c("Fe Metabolism","S phase","Gluthathione","EMT","Cell Stress (ATF3, TK1)","G2/M","TGFb","Cell Stress (ATF3, SGK1)","RNA splicing","Mitochondria",
                          "Ribosome","Cell membrane","p53 pathway","Interferon","Ribosome PLP2+")

sigResPerturbNormGA <- sigResPerturbNormRNA

for (i in 1:nrow(sigResPerturbNormRNA)){
  for (j in 1:ncol(sigResPerturbNormRNA)){
    #goi = pull(top50[top50$cluster == j-1,"gene"][1:10,]) # top n marker genes for clinical samples
    goi = top50[top50$cluster == j-1,"gene"][1:10]
    goi = goi[grep("AC092069.1",invert = TRUE,goi)]
    cells1 = paste(row.names(sigResPerturbNormRNA)[i],"NoRT",sep='_')
    cells1 = colnames(subset(dat.keep,idents=cells1))
    #cells1 = gsub("-",".",cells1)
    
    cells2 = paste(row.names(sigResPerturbNormRNA)[i],"5x1.8",sep='_')
    cells2 = colnames(subset(dat.keep,idents=cells2))
    #cells2 = gsub("-",".",cells2)
    
    pseudobulk1 = apply(df[,cells1],1,mean) # RNA normalization perturbation RT/noRT
    pseudobulk2 = apply(df[,cells2],1,mean) 
    
    score = (mean(pseudobulk2[goi],na.rm = TRUE)+0.01)/(mean(pseudobulk1[goi],na.rm = TRUE)+0.01) # signature using mean 
    sigResPerturbNormRNA[i,j] <- score
    
    pseudobulk1 = apply(gene.activities[,cells1],1,mean) # ATAC GA normalization perturbation RT/noRT
    pseudobulk2 = apply(gene.activities[,cells2],1,mean) 
    
    score = (mean(pseudobulk2[goi],na.rm = TRUE)+0.01)/(mean(pseudobulk1[goi],na.rm = TRUE)+0.01) # signature using mean 
    sigResPerturbNormGA[i,j] <- score
    
  }
}

sigResPerturbNormRNA <- sigResPerturbNormRNA[complete.cases(sigResPerturbNormRNA),]
sigResPerturbNormRNA <- sigResPerturbNormRNA[order(row.names(sigResPerturbNormRNA)),]

sigResPerturbNormGA <- sigResPerturbNormGA[complete.cases(sigResPerturbNormGA),]
sigResPerturbNormGA <- sigResPerturbNormGA[order(row.names(sigResPerturbNormGA)),]

#sigResPerturbNormRNAClinical <- sigResPerturbNormRNA[,grep("schwannoma",colnames(sigResPerturbNormRNA))] # pull out schwannoma cell states only from clinical data 
#sigResPerturbNormGAClinical <- sigResPerturbNormGA[,grep("schwannoma",colnames(sigResPerturbNormGA))] # pull out schwannoma cell states only from clinical data 

sigResPerturbNormRNA <- cbind(sigResPerturbNormRNA,sigResPerturbNormRNAClinical)
sigResPerturbNormGA <- cbind(sigResPerturbNormGA,sigResPerturbNormGAClinical)

sgRNAClassRT <- data.frame(row.names = row.names(sigResPerturbNormRNA), genecond = row.names(sigResPerturbNormRNA))
statesClass <- data.frame(row.names = colnames(sigResPerturbNormRNA), cluster = colnames(sigResPerturbNormRNA),type = c(rep("states",15),rep("types",7)))

# combine both
colnames(sigResPerturbNormRNA) <- paste(colnames(sigResPerturbNormRNA), "RNA",sep=' ') 
colnames(sigResPerturbNormGA) <- paste(colnames(sigResPerturbNormGA), "ATAC",sep=' ') 
sigResPerturbBoth <-cbind(sigResPerturbNormGA, sigResPerturbNormRNA)
statesClassBoth <- data.frame(row.names = colnames(sigResPerturbBoth), cluster = colnames(sigResPerturbBoth),type = c(rep("states - ATAC",15),rep("types - ATAC",7),rep("states - RNA",15),rep("types - RNA",7)))

# new metadata for perturbnorm data
metaDataChromVarNoCond <- metaDataChromVar
metaDataChromVarNoCond <- metaDataChromVarNoCond[grep("5x1.8",row.names(metaDataChromVar)),]
row.names(metaDataChromVarNoCond) <- gsub("_5x1.8","",row.names(metaDataChromVarNoCond))
  
pdf("vs_snarc_RNA_GeneActivity_modules_res03_and_clinical_PerturbNorm.pdf",width = 16, height=12)
Heatmap(as.matrix(t(log(sigResPerturbBoth[complete.cases(sigResPerturbBoth),],2))),
        top_annotation = HeatmapAnnotation(GO = metaDataChromVarNoCond$GO, 
                                           col = list(cond = c("NoRT" = "#1B9E77", "5x1.8" = "#D95F02"),
                                                      GO = c("Histone demethylase" = "#9E0142", "Histone acetyltransferase" = "#D53E4F",
                                                             "Histone deacetylase" = "#F46D43", "SAGA complex" = "#FDAE61",
                                                             "HBO1 complex" = "#FEE08B", "Histone methyltransferase" = "#FFFFBF",
                                                             "NuA4 complex" = "#E6F598", "NuRD complex" = "#ABDDA4","PRC2 complex" = "#66C2A5",
                                                             "SWI SNF" = "#3288BD", "Glycosylase" = "#5E4FA2", "non-targeting" = "gray50")),
                                           gamma = anno_barplot(metaDataChromVarNoCond$gamma, bar_width = 0.9,ylim = c(-1,1)), 
                                           rho = anno_barplot(metaDataChromVarNoCond$rho, bar_width = 0.9)),
        row_split = factor(statesClassBoth$type, levels = c("states - ATAC", "types - ATAC", "states - RNA", "types - RNA")),
        cluster_column_slice = FALSE,
        width = ncol(t(sigResPerturbBoth))*unit(4, "mm"),
        height = nrow(t(sigResPerturbBoth))*unit(4, "mm"),
        col = colorRamp2(seq(-3, 3, length = 3), c("#2268AD", "#EEEEEE", "#B31B2C")),
        heatmap_legend_param = list(title = "RT vs no RT"))
dev.off()

# Plot GA and RNA separately
pdf("vs_snarc_RNA_modules_res03_and_clinical_PerturbNorm.pdf",width = 13, height=8)
Heatmap(as.matrix(t(log(sigResPerturbNormRNA[complete.cases(sigResPerturbNormRNA),],2))),
        top_annotation = HeatmapAnnotation(GO = metaDataChromVarNoCond$GO, 
                                           col = list(cond = c("NoRT" = "#1B9E77", "5x1.8" = "#D95F02"),
                                                      GO = c("Histone demethylase" = "#9E0142", "Histone acetyltransferase" = "#D53E4F",
                                                             "Histone deacetylase" = "#F46D43", "SAGA complex" = "#FDAE61",
                                                             "HBO1 complex" = "#FEE08B", "Histone methyltransferase" = "#FFFFBF",
                                                             "NuA4 complex" = "#E6F598", "NuRD complex" = "#ABDDA4","PRC2 complex" = "#66C2A5",
                                                             "SWI SNF" = "#3288BD", "Glycosylase" = "#5E4FA2", "non-targeting" = "gray50")),
                                           gamma = anno_barplot(metaDataChromVarNoCond$gamma, bar_width = 0.9,ylim = c(-1,1)), 
                                           rho = anno_barplot(metaDataChromVarNoCond$rho, bar_width = 0.9)),
        row_split = factor(statesClass$type, levels = c("states", "types")),
        cluster_column_slice = FALSE,
        width = ncol(t(sigResPerturbNormRNA))*unit(4, "mm"),
        height = nrow(t(sigResPerturbNormRNA))*unit(4, "mm"),
        col = colorRamp2(seq(-3, 3, length = 3), c("#2268AD", "#EEEEEE", "#B31B2C")),
        heatmap_legend_param = list(title = "RT vs no RT"))
dev.off()

# Gene activity only
pdf("vs_snarc_GeneActivity_modules_res03_and_clinical_PerturbNorm.pdf",width = 13, height=8)
Heatmap(as.matrix(t(log(sigResPerturbNormGA[complete.cases(sigResPerturbNormGA),],2))),
        top_annotation = HeatmapAnnotation(GO = metaDataChromVarNoCond$GO, 
                                           col = list(cond = c("NoRT" = "#1B9E77", "5x1.8" = "#D95F02"),
                                                      GO = c("Histone demethylase" = "#9E0142", "Histone acetyltransferase" = "#D53E4F",
                                                             "Histone deacetylase" = "#F46D43", "SAGA complex" = "#FDAE61",
                                                             "HBO1 complex" = "#FEE08B", "Histone methyltransferase" = "#FFFFBF",
                                                             "NuA4 complex" = "#E6F598", "NuRD complex" = "#ABDDA4","PRC2 complex" = "#66C2A5",
                                                             "SWI SNF" = "#3288BD", "Glycosylase" = "#5E4FA2", "non-targeting" = "gray50")),
                                           gamma = anno_barplot(metaDataChromVarNoCond$gamma, bar_width = 0.9,ylim = c(-1,1)), 
                                           rho = anno_barplot(metaDataChromVarNoCond$rho, bar_width = 0.9)),
        row_split = factor(statesClass$type, levels = c("states", "types")),
        cluster_column_slice = FALSE,
        width = ncol(t(sigResPerturbNormGA))*unit(4, "mm"),
        height = nrow(t(sigResPerturbNormGA))*unit(4, "mm"),
        col = colorRamp2(seq(-3, 3, length = 3), c("#2268AD", "#EEEEEE", "#B31B2C")),
        heatmap_legend_param = list(title = "RT vs no RT"))
dev.off()

####################################
# Generate UMAP plots for figure 4 #
DefaultAssay(dat.keep) <- "peaks"
dat.keep <- RunTFIDF(dat.keep)
dat.keep <- FindTopFeatures(dat.keep, min.cutoff = 'q0')
dat.keep <- RunSVD(object = dat.keep)
DepthCor(dat.keep)

dat.keep <- RunUMAP(object = dat.keep, reduction = 'lsi', dims = 2:30)
dat.keep <- FindNeighbors(object = dat.keep,reduction = 'lsi',dims = 2:30)

res = 0.5
dat.keep <- FindClusters( object = dat.keep,  algorithm = 3,  resolution = res,  verbose = FALSE)

# embeddings for further exploration
embeddings <- Embeddings(dat.keep, reduction="umap")
embeddings <- data.frame(embeddings, orig.ident=as.character((dat.keep$orig.ident[row.names(embeddings)])),
                         cluster=as.character((dat.keep$peaks_snn_res.0.5[row.names(embeddings)])),
                         sgRNACond=as.character((dat.keep$sgRNACond[row.names(embeddings)])),
                         GOcond=as.character((dat.keep$GOcond[row.names(embeddings)])))
write.table(embeddings,'vs_snarc_peaks_dimplot_embeddings.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)

# dimplots
DimPlot(object = dat.keep, label = FALSE, group.by = "orig.ident")
ggsave(filename = "vs_snarc_peaks_dimplot.pdf",width = 6, height=4)

DimPlot(object = dat.keep, label = FALSE, group.by = "peaks_snn_res.0.5")
ggsave(filename = "vs_snarc_peaks_dimplot_peaks_snn_res.0.5.pdf",width = 5, height=4)

SCpubr::do_DimPlot(dat.keep, group.by = "peaks_snn_res.0.5", legend.position = "none",font.size = 4)

SCpubr::do_DimPlot(dat.keep, split.by = "sgRNACond", ncol = 8, legend.position = "none",font.size = 4)
ggsave(filename = "vs_snarc_peaks_dimplot_sgRNACond_array.pdf",width = 9, height=10)

SCpubr::do_DimPlot(dat.keep, split.by = "GOcond", ncol = 8, legend.position = "none",font.size = 4)
ggsave(filename = "vs_snarc_peaks_dimplot_GOcond_array.pdf",width = 9, height=4)

################
# UMAP for RNA #
DefaultAssay(dat.keep) <- "RNA"
dat.keep <- SCTransform(dat.keep)
dat.keep <- RunPCA(dat.keep)
ElbowPlot(object = dat.keep,ndims = 30) 
dat.keep <- RunUMAP(dat.keep,dims = 1:20, min.dist = 0.5)
dat.keep <- FindNeighbors(object = dat.keep,reduction='pca', dims = 1:20)
res = 0.5
dat.keep <- FindClusters(object = dat.keep, resolution = res,  verbose = FALSE)

#dimplots
DimPlot(object = dat.keep, label = FALSE, group.by = "orig.ident")
ggsave(filename = "vs_snarc_RNA_SCT_dimplot.pdf",width = 6, height=4)

DimPlot(object = dat.keep, label = FALSE, group.by = "SCT_snn_res.0.5")
ggsave(filename = "vs_snarc_RNA_SCT_dimplot_SCT_snn_res.0.5.pdf",width = 5, height=4)

SCpubr::do_DimPlot(dat.keep, split.by = "sgRNACond", ncol = 8, legend.position = "none",font.size = 4)
ggsave(filename = "vs_snarc_RNA_SCT_dimplot_sgRNACond_array.pdf",width = 9, height=10)

SCpubr::do_DimPlot(dat.keep, split.by = "GOcond", ncol = 8, legend.position = "none",font.size = 4)
ggsave(filename = "vs_snarc_RNA_SCT_dimplot_GOcond_array.pdf",width = 9, height=4)

# embeddings for further exploration
embeddings <- Embeddings(dat.keep, reduction="umap")
embeddings <- data.frame(embeddings, orig.ident=as.character((dat.keep$orig.ident[row.names(embeddings)])),
                         cluster=as.character((dat.keep$peaks_snn_res.0.5[row.names(embeddings)])),
                         sgRNACond=as.character((dat.keep$sgRNACond[row.names(embeddings)])),
                         GOcond=as.character((dat.keep$GOcond[row.names(embeddings)])))
write.table(embeddings,'vs_snarc_RNA_SCT_dimplot_embeddings.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)

# cell cycle scoring
dat.keep <- CellCycleScoring(dat.keep, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

SCpubr::do_FeaturePlot(dat.keep, features = c("S.Score","G2M.Score"),order = TRUE)
ggsave(filename = "vs_snarc_RNA_SCT_cellcycle.pdf",width = 7, height=4)

# Marker identification
Idents(dat.keep) <- "SCT_snn_res.0.5"
markers <- FindAllMarkers(object = dat.keep, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, assay="SCT",slot='data',test.use='wilcox',pseudocount.use = 1)
markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)

# marker heatmap
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top50 <- markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)
write.table(top50,paste("vs_snarc_RNA_SCT_top50markers_res",res,".txt",sep=''),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

plot <- DoHeatmap(subset(dat.keep, downsample = 200), features = top10$gene, size = 3) + NoLegend()
ggsave(filename = paste("vs_snarc_RNA_SCT_top10_heatmap_res",res,".pdf", sep=''),plot, width = 4, height = 6)

####################################
# Genome tracks for SETDB1 targets #
DefaultAssay(dat.keep) <- "ATAC"
Idents(dat.keep) <- "sgRNACond"

# analysis of ChIP-seq from ENCODE 
gset <- read.table('../KLF13_chip_bedidr_threshold_peakENCFF453MMH.bed',sep='\t',header=FALSE)
#gset <- read.table('../TCF3_IDR_threshold_peaks.csv',sep=',',header=TRUE)

gset=apply(gset[,1:3],1,function(x) paste(x,collapse = "-",sep=''))
gset <- gsub(' ','',gset)
gset=gset[grep("chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY",gset)]
gset=gset[grep("random",gset,invert = TRUE)]
gset <- unique(annot[nearest(StringToGRanges(gset), annot, ignore.strand=TRUE)]$gene_name)

# process plots with motif and perturbation of interest - run separately for up and down regulated regions
motif <- "KLF13"
perturb <- "KDM5C"

smp <- paste(perturb,"_NoRT",sep='')
smp2 <- paste(perturb,"_5x1.8",sep='')
dars <- resListAll[[smp2]]
dars <- dars[dars$p_val < 0.05 & dars$avg_log2FC > 0,] # all detected sites in a given perturbation
dars$gene.1 <- DARSAll.Dist[row.names(dars),"gene"]
dars <- dars[match(gset,dars$gene.1),] # Match ChIP-seq data with DARS
dars <- dars[!is.na(dars$gene),]
goi <- dars$gene.1

dat.subPlot <- subset(dat.keep,idents = c("non-targeting_NoRT","non-targeting_5x1.8",smp,smp2))
dat.subPlot <- LinkPeaks(object = dat.subPlot, peak.assay = "peaks", expression.assay = "SCT", genes.use = goi)

for (i in 1:length(goi)){
  skip_to_next <- FALSE
  tryCatch( plot <- CoveragePlot(dat.subPlot, region = goi[i],features = goi[i], extend.upstream = 10000,extend.downstream = 10000,max.downsample=30000,downsample.rate=0.1,window = 1000),
            error = function(e) { skip_to_next <<- TRUE},
            ggsave(filename = paste("coverage_plots/",perturb,"_",motif,"_DN_inRT_",goi[i],".pdf",sep=""),plot,width = 9,height = 5)
            )
  if(skip_to_next) { next }     
}


# Mass spec genes
secretedGOI = c("CUTA", "HP", "TIMP1", "INS", "APOA1", "APOA2", "TF", "HPX", "RPLP2", "SERPINE2", "CTSD", "PSAP", "CTSB", "TIMP2", "NPC2", "TRIM28",
                "ITIH4", "ECM1", "QSOX1", "CYR61", "HNRNPDL", "CALU", "PLIN3", "STC2", "LDHA", "PNP", "PLAU", "CFB", "CST3", "FGG", "FN1", "ALDOA",
                "CAPNS1", "SERPINE1", "CFI", "RPLP1", "P4HB", "PFN1", "COL1A2", "SPARC", "UCHL1", "HSPA5", "LAMC1", "MTHFD1", "COL6A2", "ACTN1",
                "PDIA4", "PRKCSH", "HSP90B1", "NME1", "PPIB", "PTX3", "CALR", "BLVRB", "PRDX5", "PDIA3", "PPP2R1A", "SRP14", "TAGLN2", "RPL3", "GARS",
                "CXCL5", "ALDH9A1", "SERPINH1", "RAD23B", "PLTP", "B2M", "NUCB1", "C1QBP", "LGALS3BP", "FSTL1", "NID2", "GANAB", "PDIA6", "PCOLCE", "PLEC",
                "RCN1", "TGFBI", "KDM3B", "TXNDC5", "VPS35", "NAA15", "RRBP1", "ADAMTS9", "SEPT9", "ANGPTL2")


dat.subPlot <- subset(dat.keep,idents = c("non-targeting_NoRT","non-targeting_5x1.8",smp,smp2))
dat.subPlot <- LinkPeaks(object = dat.subPlot, peak.assay = "peaks", expression.assay = "SCT", genes.use = secretedGOI)

for (i in 1:length(secretedGOI)){
  skip_to_next <- FALSE
  tryCatch( plot <- CoveragePlot(dat.subPlot, region = secretedGOI[i],features = secretedGOI[i], extend.upstream = 10000,extend.downstream = 10000,max.downsample=30000,downsample.rate=0.1,window = 1000),
            error = function(e) { skip_to_next <<- TRUE},
            ggsave(filename = paste("coverage_plots/",perturb,"_",motif,"_",secretedGOI[i],".pdf",sep=""),plot,width = 9,height = 5)
  )
  if(skip_to_next) { next }     
}

## TSS Enrichment by condition 
Idents(dat.keep) <- "orig.ident"
dat.keep <- TSSEnrichment(dat.keep, fast=FALSE, assay = 'ATAC',verbose=TRUE)
TSSPlot(dat.keep, group.by = 'orig.ident')
ggsave("vs_snarc_QC_TSSplot.pdf", width = 4, height = 3)

SCpubr::do_ViolinPlot(sample = dat.keep, features = c("nCount_RNA", "nFeature_RNA"),line_width =0.25,boxplot_width = 0.1,font.size = 6)
ggsave(filename = "vs_snarc_keep_QC_qc_violin.pdf",width = 4, height=2)

df <- data.frame(dat.keep$orig.ident,dat.keep$nCount_RNA,dat.keep$nFeature_RNA)
write.table(df,"vs_snarc_keep_QC_qc_violin.txt",quote=FALSE,sep='\t')

# get summary data - only those also with ATAC data 
df <- guideCallsTbl[intersect(colnames(dat.keep),colnames(dat.keep)),]
df <- melt(table(df[,c("gene","cond")]))

plot <- ggplot(df) +
  geom_bar(stat = 'identity',aes(x=gene,y=value,fill=cond),position='dodge') +
  labs(title = "Number of cells with sgRNAs detected", x="", y="") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("vs_snarc_keep_sgRNAs_hist.pdf", plot, width =3.5, height = 2)

write.table(df,"vs_snarc_keep_sgRNAs_hist.txt",quote=FALSE,sep='\t')
