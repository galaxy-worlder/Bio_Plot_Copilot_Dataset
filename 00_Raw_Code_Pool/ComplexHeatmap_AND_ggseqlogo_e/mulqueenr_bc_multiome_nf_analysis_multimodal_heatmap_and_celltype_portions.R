#module load singularity
#sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
#singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome $sif

library(Signac)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(grid)
library(dplyr) 
library(ggplot2)
library(ggrepel)
library(patchwork)
library(presto)
library(seriation)
library(optparse)

option_list = list(
  make_option(c("-i", "--object_input"), type="character", default=NULL, 
              help="List of sample RDS files", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("/home/groups/CEDAR/mulqueen/bc_multiome/nf_analysis_round3/seurat_objects")
opt$object_input="merged.geneactivity.SeuratObject.rds"
dat=readRDS(opt$object_input)
table(dat$sample)

#Grab top overlapping TFs
topTFs <- function(markers_list,celltype, padj.cutoff = 1e-2,rna=NA,ga=NA,motifs=NA) {
  ctmarkers_rna <- dplyr::filter(rna, RNA.group == celltype) %>% 
    arrange(-RNA.auc)

    if(is.data.frame(motifs)) {
    ctmarkers_motif <- dplyr::filter(motifs, chromvar.group == celltype) %>% 
      arrange(-chromvar.auc)
    }

    if(is.data.frame(ga)) {
    ctmarkers_ga<- dplyr::filter(ga, GeneActivity.group == celltype) %>% 
      arrange(-GeneActivity.auc)
    }

    if(is.data.frame(motifs) && is.data.frame(ga)){    
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
      )
      top_tfs <- inner_join(
        x = top_tfs ,
        y = ctmarkers_ga [,c(2, 11, 6, 7)], by = "gene"
      )
    }else if(is.data.frame(motifs)) {
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
      )
    } else if (is.data.frame(ga)) {
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_ga[,c(2, 11, 6, 7)], by = "gene"
      )
    } 
  auc_colnames<-grep(".auc$",colnames(top_tfs))
  top_tfs$avg_auc <-  rowMeans(top_tfs[auc_colnames])
  top_tfs <- arrange(top_tfs, -avg_auc)
  top_tfs$celltype<-celltype
  return(top_tfs)
}


#Identify top markers
Identify_Marker_TFs<-function(x,group_by.="predicted.id",assay.="RNA",prefix.,pval_filt=1){
    #x[[assay.]]<-as(object = x[[assay.]], Class = "Assay")
    markers <- presto:::wilcoxauc.Seurat(X = x, group_by = group_by., 
      groups_use=unname(unlist(unique(x@meta.data[group_by.]))),
      y=unname(unlist(unique(x@meta.data[group_by.]))), 
      assay = 'data', seurat_assay = assay.)
    markers<-markers[markers$padj<=pval_filt,]
    write.table(markers,file=paste0(prefix.,"_",assay.,"_DE_table.tsv"),sep="\t",row.names=F,col.names=T,quote=F)
    colnames(markers) <- paste(assay., colnames(markers),sep=".")
    if (assay. == "chromvar") {
      motif.names <- markers[,paste0(assay.,".feature")]
      markers$gene <- ConvertMotifID(x, id = motif.names,assay="peaks")
    } else {
    markers$gene <- markers[,paste0(assay.,".feature")]
    }
    return(markers) 
}

#Average markers across groups
average_features<-function(x=hg38_atac,features=da_tf_markers$motif.feature,assay="chromvar",group_by.="predicted.id"){
    #Get gene activity scores data frame to summarize over subclusters (limit to handful of marker genes)
    dat_motif<-x[[assay]]@data[features,]
    dat_motif<-as.data.frame(t(as.data.frame(dat_motif)))
    sum_motif<-split(dat_motif,x@meta.data[,group_by.]) #group by rows to seurat clusters
    sum_motif<-lapply(sum_motif,function(x) apply(x,2,mean,na.rm=T)) #take average across group
    sum_motif<-do.call("rbind",sum_motif) #condense to smaller data frame

    sum_motif<-t(scale(sum_motif))
    sum_motif<-sum_motif[row.names(sum_motif)%in%features,]
    sum_motif<-sum_motif[complete.cases(sum_motif),]
    return(sum_motif)
}

#Make a heatmap of aligned multiple modalities
plot_top_TFs<-function(x=stromal,tf_markers=da_tf_markers,prefix="stromal",group_by.="predicted.id",CHROMVAR=TRUE,GA=TRUE,height.){
    tf_rna<-average_features(x=x,features=tf_markers$gene,assay="RNA",group_by.=group_by.)
    tf_rna<-tf_rna[row.names(tf_rna) %in% tf_markers$gene,]

    tf_motif<-average_features(x=x,features=tf_markers$chromvar.feature,assay="chromvar",group_by.=group_by.)
    tf_motif<-tf_motif[row.names(tf_motif) %in% tf_markers$chromvar.feature,]
    row.names(tf_motif)<-tf_markers[tf_markers$chromvar.feature %in% row.names(tf_motif),]$gene
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_motif)))
    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]

    tf_ga<-average_features(x=x,features=tf_markers$gene,assay="GeneActivity",group_by.=group_by.)
    tf_ga<-tf_ga[row.names(tf_ga) %in% tf_markers$gene,]
    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_ga)))
    tf_rna<-tf_rna[markers_list,]
    tf_ga<-tf_ga[markers_list,]

    markers_list<-Reduce(intersect, list(row.names(tf_rna),row.names(tf_motif),row.names(tf_ga)))
    tf_rna<-tf_rna[markers_list,]
    tf_motif<-tf_motif[markers_list,]
    tf_ga<-tf_ga[markers_list,]

    #set up heatmap seriation and order by RNA
    o = seriate(max(tf_rna) - tf_rna, method = "BEA_TSP")
    saveRDS(o,file=paste0(prefix,".geneactivity.dend.rds")) 
    side_ha_rna<-data.frame(ga_motif=tf_markers[get_order(o,1),]$RNA.auc)
    colfun_rna=colorRamp2(quantile(unlist(tf_rna), probs=c(0.5,0.80,0.95)),plasma(3))

    side_ha_motif<-data.frame(chromvar_motif=tf_markers[get_order(o,1),]$chromvar.auc)
    colfun_motif=colorRamp2(quantile(unlist(tf_motif), probs=c(0.5,0.80,0.95)),cividis(3))
    #Plot motifs alongside chromvar plot, to be added to the side with illustrator later
    motif_list<-tf_markers[tf_markers$gene %in% row.names(tf_motif),]$chromvar.feature
    plt<-MotifPlot(object = x,assay="peaks",motifs = motif_list[get_order(o,1)],ncol=1)+theme_void()+theme(strip.text = element_blank())
    ggsave(plt,file=paste0(prefix,".tf.heatmap.motif.pdf"),height=100,width=2,limitsize=F)

    side_ha_ga<-data.frame(ga_auc=tf_markers[get_order(o,1),]$GeneActivity.auc)
    colfun_ga=colorRamp2(quantile(unlist(tf_ga), probs=c(0.5,0.80,0.95)),magma(3))


    side_ha_col<-colorRamp2(c(0,1),c("white","black"))
    gene_ha = rowAnnotation(foo = anno_mark(at = c(1:nrow(tf_rna)), labels =row.names(tf_rna),labels_gp=gpar(fontsize=6)))


    rna_auc<-Heatmap(side_ha_rna,
        row_order = get_order(o,1),
        col=side_ha_col,
        show_column_names=FALSE,
        row_names_gp=gpar(fontsize=7))

    rna_plot<-Heatmap(tf_rna,
        row_order = get_order(o,1),
        column_order = get_order(o,2),
        name="RNA",
        column_title="RNA",
        col=colfun_rna,
        column_names_gp = gpar(fontsize = 8),
        show_row_names=FALSE,
        column_names_rot=90)

      ga_auc<-Heatmap(side_ha_ga,
          row_order = get_order(o,1),
          col=side_ha_col,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      ga_plot<-Heatmap(tf_ga,
          row_order = get_order(o,1),
          column_order = get_order(o,2),
          name="Gene Activity",
          column_title="Gene Activity",
          col=colfun_ga,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90)

      motif_auc<-Heatmap(side_ha_motif,
          row_order = get_order(o,1),
          col=side_ha_col,
          show_row_names=FALSE,
          show_column_names=FALSE,
          row_names_gp=gpar(fontsize=7))

      motif_plot<-Heatmap(tf_motif,
          row_order = get_order(o,1),
          column_order = get_order(o,2),
          name="TF Motif",
          column_title="TF Motif",
          col=colfun_motif,
          #top_annotation=top_ha,
          column_names_gp = gpar(fontsize = 8),
          show_row_names=FALSE,
          column_names_rot=90,
          right_annotation=gene_ha)
      
      plt1<-draw(ga_auc+ga_plot+rna_auc+rna_plot+motif_auc+motif_plot)


    pdf(paste0(prefix,".tf.heatmap.pdf"),height=height.)
    print(plt1)
    dev.off()
}

#Final wrapper function
run_top_TFs<-function(obj=stromal,prefix="stromal",i="predicted.id",n_markers=5,plot_height=10){
  markers<-lapply(c("RNA","GeneActivity","chromvar"),function(assay) Identify_Marker_TFs(x=obj,group_by.=i,assay.=assay,prefix.=prefix))
  names(markers)<-c("RNA","GeneActivity","chromvar")
  markers_out<-do.call("rbind",lapply(unique(obj@meta.data[,i]),function(x) head(topTFs(markers_list=markers,celltype=x,rna=markers$RNA,ga=markers$GeneActivity,motifs=markers$chromvar),n=n_markers))) #grab top 5 TF markers per celltype
  dim(markers_out)
  markers_out<-markers_out[!duplicated(markers_out$gene),]
  dim(markers_out)
  saveRDS(markers_out,file=paste0(prefix,"_celltype_TF_markers.RDS"))
  da_tf_markers<-readRDS(paste0(prefix,"_celltype_TF_markers.RDS"))
  plot_top_TFs(x=obj,tf_markers=da_tf_markers,prefix=prefix,group_by.=i)
 
}


#devtools::install_github("immunogenomics/presto")
#install.packages("seriation")
#install.packages('ggseqlogo')
x="test"
dat$cell_diag<-paste(dat$HBCA_predicted.id,dat$Diagnosis,dat$Mol_Diagnosis,sep="|")
DefaultAssay(dat)<-"RNA"
dat<-NormalizeData(dat)
dat <- FindVariableFeatures(dat)
dat <- ScaleData(dat)
dat [["RNA"]] <- as(object = dat[["RNA"]], Class = "Assay")

run_top_TFs(obj=dat,prefix=x,i="HBCA_predicted.id",n_markers=5,plot_height=15)

#plot sample composition heatmap

pdf("test_sample_comp.pdf")
Heatmap(t(proportions(t(table(dat$sample,dat$HBCA_predicted.id)),2)*100),col=colorRamp2(c(0,25,100),c("white","grey","black")))
dev.off()

dat_epi<-subset(dat,HBCA_predicted.id %in% c("luminal epithelial cell of mammary gland","basal cell"))
run_top_TFs(obj=dat_epi,prefix=x,i="cell_diag",n_markers=3,plot_height=15)
