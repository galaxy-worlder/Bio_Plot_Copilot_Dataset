#library(biomaRt)
library(GenomicRanges)
#library(httr)
#library(AnnotationHub)
library(genomation)
library(ggplot2)
library(rtracklayer)
library("readr")
library("MASS")
#logo
library("ggseqlogo")

library(grid)
library(cowplot)
library(ggplotify)
library("ggpubr")
library("dplyr")
library("corrplot")
library(RColorBrewer)
library("tidyverse")
library("ComplexHeatmap")
############# Chromosome lengths #################################

.libPaths("/projappl/project_2007567/project_rpackages_4.3.0")

representatives=read.table("/projappl/project_2006203/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", sep="\t", header=TRUE)

spacing_Yimeng=readRDS( file="/projappl/project_2007567/RProject/RData/spacing_Yimeng.Rds") #184

source("/projappl/project_2007567/code/functions_find_occurrences_of_kmers_in_spacing_motifs.R", echo=FALSE)

library(readxl)
bHLH_homeo_spacing_kmerfreq <- read_excel("/projappl/project_2007567/from_Yimeng/bHLH_homeodomain_new_spacing_all.kmerfreq12.xlsx")

bHLH_homeo_spacing_kmerfreq  <- bHLH_homeo_spacing_kmerfreq  %>%
  rename(
    "symbol"= "...1" ,
    "ligand"="...2",
    "info"= "...3" 
  )

#Draw heatmap of BARHL2_TCF4 k-mer frequencies

symbol="NHLH1_EVX2"

BARHL2_TCF4=bHLH_homeo_spacing_kmerfreq[which(bHLH_homeo_spacing_kmerfreq$symbol==symbol),]
BARHL2_TCF4_matrix=as.matrix(BARHL2_TCF4[, 5:33])
BARHL2_TCF4_matrix=BARHL2_TCF4_matrix/sum(BARHL2_TCF4_matrix)
BARHL2_TCF4_matrix=BARHL2_TCF4_matrix[, 1:10]

library("ComplexHeatmap")
Heatmap(BARHL2_TCF4_matrix, cluster_columns = FALSE, cluster_rows = FALSE, show_row_dend = FALSE,row_names_side = "left", 
        name="frequency",
        heatmap_legend_param=list(legend_direction="vertical"), 
        height = nrow(BARHL2_TCF4_matrix)*unit(1,"cm"), 
        width = ncol(BARHL2_TCF4_matrix)*unit(1, "cm"),)

Heatmap(t(BARHL2_TCF4_matrix), cluster_columns = FALSE, cluster_rows = FALSE, show_row_dend = FALSE,row_names_side = "left", 
        name="frequency",
        heatmap_legend_param=list(legend_direction="vertical"), 
        height = 1*unit(1,"cm"), 
        width = length(BARHL2_TCF4_matrix)*unit(1, "cm"),)


all_kmers=do.call(rbind, strsplit(bHLH_homeo_spacing_kmerfreq$`Spacing (bp)`, "\\."))

all_kmers=unique(c(all_kmers[,1],all_kmers[,2]))

homeodomain_kmers=c("TAATTA","TAAACG","CGTAAA", "CGTTTA","CACTTA")
bHLH_kmers=c("CAGCTG","CACCTG","CTGTCA", "CATATG",  "CACGTG", "CAGGTG")
library(Biostrings)
setwd("/projappl/project_2006203/TFBS/RProjects/TFBS/")
library("stringdist")

# motif_order=c("BHLHA15_EN2_TTTATG40NTAA_YRIII_TAATTANNNCATATG_m1_c3b0", 
#         "ATOH1_HOXB5_TCCATG40NCAC_YVIII_NTAATTANNNNCAKMTGN_m1_c3b0",  
#         "BARHL1_TCF12_TTGCGA40NTACG_YADIII_NTAATTGNNNNCACCTGN_m1_c3b0", 
#         "NHLH1_EVX2_TGCGCG40NTAG_YZIIII_NMNCAKSTGNNNNNTAATTAN_m1_c3b0", 
#         "OLIG2_EVX2_TGGCCG40NCTTT_YZIIII_NTAATTANNNNNCATATGN_m1_c3b0",
#         "OLIG2_HOXB5_TCTAAA40NCGC_YVIII_NCATATGNNNNNNYMATTAN_m1_c3b0",  
#         "ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNNCAKMTGN_m1_c3b0",    
#         "NHLH1_EVX2_TGCGCG40NTAG_YZIIII_NMNCAKSTGNNNNNNTAATTAN_m1_c3b0",
#         "CLOCK_ISL2_TCTGAA40NCAT_YUIII_NCAMTTANNNNNNCACGTGN_m1_c3b0",   
#         "NHLH1_HOXD8_TGAGAA40NGCG_YZIII_NCASCTGNNNNNNNNYAATTRN_m1_c3b0")

#row_split=factor(as.character(c(3, 4, 4, 5, 5, 6, 6, 6, 6, 8)), levels=c("3", "4", "5", "6", "8"))

motif_order=c("BHLHA15_EN2_TTTATG40NTAA_YRIII_TAATTANNNCATATG_m1_c3b0", 
        "BHLHE22_EVX2_TACGAT40NCCG_YZIIII_NRNCATATGNYNNTAATTAN_m1_c3b0", #coordinator
        "ATOH1_HOXB5_TCCATG40NCAC_YVIII_NTAATTANNNNCAKMTGN_m1_c3b0",  #this represents the coordinator
        "BARHL1_TCF12_TTGCGA40NTACG_YADIII_NTAATTGNNNNCACCTGN_m1_c3b0", 
        "NHLH1_EVX2_TGCGCG40NTAG_YZIIII_NMNCAKSTGNNNNNTAATTAN_m1_c3b0", 
        "NEUROG2_EVX2_TCACAT40NATG_YZIIII_NYMATTANNNNNCATATGN_m1_c3b0", #coordinator + 1
        "OLIG2_EVX2_TGGCCG40NCTTT_YZIIII_NTAATTANNNNNCATATGN_m1_c3b0", #this represents the coordinator + 1
        "OLIG2_HOXB5_TCTAAA40NCGC_YVIII_NCATATGNNNNNNYMATTAN_m1_c3b0",  
        "ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNNCAKMTGN_m1_c3b0",    
        "NHLH1_EVX2_TGCGCG40NTAG_YZIIII_NMNCAKSTGNNNNNNTAATTAN_m1_c3b0",
        "CLOCK_ISL2_TCTGAA40NCAT_YUIII_NCAMTTANNNNNNCACGTGN_m1_c3b0",   
        "NHLH1_HOXD8_TGAGAA40NGCG_YZIII_NCASCTGNNNNNNNNYAATTRN_m1_c3b0")

reps=c(1,0,1,1,1,0,1,1,1,1,1,1)

row_split=factor(as.character(c(3, 4,4, 4, 5,5, 5, 6, 6, 6, 6, 8)), levels=c("3", "4", "5", "6", "8"))



capselex_spacing=list()


for(motif in motif_order){
   #motif=motif_order[3]
   print(motif)

  pssm_dimer=get_pssm(filename=spacing_Yimeng %>% filter(ID==motif) %>% select(filename) %>% pull(filename))
  
  #Read k-mers
  
  kmers=bHLH_homeo_spacing_kmerfreq %>% #This may give multiple k-mers
    filter(symbol==representatives %>% filter(ID==motif) %>% select(symbol) %>% pull(symbol) ) %>% 
    select(`Spacing (bp)`) %>% 
    pull(`Spacing (bp)`)
  
  min_index=NA
  
  if( length(kmers)==1 ){
    kmers=strsplit(kmers, "\\.")[[1]]
  }else{
    #multiple k-mer pairs for the same symbol (TF-pair)
    #Extract the seeds from the original motif
    seed=representatives %>% filter(ID==motif) %>% select(seed) %>% pull(seed)
    
    motif_kmers <- strsplit(seed, "N+")[[1]][2:3]
    paste0(motif_kmers[1], ".", motif_kmers[2])
    
    distances <- stringdist::stringdist(paste0(motif_kmers[1], ".", motif_kmers[2]), kmers, method = "hamming")
    min_index <- which.min(distances)
    kmers=kmers[min_index]
    kmers=strsplit(kmers, "\\.")[[1]]
  }
  
  
  # Define the k-mer
  kmer1 <- kmers[1] #TAATTA is homeodomain
  kmer2 <- kmers[2] #CATATG is bHLH
  
  #Monomers
  
  kmer1_type=NA
  if(kmer1 %in% homeodomain_kmers){
    kmer1_type="Homeodomain"
  }else{
    kmer1_type="bHLH"
  }
  
  kmer2_type=NA
  if(kmer2 %in% homeodomain_kmers){
    kmer2_type="Homeodomain"
  }else{
    kmer2_type="bHLH"
  }
  
  
  
  info=spacing_Yimeng %>% filter(ID==motif) 
  
  monomer_families=strsplit(info$Lambert2018_families, "_")[[1]]
  if(info$switched_order==TRUE){
    monomer_families=rev(monomer_families)
  }
    
  if(monomer_families[1]==kmer1_type & monomer_families[2]==kmer2_type){
    print("Good!")
  }else{
    print("Bad!")
  }
  
  
  pssm_monomer1=get_pssm(filename=representatives %>% filter(ID==info$true.first.monomer) %>% select(filename) %>% pull(filename))
  pssm_monomer2=get_pssm(filename=representatives %>% filter(ID==info$true.second.monomer) %>% select(filename) %>% pull(filename))
  
  
  
  #Which motif is homeodomain, which one is bHLH?
  
  start1=find_best_match_position(kmer=kmer1, pssm_monomer1)$best_position
  end1=start1+nchar(kmer1)-1
  
  if(info$true.first.orientation=="-"){
    start1=ncol(pssm_monomer1)-start1+1
    end1=ncol(pssm_monomer1)-end1+1
  }
  
  
  start2=find_best_match_position(kmer=kmer2, pssm_monomer2)$best_position
  end2=start2+nchar(kmer2)-1
  
  if(info$true.second.orientation=="-"){
    start2=ncol(pssm_monomer2)-start2+1
    end2=ncol(pssm_monomer2)-end2+1
  }
  
  #Need to also find the corresponding positions in the heterodimer
  
  start_kmer1=find_best_match_position(kmer=kmer1, pssm_dimer)$best_position
  end_kmer1=start_kmer1+nchar(kmer1)-1
  
  start_kmer2=find_best_match_position(kmer=kmer2, pssm_dimer)$best_position
  end_kmer2=start_kmer2+nchar(kmer2)-1
  
  if(monomer_families[1]==kmer1_type & monomer_families[2]==kmer2_type){
    print("Ok!")
    #What is the spacing between kmers
    
    kmer_spacing=start_kmer2-end_kmer1-1
    monomer_spacing=info$gap_length
    if(is.na(monomer_spacing)){
      monomer_spacing=info$overlap_length
    }
    
  }else{
    print("Reverse order")
    #What is the spacing between kmers
    
    kmer_spacing=start_kmer1-end_kmer2-1
    monomer_spacing=info$gap_length
    
    if(is.na(monomer_spacing)){
      monomer_spacing=info$overlap_length
    }
    
  }
  
  #Shifted counts
  spacing_counts=as.data.frame(bHLH_homeo_spacing_kmerfreq %>% 
    filter(symbol==representatives %>% filter(ID==motif) %>% select(symbol) %>% pull(symbol) ))
  
  #Remove rowSum
  spacing_counts=spacing_counts[,-ncol(spacing_counts)]
  
  #Convert the counts to frequencies
  
  spacing_counts_matrix=as.matrix(spacing_counts[, 5:33])
  spacing_counts_matrix=spacing_counts_matrix/sum(spacing_counts_matrix)
  
  #There can be multiple rows, which one to choose. Choose the one with min distance to the motif seed
  if(nrow(spacing_counts)>1){
    spacing_counts=spacing_counts[min_index,]
    spacing_counts_matrix=spacing_counts_matrix[min_index,]
  }
  
  spacing_counts=spacing_counts[,5:33]
  
  spacing_counts_names=as.numeric(names(spacing_counts))
  
  spacing_counts_matrix=as.numeric(spacing_counts_matrix)
  
  names(spacing_counts_matrix)=as.character(spacing_counts_names+monomer_spacing-kmer_spacing)

  capselex_spacing[[motif]]=spacing_counts_matrix
  
}



start_min=min(sapply(capselex_spacing, function(x) min(as.numeric(names(x))))) #-7, lets make all start from this

capselex_spacing_extended=list()

for(i in 1:length(capselex_spacing)){
  cp=capselex_spacing[[i]]
  
  if(as.numeric(names(cp)[1]) > start_min){ #Put this later
    
    spacings=as.character(seq(start_min, as.numeric(names(cp)[1])-1,1 ))
    tmp=rep(NA, length(spacings))
    names(tmp)=spacings
    cp=c(tmp, cp)
  }
  
  capselex_spacing_extended[[names(capselex_spacing)[i]]]=cp
  
}

end_max=max(sapply(capselex_spacing_extended, function(x) max(as.numeric(names(x))))) #26

capselex_spacing_extended2=list()

for(i in 1:length(capselex_spacing_extended)){
  cp=capselex_spacing_extended[[i]]
  
  if(as.numeric(names(cp)[length(cp)]) < end_max){ #Put this later
    
    spacings=as.character(seq(as.numeric(names(cp)[length(cp)])+1, end_max,1 ))
    tmp=rep(NA, length(spacings))
    names(tmp)=spacings
    cp=c(cp,tmp)
  }
  
  capselex_spacing_extended2[[names(capselex_spacing_extended)[i]]]=cp
  
}

#Are all now of same size, all 34
sapply(capselex_spacing_extended2, length)

capselex_spacing_table=do.call(rbind, capselex_spacing_extended2 )
#Convert to frequencies




capselex_spacing_table_freqs=capselex_spacing_table[,as.character(seq(-7,11,1))]

match_ind=match( rownames(capselex_spacing_table_freqs), representatives$ID)

shorter_names=paste(representatives$symbol[match_ind], representatives$seed[match_ind], sep="_") 

result=as.data.frame(do.call(rbind, strsplit(representatives$Lambert2018_families[match_ind],"_")))
names(result)=c("X1", "X2")

match_ind=match( rownames(capselex_spacing_table_freqs), spacing_Yimeng$ID)
capselex_spacing_metadata=spacing_Yimeng[match_ind,]

rownames(capselex_spacing_table_freqs)=shorter_names


draw(Heatmap(capselex_spacing_table_freqs, cluster_columns = FALSE, cluster_rows = FALSE, show_row_dend = FALSE,row_names_side = "left", 
             name="frequency",
             heatmap_legend_param=list(legend_direction="horizontal")), heatmap_legend_side = "bot")

#How many NA entries there are in the beginning of the rows

na_maxind=rowSums(is.na(capselex_spacing_table_freqs))

for(i in 1:nrow(capselex_spacing_table_freqs)){
  if(na_maxind[i]!=0){
    capselex_spacing_table_freqs[i,]=t(c(capselex_spacing_table_freqs[i,(na_maxind[i]+1):ncol(capselex_spacing_table_freqs)], rep(NA, na_maxind[i])))
  }
}

saveRDS(na_maxind,file="/projappl/project_2007567/RProject/RData/na_maxind.Rds")

colnames(capselex_spacing_table_freqs)=as.character(seq(0,ncol(capselex_spacing_table_freqs)-1,1))

Protein_family_colors=c("#9d4977",
                        "#45bf52",
                        "#8d50ca",
                        "#74bf3e",
                        "#d975e2",
                        "#a2b937",
                        "#4e63d1",
                        "#c6b13d",
                        "#a782e8",
                        "#519531",
                        "#d247aa",
                        "#38c481",
                        "#e53580",
                        "#47cebe",
                        "#d03f36",
                        "#3bbac6",
                        "#e1733a",
                        "#5588e5",
                        "#df9a37",
                        "#a541a0",
                        "#86c170",
                        "#8c4fa1",
                        "#3b772c",
                        "#d25189",
                        "#4e995b",
                        "#ce435e",
                        "#6ebf92",
                        "#6e5397",
                        "#6a7a1d",
                        "#9a79bf",
                        "#9e8027",
                        "#4a71ac",
                        "#a64c28",
                        "#4facdc",
                        "#996531",
                        "#98a2e4",
                        "#6c6c2e",
                        "#db8bc0",
                        "#327243",
                        "#e6867a",
                        "#379577",
                        "#b86065",
                        "#206e54",
                        "#904455",
                        "#a1ac66",
                        "#d3a36c")


names(Protein_family_colors)=c("AP-2","TFAP","bHLH","Grainyhead","bZIP","EBF1",                   
                               "CENPB","DM",  "Ets; AT hook", "Ets", "E2F", "Forkhead",               
                               "GATA","GCM", "Nuclear receptor",        "HMG", "HMG/Sox",                 "Homeodomain",            
                               "CUT; Homeodomain",        "Homeodomain; Paired box", "Homeodomain; POU",        "POU", "HSF", "IRF",
                               "SMAD","Rel", "MADS box",                "Paired box" ,             "Prospero" ,               "Myb/SANT"   ,            
                               "MEIS","CSD", "NRF", "p53", "XPA", "RFX",
                               "RRM", "Runt","SAND","T-box","TEA", "C2H2 ZF",                
                               "C2H2 ZF; AT hook"  ,      "CCCH ZF"      ,           "Znf", "BED ZF")




#All are representatives
#groups=factor(spacing_Yimeng$new_representative, levels=c("YES", "NO") )
#library(RColorBrewer)
#group_colors <- c("red", "black") #"#E41A1C"             "#377EB8"             "#4DAF4A" 
#names(group_colors) <- levels(groups)

fontsize=12

library("ComplexHeatmap")



pdf("/projappl/project_2007567/capselex_spacings.pdf", width=15, height=10)
draw(Heatmap(capselex_spacing_table_freqs, cluster_columns = FALSE, cluster_rows=FALSE, 
             show_row_dend = FALSE,row_names_side = "left", name="frequency",
             #row_names_gp = gpar(col = group_colors[groups], font = 2),
        heatmap_legend_param=list(legend_direction="horizontal"),
        #height = nrow(p_values_all)*unit(0.618,"cm"), 
        height = nrow(capselex_spacing_table_freqs)*unit(1,"cm"), 
        width = ncol(capselex_spacing_table_freqs)*unit(1, "cm"),

        ), heatmap_legend_side = "bot")

#Highlight the observed spacing



decorate_heatmap_body("frequency", {
  
  for(i in 1:nrow(capselex_spacing_metadata)){
    print(i)
    site=capselex_spacing_metadata$gap_length[i]+7+1-na_maxind[i] 
    if(is.na(site)){
      print("composite")
      site=capselex_spacing_metadata$overlap_length[i]+7+1-na_maxind[i]
    }
    print(site)
    grid.rect(x = mean(c( (site-0.5)/ncol(capselex_spacing_table_freqs), (site-0.5)/ncol(capselex_spacing_table_freqs) )), #which(cluster_matrix["BYSL",]!=0)
                y = 1-mean(c((i-0.5)/nrow(capselex_spacing_table_freqs), (i-0.5)/nrow(capselex_spacing_table_freqs))), #which(rownames(cluster_matrix)=="BYSL")
                width = (abs((site) - (site)+1) )/ncol(capselex_spacing_table_freqs),
                height = (abs(i - i)+1 )/nrow(capselex_spacing_table_freqs),
                gp = gpar(col = "black", lwd = 2, fill = NA),
                default.units = "native")

   
  } 
})


dev.off()


#Add also coordinator + one more spacing conserved motif

grep("BHLH22_EVX2", spacing_Yimeng$symbol)


order=c("BHLHA15_EN2_TTTATG40NTAA_YRIII_TAATTANNNCATATG_m1_c3b0", 
        "ATOH1_HOXB5_TCCATG40NCAC_YVIII_NTAATTANNNNCAKMTGN_m1_c3b0",  
        "BARHL1_TCF12_TTGCGA40NTACG_YADIII_NTAATTGNNNNCACCTGN_m1_c3b0", 
        "NHLH1_EVX2_TGCGCG40NTAG_YZIIII_NMNCAKSTGNNNNNTAATTAN_m1_c3b0", 
        "OLIG2_EVX2_TGGCCG40NCTTT_YZIIII_NTAATTANNNNNCATATGN_m1_c3b0",
        "OLIG2_HOXB5_TCTAAA40NCGC_YVIII_NCATATGNNNNNNYMATTAN_m1_c3b0",  
        "ATOH1_EVX2_TTCGGG40NGAG_YZIIII_NTAATTANNNNNNCAKMTGN_m1_c3b0",    
        "NHLH1_EVX2_TGCGCG40NTAG_YZIIII_NMNCAKSTGNNNNNNTAATTAN_m1_c3b0",
        "CLOCK_ISL2_TCTGAA40NCAT_YUIII_NCAMTTANNNNNNCACGTGN_m1_c3b0",   
        "NHLH1_HOXD8_TGAGAA40NGCG_YZIII_NCASCTGNNNNNNNNYAATTRN_m1_c3b0")


row_split=factor(as.character(c(3, 4, 4, 5, 5, 6, 6, 6, 6, 8)), levels=c("3", "4", "5", "6", "8"))

#Only new representative bHLH-Homeodomain motifs 12
bHLH_homeo=spacing_Yimeng %>% filter(ID %in% order)

p_values_bHLH_homeo=p_values_all_orig[which(rownames(p_values_all_orig) %in% bHLH_homeo$ID),]

#Order bHLH_homeo and p_values_bHLH according to order

match_ind=match( order, rownames(p_values_bHLH_homeo))
p_values_bHLH_homeo=p_values_bHLH_homeo[match_ind,]
rownames(p_values_bHLH_homeo)

match_ind=match( order, bHLH_homeo$ID)
bHLH_homeo=bHLH_homeo[match_ind,]
bHLH_homeo$ID


shorter_names=paste(bHLH_homeo$symbol, bHLH_homeo$seed, sep="_")

long_names=rownames(p_values_bHLH_homeo)
rownames(p_values_bHLH_homeo)=shorter_names


library("ComplexHeatmap")

#Draw also reverse complements
#image_svg = paste0("/scratch/project_2006203/spacek/Figures_version2.2_Logos/svg/", long_names, ".svg")

alignment_svg = paste0("/scratch/project_2006203/spacek/Figures_version2.2/", long_names, ".svg")

# we only draw the image annotation for PNG images, while the others are the same
#install.packages("grImport2")
#install.packages("rsvg")

ha = HeatmapAnnotation(Logo = anno_image(alignment_svg, height= unit(8, "cm"), 
                                         border = FALSE,
                                         width= unit(8, "cm"),
                                         space = unit(0, "mm")), which="row") #, height=unit(0.618,"cm")) #, 


pdf("bHLH_homeo.pdf", width=15, height=7)
ht_list=Heatmap(p_values_bHLH_homeo, cluster_columns = FALSE, cluster_rows=FALSE, 
             show_row_dend = FALSE,row_names_side = "left", name="p-value",
             #row_split=row_split,
             #row_names_gp = gpar(col = group_colors[groups], font = 2),
             heatmap_legend_param=list(legend_direction="horizontal"),
             # left_annotation=rowAnnotation(`Protein family 1`=result$X1, `Protein family 2`= result$X2,  na_col = "white",
             #                               col=list(`Protein family 1`=Protein_family_colors,`Protein family 2`=Protein_family_colors ),
             #                               annotation_legend_param=list(#at=names(sort(table(c(anno$X1, anno$X2)), decreasing = TRUE)),
             #                                 gp=gpar(fontsize=fontsize)),
             #                               annotation_name_side="top",
             #                               #annotation_name_rot=0,
             #                               annotation_label=c("Protein family 1","Protein family 2"),
             #                               #annotation_label="Protein\nfamily",
             #                               gp=gpar(fontsize=fontsize),
             #                               width = unit(10, "cm")),
             #height = nrow(p_values_bHLH_homeo)*unit(0.618,"cm"), 
             height = nrow(p_values_bHLH_homeo)*unit(1,"cm"), 
             width = ncol(p_values_bHLH_homeo)*unit(1, "cm"),
             
)

draw(ht_list+ha, #row_title="spacing", 
     row_title_gp=gpar(fontface="bold"), 
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

#Highlight the observed spacing

decorate_heatmap_body("p-value", {
  
  for(i in 1:nrow(bHLH_homeo)){
    
    site=bHLH_homeo$gap_length[i]+7+1
    if(is.na(site)){
      site=bHLH_homeo$overlap_length[i]+7+1
    }
    
    grid.rect(x = mean(c( (site-0.5)/ncol(p_values_bHLH_homeo), (site-0.5)/ncol(p_values_bHLH_homeo) )), #which(cluster_matrix["BYSL",]!=0)
              y = 1-mean(c((i-0.5)/nrow(p_values_bHLH_homeo), (i-0.5)/nrow(p_values_bHLH_homeo))), #which(rownames(cluster_matrix)=="BYSL")
              width = (abs((site) - (site)+1) )/ncol(p_values_bHLH_homeo),
              height = (abs(i - i)+1 )/nrow(p_values_bHLH_homeo),
              gp = gpar(col = "black", lwd = 2, fill = NA),
              default.units = "native")
    
    
    for(j in 1:ncol(p_values_bHLH_homeo)){
      #print(j)
      if(!is.na(p_values_bHLH_homeo[i,j]) & p_values_bHLH_homeo[i,j] < 0.05){
        print(i)
        print(j)
        site=j
        grid.text("*",
                  x = mean(c( (site-0.5)/ncol(p_values_bHLH_homeo), (site-0.5)/ncol(p_values_bHLH_homeo) )), #which(cluster_matrix["BYSL",]!=0)
                  y = 1-mean(c((i-0.5)/nrow(p_values_bHLH_homeo), (i-0.5)/nrow(p_values_bHLH_homeo))))
        
      }
      
    }
    
  } 
})




dev.off()
