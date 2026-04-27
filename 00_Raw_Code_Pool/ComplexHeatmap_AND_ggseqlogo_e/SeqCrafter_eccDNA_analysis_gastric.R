# samplenames -------------------------------------------------------------
samples<-c("JXM-1", "JXM-2", "JXM-3", "JXM-4", "JXM-5", "JXM-6",
           "JXM-7", "JXM-8", "JXM-9", "JXM-10", "JXM-11", "JXM-12",
           "JXM-13", "JXM-14", "JXM-15", "JXM-16", "JXM-17",
           "JXM-18", "JXM-19", "JXM-20", "JXM-21", "JXM-22", "JXM-23",
           "JXM-24", "JXM-25", "JXM-26","JXM-27")
cases<-c("JXM-1", "JXM-2", "JXM-3", "JXM-4", "JXM-5", "JXM-6",
         "JXM-7", "JXM-8", "JXM-9", "JXM-10", "JXM-11", "JXM-12",
         "JXM-13")
normals<-c("JXM-14", "JXM-15", "JXM-16", "JXM-17",
           "JXM-18", "JXM-19", "JXM-20", "JXM-21", "JXM-22", "JXM-23",
           "JXM-24", "JXM-25", "JXM-26","JXM-27")
randoms<-c("random_dataset.num1","random_dataset.num2","random_dataset.num3")

reference<-c(paste0("chr",seq(1,22)),"chrX","chrY","chrM")
getmapping<-function(x){
  read_tsv(stringr::str_glue("mapStat/{x}.mapping.stats.tsv"))%>%
    select(contig,mappedRatio)%>%
    setNames(c("refs",x))
}
totalM<-Reduce(\(x,y) full_join(x,y,by="refs") , lapply(samples, getmapping))

getmapping2<-function(x){
  read_tsv(stringr::str_glue("mapStat/{x}.align2plasmid.stat.tsv"))%>%
    select(contig,mappedRatio)%>%
    setNames(c("refs",x))
}
totalM2<-Reduce(\(x,y) full_join(x,y,by="refs") , lapply(samples, getmapping2))
totalM<-bind_rows(totalM,totalM2)

totalM$refs<-factor(totalM$refs,levels = c("P895","P1035",reference))

totalM<-totalM%>%
  arrange(refs)

openxlsx::write.xlsx(totalM,"figures/totalMapping.xlsx")
plotData<-totalM%>%
  tidyr::gather(sample,value,-refs)
plotData$sample<-factor(plotData$sample,levels = samples)
simplevis::gg_hbar_facet(data = plotData,
              y_var = refs,
              x_var = value,
              facet_var = sample,
              y_labels = \(x) x,
              facet_labels = \(x) x)
ggplot2::ggsave("figures/maping.pdf",width =9.24 ,height = 30)

# eccDNA from chrom -------------------------------------------------------


get_eccDNA_chrom<-function(x){
  reference<-c(paste0("chr",seq(1,22)),"chrX","chrY")
  chromSize<-read_tsv("../qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size",col_names = F)
  names(chromSize)<-c("Chromosome","Size")
  test<-read_tsv(stringr::str_glue("bedFile/{x}.ecc.bed"),col_names = F)%>%
    group_by(X1)%>%
    summarise(count=n())%>%
    filter(X1%in%reference)%>%
    rename(Chromosome=X1)%>%
    left_join(chromSize,by="Chromosome")%>%
    mutate(normalized_count=count/Size,
           normalized_percentage=normalized_count/sum(normalized_count))
  test$sample<-x
  test
}

get_eccDNA_chrom<-function(x){
  reference<-c(paste0("chr",seq(1,22)),"chrX","chrY")
  test<-read_tsv(stringr::str_glue("bedFile/{x}.ecc.bed"),col_names = F)%>%
    group_by(X1)%>%
    summarise(count=n())%>%
    filter(X1%in%reference)%>%
    rename(Chromosome=X1)
  test$sample<-x
  test
}

fin<-do.call("bind_rows",lapply(samples, get_eccDNA_chrom))
ct<-fin%>%
  select(Chromosome,count,sample)%>%
  tidyr::pivot_wider(names_from = "sample",values_from = "count")
cp<-fin%>%
  select(Chromosome,normalized_percentage,sample)%>%
  tidyr::pivot_wider(names_from = "sample",values_from = "normalized_percentage")

openxlsx::write.xlsx(list(count=ct,normalized=cp),"figures/eccDNA_derived_chrom.xlsx")
fin$sample<-factor(fin$sample,levels = samples)
fin$Chromosome<-factor(fin$Chromosome,levels = reference)
simplevis::gg_hbar_facet(data = fin,
              y_var = Chromosome,
              x_var = normalized_percentage,
              facet_var = sample,
              y_labels = \(x) x,
              facet_labels = \(x) x)
ggplot2::ggsave("figures/eccDNA_derived_chrom.pdf",width =9.24 ,height = 30)

# eccDNA numbers ----------------------------------------------------------


get_eccDNA_chrom<-function(x){
  test<-read_tsv(stringr::str_glue("bedFile/{x}.ecc.bed"),col_names = F)%>%
    filter(X1!="chrM")
  nrow(test)
}

fin<-tibble(sample=samples,eccnumber=sapply(samples,get_eccDNA_chrom))
fin$sample<-factor(fin$sample,levels = samples)
simplevis::gg_bar(data=fin,
       x_var = sample,
       y_var = eccnumber,
       x_title = "",
       y_title = "Total eccDNA number from Circle-Map",
       y_title_wrap = 30,
       x_labels = function(x) x)+
  ggplot2::theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1))
ggsave("figures/eccDNA_numbers.pdf",width = 6.99,height = 3.06)
openxlsx::write.xlsx(fin,"figures/eccDNA_numbers.xlsx")


# length ------------------------------------------------------------------

got_length<-function(x){
  test<-read_tsv(stringr::str_glue("bedFile/{x}.ecc.bed"),col_names = F)%>%
    mutate(length=X3-X2)%>%
    select(length)
  test$sample<-x
  test
}

test<-do.call('rbind',lapply(samples,got_length))

test<-test%>%
  mutate(gp=if_else(sample %in% cases,"GCT","NAT"))
plotData<-test%>%
  filter(length<=2000)

ggplot(plotData,aes(x=length))+
  stat_ecdf(aes(color=gp))+
  scale_color_manual(values = c("GCT"="#DC8F8E","NAT"="#689EF4"))+
  xlab("eccDNA length")+
  ylab("length accumulation")+
  theme_pubr()

ggsave("length_accu.pdf",width = 4.64,height =3.36 )

plotData2<-test%>%
  mutate(catergory=cut(length,
                       breaks = c(-1,200,500,1000,2000,Inf),
                       labels = c("<200","200~500","500~1K","1K~2K","2K~")))%>%
  group_by(sample,catergory)%>%
  summarise(number=n())%>%
  mutate(percentage=number/sum(number))

plotData$sample<-factor(plotData$sample,levels = samples)
gg_hbar_col(data=plotData,
            x_var = percentage,
            y_var = sample,
            col_var = catergory,
            stack = T,
            y_labels = function(x) x,
            col_labels = function(x) x,
            pal = pals::brewer.set3(4))
ggsave("figures/eccDNA_length_cat.pdf",width =3.88 ,height =5.48)
shuju<-plotData2%>%select(sample,catergory,percentage)%>%tidyr::pivot_wider(names_from = "catergory",values_from = "percentage",values_fill = 0)
openxlsx::write.xlsx(shuju,"figures/eccDNA_length_cat.xlsx")

#Gcs ------------------------------------------------------------------
getGCcontents<-function(x,sample=NULL){
  df<-read_tsv(x,col_names=F)%>%
    select(1:3)%>%
    setNames(c("seqname","start","end"))%>%
    mutate(start=start+1)%>%
    mutate(eccName=paste0(sample,"_",1:nrow(.)))
  gr <- makeGRangesFromDataFrame(df,keep.extra.columns = T,ignore.strand = T)
  gr <- keepStandardChromosomes(gr,species = "Homo_sapiens",pruning.mode = "coarse")
  tmpinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  tmpinfo <- keepStandardChromosomes(tmpinfo,species = "Homo_sapiens")
  seqinfo(gr)<-tmpinfo
  upstreams <- flank(gr, width(gr))
  downstrams <- flank(gr,width(gr),start = F)
  ecGC<-GCcontent(BSgenome.Hsapiens.UCSC.hg38, gr)
  upGC<-GCcontent(BSgenome.Hsapiens.UCSC.hg38, upstreams)
  downGC<-GCcontent(BSgenome.Hsapiens.UCSC.hg38, downstrams)
  tibble(ecName=gr$eccName,self=as.numeric(ecGC),upstream=as.numeric(upGC),downstream=as.numeric(downGC))
}

test1.gc<-getGCcontents("bedFile/JXM-1.ecc.bed",sample = "JXM-1")


readGC<-function(x){
  df<-read_tsv(stringr::str_glue("GCs/{x}.gcContents.txt"))%>%
    tidyr::gather(type,value,-ecc)%>%
    mutate(type2=case_when(type=="self" ~ "eccDNA",
                           type=="downstream" ~ "Downstream",
                           type=="upstream" ~ "Upstream"))
  df$sample<-x
  df
}


fin<-do.call('bind_rows',lapply(samples,readGC))
fin2<-do.call("bind_rows",lapply(randoms, readGC))
fin<-bind_rows(fin,fin2)
plotData2<-fin%>%
  mutate(type=case_when(type2=="Downstream" & (sample %in%samples) ~ "Downstream ecc",
                        type2=="Upstream"& (sample %in%samples) ~ "Upstream ecc",
                        type2=="eccDNA"& (sample %in%cases) ~ "GCT-eccDNA",
                        type2=="eccDNA"& (sample %in%normals) ~ "NAT-eccDNA",
                        type2=="Downstream" & (sample %in%randoms) ~ "In silico downstream ecc",
                        type2=="Upstream"& (sample %in%randoms) ~ "In silico upstream ecc",
                        type2=="eccDNA"& (sample %in%randoms) ~ "In silico eccDNA"))

plotData2$type<-factor(plotData2$type,levels = rev(c("GCT-eccDNA","NAT-eccDNA","Upstream ecc",
                                                     "Downstream ecc","In silico eccDNA","In silico upstream ecc",
                                                     "In silico downstream ecc")))
p<-ggplot(plotData2,aes(x=value,color=type))+
  geom_density(alpha=0.5)+
  cowplot::theme_cowplot()+
  xlab("GC Content(%)")+
  ylab("Density")+
  scale_x_continuous(labels = c("0","25","50","75","100"),breaks = c(0,0.25,0.5,0.75,1))+
  geom_vline(xintercept = 0.39,linetype="dashed")+
  theme(legend.position = "none")+
  scale_fill_brewer(palette = "Set1")
ggsave("figures/GCs.pdf",width = 5.63,height = 4.90)

# heatmap -----------------------------------------------------------------

readGC<-function(x){
  df<-read_tsv(stringr::str_glue("GCs/{x}.gcContents.txt"))%>%
    setNames(c("ecc","Upstream","eccDNA","Downstream"))
  df$sample<-x
  df
}

fin<-do.call('bind_rows',lapply(c(cases,normals),readGC))
fin<-fin%>%mutate(type=if_else(sample%in%cases,"GCT","NAT"))
fin$label<-seq(1,nrow(fin))
fin2<-fin%>%
  group_by(type)%>%
  slice_sample(prop = 0.001,replace = F)
mts<-fin2%>%
  ungroup()%>%
  select(label,Upstream,eccDNA,Downstream)%>%
  tibble::column_to_rownames(var="label")%>%
  as.matrix()
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.4, 0.8), c("#053061", "#f7f7f7", "#67001f"))
ht<-Heatmap(mts,
            cluster_rows = T,
            cluster_columns = F,
            show_row_names = F,
            show_row_dend = F,
            row_split = c(rep("GCT",3673),rep("NAT",506)),
            col = col_fun,
            name = "GC content(%)")
pdf("figures/Fig1B.pdf",width =3.98,height =4.82 )
draw(ht)
dev.off()

# length dist -------------------------------------------------------------


test<-test%>%
  mutate(gp=if_else(sample%in%cases,"GCT","NAT"))

plotData<-test%>%
  filter(length<=2000)

p<-ggplot(plotData,aes(x=length,fill=gp))+
  geom_density(alpha=0.6)+
  scale_x_continuous(limits = c(0,2000))+
  scale_fill_manual(values = c("GCT"="#b2182b","NAT"="#5581B0"))+
  cowplot::theme_cowplot()+
  xlab("Fragment length(bp)")+
  ylab("Density")
ggsave("figures/length_dst.pdf",width = 6.10,height =2.93 )

plotly::ggplotly(p)
##188,364,568,751,955,1135

# accumulate --------------------------------------------------------------

ggplot(plotData,aes(x=length,fill=gp))+
  stat_ecdf(geom="density",alpha=0.5)+
  scale_x_continuous(limits = c(0,2000))+
  scale_fill_manual(values = c("GCT"="#b2182b","NAT"="#5581B0"))+
  geom_vline(xintercept = 358,linetype="dashed")+
  cowplot::theme_cowplot()+
  xlab("Fragment length(bp)")+
  ylab("Cumulative Frequency")
ggsave("figures/length_cumulate.pdf",width = 6.10,height =3.28 )

### circos

getIN<-function(x){
  df<-read_tsv(stringr::str_glue("bedFile/{x}.ecc.bed"),col_names = F)
  df<-df%>%
    select(1:3)%>%
    mutate(sample=x)
  df
}

fin<-do.call("bind_rows",lapply(c(normals,cases), getIN))
fin<-fin%>%
  mutate(type=if_else(sample %in%cases,"GCT","NAT"))%>%
  mutate(Y_random=runif(nrow(fin),min = -1))
fin2<-fin%>%
  group_by(type)%>%
  slice_sample(prop = 0.01,replace = F)
library(circlize)
pdf("figures/all_samples.cicos.pdf",width = 9.6,height =6.68 )
circos.clear()
circos.par("start.degree"=87,
           gap.after=c(rep(2,23),10),
           track.height=0.1)
circos.initializeWithIdeogram(species = "hg38",ideogram.height=0.04)
circos.genomicTrack(
  fin2%>%filter(type=="NAT"),
  numeric.colum=6,
  bg.col="white",
  bg.border=F,
  panel.fun=function(region,value,...){
    circos.genomicPoints(region,value,pch=16,cex=0.25,col="#08519c",...)
  }
)

circos.genomicTrack(
  fin2%>%filter(type=="GCT"),
  numeric.colum=6,
  bg.col="white",
  bg.border=F,
  panel.fun=function(region,value,...){
    circos.genomicPoints(region,value,pch=16,cex=0.25,col="#ef3b2c",...)
  }
)
dev.off()

#  EPM case control -------------------------------------------------------


allmappings<-openxlsx::read.xlsx("mapStat/total_mapping_reads.xlsx")

get_totals<-function(x){
  df<-read_tsv(stringr::str_glue("bedFile/{x}.ecc.bed"),col_names = F)
  tibble(sample=x,numbers=nrow(df))
}

fin<-do.call("bind_rows",lapply(samples,get_totals))
names(allmappings)<-c("sample","mappings")
fin<-fin%>%left_join(allmappings,by="sample")
fin$mappings<-as.numeric(fin$mappings)
fin<-fin%>%
  mutate(EPM=numbers/mappings*(10^6))
fin<-fin%>%
  mutate(gp=if_else(sample %in%cases,"GCT","NAT"))

df_pvalue<-fin%>%
  rstatix::t_test(EPM~gp)%>%
  rstatix::add_xy_position()
ggplot(fin,aes(x=gp,y=EPM))+
  geom_boxplot(aes(fill=gp),width=0.8,outlier.alpha = 0)+
  geom_jitter()+
  ggprism::theme_prism(base_line_size = 0.4,
                       base_size = 11)+
  scale_fill_manual(values = c("GCT"="#990000","NAT"="#084594"))+
  ggprism::add_pvalue(df_pvalue)+
  xlab("")+
  ylab("eccDNA per-million \nmapping reads")

ggsave("figures/EPM.pdf",width = 3.75,height = 3.34)

# ecc per chrom -----------------------------------------------------------

stat_chromo_num<-function(x) {
  data<-read_tsv(stringr::str_glue("bedFile/{x}.ecc.bed"),col_names = F)
  names(data)[1]<-"Chromosome"
  result<-data%>%
    group_by(Chromosome)%>%
    dplyr::summarise(num=n())
  result$samples<-x
  result
}
chromos<-c(paste0("chr",seq(1,22)),"chrX","chrY")
fin<-do.call('rbind',lapply(samples,stat_chromo_num))
fin<-fin%>%
  filter(Chromosome %in% chromos)
chromosome_size<-read_tsv("../qd-ECC4/S/CircleSeq/ECC_report/DataBase/hg38.chromo.size",col_names = F)
names(chromosome_size)<-c("Chromosome","size")
need_size<-chromosome_size%>%filter(Chromosome%in%chromos)
plotData<-fin%>%left_join(need_size,by="Chromosome")%>%
  tidyr::replace_na(replace = list(num=0))%>%
  mutate(type=if_else(samples%in%cases,"GCT","NAT"))


plotData<-plotData%>%
  mutate(countP=num/size*10^6)%>%
  group_by(samples)%>%
  mutate(countP=countP/sum(countP)*100)%>%
  ungroup()



plotData$Chromosome<-factor(plotData$Chromosome,levels = chromos)

male<-plotData%>%filter(type=="GCT")
female<-plotData%>%filter(type=="NAT")

p1<-ggplot(female,aes(x=Chromosome,y=countP))+
  geom_boxplot(fill='#377eb8',width=0.45,color="black",outlier.size = 0.1)+
  geom_hline(yintercept = 4.3,linetype="dashed")+
  ggprism::theme_prism(base_size = 11,
                       axis_text_angle = 45,
                       border = T,
                       base_line_size = 0.5,
                       base_fontface = "plain")+
  scale_y_continuous(limits = c(0,8))+
  xlab("")+
  ylab("Percent eccDNA Per Mb")

p2<-ggplot(male,aes(x=Chromosome,y=countP))+
  geom_boxplot(fill='#b2182b',width=0.45,color="black",outlier.size = 0.1)+
  geom_hline(yintercept = 4.3,linetype="dashed")+
  ggprism::theme_prism(base_size = 11,
                       axis_text_angle = 45,
                       border = T,
                       base_line_size = 0.5,
                       base_fontface = "plain")+
  scale_y_continuous(limits = c(0,8))+
  xlab("")+
  ylab("Percent eccDNA Per Mb")
p2/p1
ggsave("figures/ecc_chrom_dist.pdf",width =8.14 ,height =5.87)


# element -----------------------------------------------------------------


## must remove chrM only
dbs<-read_tsv("../CKDs/hg38.genetic_element.sort.bed.stat")
genome<-read_tsv("../qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size",col_names = F)

##########
get_element<-function(x,df2){
  df<-read_tsv(stringr::str_glue("elementAnno/{x}.startAnno.bed"),col_names = F)
  lishu<-ncol(df)
  names(df)[c(1:3,lishu)]<-c("chrm","st","ed","elements")
  df<-df%>%
    select(chrm,st,ed,elements)%>%
    filter(chrm!="chrM")%>%
    mutate(ecc=paste(paste(chrm,st,sep=":"),ed,sep="-"))%>%
    tidyr::separate(elements,into=c("gene","elements"),sep="\\|")%>%
    group_by(elements)%>%
    summarise(n=n_distinct(ecc))%>%
    left_join(df2,by="elements")%>%
    mutate(ratio=n/sum(n))%>%
    mutate(ratio2=ratio/total_p)%>%
    select(elements,ratio2)
  df$sample<-x
  df
}
fin<-do.call("bind_rows",lapply(c(normals,cases),get_element,df2=dbs))

fin<-fin%>%
  mutate(type=if_else(sample%in%cases,"GCT","NAT"))
df_p_val <- fin %>%
  rstatix::group_by(elements) %>%
  rstatix::wilcox_test(ratio2 ~ type) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "elements")

ggplot(fin,aes(x=elements,y=ratio2))+
  geom_boxplot(aes(fill=type),width=0.8,color="black",outlier.size = 0.1)+
  ggprism::theme_prism(base_size = 11,
                       axis_text_angle = 45,
                       border = T,
                       base_line_size = 0.5,
                       base_fontface = "plain")+
  scale_fill_manual(values = c("GCT"="#b2182b","NAT"="#377eb8"))+
  xlab("")+
  ylab("Normalized ratio of eccDNA\n in different elements")+
  geom_hline(yintercept = 1.21,linetype="dashed")

ggsave("figures/elements.anno.pdf",width = 6.78,height = 2.82)

# numbers -----------------------------------------------------------------


get_totals<-function(x){
  df<-read_tsv(stringr::str_glue("bedFile/{x}.ecc.bed"),col_names = F)
  tibble(sample=x,numbers=nrow(df))
}

fin<-do.call("bind_rows",lapply(c(cases,normals),get_totals))
write_tsv(fin,"figures/TableS1.tsv")

#  coding vs. gentic ------------------------------------------------------

randoms_nums<-tibble(sample=c("random_dataset.num1","random_dataset.num2","random_dataset.num3"),numbers=c(10000,10000,10000))
eccnumbers<-read_tsv("figures/TableS1.tsv")
eccnumbers<-bind_rows(eccnumbers,randoms_nums)
read_gentic<-function(x){
  df<-read_tsv(stringr::str_glue("coding/{x}.genetic.bed"),col_names = F)
  shuliang<-df%>%
    select(1:3)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep = "-"))%>%
    distinct(ecc)%>%
    nrow()
  tibble(sample=x,nums=shuliang)
}
total_map2gene<-do.call("bind_rows",lapply(c(cases,normals,randoms),read_gentic))

read_coding<-function(x){
  df<-read_tsv(stringr::str_glue("coding/{x}.coding.bed"),col_names = F)
  shuliang<-df%>%
    select(1:3)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep = "-"))%>%
    distinct(ecc)%>%
    nrow()
  tibble(sample=x,coding=shuliang)
}
total_map2coding<-do.call("bind_rows",lapply(c(cases,normals,randoms),read_coding))
fin<-full_join(eccnumbers,total_map2coding,by="sample")%>%
  full_join(total_map2gene,by="sample")
plotData<-fin%>%
  mutate(gentic_ratio=nums/numbers*100,
         coding_ratio=coding/numbers*100)%>%
  select(1,5,6)%>%
  tidyr::gather(gp2,ratio,-sample)%>%
  mutate(gp=case_when(sample %in% cases ~ "GCT",
                      sample %in% normals ~ "NAT",
                      sample %in% randoms ~ "in silico"))
df_p_val <- plotData %>%
  rstatix::group_by(gp2) %>%
  rstatix::t_test(ratio ~ gp) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "gp2")
ggplot(plotData,aes(x=gp2,y=ratio))+
  geom_boxplot(aes(fill=gp),width=0.8)+
  ggprism::theme_prism(base_line_size = 0.4,
                       base_size = 11)+
  scale_fill_manual(values = c("GCT"="#990000","NAT"="#084594","in silico"="#bdbdbd"))+
  xlab("")+
  ylab("eccDNA percentage(%)")+
  ggprism::add_pvalue(df_p_val, 
                      xmin = "xmin", 
                      xmax = "xmax",
                      label = "p.adj.signif",
                      tip.length = 0.01)
ggsave("figures/gentic_regions.pdf",width = 3.78,height = 3.67)


# repeat ------------------------------------------------------------------

dbs<-read_tsv("../CKDs/repeats/UCSC_repeatMarsker.region.bed",col_names = F)
stst<-dbs%>%
  filter(X1 !="chrM")%>%
  mutate(length=X3-X2)%>%
  group_by(X4)%>%
  summarise(total_l=sum(length))

chrmSize<-read_tsv("../甲状腺癌/old/plasma/dbs/hg38.chromo.size",col_names = F)
chrmSize<-chrmSize%>%
  filter(X1!="chrM")%>%
  summarise(total_l=sum(X2))

##normalized repeats
stst<-stst%>%
  mutate(ratio=total_l/(chrmSize$total_l))
names(stst)[1]<-"name"
stst<-stst%>%
  select(1,3)%>%
  rename(lgthR=ratio)

getReps<-function(x,dd=stst){
  df<-read_tsv(stringr::str_glue("repeatStat/repStats.{x}.tsv"))
  df<-df%>%
    left_join(dd,by="name")%>%
    mutate(fin=ratio/lgthR)%>%
    select(name,fin)%>%
    setNames(c("name",x))
  df
}
fin<-Reduce(\(x,y){full_join(x,y,by="name")},lapply(c(samples), getReps))
write_tsv(fin,"figures/repeats.tsv")

# gene drop ---------------------------------------------------------------

genes<-read_tsv("../乳腺癌/bedFile/dbs/hg38.coding.bed",col_names = F)
names(genes)<-c("chr","Start","End","gene")
genes<-genes%>%
  mutate(length=End-Start)
get_genes<-function(x,dbs=genes){
  df<-read_tsv(stringr::str_glue("genedrops/{x}.startAnno.bed"),col_names = F)
  df<-df%>%
    select(1:4,8)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
    group_by(X8)%>%
    distinct(ecc,.keep_all = T)%>%
    summarise(count=n())%>%
    arrange(desc(count))%>%
    rename(gene=X8)%>%
    left_join(dbs,by="gene")%>%
    mutate(pct=count/length)%>%
    mutate(pct2=pct/sum(pct)*(10^6))%>%
    select(1,8)
  names(df)[2]<-x
  df
}
hebing<-function(x,y){full_join(x,y,by="gene")}
fin<-Reduce(hebing,lapply(c(cases,normals), get_genes))
fin<-fin%>%tidyr::gather(sample,value,-gene)%>%tidyr::replace_na(replace = list(value=0))%>%tidyr::pivot_wider(names_from = "sample",values_from = "value")
saveRDS(fin,"gene_drop.RDS")

dffs<-fin%>%
  tidyr::gather(sample,TPM,-gene)%>%
  tidyr::replace_na(replace = list(TPM=0))%>%
  mutate(gp=if_else(sample %in% cases,"GCT","NAT"))%>%
  group_by(gene)%>%
  rstatix::pairwise_wilcox_test(TPM ~ gp,p.adjust.method = "BH")

saveRDS(dffs,"gene_drop_diffs.RDS")


fin<-readRDS("gene_drop.RDS")
dffs<-readRDS("gene_drop_diffs.RDS")
plotData<-fin%>%tibble::column_to_rownames(var="gene")%>%as.matrix()
nima<-dffs%>%filter(p.adj.signif%in%c("***","****"))
plt<-plotData[nima$gene,]
plt[is.na(plt)]<-0
plt<-log2(plt+1)


ht<-Heatmap(plt,
            #col = col_fun,
            show_row_names = F,
            show_column_names = F,
            cluster_rows = T,
            cluster_columns = F,
            top_annotation = HeatmapAnnotation(group=c(rep("GCT",13),rep("NAT",13)),
                                               col = list(group=c("GCT"="#e41a1c","NAT"="#377eb8"))),
            name = "Log normalized \neccDNA count")


kygenes<-tibble(genes=rownames(plt)[row_order(ht)])
openxlsx::write.xlsx(kygenes,"figures/enrich_genes.xlsx")
pdf("figures/gene_diffs.pdf",width = 6.22,height = 6.12)
draw(ht)
dev.off()

# miRNA  ------------------------------------------------------------------

miRNAs<-read_tsv("ecc-enhancer.bed",col_names = F)%>%
  mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
  tidyr::separate(X4,into = c("sample","smalllabel"),sep="_")%>%
  select(1:5,9,10)

exp_mat<-miRNAs%>%
  select(X8,sample,ecc)%>%
  group_by(sample,X8)%>%
  summarise(n=n())
exp_mat$n<-1
exp_mat<-exp_mat%>%tidyr::pivot_wider(names_from = "sample",values_from = n)
#exp_mat$`JXM-20`<-NA
#exp_mat$`JXM-21`<-NA
#exp_mat$`JXM-25`<-NA
wori<-exp_mat%>%tibble::column_to_rownames(var='X8')%>%as.matrix()
wori[is.na(wori)]<-0
case_mat<-wori[,cases]
control_mat<-wori[,normals]
case_shunxu<-tibble(miRNA=rownames(case_mat),numbers=apply(case_mat,1,function(x){sum(x)}))
control_shunxu<-tibble(miRNA=rownames(control_mat),numbers=apply(control_mat,1,function(x){sum(x)}))

case_shunxu<-case_shunxu%>%
  mutate(nubers2=13-numbers)
control_shunxu<-control_shunxu%>%
  mutate(nubers2=13-numbers)%>%
  dplyr::rename(wocao1=numbers,wocao2=nubers2)

haha<-full_join(case_shunxu,control_shunxu,by="miRNA")

haha<-haha%>%
  rowwise()%>%
  mutate(pvalues=fisher.test(matrix(c(numbers,nubers2,wocao1,wocao2),
                                    nrow = 2,
                                    byrow = T),
                             alternative = "two.sided")$p)
haha<-haha%>%arrange(pvalues,desc(numbers))
needs<-haha[1:108,]
need2<-miRNAs%>%filter(X8%in%(needs$miRNA))
rst<-list(dff=needs,raw=need2)
openxlsx::write.xlsx(rst,"figures/miRNA.data.xlsx")
#genes<-read_tsv("figures/fullgene.data.tsv")
#genes<-genes[1:19,]$miRNA
#genes<-miRNAData[1:7,]$miRNA
genes<-haha%>%filter(pvalues<0.05)%>%pull(miRNA)
plotData<-wori[genes,c(cases,normals)]
lieshu<-apply(plotData,2,function(x) sum(x,na.rm=T))
col_fun<-colorRamp2(c(0,3),c("white", "firebrick3"))
ht<-Heatmap(plotData,
            show_column_names = F,
            show_row_names = T,
            cluster_columns = F,
            cluster_rows = F,
            show_column_dend = F,
            show_row_dend = F,
            col = col_fun,
            na_col = "#d9d9d9",
            rect_gp = gpar(col="white",lwd=1.5),
            row_names_gp = gpar(fontface="italic",cex=0.7),
            show_heatmap_legend = T,
            name = "eccDNA Frequency",
            top_annotation = HeatmapAnnotation(genes=anno_barplot(lieshu,
                                                                  border=F,
                                                                  axis=T,
                                                                  gp=gpar(fill="firebrick3",
                                                                          col="transparent"))),
            bottom_annotation = HeatmapAnnotation(group=c(rep("case",13),rep("control",13)),
                                                  col = list(group=c("case"="#b2182b","control"="#2166ac")),
                                                  simple_anno_size = unit(0.2, "cm"),show_legend = F))

plotData[is.na(plotData)]<-0

pheatmap::pheatmap(plotData,
         show_rownames = F,
         color = colorRampPalette(c("white", "firebrick3"))(100),
         cluster_rows = F,
         cluster_cols = F,
         border_color = "white",
         filename = "figures/miRNA-sup.pdf",width =5.44 ,height =5.42 )

pdf("figures/enhancer-ecc.pdf",width =5.29 ,height =5.15)
draw(ht)
dev.off()

# gene60  ------------------------------------------------------------------

miRNAs<-read_tsv("ecc-gene60.bed",col_names = F)%>%
  filter(X9>=60)%>%
  mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
  tidyr::separate(X4,into = c("sample","smalllabel"),sep="_")%>%
  select(1:5,9,11)

exp_mat<-miRNAs%>%
  select(X8,sample,ecc)%>%
  group_by(sample,X8)%>%
  summarise(n=n())
exp_mat$n<-1
exp_mat<-exp_mat%>%tidyr::pivot_wider(names_from = "sample",values_from = n)

wori<-exp_mat%>%tibble::column_to_rownames(var='X8')%>%as.matrix()
wori[is.na(wori)]<-0
##need-genes

case_mat<-wori[,cases]
control_mat<-wori[,normals]
case_shunxu<-tibble(miRNA=rownames(case_mat),numbers=apply(case_mat,1,function(x){sum(x)}))
control_shunxu<-tibble(miRNA=rownames(control_mat),numbers=apply(control_mat,1,function(x){sum(x)}))

case_shunxu<-case_shunxu%>%
  mutate(nubers2=13-numbers)
control_shunxu<-control_shunxu%>%
  mutate(nubers2=13-numbers)%>%
  dplyr::rename(wocao1=numbers,wocao2=nubers2)

haha<-full_join(case_shunxu,control_shunxu,by="miRNA")

haha<-haha%>%
  rowwise()%>%
  mutate(pvalues=fisher.test(matrix(c(numbers,nubers2,wocao1,wocao2),
                                    nrow = 2,
                                    byrow = T),
                             alternative = "two.sided")$p)
haha<-haha%>%arrange(pvalues,desc(numbers))
needs<-haha
need2<-miRNAs
rst<-list(dff=needs,raw=need2)
openxlsx::write.xlsx(rst,"figures/gene-60.data.xlsx")

#genes<-read_tsv("gene60.need.txt",col_names = F)$X1

plotData<-wori[genes,c(cases,normals)]
tmp<-plotData[1:44,]
lieshu<-apply(tmp,2,function(x) sum(x,na.rm=T))
col_fun<-colorRamp2(c(0,10),c("white", "firebrick3"))
ht<-Heatmap(tmp,
            show_column_names = F,
            show_row_names = T,
            cluster_columns = F,
            cluster_rows = F,
            show_column_dend = F,
            show_row_dend = F,
            col = col_fun,
            na_col = "#d9d9d9",
            rect_gp = gpar(col="white",lwd=1.5),
            row_names_gp = gpar(fontface="italic",cex=0.7),
            show_heatmap_legend = T,
            name = "eccDNA Frequency",
            top_annotation = HeatmapAnnotation(genes=anno_barplot(lieshu,
                                                                  border=F,
                                                                  axis=T,
                                                                  gp=gpar(fill="firebrick3",
                                                                          col="transparent"))),
            bottom_annotation = HeatmapAnnotation(group=c(rep("case",13),rep("control",13)),
                                                  col = list(group=c("case"="#b2182b","control"="#2166ac")),
                                                  simple_anno_size = unit(0.2, "cm"),show_legend = F))






pdf("figures/gene60-key2.pdf",width =5.48 ,height =6.14)
draw(ht)
dev.off()


# junction motif ----------------------------------------------------------

up<-read_tsv("control_ecc.up_15bp.fa",col_names = F)
up$X2<-stringr::str_to_upper(up$X2)
down<-read_tsv("control_ecc.down_15bp.fa",col_names = F)
down$X2<-stringr::str_to_upper(down$X2)
up<-up%>%filter(stringr::str_length(X2)==30)
down<-down%>%filter(stringr::str_length(X2)==30)
upS<-up%>%slice_sample(n=10000,replace = F)
downS<-down%>%slice_sample(n=10000,replace = F)
p1<-ggseqlogo(upS$X2,seq_type="dna")+scale_y_continuous(limits = c(0,0.15))
p2<-ggseqlogo(downS$X2,seq_type="dna")+scale_y_continuous(limits = c(0,0.15))
p1
ggsave("figures/control.upstream.motif.pdf",width =6.52,height =1.83)
p2
ggsave("figures/control.downstream.motif.pdf",width =6.52,height =1.83)


# cancer gene bed ---------------------------------------------------------

df <- read_tsv("cancerGeneList.tsv")
db <- read_tsv("../乳腺癌/bedFile/dbs/hg38.coding.bed", col_names = F)
names(db) <- c("chr", "Start", "End", "gene")

db_new <- db %>%
  filter(gene %in% (df$`Hugo Symbol`))

write_tsv(db_new, "cancer_gene_bedintervals.bed", col_names = F)

readin <- function(x) {
  df <- read_tsv(stringr::str_glue("allgenes.{x}.bed"), col_names = F)
  df$sample <- x
  df
}

fin <- do.call("rbind", lapply(c("JXM-11", "JXM-13", "JXM-16", "JXM-1", "JXM-23", "JXM-3", "JXM-4", "JXM-5", "JXM-7", "JXM-8", "JXM-9"), readin))
fin <- fin %>%
  select(9, 1, 2, 3, 8) %>%
  setNames(c("sample", "chrom", "start", "end", "oncogenes"))
write_tsv(fin, "all_cancerGenes.eccDNA.bed")

####################################
preprocess<-function(x){
  df<-read_tsv(stringr::str_glue("WGS/CNVkit/{x}.cs.rmdup.sort.anno.cns"))
  df%>%
    select(gene,log2)%>%
    tidyr::separate_rows(gene,sep = ",")%>%
    filter(gene!="-")%>%
    arrange(gene,desc(log2))%>%
    setNames(c("Geneid",paste0(x,"T")))%>%
    group_by(Geneid)%>%
    slice_head(n=1)%>%
    distinct(Geneid,.keep_all = T)
}



CNVs<-Reduce(function(x,y){full_join(x,y,by="Geneid")},lapply(c("J1","J3","J12"),preprocess))
saveRDS(CNVs,"figures/CNV_mat.RDS")

mat_ecc<-ecc_mat%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.data.frame()


mat_cnv<-CNVs%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.data.frame()


sourceOmics<-mat_cnv
targetOmics<-mat_ecc

#####process chrome location###


chrome_sourceOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
                        "12","13","14","15","16","17","18","19","20",
                        "21","22","X","Y")

chrome_targetOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
                        "12","13","14","15","16","17","18","19","20",
                        "21","22","X","Y")

#######Extract sub list#########
genelocate<-read_tsv("../Bladder/dbs/hg38.coding.bed",col_names = F)
genelocate<-genelocate%>%
  select(4,1,2,3)%>%
  mutate(X1=stringr::str_remove(X1,"chr"))%>%
  setNames(c("Symbol","chrom","start","end"))%>%
  as.data.frame()

genelocate_sourceOmics <- genelocate[genelocate[,2] %in%
                                       chrome_sourceOmics,]
genelocate_targetOmics <- genelocate[genelocate[,2] %in%
                                       chrome_targetOmics,]

intG <- intersect(rownames(targetOmics),genelocate_targetOmics[,1])

targetOmics <- targetOmics[intG,]

source_gene <- rownames(sourceOmics)
source_gene_locate <-intersect(unique(genelocate_sourceOmics[,1]),source_gene)
source_gene <- sourceOmics[source_gene_locate,]
genelocate_sourceOmics <- genelocate_sourceOmics[genelocate_sourceOmics[,1] %in% source_gene_locate,]
genelocate_targetOmics <- genelocate_targetOmics[genelocate_targetOmics[,1] %in% intG,]
###Calculate the correlation between cna and other omics data######
corrArray <-calculateCorForTwoMatrices(source_gene,targetOmics,0.01)
saveRDS(corrArray,"final.ecc2cna.RDS")

## functions#############
calculateChromLength <- function(chromLength,selectedChrom,genelocate){
  chromLength <- chromLength[chromLength[,1] %in% selectedChrom,,drop=FALSE]
  
  if(length(selectedChrom)==1){
    x <- 0
  }else{
    x <- c(0,chromLength[1:(nrow(chromLength)-1),2])
  }
  chromLength[,3] <- cumsum(as.numeric(x))
  chromLength[,4] <- cumsum(as.numeric(chromLength[,2]))
  
  genelocate <- cbind(genelocate,0,0)
  
  colnames(genelocate)[5:6] <- c("finalstart","finalend")
  
  for(i in c(1:nrow(genelocate))){
    chr <- genelocate[i,2]
    s <- genelocate[i,3]
    e <- genelocate[i,4]
    cs <- chromLength[chromLength[,1]==chr,3]
    genelocate[i,5] <- s+cs
    genelocate[i,6] <- e+cs
  }
  re <- list(chromLength=chromLength,genelocate=genelocate)
  return(re)
}
plotHeatMap <- function(corrArray,genelocate_sourceOmics,
                        chromLength_sourceOmics,genelocate_targetOmics,chromLength_targetOmics,
                        sourceOmicsName,targetOmicsName,dim=1){
  
  allChromlen_sourceOmics <-
    chromLength_sourceOmics[nrow(chromLength_sourceOmics),4]
  allChromlen_targetOmics <-
    chromLength_targetOmics[nrow(chromLength_targetOmics),4]
  
  if(dim==1){
    par(mar=c(4,4,4,0))
  }else{
    par(mar=c(0,4,4,0))
  }
  
  p <- which(corrArray!=0,arr.ind=TRUE)
  allcnagene <- rownames(corrArray)
  allovgene <- colnames(corrArray)
  
  la <- 1
  for(i in c(1:nrow(p))){
    
    cnag <- allcnagene[p[i,1]]
    ovg <- allovgene[p[i,2]]
    cnagp <- genelocate_sourceOmics[genelocate_sourceOmics[,1]==cnag,5]
    ovgp <- genelocate_targetOmics[genelocate_targetOmics[,1]==ovg,5]
    
    if(length(cnagp)==0 || length(ovgp)==0){
      next
    }
    
    cov <- corrArray[cnag,ovg]
    color <- ifelse(cov>0,"#E63126","#0932E3")
    
    if(la==1){
      if(dim==1){
        plot(cnagp,ovgp,main=paste(sourceOmicsName,"-", targetOmicsName,
                                   " correlation",sep=""), xlim=c(0,allChromlen_sourceOmics),
             ylim=c(0,allChromlen_targetOmics),xaxt="n",yaxt="n",
             frame.plot=FALSE,xlab=paste(sourceOmicsName,
                                         " chromosomal location",sep=""),ylab=paste(targetOmicsName,
                                                                                    " chromosomal location",sep=""),pch=20,col=color,cex=0.2)
        axis(side=1,at=(chromLength_sourceOmics[,4]-
                          chromLength_sourceOmics[,2]/2),
             labels=chromLength_sourceOmics[,1])
      }else{
        plot(cnagp,ovgp,main=paste(sourceOmicsName,"-", 
                                   targetOmicsName," correlation",sep=""),
             xlim=c(0,allChromlen_sourceOmics),
             ylim=c(0,allChromlen_targetOmics),xaxt="n",yaxt="n",
             frame.plot=FALSE,ylab=paste(targetOmicsName," chromosomal
location",sep=""),
             xlab="",pch=20,col=color,cex=0.2)
      }
      axis(side=2,at=(chromLength_targetOmics[,4]-
                        chromLength_targetOmics[,2]/2),
           labels=chromLength_targetOmics[,1])
      
      #abline(h=c(0,chromLength_targetOmics[,4]),v=c(0,chromLength_sourceOmics[,4]),
      #       col="gray",lty=3)
      la <- la+1
    }else{
      for(u in seq_len(length(cnagp))){
        for(v in seq_len(length(ovgp))){
          points(cnagp[u],ovgp[v],pch=20,col=color,cex=0.2)
        }
      }
    }
  }
}
plotSummaryBar <- function(corrArray,chromLength_sourceOmics,
                           genelocate_sourceOmics,sourceOmicsName){
  allChromlen <- chromLength_sourceOmics[nrow(chromLength_sourceOmics),4]
  haha<-tibble(gene=rownames(corrArray),
               numP=apply(corrArray,1,function(x){sum(x!=0)}))
  raw_gene_coord<-genelocate_sourceOmics%>%
    as_tibble()%>%
    select(1,5)%>%
    rename(gene=Symbol)
  haha<-haha%>%
    left_join(raw_gene_coord,by="gene")
  par(mar=c(4,4,0,0))
  plot(0,0,xlim=c(0,allChromlen),ylim=c(0,2000),type="n",xaxt="n", ##ylim=c(-1583,503)
       frame.plot=FALSE,xlab=paste(sourceOmicsName," chromosomal
      location",sep=""),
       ylab="Number of significant \ncorrelations")
  axis(side=1,at=(chromLength_sourceOmics[,4]-
                    chromLength_sourceOmics[,2]/2),
       labels=chromLength_sourceOmics[,1])
  axis(side=2,at=NULL,labels=TRUE)
  abline(v=c(0,chromLength_sourceOmics[,4]),col="gray",lty=3)
  points(haha$finalstart,haha$numP,cex=0.2,type="h",col="#E63126")
  #points(haha$finalstart,haha$numN,cex=0.2,type="h",col="#0932E3")
}
#####################
##Calculate the location of genes in the heatmap
chromLength <- read_tsv("../../qd-ECC4/S/CircleSeq/ECC_report/DataBase/hg38.chromo.size",col_names = F)
chromLength<-chromLength%>%
  filter(X1%in%c(paste0("chr",seq(1,22)),"chrX","chrY"))%>%
  mutate(X1=stringr::str_remove(X1,"chr"))%>%
  setNames(c("V1","V2"))
chromLength$V1<-factor(chromLength$V1,levels = c(seq(1,22),"X","Y"))
chromLength<-chromLength%>%
  arrange(V1)%>%
  as.data.frame()

re <-calculateChromLength(chromLength,chrome_sourceOmics,genelocate_sourceOmics)
genelocate_sourceOmics <- re$genelocate
chromLength_sourceOmics <- re$chromLength
re <-calculateChromLength(chromLength,chrome_targetOmics,genelocate_targetOmics)
genelocate_targetOmics <- re$genelocate
chromLength_targetOmics <- re$chromLength
###plot
png("ecc_CNA_heatmap.png",width = 3.97,height =4.22,res=300,units = "in",type = "cairo")
plotHeatMap(corrArray,genelocate_sourceOmics,chromLength_sourceOmics,
            genelocate_targetOmics,chromLength_targetOmics,"eccDNA",
            "CNA",dim=1)
dev.off()

library(circlize)
## data

df<-read_tsv("WGS/all_SVs/format.J1.sv.delly.txt")
df_f<-df%>%
  filter(filter=="PASS")
df2<-read_tsv("WGS/all_SVs/format.J14.sv.delly.txt")
df2_f<-df2%>%
  filter(filter=="PASS")

get_numbers<-function(x){
  df<-read_tsv(stringr::str_glue("WGS/all_SVs/format.{x}.sv.delly.txt"))%>%
    nrow()
  #df%>%
  #  filter(filter=="PASS")%>%
  #  nrow()
}

get_numbers2<- function(x) {
  df<-read_tsv(stringr::str_glue("WGS/all_SVs/format.{x}.sv.delly.txt"))
  df<-df%>%
    filter(filter=="PASS")
  part1<-df%>%
    select(1:3,7)%>%
    setNames(c("chrom","start","end","id"))
  part2<-df%>%
    select(4:7)%>%
    setNames(c("chrom","start","end","id"))
  tmp<-bind_rows(part1,part2)%>%
    group_by(chrom)%>%
    summarise(num=n_distinct(id))
  tmp$sample<-x
  tmp
}


fin<-tibble(sample=c("J1","J3","J12","J14","J16","J25"),numbers = sapply(c("J1","J3","J12","J14","J16","J25"),get_numbers))
fin<-do.call("bind_rows",lapply(c("J1","J3","J12","J14","J16","J25"),get_numbers2))
bed1 <- df%>%select(1:3)
bed2 <- df%>%select(4:6)
circos.initializeWithIdeogram(species = "hg38")
circos.genomicLink(bed1, bed2)

### split reads
get_split_ratio<-function(x){
  df1<-read_tsv(stringr::str_glue("split_reads/{x}.uniq_reads.txt"),col_names = F)
  total_reads<-sum(df1$X2)
  df2<-read_tsv(stringr::str_glue("split_reads/{x}.uniq__split_reads.txt"),col_names = F)
  total_splits<-sum(df2$X2)
  tibble(sample=x,number=total_splits/total_reads)
}

fin<-do.call("bind_rows",lapply(c("J1","J3","J12","J14","J16","J25"),get_split_ratio))

## chrom
get_split_ratio<-function(x){
  df1<-read_tsv(stringr::str_glue("split_reads/{x}.uniq_reads.txt"),col_names = F)
  total_reads<-sum(df1$X2)
  df2<-read_tsv(stringr::str_glue("split_reads/{x}.uniq__split_reads.txt"),col_names = F)
  tmp<-df2%>%
    mutate(split_ratio=X2/total_reads)
  tmp$sample<-x
  tmp
}
fin2<-do.call("bind_rows",lapply(c("J1","J3","J12","J14","J16","J25"),get_split_ratio))
fin2<-fin2%>%select(-X2)%>%tidyr::pivot_wider(names_from = "sample",values_from = "split_ratio")
write_tsv(fin2,"split_reads/chrom_split_ratio.tsv")

###CNV
get_numbers<-function(x){
  df<-read_tsv(stringr::str_glue("WGS/CNVkit/{x}.cs.rmdup.sort.call.cns"))
  df<-df%>%
    filter(cn!=2)%>%
    group_by(chromosome)%>%
    summarise(count=n())
  df$sample<-x
  df
}

fin<-do.call("bind_rows",lapply(c("J1","J3","J12"),get_numbers))

### gistic2
scores<-read_tsv("WGS/CNVkit/scores.gistic")
scores<-scores%>%
  select(Chromosome,Start,End,`G-score`,Type)%>%
  mutate(`G-score`=if_else(Type=="Amp",`G-score`,-`G-score`))%>%
  setNames(c("seqnames","start","end","score","type"))

plot<-makeGRangesFromDataFrame(scores,keep.extra.columns = T)
seqlevelsStyle(plot)<-"UCSC"
sqinfos<-seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
sqinfos<-keepSeqlevels(sqinfos,paste0("chr",1:22))
seqinfo(plot)<-sqinfos

plot2<-transformToGenome(plot,space.skip=0)
bks<-getXScale(plot2)

png("figures/cna.plot.png",width = 7.31,height = 3.36,res = 300,units = "in")
ggplot(plot2)+
  geom_area(aes(x=.start,y=score,fill=type))+
  scale_x_continuous(breaks = bks$breaks)+
  scale_fill_manual(values = c("Amp"="#C5362C","Del"="#4475A7"))+
  theme_clear(axis.ticks.x = TRUE,axis.line.color = "black")+
  theme(axis.text.x = element_text(angle=90))
dev.off()

