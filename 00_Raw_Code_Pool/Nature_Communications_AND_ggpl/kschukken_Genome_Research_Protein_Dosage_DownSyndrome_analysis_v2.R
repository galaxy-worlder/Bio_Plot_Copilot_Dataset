##### Down syndrome gene expression differences, relative to human cancer ####
## 210507 
## Liu et al. 2017 aneuploid down syndrome patient proteomics data
## Author: Klaske Schukken
## Compare chromosome loss/gain difference with protein expression changes. 
## Data from Liu et al. 2017, Nature Communications
## Systematic proteome and proteostasis profiling in human Trisomy 21 fibroblast cells

## Get proteomics data


library('ggplot2')
library('tidyverse')
library(xlsx)
library(readxl)
library(reshape2)
library('BBmisc')

###
#plan: 
## Step 1: get Proteomics  data 
##      down syndrome: several different aneuploid down syndrome fibroblast data relative to control group
##      Data type: Proteomics between 11 down syndrome patient fibroblasts and matched controls
##      log2 transformed difference between aneuploid and controls cells
##      also RNA data from Letourneau et al. between matched twins discordant for chrm 21
##      Protein in aneuploid chromosomes (Depmap)
## Step 2: 2d density plot Protein difference in down syndrome vs cancer. 
## Step 3: 2d density plot RNA difference in down syndrome vs cancer. 



##### Step 1: get Proteomics  data ####
## Compare chromosome loss/gain difference with protein expression changes. 
## Data from Liu et al. 2017, Nature Communications
## Systematic proteome and proteostasis profiling in human Trisomy 21 fibroblast cells

##
## Lui et al data was downloaded and correlated with CCLE human cancer data in excel 
# Protein data
setwd()
Ts21.Prot<- read_xlsx("trisomy 21 comparison.xlsx")

# RNA data
setwd()
# from Letourneau et al. supplementary figure
# RNA log2 transformed transcriptomic data comparing twins discordant for chrm 21
Ts21.RNA<- read_xlsx("41586_2014_BFnature13200_MOESM122_ESM.xlsx")



##### Step 2: 2d density plot Protein difference in down syndrome vs cancer. ####

###Plot difference upon gain in Down Syndrome vs human cancer cell lines
## Plot all aneuploid gene difference: Liu vs. depmap

#scatterplot
ggplot(Ts21.Prot, aes(y=`Ts vs. Ds`, x=`Ts21 vs. WT...4`))+ #dotplot
  geom_point(color="black")+ 
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  #coord_cartesian(ylim=c(-0.5,1), xlim=c(-0.5,1))+
  geom_smooth(method="lm", color="red")+
  ylab("Depmap: Protein difference\nupon chrm 21 gain in cancer")+
  xlab("Liu: Protein difference\nin down syndrome patients")
#4x4
# plot.Liu.Depmap.Chrm21MeanDiff_v2

#density plot
ggplot(Ts21.Prot, aes(y=`Ts vs. Ds`, x=`Ts21 vs. WT...4`))+ #density
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  #coord_cartesian(ylim=c(-0.5,1), xlim=c(-0.5,1))+
  geom_smooth(method="lm", color="red")+
  ylab("Depmap: Protein difference\nupon chrm 21 gain in cancer")+
  xlab("Liu: Protein difference\nin down syndrome patients")
#4x4
# plot.Liu.Depmap.Chrm21MeanDiff_density

#get pearson correlation coefficient and p-value
Ts21Corr<-cor.test(Ts21.Prot$`Ts vs. Ds`, Ts21.Prot$`Ts21 vs. WT...4`, method="pearson")
# Pearson corr= 0.538 
# p-value = 3E-04 
# n=42



##### Step 3: 2d density plot RNA difference in down syndrome vs cancer. #####

Ts21.RNA<- subset(Ts21.RNA, chr=="chr21") #only get RNAs on chrm 21

Ts21.RNA.depmap<- merge(Ts21.RNA, y= CN.Diff.xRNA.yProt.ThreeGroups, 
                        by.x="gene_name" , by.y="RNA_Name")

###Plot difference upon gain in yeast vs human

ggplot(Ts21.RNA.depmap, aes(y=RNA.Diff.Gain, x=logFC))+ #density
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  coord_cartesian(ylim=c(-0.1,0.7), xlim=c(-1,2))+
  geom_smooth(method="lm", color="red")+
  ylab("CCLE: RNA difference\nupon chrm 21 gain in cancer")+
  xlab("Letourneau: RNA difference\nin down syndrome fibroblasts")
#4x4
# plot.Letourneau.Depmap.Chrm21MeanDiff.RNA_density

# paerson correlation coefficient and p-value for RNA part
Ts21Corr<-cor.test(Ts21.RNA.depmap$logFC, Ts21.RNA.depmap$RNA.Diff.Gain, method="pearson")
Ts21Corr
# Pearson corr= 0.294 
# p-value = 4E-03
# n=93

