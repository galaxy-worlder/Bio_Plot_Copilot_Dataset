#  Loading libraries

suppressPackageStartupMessages({
  library(MAST);  library(immunarch);  library(plyr); library(data.table);  library(reshape);  library(scater); library(ComplexHeatmap);
  library(edgeR);  library(ROCR);  library(scRNA.seq.funcs);  library(gridExtra);  library(NMF);  library(RColorBrewer);  library(ggplot2);
  library(dplyr);  library(MASS);  library(pscl);  library(ape);  library(factoextra);  library(NbClust);  library(igraph);  library(stringdist);
  library(tcR);  library(purrr);  library(tidyr);  library(dendextend);  library(doParallel);  library(foreach);  library(dtplyr);  library(tidyverse);
  library(ggpubr);  library(corrplot);  library(e1071);  library(bootstrap);  library(DAAG);  library(RColorBrewer);  library(ppcor);  library(factoextra);
  library(mixOmics);  library(glmnet);  library(psych);  library(caret);  library(randomForest);  library(PRROC);  library(viridis); library(grid);
  library(ggplotify);  library(ROCR);   library(jpeg);   library(lme4);   library(glmmTMB);   library(reshape);   library(protr);  library(msa);  library(Biostrings);  library(ggseqlogo)
})

# Analyses

df.is<-readRDS("Data2analyze/IsoPerc.rds") # This file can be downloaded through the website of the Rheumatology Research Group.
isotypes<-as.character(unique(df.is$Isotype))
df.pvals<-data.frame("Isotype"=rep("NA",length(isotypes)), "Beta"=rep("NA",length(isotypes)), "Pval"=rep("NA",length(isotypes)), stringsAsFactors=F)
count=0
for (i in isotypes) {
  count=count+1
  df2use<-df.is[df.is$Isotype==i,]
  df2use$IMID.num<-as.numeric(df2use$IMID);  df2use$IMID.num[df2use$IMID.num==1]<-0;  df2use$IMID.num[df2use$IMID.num==2]<-1
  df2use$Reads.Chain<-as.numeric(df2use$Reads.Chain)
  model<-glm(IMID.num~Perc+Age+Gender, data=df2use, family="binomial")
  beta<-coef(summary(model))[2,'Estimate']
  pval<-coef(summary(model))[2,'Pr(>|z|)']
  df.pvals[count,]<-c(i,beta,pval)
}
df.pvals.final<-df.pvals[order(df.pvals$Pval),]
df.pvals.final$FDR<-p.adjust(df.pvals.final$Pval,method='fdr')
p1<-ggboxplot(df.is, x="Isotype", y="Perc", color="IMID", add="jitter", fill="IMID") +
  labs(x="", y="\nIsotype Use (%)\n") + scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=c("#D8EBD8","#F69191")) + theme(legend.key.size=unit(1.5, "cm")) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  annotate(geom="text", x=1, y=max(df.is$Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IgA","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=2, y=max(df.is$Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IgD","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=3, y=max(df.is$Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IgE","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=4, y=max(df.is$Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IgG","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=5, y=max(df.is$Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IgM","Pval"]), format="e", digits=2)), color="black", size=4)
jpeg("Output_Data/IsotypePercentage_RAvsCTRL.jpeg", width=4500, height=3000, res=300)
plot(p1)
dev.off()

phenos<-c("RespGM","ActB","ACPA","RF")
isotypes<-as.character(unique(df.is$Isotype))
df.pvals<-data.frame("Phenotype"=rep("NA",length(isotypes)*length(phenos)), "Isotype"=rep("NA",length(isotypes)), "Beta"=rep("NA",length(isotypes)), "Pval"=rep("NA",length(isotypes)), stringsAsFactors=F)
count=0
for (p in phenos) {
  df1<-df.is[df.is$IMID=="RA",]
  for (i in isotypes) {
    count=count+1
    df2<-df1[df1$Isotype==i,]
    meta<-read.csv("Data2analyze/Metadata.csv", sep="\t");   rownames(meta)<-meta$Sample;
    df3<-merge(df2,meta,by="Sample")
    if (p=="RespGM") {pheno2keep="Resp"} else {pheno2keep=p}
    df2use<-df3[,c("Sample","Reads.Chain","Age.x","Gender.x",pheno2keep,"Perc")]
    colnames(df2use)[5]<-"Pheno"
    if (p=="RespGM") {df2use$Pheno[df2use$Pheno=="MOD"]<-"GOOD"} else {df2use<-df2use}
    df2use<-df2use[!(is.na(df2use$Pheno)),]
    df2use$Pheno<-as.character(df2use$Pheno)
    df2use$Pheno[df2use$Pheno==unique(df2use$Pheno)[1]]<-"0"
    df2use$Pheno[df2use$Pheno==unique(df2use$Pheno)[2]]<-"1"
    df2use$Pheno<-as.numeric(df2use$Pheno)
    model<-glm(Pheno~Perc+Age.x+Gender.x, data=df2use, family="binomial")
    beta<-coef(summary(model))[2,'Estimate']
    pval<-coef(summary(model))[2,'Pr(>|z|)']
    df.pvals[count,]<-c(p,i,beta,pval)
  }
}
df.pvals.final<-df.pvals[order(df.pvals$Pval),]
df.pvals.final$FDR<-p.adjust(df.pvals.final$Pval, method='fdr')
write.csv(df.pvals.final[,c(1:4)],"Output_Data/IsotypePercentage_RApheno.csv", row.names=F)

df.csr<-readRDS("Data2analyze/Class_Switching.rds") # This file can be downloaded through the website of the Rheumatology Research Group.
csrs<-as.character(unique(df.csr$CSR.Name))
df.pvals<-data.frame("CSR"=rep("NA",length(csrs)), "Estimate"=rep("NA",length(csrs)), "Pval"=rep("NA",length(csrs)), stringsAsFactors=F)
count=0
for (i in csrs) {
  count=count+1
  df2use<-df.csr[df.csr$CSR.Name==i,]
  df2use$IMID.num<-as.numeric(df2use$IMID);  #df2use$IMID.num[df2use$IMID.num==1]<-0;  df2use$IMID.num[df2use$IMID.num==2]<-1
  df2use$Reads.Chain<-as.numeric(df2use$Reads.Chain)
  model<-glm(IMID.num~CSR.Index+Reads.Chain+Age+Gender, data=df2use, family="Gamma")
  beta<-coef(summary(model))[2,'Estimate']
  pval<-coef(summary(model))[2,'Pr(>|t|)']
  df.pvals[count,]<-c(i,beta,pval)
}
df.pvals.final<-df.pvals[order(df.pvals$Pval),]
df.pvals.final$FDR<-p.adjust(df.pvals.final$Pval,method='fdr')
p1<-ggboxplot(df.csr, x="CSR.Name", y="CSR.Index", color="IMID", add="jitter", fill="IMID") +
  labs(x="", y="\nCSR Index\n") + scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=c("#D8EBD8","#F69191")) + theme(legend.key.size=unit(1.5, "cm")) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  annotate(geom="text", x=1, y=max(df.csr$CSR.Index)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$CSR=="IgA-IgDM","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=2, y=max(df.csr$CSR.Index)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$CSR=="IgA-IgE","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=3, y=max(df.csr$CSR.Index)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$CSR=="IgA-IgG","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=4, y=max(df.csr$CSR.Index)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$CSR=="IgDM-IgE","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=5, y=max(df.csr$CSR.Index)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$CSR=="IgDM-IgG","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=6, y=max(df.csr$CSR.Index)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$CSR=="IgE-IgG","Pval"]), format="e", digits=2)), color="black", size=4)
jpeg("Output_Data/CSR_RAvsCTRL.jpeg", width=5000, height=3000, res=300)
plot(p1)
dev.off()

phenos<-c("RespGM","ActB","ACPA","RF")
csrs<-as.character(unique(df.csr$CSR.Name))
df.pvals<-data.frame("Phenotype"=rep("NA",length(csrs)*length(phenos)), "CSR"=rep("NA",length(csrs)*length(phenos)), "Estimate"=rep("NA",length(csrs)*length(phenos)), "Pval"=rep("NA",length(csrs)*length(phenos)), stringsAsFactors=F)
count=0
for (p in phenos) {
  df1<-df.csr[df.csr$IMID=="RA",]
  for (i in csrs) {
    count=count+1
    df2<-df1[df1$CSR.Name==i,]
    meta<-read.csv("Data2analyze/Metadata.csv", sep="\t");   rownames(meta)<-meta$Sample;
    df3<-merge(df2,meta,by="Sample")
    if (p=="RespGM") {pheno2keep="Resp"} else {pheno2keep=p}
    df2use<-df3[,c("Sample","Reads.Chain",pheno2keep,"Age.x","Gender.x","CSR.Index")]
    colnames(df2use)[3]<-"Pheno"
    if (p=="RespGM") {df2use$Pheno[df2use$Pheno=="MOD"]<-"GOOD"} else {df2use<-df2use}
    df2use<-df2use[!(is.na(df2use$Pheno)),]
    df2use$Pheno<-as.character(df2use$Pheno)
    df2use$Pheno[df2use$Pheno==unique(df2use$Pheno)[1]]<-"1"
    df2use$Pheno[df2use$Pheno==unique(df2use$Pheno)[2]]<-"2"
    df2use$Pheno<-as.numeric(df2use$Pheno)
    model<-glm(Pheno~CSR.Index+Reads.Chain+Age.x+Gender.x, data=df2use, family="Gamma")
    beta<-coef(summary(model))[2,'Estimate']
    pval<-coef(summary(model))[2,'Pr(>|t|)']
    df.pvals[count,]<-c(p,i,beta,pval)
  }
}
df.pvals.final<-df.pvals[order(df.pvals$Pval),]
df.pvals.final$FDR<-p.adjust(df.pvals.final$Pval, method='fdr')
write.csv(df.pvals.final[,c(1:4)],"Output_Data/CSR_RApheno.csv", row.names=F)


df.shm<-readRDS("Data2analyze/SHM.rds") # This file can be downloaded through the website of the Rheumatology Research Group.
shms<-as.character(unique(df.shm$Isotype))
df.pvals<-data.frame("Isotype"=rep("NA",length(shms)), "Beta"=rep("NA",length(shms)), "Pval"=rep("NA",length(shms)), stringsAsFactors=F)
count=0
for (i in shms) {
  count=count+1
  if (i=="IgE" | i=="IgD" | i=="IgG34" | i=="IgA" | i=="IgG12" | i=="IgM") {df2use<-df.shm[(df.shm$Isotype==i & df.shm$Chain=="IGH"),]}
  if (i=="IGL") {df2use<-df.shm[(df.shm$Isotype==i & df.shm$Chain=="IGL"),]}
  if (i=="IGK") {df2use<-df.shm[(df.shm$Isotype==i & df.shm$Chain=="IGK"),]}
  df2use$IMID.num<-as.numeric(df2use$IMID);  df2use$IMID.num[df2use$IMID.num==1]<-0;  df2use$IMID.num[df2use$IMID.num==2]<-1
  model<-glm(IMID.num~SHM.Perc+Age+Gender, data=df2use, family="binomial")
  beta<-coef(summary(model))[2,'Estimate']
  pval<-coef(summary(model))[2,'Pr(>|z|)']
  df.pvals[count,]<-c(i,beta,pval)
}
df.pvals.final<-df.pvals[order(df.pvals$Pval),]
df.pvals.final$FDR<-p.adjust(df.pvals.final$Pval,method='fdr')

df.shm.igh<-df.shm[(df.shm$Isotype!="IGL" & df.shm$Isotype!="IGK" & df.shm$Chain=="IGH"),]
df.shm.igl<-df.shm[(df.shm$Isotype=="IGL" & df.shm$Chain=="IGL"),]
df.shm.igk<-df.shm[(df.shm$Isotype=="IGK" & df.shm$Chain=="IGK"),]
df2plot<-rbind(df.shm.igh,df.shm.igl,df.shm.igk)
p1<-ggboxplot(df2plot, x="Isotype", y="SHM.Perc", color="IMID", add="jitter", fill="IMID") +
  labs(x="", y="\nSHM (%)\n") + scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=c("#D8EBD8","#F69191")) + theme(legend.key.size=unit(1.5, "cm")) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  annotate(geom="text", x=1, y=max(df2plot$SHM.Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IgA","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=2, y=max(df2plot$SHM.Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IgM","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=3, y=max(df2plot$SHM.Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IgD","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=4, y=max(df2plot$SHM.Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IgG34","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=5, y=max(df2plot$SHM.Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IgG12","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=6, y=max(df2plot$SHM.Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IgE","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=7, y=max(df2plot$SHM.Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IGK","Pval"]), format="e", digits=2)), color="black", size=4) +
  annotate(geom="text", x=8, y=max(df2plot$SHM.Perc)*1.1, label=paste0("P=",formatC(as.numeric(df.pvals[df.pvals$Isotype=="IGL","Pval"]), format="e", digits=2)), color="black", size=4)
jpeg("Output_Data/SHM_RAvsCTRL.jpeg", width=5500, height=3000, res=300)
plot(p1)
dev.off()

df.mut<-readRDS("Data2analyze/MutatedIgDM.rds") # This file can be downloaded through the website of the Rheumatology Research Group.
model<-glm(IMID~Perc.IgDM.Mut+Age+Gender, data=df.mut, family="binomial")
beta<-coef(summary(model))[2,'Estimate']
pval<-coef(summary(model))[2,'Pr(>|z|)']
df.pvals.final<-df.pvals[order(df.pvals$Pval),]
df.pvals.final$FDR<-p.adjust(df.pvals.final$Pval,method='fdr')
p1<-ggboxplot(df.mut, x="IMID", y="Perc.IgDM.Mut", color="IMID", add="jitter", fill="IMID") +
  labs(x="", y="\nIgD/M Mutated (%)\n") + scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=c("#D8EBD8","#F69191")) + theme(legend.key.size=unit(1.5, "cm")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  annotate(geom="text", x=1.5, y=101, label="P=4.31e-03", color="black", size=6)
jpeg("Output_Data/MutatedIgDM_RAvsCTRL.jpeg", width=2500, height=3000, res=300)
plot(p1)
dev.off()

phenos<-c("RespGM","ActB","ACPA","RF")
isotypes<-as.character(unique(df.shm$Isotype))
df.pvals<-data.frame("Phenotype"=rep("NA",length(isotypes)*length(phenos)), "Isotype"=rep("NA",length(isotypes)), "Beta"=rep("NA",length(isotypes)), "Pval"=rep("NA",length(isotypes)), stringsAsFactors=F)
count=0
for (p in phenos) {
  df1<-df.shm[df.shm$IMID=="RA",]
  for (i in isotypes) {
    count=count+1
    if (i=="IgE" | i=="IgD" | i=="IgG34" | i=="IgA" | i=="IgG12" | i=="IgM") {df2<-df1[(df1$Isotype==i & df1$Chain=="IGH"),]}
    if (i=="IGL") {df2<-df1[(df1$Isotype==i & df1$Chain=="IGL"),]}
    if (i=="IGK") {df2<-df1[(df1$Isotype==i & df1$Chain=="IGK"),]}
    meta<-read.csv("Data2analyze/Metadata.csv", sep="\t");   rownames(meta)<-meta$Sample;
    df3<-merge(df2,meta,by="Sample")
    if (p=="RespGM") {pheno2keep="Resp"} else {pheno2keep=p}
    df2use<-df3[,c("Sample","Reads.Chain",pheno2keep,"Age.x","Gender.x","SHM.Perc")]
    colnames(df2use)[3]<-"Pheno"
    if (p=="RespGM") {df2use$Pheno[df2use$Pheno=="MOD"]<-"GOOD"} else {df2use<-df2use}
    df2use<-df2use[!(is.na(df2use$Pheno)),]
    df2use$Pheno<-as.character(df2use$Pheno)
    df2use$Pheno[df2use$Pheno==unique(df2use$Pheno)[1]]<-"0"
    df2use$Pheno[df2use$Pheno==unique(df2use$Pheno)[2]]<-"1"
    df2use$Pheno<-as.numeric(df2use$Pheno)
    model<-glm(Pheno~SHM.Perc+Age.x+Gender.x, data=df2use, family="binomial")
    beta<-coef(summary(model))[2,'Estimate']
    pval<-coef(summary(model))[2,'Pr(>|z|)']
    df.pvals[count,]<-c(p,i,beta,pval)
  }
}
df.pvals.final<-df.pvals[order(df.pvals$Pval),]
df.pvals.final$FDR<-p.adjust(df.pvals.final$Pval, method='fdr')
write.csv(df.pvals.final[,c(1:4)],"Output_Data/SHM_RApheno.csv", row.names=F)

phenos<-c("RespGM","ActB","ACPA","RF")
df.pvals<-data.frame("Phenotype"=rep("NA",length(phenos)), "Beta"=rep("NA",length(phenos)), "Pval"=rep("NA",length(phenos)), stringsAsFactors=F)
count=0
for (p in phenos) {
  count=count+1
  meta<-read.csv("Data2analyze/Metadata.csv", sep="\t");   rownames(meta)<-meta$Sample;
  df.mut.ra<-df.mut[df.mut$IMID=="RA",]
  df3<-merge(df.mut.ra,meta,by="Sample")
  if (p=="RespGM") {pheno2keep="Resp"} else {pheno2keep=p}
  df2use<-df3[,c("Sample","Reads.Chain",pheno2keep,"Age.x","Gender.x","Perc.IgDM.Mut")]
  colnames(df2use)[3]<-"Pheno"
  if (p=="RespGM") {df2use$Pheno[df2use$Pheno=="MOD"]<-"GOOD"} else {df2use<-df2use}
  df2use<-df2use[!(is.na(df2use$Pheno)),]
  df2use$Pheno<-as.character(df2use$Pheno)
  df2use$Pheno[df2use$Pheno==unique(df2use$Pheno)[1]]<-"0"
  df2use$Pheno[df2use$Pheno==unique(df2use$Pheno)[2]]<-"1"
  df2use$Pheno<-as.numeric(df2use$Pheno)
  model<-glm(Pheno~Perc.IgDM.Mut+Age.x+Gender.x, data=df2use, family="binomial")
  beta<-coef(summary(model))[2,'Estimate']
  pval<-coef(summary(model))[2,'Pr(>|z|)']
  df.pvals[count,]<-c(p,beta,pval)
}
df.pvals.final<-df.pvals[order(df.pvals$Pval),]
df.pvals.final$FDR<-p.adjust(df.pvals.final$Pval, method='fdr')
write.csv(df.pvals.final[,c(1:3)],"Output_Data/MutIgDM_RApheno.csv", row.names=F)

# Note: the generated figures has been reshaped and merged for publication