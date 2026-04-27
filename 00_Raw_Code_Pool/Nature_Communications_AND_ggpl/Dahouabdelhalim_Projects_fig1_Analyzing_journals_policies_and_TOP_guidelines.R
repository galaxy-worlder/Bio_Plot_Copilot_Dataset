# Analysis of the journal's TOP guidelines
# INPUT DATA downloaded from:  https://osf.io/qatkz/  (Jan 2022)
# Results embeded in https://zenodo.org/record/5871968


library(ggplot2)
library(reshape)
library(stringr)

# IMPORT
top=read.csv("C:/Users/ipatarc/Desktop/PR_TOP/top-factor -v33.csv")

# # Policies of interest 
# 


    policies=which(str_detect(colnames(top),"score"))
    
    filtered.j=top[,policies]
    per_score_policy=rbind(apply(filtered.j,2,function(x){length(which(x==1))}),
                          apply(filtered.j,2,function(x){length(which(x==2))}),
                          apply(filtered.j,2,function(x){length(which(x==3))}))
    
    
    # Data.citation.score Data.transparency.score Analysis.code.transparency.score
    # [1,]                 734                     724                              427
    # [2,]                 446                     158                              120
    # [3,]                  12                      24                               36
    # Materials.transparency.score Design...analysis.reporting.guidelines.score
    # [1,]                          181                                          545
    # [2,]                           84                                          142
    # [3,]                            9                                           28
    # Study.preregistration.score Analysis.plan.preregistration.score Replication.score
    # [1,]                         165                                 150               183
    # [2,]                          22                                  22                 2
    # [3,]                           6                                   6               110
    # Registered.reports...publication.bias.score Open.science.badges.score
    # [1,]                                          20                        15
    # [2,]                                           1                       108
    # [3,]                                         181                         0



# Results slide 17 in in https://zenodo.org/record/5871968 - general journal's staistics

#################
#adjusting data for plotting

    df.policies=melt((per_score_policy))
    colnames(df.policies)=c("TOP_Level","TOP_Score","N_journals")
    df.policies$TOP_Score=str_replace(paste0(rep(LETTERS[10:1],each=3),".",
                                  str_replace_all(str_replace(df.policies$TOP_Score,
                                             ".score",""),"\\\\."," ")),"   ","/")
    
# df.policies=df.policies[24:1,]
png("C:/Users/ipatarc/Desktop/PR_TOP/Fig1c_guidelines_across_journals_barplot.png", width = 800,height = 340)
ggplot(data = df.policies)+
  geom_bar(aes(x = N_journals,y=as.factor(TOP_Score),fill = as.factor(TOP_Level)), stat="identity")  +
  scale_fill_discrete(breaks=c("1","2","3"))+
  theme_bw()+labs(y="TOP standard",x="Number of journals") + 
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold"),
                    plot.title=element_text(size=14, face="bold"))+
  
  labs(fill='TOP policy Level')
dev.off()

#getwd()

################################################





# 
# top[top$Journal=="Science",policies]
# top[top$Journal=="Nature",policies]
# top[top$Journal=="Nature Methods",policies]
# top[top$Journal=="eLife",policies]
# top[top$Journal=="Nature Communications",policies]
# top[top$Journal=="PLOS Computational Biology",policies]
# top[top$Journal=="PLOS Biology",policies]
# top[top$Journal=="Proceedings of the National Academy of Sciences",policies]
# top[top$Journal=="Biophysical Journal",policies]
# top[top$Journal=="Cell",policies]
# top[top$Journal=="Cell Reports",policies]
# top[top$Journal=="Genome Research",policies]
# top[top$Journal=="EMBO Journal",policies]
# 
# 
# # Checking journals with specific policy levels
# top$Journal[which(top$Data.citation.score==3)]
# top$Journal[which(top$Data.transparency.score==3)]
# top$Journal[which(top$Analysis.code.transparency.score==3)]
# top$Journal[which(top$Materials.transparency.score==3)]
# top$Journal[which(top$Design...analysis.reporting.guidelines.score==3)]
# top$Journal[which(top$Study.preregistration.score==3)]
# top$Journal[which(top$Analysis.plan.preregistration.score==3)]
# top$Journal[which(top$Replication.score==3)]


# Open.science.badges.score 
# 123 
# Analysis.plan.preregistration.score 
# 178 
# Study.preregistration.score 
# 193 
# Registered.reports...publication.bias.score 
# 202 
# Materials.transparency.score 
# 274 
# Replication.score 
# 295 
# Analysis.code.transparency.score 
# 583 
# Design...analysis.reporting.guidelines.score 
# 715 
# Data.transparency.score 
# 906 
# Data.citation.score 
# 1192 


#getting information how often each policy is present
#idea get all scores as a vector and do summary staistics on the vector

all_policies= (unlist(filtered.j))

summary(all_policies[all_policies!=0])
table(all_policies[all_policies!=0])

# total number of policies
# 1    2    3 
# 3144 1105  412 



# how does statistics look like when different stringency levels compared across 
# guidelines

round(cbind(100*per_score_policy[1,]/(colSums(per_score_policy)),
100*per_score_policy[2,]/(colSums(per_score_policy)),
100*per_score_policy[3,]/(colSums(per_score_policy))),1)
