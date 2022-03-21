
load("TARGET_AML_HTSeq-Count.Normalized.RData")

library(data.table)
library(dplyr)


library(readxl)
normal_expression <- read_xlsx("/Users/yinjie/Project/GTEx/All.Normal.CDK2vsGAPDH.xlsx")
blood_sample <- read_xlsx("/Users/yinjie/Project/GTEx/GTEx_Blood_Normal.xlsx")

blood_normal_expression <- filter(normal_expression, Description %in% blood_sample$SAMPID)

blood_normal_expression_merge <- bind_cols(blood_normal_expression,blood_sample)
blood_normal_expression_final <- select(blood_normal_expression_merge, c("SAMPID","GAPDH","CDK2","SMTSD"))

rm(blood_normal_expression,blood_sample)
save.image("GTEx_BloodNormal.RData")

AML <- read_xlsx("AML_CDK2vsGAPDH.xlsx")
All_AML_Normal <- bind_rows(blood_normal_expression_final,AML)


#ggplot

library(ggpubr)

ggboxplot(All_AML_Normal, x = "SMTSD", y = "CDK2", width = 0.8)

All_AML_Blood <-  filter(All_AML_Normal,SMTSD %in% c("Whole Blood","AML"))



p <- ggboxplot(All_AML_Blood, x = "SMTSD", y = "CDK2",
               color = "SMTSD",
               add = "jitter")
p + stat_compare_means()

All_AML_Blood_ratio <- mutate(All_AML_Blood,Ratio=log2(CDK2/GAPDH))

my_comparison <- list(c("AML","Whole Blood"))
relevel
#reorder(SMTSD, CDK2, FUN = median)
p_tpm <- ggplot(All_AML_Blood,aes(x=reorder(SMTSD, CDK2, FUN = median),y=log2(CDK2+1),fill=SMTSD)) +
  geom_violin() +
  scale_fill_manual(values=c('lightsalmon','steelblue')) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
  geom_boxplot(width=0.1, fill="white") +
  theme(axis.text.x = element_blank()) +
  theme_pubr() +
  theme(legend.position="none") +
  xlab("AML vs GTEx") + ylab("CDK2 log2(TPM+1)") +
  stat_compare_means(aes(label=..p.format..),comparisons = my_comparison)
p_tpm

shapiro.test(All_AML_Blood$CDK2)
shapiro.test(All_AML_Normal$CDK2)
shapiro.test()

re <- filter(recurred,Type=="Recurrent")
shapiro.test(re$CDK2)
pr <- filter(recurred,Type!="Recurrent")
shapiro.test(pr$CDK2)
p_ratio <- ggplot(All_AML_Blood_ratio,aes(x=reorder(SMTSD, CDK2, FUN = median),y=Ratio,fill=SMTSD)) +
  geom_violin() +
  scale_fill_manual(values=c('lightsalmon','steelblue')) +
  geom_boxplot(width=0.1, fill="white") +
  theme(axis.text.x = element_blank()) +
  theme_pubr() +
  theme(legend.position="none") +
  xlab("AML vs GTEx") + ylab("log2(CDK2/GAPDH)") +
  stat_compare_means(aes(label=..p.format..),comparisons = my_comparison)
  #stat_compare_means(aes(label=..p.signif..),comparisons = my_comparison)
p_ratio


save.image("AML_Plot.RData")
#load clinical information


clinical <- read.csv("/Users/yinjie/Project/TARGET_AML/clinical.cases_selection.2020-01-14/clinical.tsv",sep="\t",header=T)
clinical_sub <- select(clinical,c("case_id","days_to_last_follow_up","gender","vital_status","age_at_diagnosis"))
write.table(clinical_sub,"Clinical.information.tsv",quote = FALSE,col.names = T,sep="\t")



All_AML_Blood_Clincal <- filter(All_AML_Blood_ratio, SAMPID %in% clinical_info$case_id)

fileid_caseid <- read_xlsx("FileID_CaseID.xlsx")
caseid_uuid <- read_xlsx("CaseID_UUID.xlsx")

merge <- full_join(fileid_caseid,caseid_uuid)

colnames(clinical_sub) <- c("UUID","Time","Gender","Censored","Age")

merge_all <- full_join(merge,clinical_sub)



merge_data <- merge_all
merge_data$`File Name` <- gsub("(.*).htseq.counts.gz","\\1",merge_all$`File Name`)

merge_data$`SAMPID` <- merge_data$`File Name`
data_all <- full_join(merge_data,AML)
write.table(data_all,"AML_with_all_clinicalinformation.tsv",quote=F,col.names = T,row.names = F,sep="\t")
data4clin <- select(data_all,c("Case ID","SAMPID","Gender","Age","Censored","Time","Sample Type","CDK2","GAPDH"))
cutoff <- median(data4clin$CDK2)
data3clin <- mutate(data4clin,Label=ifelse(data4clin$CDK2 >= cutoff,"High","Low"))
write.table(data3clin,"AML_with_all_needinformation.tsv",quote=F,col.names = T,row.names = F,sep="\t")


library(survival)
library(survminer)



data3clin <- mutate(data3clin,Status=ifelse(data3clin$Censored=="Dead",TRUE,FALSE))
data3clin.df <- as.data.frame(data3clin)

data_sur <- select(data3clin,c("Case ID","Time","Status","Label")) %>% unique()

write.table(data_sur,"Survival.All.tsv",sep="\t",quote = F)

library(readxl)


data_sur2 <- read_xlsx("Survial_Dedup.xlsx")


fit <- survfit(Surv(Time/30, Status) ~ Label, data = data_sur2)

#sur <- ggsurvplot(fit,data=data3clin.df,size = 1,pval = TRUE,risk.table = TRUE,ggtheme = theme_pubr())
sur <- ggsurvplot(fit,data=data_sur2,size = 1,pval = TRUE,risk.table = TRUE,ggtheme = theme_pubr(),xlab = "Time (Months)")

sur

recurred <- read_xlsx("RecurredSample.xlsx")
my_comparison2 <- list(c("Primary","Recurrent"))
p_recur_tpm <- ggplot(recurred, aes(x= Type,y = log2(CDK2+1),fill=Type)) +
  geom_violin() +
  scale_fill_manual(values=c('steelblue','lightsalmon')) +
  geom_boxplot(width=0.1, fill="white") +
  #geom_boxplot(fill = c('steelblue','lightsalmon')) +
  geom_point(size = 1) +
  theme(legend.position="none") +
  geom_line(aes(group = `Case ID`),colour="gray") +
  theme(axis.text.x = element_blank()) +
  theme_pubr() +xlab("AML Primary vs. Recurrent ( n = 31)") + ylab("CDK2 log2(TPM+1)") + 
  theme(legend.position="none") +
  #ylim(2.45,3.45) +
  stat_compare_means(aes(label = ..p.format..),label.x = 1.5, label.y = 6,paired = TRUE,comparisons = my_comparison2)
p_recur_tpm


p_recur_ratio <- ggplot(recurred, aes(x= Type,y = log(CDK2/GAPDH,base=2),fill=Type)) +
  geom_violin() +
  #geom_boxplot(fill = c('steelblue','lightsalmon')) +
  scale_fill_manual(values=c('steelblue','lightsalmon')) +
  geom_boxplot(width=0.1, fill="white") +
  geom_point(size = 1) +
  geom_line(aes(group = `Case ID`),colour="gray") +
  theme(axis.text.x = element_blank()) +
  theme_pubr() +xlab("AML Primary vs. Recurrent ( n = 31)") + ylab("log2(CDK2/GAPDH)") + 
  #ylim(-2.5,-11) +
  theme(legend.position="none") +
  stat_compare_means(aes(label = ..p.format..),label.x = 1.5, label.y = -4,paired = TRUE,comparisons = my_comparison2)
p_recur_ratio

fig<-ggarrange(p_tpm,p_ratio,p_recur_tpm,p_recur_ratio,sur$plot,sur$table,labels = c("A","B","C","D","E","F"),ncol=2,nrow=3)


ggsave("Expression_Figure4.pdf", plot = fig,device = "pdf", 
       scale = 1, width = 8, height = 10, units = c("in"),
       dpi = 300 )


library("gridExtra")
library("cowplot")

save.image("AML_Plot.RData")


load("AML_Plot.RData")


