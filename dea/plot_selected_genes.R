genes = c(
  "Ddit3",
  "Atf4",
  "Atf5",
  "Fgf21",
  "Asns",
  "Hspa9",
  "Hspd1",
  "Ehmt1",
  "Baz2b",
  "Lonp1",
  "Ubl5",
  "Abcb10",
  "Defa5",
  "Clpp",
  "Satb1",
  "Satb2",
  "Kdm7a",
  "Kdm6b"
)

library(MotrpacRatTraining6moData)
library(dplyr)
library(reshape2)
FEATURE_TO_GENE = FEATURE_TO_GENE[!grepl("cluster",FEATURE_TO_GENE$feature_ID),]
FEATURE_TO_GENE = FEATURE_TO_GENE[!grepl("^chr\\d",FEATURE_TO_GENE$feature_ID,perl=T),]
FEATURE_TO_GENE = FEATURE_TO_GENE[FEATURE_TO_GENE$gene_symbol %in% genes,]
length(unique(FEATURE_TO_GENE$gene_symbol))

datasets = list(
  "BLOOD" = TRNSCRPT_BLOOD_DA,
  "ADRNL" = TRNSCRPT_ADRNL_DA,
  "COLON" = TRNSCRPT_COLON_DA,
  "LIVER" = TRNSCRPT_LIVER_DA,
  "HEART" = TRNSCRPT_HEART_DA,
  "SKMGN" = TRNSCRPT_SKMGN_DA,
  "SKMVL" = TRNSCRPT_SKMVL_DA,
  "BAT" = TRNSCRPT_BAT_DA,
  "WATSC" = TRNSCRPT_WATSC_DA
)

datasets = lapply(
  datasets,
  function(x)x[x$feature_ID %in% FEATURE_TO_GENE$feature_ID,]
)

sapply(datasets,dim)

df = do.call(rbind,datasets)
df = merge(df,FEATURE_TO_GENE,by="feature_ID")
df = df[,c("gene_symbol","tissue","sex","logFC","comparison_group")]
df = df[df$comparison_group == "1w",]
dim(df)
df$group = paste(df$tissue,df$sex,sep="_")

m_genes = unique(df$gene_symbol)
m_tissues = sort(unique(df$tissue))
m_male = matrix(NA,nrow=length(m_genes),ncol=length(m_tissues),dimnames = list(m_genes,m_tissues))
m_female = matrix(NA,nrow=length(m_genes),ncol=length(m_tissues),dimnames = list(m_genes,m_tissues))
for(tis in m_tissues){
  currdf = df[df$tissue == tis &df$sex == "male",]
  curr_fcs = tapply(currdf$logFC,currdf$gene_symbol,function(x)x[abs(x)==max(abs(x))][1])
  m_male[names(curr_fcs),tis] = curr_fcs
  
  currdf = df[df$tissue == tis &df$sex == "female",]
  curr_fcs = tapply(currdf$logFC,currdf$gene_symbol,function(x)x[abs(x)==max(abs(x))][1])
  m_female[names(curr_fcs),tis] = curr_fcs
}

library(gplots)
heatmap.2(m_male,col=bluered(100),mar=c(10,5),trace="none",scale = "none",
          key.title = "log FC",key.xlab = "")
heatmap.2(m_female,col=bluered(100),mar=c(10,5),trace="none",scale = "none",
          key.title = "log FC",key.xlab = "")


