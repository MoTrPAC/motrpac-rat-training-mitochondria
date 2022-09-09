#!/bin/R
# David Amar, Nicole Gay

library(data.table)
library(dplyr)
library(tidyverse)
library(IHW)

repo_local_dir = "~/Desktop/repos/"
source(paste0(repo_local_dir,"motrpac-mawg/pass1b-06/integrative/clustering/cluster_viz_fx.R"))
library(DO.db)
library(DOSE)

#####################################################################
# DO preprocessing
# get structure, higher level terms, and gene sets
doterms = as.list(DOTERM)
dochild = as.list(DOCHILDREN)
doancestors = as.list(DOANCESTOR)
root = "DOID:4"
max_slim_dist = 3
DO_layers = list()
DO_layers[[1]] = root
for(j in 1:max_slim_dist){
  DO_layers[[j+1]] = unique(unlist(dochild[DO_layers[[j]]]))
  if(length(DO_layers[[j+1]])==0){
    DO_layers[[j+1]] = DO_layers[[j]]
  }
  unique_set = setdiff(DO_layers[[j+1]],DO_layers[[j]])
  if(length(unique_set)>0){
    DO_layers[[j+1]] = unique_set
  }
}
do_slim_terms = sapply(DO_layers[[1+max_slim_dist]],function(x)doterms[[x]])
do_slim_terms = do_slim_terms[!sapply(do_slim_terms,is.null)]
do_slim_terms = data.frame(
  doid = names(do_slim_terms),
  term = sapply(do_slim_terms,Term)
)
rownames(do_slim_terms) = do_slim_terms$doid
classlist = sapply(
  names(doterms),
  function(x,df){
    df[intersect(doancestors[[x]],df[,1]),2]
  },
  df = do_slim_terms
)

# Add higher term DO ids that do not appear in the classlist
for(j in 1:max_slim_dist){
  curr_set = names(which(sapply(classlist[DO_layers[[j]]],length)==0))
  for(d in curr_set){
    classlist[[d]] = Term(doterms[[d]])
  }
}

#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
symb2eg = as.list(org.Hs.egSYMBOL2EG)
human_feature_to_gene = data.frame(
  feature_ID = names(unlist(symb2eg)),
  entrez_gene = as.character(unlist(symb2eg)),
  stringsAsFactors = F
)
human_feature_to_gene$gene_symbol = human_feature_to_gene$feature_ID
human_feature_to_gene$ensembl_gene = NA
human_feature_to_gene$kegg_id = NA

# map genes to DO terms:
data(DO2ALLEG)
do2gene = DO2ALLEG

#####################################################################
# Get tissue-specific DO terms
library(TissueEnrich)
load(file = system.file("extdata", "combine-expression.rda", package = "TissueEnrich"))
ens2entrez = as.list(org.Hs.egENSEMBL2EG)

gtex_tissue_spec = dataset$`GTEx-Combine`$tissueSpecificGenes
table(gtex_tissue_spec$Tissue)

gtex_bg = gtex_tissue_spec[,"Gene"]
gtex_bg = unique(unname(unlist(ens2entrez[gtex_bg])))

simple_hyper_test<-function(s1,s2,bg){
  x1 = rep(F,length(bg))
  x2 = rep(F,length(bg))
  names(x1) = bg
  names(x2) = bg
  x1[intersect(bg,s1)] = T
  x2[intersect(bg,s2)] = T
  if(sum(x1)==0 || sum(x2)==0){return(1)}
  tb = table(x1,x2)
  return(fisher.test(tb,alternative = "greater")$p.value)
}

gtex2motrpac_tissue = c(
  "Heart" = "HEART",
  "Adrenal.Gland" = "ADRNL",
  "Muscle" = "SKM-GN",
  "Liver"  = "LIVER",
  "Small.Intestine" = "SMLINT",
  "Brain" = "CORTEX",
  "Spleen" = "SPLEEN",
  "Lung" = "LUNG",
  "Adipose.Tissue" = "WAT-SC",
  "Kidney" = "KIDNEY"
)

tissue2do_terms = list()
for(tissue in names(gtex2motrpac_tissue)){
  motrpac_tissue = gtex2motrpac_tissue[tissue]
  tissue2do_terms[[motrpac_tissue]] = c()
  gtex_set = gtex_tissue_spec[gtex_tissue_spec$Tissue == tissue,"Gene"]
  gtex_set = unique(unname(unlist(ens2entrez[gtex_set])))
  for(dterm in names(do2gene)){
    curr_p = simple_hyper_test(do2gene[[dterm]],gtex_set,gtex_bg)
    if(curr_p < 0.001){
      tissue2do_terms[[motrpac_tissue]] = c(tissue2do_terms[[motrpac_tissue]],dterm)
    }
  }
  print(paste(motrpac_tissue,length(tissue2do_terms[[motrpac_tissue]])))
}
sapply(tissue2do_terms,length)

#####################################################################
# Add MoTrPAC data
scratch = "~/Desktop/MoTrPAC/data/pass1b_6m/"
data = load_graph_vis_data('gsutil',scratch)
# load(paste0(scratch,"graphical_analysis_results_20220126.RData"))
# load(paste0(scratch,"cluster-viz-inputs_20211116.RData"))
rdg_mapping = fread(
  sprintf("%s/gencode.v39.RGD.20201001.human.rat.gene.ids.txt",scratch),
  stringsAsFactors=F,data.table=F
)

feature_to_gene = as.data.frame(data$feature_to_gene)
rat_feature_to_gene = split(feature_to_gene$gene_symbol,feature_to_gene$feature_ID)
rat_symb_to_human_entrez = split(rdg_mapping$HUMAN_ORTHOLOG_NCBI_GENE_ID,
                                 rdg_mapping$RAT_SYMBOL)

# Add mito annotations
load(paste0(scratch,"mito_gene_annotation.RData"))
mt_symbols = mt$gene_symbol
mt_entrez_human = unlist(rat_symb_to_human_entrez[mt_symbols])

min_set_size = 20
nodes_df = c()
nodes_for_analysis = names(data$node_sets)[grepl("8w",names(data$node_sets))]
nodes_for_analysis = setdiff(nodes_for_analysis,"8w_F0_M0")
for(node in nodes_for_analysis){
  feature_set = data$node_sets[[node]]
  arrs = strsplit(feature_set,split=";")
  omes = sapply(arrs,function(x)x[1])
  tissues = sapply(arrs,function(x)x[2])
  f_ids = sapply(arrs,function(x)x[3])
  for(ome in unique(omes)){
    for(tissue in unique(tissues)){
      curr_analysis_set = f_ids[tissues==tissue & omes==ome]
      if(length(curr_analysis_set) < min_set_size){next}
      table(curr_analysis_set %in% feature_to_gene$feature_ID)
      table(curr_analysis_set %in% names(rat_feature_to_gene))
      curr_rat_genes = unlist(rat_feature_to_gene[curr_analysis_set])
      curr_rat_genes = unique(na.omit(curr_rat_genes))
      curr_human_genes = unlist(rat_symb_to_human_entrez[curr_rat_genes])
      curr_human_genes = unique(na.omit(curr_human_genes))
      if(length(curr_human_genes) < min_set_size){next}
      set_name = paste(node,ome,tissue,sep=";")
      print(set_name)
      df = data.frame("set_name" = set_name,"genes"=curr_human_genes)
      nodes_df = rbind(nodes_df,df)
    }
  }
}

# OPTIONAL: remove the ome component
newsets = strsplit(nodes_df$set_name,split=";")
newsets = sapply(newsets,function(x)paste(x[1],x[3],sep=";"))
nodes_df$set_name = newsets
dim(nodes_df)
nodes_df = unique(nodes_df)
dim(nodes_df)

# OPTIONAL: merge female and male data
# first, remove sex discordant clusters
sex_discordant_sets = grepl("F-1_M1",nodes_df$set_name) | grepl("F1_M-1",nodes_df$set_name)
nodes_df = nodes_df[!sex_discordant_sets,]
newsets = strsplit(nodes_df$set_name,split=";")
da_direction = sapply(newsets,function(x)x[1])
da_tissue = sapply(newsets,function(x)x[2])
da_direction[grepl("M1",da_direction)] = "Up"
da_direction[grepl("F1",da_direction)] = "Up"
da_direction[grepl("M-1",da_direction)] = "Down"
da_direction[grepl("F-1",da_direction)] = "Down"
nodes_df$sex_direction_set = nodes_df$set_name
nodes_df$set_name = paste(da_direction,da_tissue,sep=";")
table(nodes_df$set_name)

# Go over each set and run the enrichment analysis with the proper
# background set
all_sets = unique(nodes_df$set_name)
do_enrichment_analysis_results = c()
use_expression_bg=T
for(sname in all_sets){
  arr = strsplit(sname,split=";")[[1]]
  
  curr_tissue = arr[length(arr)]
  if(! curr_tissue %in% names(tissue2do_terms)){next}
  curr_do_terms = tissue2do_terms[[curr_tissue]]
  curr_genes = nodes_df[nodes_df$set_name==sname,2]
  
  curr_bg = data$universes_list$gene_symbol[[arr[2]]][[arr[3]]]
  curr_gene_bg = unlist(rat_symb_to_human_entrez[curr_bg])
  curr_gene_bg = unique(na.omit(curr_gene_bg))
  
  expression_bg = data$universes_list$gene_symbol$TRNSCRPT[[curr_tissue]]
  curr_expression_bg_gene_bg = unlist(rat_symb_to_human_entrez[expression_bg])
  curr_expression_bg_gene_bg = unique(na.omit(curr_expression_bg_gene_bg))
  
  # if set name does not have the "ome" component:
  if(use_expression_bg){
    curr_gene_bg = curr_expression_bg_gene_bg
  }
  
  do_term_bg_overlap = sapply(do2gene,
    function(x,y)length(intersect(x,y))/length(x),
    y = curr_expression_bg_gene_bg
  )
  
  # Limit the data to mito genes and proceed
  curr_genes = intersect(mt_entrez_human,curr_genes)
  curr_gene_bg = intersect(mt_entrez_human,curr_gene_bg)
  if(length(curr_genes)<3){next}
  
  curr_DO_res = enrichDO(
    curr_genes,
    ont = "DO",
    pvalueCutoff = 1,
    pAdjustMethod = "none",
    universe = curr_gene_bg,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 1,
    readable = T
  )
  if(is.null(curr_DO_res)){next}
  curr_DO_res = curr_DO_res@result
  curr_DO_res = curr_DO_res[curr_DO_res$ID %in% curr_do_terms,]
  if(length(curr_DO_res)==0 || nrow(curr_DO_res)==0){next}
  curr_DO_res$do_term_bg_overlap = do_term_bg_overlap[curr_DO_res$ID]
  curr_DO_res$setname = sname
  print(paste(c(sname,curr_DO_res[1,c(2,5)]),collapse="      "))
  do_enrichment_analysis_results = rbind(do_enrichment_analysis_results,curr_DO_res)
}

is_up = grepl("Up",do_enrichment_analysis_results$setname)
par(mfrow=c(2,1))
hist(do_enrichment_analysis_results$pvalue[is_up],
     main = "Up-regulated sets",col="red")
hist(do_enrichment_analysis_results$pvalue[!is_up],
     main = "Down-regulated sets",col="blue")
qqplot(
  -log10(do_enrichment_analysis_results$pvalue[is_up]),
  -log10(do_enrichment_analysis_results$pvalue[!is_up]),
  xlab = "Up-regulated sets (-log10 p-values)",
  ylab = "Down-regulated sets (-log10 p-values)",
  pch=20,col="gray",cex=1.2
);abline(0,1,col="red",lty=2)
# arrs = strsplit(do_enrichment_analysis_results$setname,split=";")
# do_enrichment_analysis_results$tp = sapply(arrs,function(x)x[1])
# do_enrichment_analysis_results$ome = sapply(arrs,function(x)x[2])
# do_enrichment_analysis_results$tissue = sapply(arrs,function(x)x[3])
# do_enrichment_analysis_results$week = sapply(
#   do_enrichment_analysis_results$tp,function(x)strsplit(x,split="_")[[1]][1])
# do_enrichment_analysis_results = do_enrichment_analysis_results[
#   do_enrichment_analysis_results$week != "0w",]

# Use IHW (requires many p-values)
# o = ihw(do_enrichment_analysis_results$pvalue,
#         as.factor(do_enrichment_analysis_results$tissue),alpha=0.05)
# do_enrichment_analysis_results$qvalue = adj_pvalues(o)
do_enrichment_analysis_results$qvalue = p.adjust(
   do_enrichment_analysis_results$pvalue,method="BH")

# Add columns to fit what enrichment_network_vis needs
do_enrichment_analysis_results$adj_p_value = do_enrichment_analysis_results$qvalue
do_enrichment_analysis_results$intersection = unname(sapply(do_enrichment_analysis_results$geneID,
  gsub,pattern='/',replacement=","))
do_enrichment_analysis_results$computed_p_value = do_enrichment_analysis_results$pvalue
do_enrichment_analysis_results$term_size = sapply(do_enrichment_analysis_results$GeneRatio,
  function(x)as.numeric(strsplit(x,split='/')[[1]][2]))
do_enrichment_analysis_results$query_size = sapply(do_enrichment_analysis_results$BgRatio,
  function(x)as.numeric(strsplit(x,split='/')[[1]][1]))
do_enrichment_analysis_results$intersection_size  = do_enrichment_analysis_results$Count
do_enrichment_analysis_results$term_name = do_enrichment_analysis_results$Description
do_enrichment_analysis_results$term_id = do_enrichment_analysis_results$ID

# write a supplementary table with all unadjusted results
supp_table_columns = c("setname","term_id","term_name","query_size","term_size","pvalue",
                       "qvalue","intersection")
do_supp_table = do_enrichment_analysis_results[,supp_table_columns]
do_supp_table = do_supp_table[order(do_supp_table$qvalue),]
write.table(do_supp_table,file="supp_table_disease_ontology_enrichment.txt",
            row.names = F,col.names = T,quote = F,sep="\t")

# select the significant and meaningful results
significant_results = do_enrichment_analysis_results[
  do_enrichment_analysis_results$qvalue < 0.2 &
    do_enrichment_analysis_results$Count > 2 ,]
significant_results[,c("Description","qvalue","setname")]
dim(significant_results)
significant_results_gene_sets = lapply(significant_results$intersection,
    function(x)strsplit(x,split=",")[[1]])

# Rank the top genes for each setname
set2top_genes = c()
for (sname in unique(significant_results$setname)){
  curr_sets = significant_results_gene_sets[
    significant_results$setname == sname
  ]
  names(curr_sets) = significant_results[
    significant_results$setname == sname, "Description"
  ]
  curr_gene_ranks = sort(table(unlist(curr_sets)),decreasing = T)
  set2top_genes[[sname]] = curr_gene_ranks
}

# Make some disease names simpler
significant_results$Description = gsub(
  "Human immunodeficiency virus",
  "HIV",
  significant_results$Description
)

significant_results$Description = gsub(
  "chronic obstructive pulmonary disease",
  "COPD",
  significant_results$Description
)

significant_results$Description = gsub(
  "non-small cell lung carcinoma",
  "NSCLC",
  significant_results$Description
)

significant_results$Description = gsub(
  "hypersensitivity reaction type IV disease",
  "DTH",
  significant_results$Description
)

significant_results$Description = gsub(
  "type 2 diabetes mellitus",
  "T2D",
  significant_results$Description
)

significant_results$Description = gsub(
  "hypersensitivity reaction type II disease",
  "hypersensitivity type II",
  significant_results$Description
)

# Create a set-do-gene network
edges = c()
nodes  = c()
for(sname in unique(significant_results$setname)){
  arr = strsplit(sname,split=";")[[1]]
  tp = arr[1]
  tp = gsub("8w_","",tp)
  tissue = arr[3]
  ome = arr[2]
  set_node = paste(tissue,tp,sep=";")
  
  curr_df = significant_results[significant_results$setname == sname,]
  rownames(curr_df)  = curr_df$Description
  curr_genes = lapply(curr_df$intersection,
    function(x)strsplit(x,split=",")[[1]])
  names(curr_genes) = rownames(curr_df)
  
  # Step 1: add all disease-set edges
  curr_disease_set_edges = data.frame(
    A = set_node,B=curr_df$Description,
    Type = ome,Score = -log10(curr_df[,"pvalue"])
  )
  
  # Step 2: add all disease-gene edges:
  curr_gene_disease_edges = c()
  for(gsname in names(curr_genes)){
    gs_df  = data.frame(
      A = gsname, B = curr_genes[[gsname]],
      Type = "Disease-Gene",Score = -log10(curr_df[gsname,"pvalue"])
    )
    curr_gene_disease_edges = rbind(
      curr_gene_disease_edges,gs_df
    )
  }
  
  # Step 3: add all set-gene edges:
  curr_gene_tb = table(unlist(curr_genes))
  curr_gene_sname_edges = data.frame(
    A = set_node,B = names(curr_gene_tb),
    Type = "Set-Gene",Score = unname(as.numeric(curr_gene_tb))
  )
  
  edges = rbind(
    edges,
    curr_disease_set_edges,
    curr_gene_disease_edges,
    curr_gene_sname_edges
  )
  
  nodes = rbind(
    nodes,
    c(set_node,"Set"),
    cbind(names(curr_gene_tb),"Gene"),
    cbind(curr_df$Description,"Disease")
  )
  nodes = unique(nodes)
}
nodes = as.data.frame(nodes,stringsAsFactors = F)
colnames(nodes) = c("Node","Type")
nodes$tp = NA
nodes$tissue = NA
setinds = nodes$Type == "Set"
nodes[setinds,"tp"] = sapply(nodes[setinds,"Node"],function(x)strsplit(x,split=";")[[1]][2])
nodes[setinds,"tissue"] = sapply(nodes[setinds,"Node"],function(x)strsplit(x,split=";")[[1]][1])

# Write text files for analysis/viz in cytoscape
write.table(edges,file="DO_res_edges.txt",sep="\t",quote = F,row.names = F,col.names = T)
write.table(nodes,file="DO_res_nodes.txt",sep="\t",quote = F,row.names = F,col.names = T)


# greedy_top_results<-function(m,jacc_thr = 0.5){
#   if(length(m)==0 || is.null(m) || is.null(dim(m))){return(NA)}
#   m = m[order(m$computed_p_value),]
#   selected_inds = c(1)
#   sets = strsplit(m$intersection,split=",")
#   for(j in 2:nrow(m)){
#     add_j = T
#     for(j2 in selected_inds){
#       curr_j = length(intersect(sets[[j]],sets[[j2]]))/length(union(sets[[j]],sets[[j2]]))
#       if(curr_j > jacc_thr){
#         add_j = F
#         break
#       }
#     }
#     if(add_j){selected_inds = c(selected_inds,j)}
#   }
#   return(m[selected_inds,])
# }
# 
# greedy_reduced_sig_results = c()
# for(sname in unique(significant_results$setname)){
#   m = significant_results[significant_results$setname==sname,]
#   if(nrow(m)>2){
#     m = greedy_top_results(m,jacc_thr = 0.2)
#   }
#   greedy_reduced_sig_results = rbind(greedy_reduced_sig_results,m)
# }

# #####################################################################################
# # Plot enrichments by state, overall tissues
# df = data.frame(counts = sort(table(significant_results$tp)),stringsAsFactors = F)
# names(df) = c("sname","count")
# df$sname = as.character(df$sname)
# df$state =  sapply(df$sname,function(x)strsplit(x,split="w_")[[1]][2])
# df$x = 1:nrow(df)
# df$week = sapply(df$sname,function(x)strsplit(x,split="_")[[1]][1])
# p<-ggplot(data=df, aes(x=state, y=count,fill=week)) +
#   geom_bar(stat="identity",position="dodge") +
#   #geom_text(aes(label=sname), vjust=1.6, color="white", size=3.5)+
#   theme_minimal()
# p
# #####################################################################################
# # show a summary matrix per tissue
# tissue = "SMLINT"
# tp = "8w_F-1_M-1"
# 
# # network of enrichment results 
# curr_df = significant_results[
#   significant_results$tissue == tissue &
#     significant_results$tp == tp,]
# xx = enrichment_network_vis(curr_df,human_feature_to_gene,classlist,
#                             include_metab = F)
# 
# curr_df = significant_results[significant_results$tissue == tissue,]
# curr_dos = unique(curr_df$Description)
# curr_omes = unique(curr_df$ome)
# curr_tps = unique(curr_df$tp)
# 
# q_df = c()
# for(ome in curr_omes){
#   for(tp in curr_tps){
#     curr_ps = rep(NA,length(curr_dos))
#     names(curr_ps) = curr_dos
#     curr_df_subset = curr_df[curr_df$tp==tp & curr_df$ome==ome,]
#     curr_ps[curr_df_subset$Description] = curr_df_subset$qvalue
#     df = data.frame("q"=curr_ps)
#     if(length(q_df)==0){
#       q_df = df
#     }
#     else{
#       q_df = cbind(q_df,df)
#     }
#     colnames(q_df)[ncol(q_df)] = paste(ome,tp,sep=",")
#   }
# }
# 
# to_rem = rowSums(q_df<0.05,na.rm=T) == 0
# q_df = q_df[!to_rem,]
# to_rem2 = colSums(q_df<0.05,na.rm=T) == 0
# q_df = q_df[,!to_rem2]
# dim(q_df)
# q_df[is.na(q_df)] = 1
# library(gplots)
# heatmap.2(t(-log10(as.matrix(q_df))),scale="none",trace = "none",mar=c(18,20),
#           key.title = "-log10(q-value)",key.xlab = "",col=bluered(100),
#           main=tissue,cexRow = 1.5,cexCol = 0.5)
# 
# #####################################################################################
# tp = "8w_F-1_M-1"
# curr_df = significant_results[significant_results$tp == tp,]
# curr_dos = unique(curr_df$Description)
# curr_omes = unique(curr_df$ome)
# curr_tissues = unique(curr_df$tissue)
# q_df = c()
# for(ome in curr_omes){
#   for(tissue in curr_tissues){
#     curr_ps = rep(NA,length(curr_dos))
#     names(curr_ps) = curr_dos
#     curr_df_subset = curr_df[curr_df$tissue==tissue & curr_df$ome==ome,]
#     curr_ps[curr_df_subset$Description] = curr_df_subset$qvalue
#     df = data.frame("q"=curr_ps)
#     if(length(q_df)==0){
#       q_df = df
#     }
#     else{
#       q_df = cbind(q_df,df)
#     }
#     colnames(q_df)[ncol(q_df)] = paste(tissue,ome,sep=",")
#   }
# }
# 
# to_rem = rowSums(q_df<0.05,na.rm=T) == 0
# q_df = q_df[!to_rem,]
# to_rem2 = colSums(q_df<0.05,na.rm=T) == 0
# q_df = q_df[,!to_rem2]
# dim(q_df)
# q_df[is.na(q_df)] = 1
# library(gplots)
# heatmap.2(t(-log10(as.matrix(q_df))),scale="none",trace = "none",mar=c(15,10),
#           key.title = "-log10(q-value)",key.xlab = "",col=bluered(100),
#           main=tp,cexRow = 0.6)


# ####################################################################
# # Play around with visualization of enrichment results 
# library(ggplot2)
# library(data.table)
# 
# head(significant_results)
# sigresults_dt = data.table(significant_results)
# sigresults_dt[,parent := classlist[term_id]] 
# sigresults_dt[,single_parent := sapply(parent, function(x) x[1])] # just take the first parent category in the list for now
# table(sigresults_dt[,single_parent])
# table(sigresults_dt[,tp])
# 
# # what if we just plot everything in the "8w_F-1_M-1"?
# ig = enrichment_network_vis(sigresults_dt[tp=="8w_F-1_M-1"],
#                       human_feature_to_gene,
#                       classlist,
#                       corr_thresh = 0.375,
#                       adj_pval_cutoff = 0.05,
#                       title = "8w_F-1_M-1",
#                       add_group_label_nodes = FALSE,
#                       intersection_id_type='gene_symbol',
#                       return_graph_for_cytoscape = TRUE)
# 
# # export to look at it in Cytoscape
# save(ig, file="~/DO_8w_F-1_M-1_igraph.RData")
# 
# # what if we plot all sex-consistent enrichments?
# enrichment_network_vis(sigresults_dt[tp %in% c("8w_F1_M1","8w_F-1_M-1")],
#                        human_feature_to_gene,
#                        classlist,
#                        corr_thresh = 0.375,
#                        adj_pval_cutoff = 0.05,
#                        title = "Sex-consistent",
#                        add_group_label_nodes = FALSE,
#                        intersection_id_type='gene_symbol')
# 
# ####################################################################
