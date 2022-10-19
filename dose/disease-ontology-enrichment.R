#!/bin/R
# David Amar

library(data.table)
library(dplyr)
library(tidyverse)
library(IHW)

repo_local_dir = "~/Desktop/repos/"
source(paste0(repo_local_dir,"motrpac-mawg/pass1b-06/integrative/clustering/cluster_viz_fx.R"))
library(DO.db)
library(DOSE)
library(MotrpacRatTraining6moData)
data_dir = "~/Desktop/repos/motrpac-rat-training-mitochondria/data/"

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
    if(curr_p < 0.01){
      tissue2do_terms[[motrpac_tissue]] = c(tissue2do_terms[[motrpac_tissue]],dterm)
    }
  }
  print(paste(motrpac_tissue,length(tissue2do_terms[[motrpac_tissue]])))
}
sapply(tissue2do_terms,length)

#####################################################################
# Add MoTrPAC data
feature2gene = FEATURE_TO_GENE[!grepl("chr\\d",FEATURE_TO_GENE$feature_ID),]
feature2gene = feature2gene[!grepl("cluster\\d",feature2gene$feature_ID),]
feature2gene = feature2gene[!grepl("chrX",feature2gene$feature_ID),]
feature2gene = merge(feature2gene,RAT_TO_HUMAN_GENE,by.x = "gene_symbol",by.y = "RAT_SYMBOL")
feature2human_gene = split(feature2gene$HUMAN_ORTHOLOG_SYMBOL,feature2gene$feature_ID)
feature2human_gene = lapply(feature2human_gene, unique)

# Add mito annotations
load(paste0(data_dir,"mito_gene_annotation.RData"))
# Human MitoCarta
mt_pw = fread(paste0(data_dir,"Human.MitoCarta3.0.txt"),header = T,
              stringsAsFactors = F,data.table = F)
mt_genes = unique(unlist(strsplit(mt_pw$Genes,split=",\\s+")))
mitocarta_pathways = strsplit(mt_pw$Genes,split=",\\s+")
names(mitocarta_pathways) = mt_pw$MitoPathway
symb2enter = as.list(org.Hs.egALIAS2EG)
mt_entrez_gene = unique(unlist(symb2enter[mt_genes]))
graph_nodes = merge(GRAPH_STATES,feature2gene,by.x="feature_ID",by.y="feature_ID")
graph_nodes = graph_nodes[!is.na(graph_nodes$state_8w) & graph_nodes$state_8w!="F0_M0",]
graph_nodes = graph_nodes[graph_nodes$tissue %in% gtex2motrpac_tissue,]

# # OPTIONAL: merge female and male data
# # first, remove sex discordant clusters
# sex_discordant_sets = grepl("F-1_M1",nodes_df$set_name) | grepl("F1_M-1",nodes_df$set_name)
# nodes_df = nodes_df[!sex_discordant_sets,]
# newsets = strsplit(nodes_df$set_name,split=";")
# da_direction = sapply(newsets,function(x)x[1])
# da_tissue = sapply(newsets,function(x)x[2])
# da_direction[grepl("M1",da_direction)] = "Up"
# da_direction[grepl("F1",da_direction)] = "Up"
# da_direction[grepl("M-1",da_direction)] = "Down"
# da_direction[grepl("F-1",da_direction)] = "Down"
# nodes_df$sex_direction_set = nodes_df$set_name
# nodes_df$set_name = paste(da_direction,da_tissue,sep=";")
# table(nodes_df$set_name)


# Go over each set and run the enrichment analysis with the proper
# background set
min_set_size = 10
do_enrichment_analysis_results = c()
for(sname in unique(graph_nodes$state_8w)){
  curr_df = graph_nodes[graph_nodes$state_8w==sname,]
  for(tissue in unique(curr_df$tissue)){
    curr_do_terms = tissue2do_terms[[tissue]]
    cross_ome_genes = c()
    cross_ome_bg = c()
    for(ome in unique(curr_df$ome)){
      curr_genes = curr_df[curr_df$tissue==tissue & curr_df$ome==ome,
                           "HUMAN_ORTHOLOG_NCBI_GENE_ID"]
      if(length(curr_genes)<min_set_size){next}
      print(paste(sname,tissue,ome))
      curr_gene_bg = GENE_UNIVERSES$entrez_gene[[ome]][[tissue]]
      curr_gene_bg = unique(feature2gene[
        feature2gene$entrez_gene %in% curr_gene_bg,"HUMAN_ORTHOLOG_NCBI_GENE_ID"])
      do_term_bg_overlap = sapply(do2gene,
          function(x,y)length(intersect(x,y))/length(x),
          y = curr_gene_bg)
      
      curr_DO_res = enrichDO(
        curr_genes,
        ont = "DO",
        pvalueCutoff = 1,
        pAdjustMethod = "none",
        universe = curr_gene_bg,
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 1,
        readable = T)
      if(is.null(curr_DO_res)){next}
      curr_DO_res = curr_DO_res@result
      curr_DO_res = curr_DO_res[curr_DO_res$ID %in% curr_do_terms,]
      if(length(curr_DO_res)==0 || nrow(curr_DO_res)==0){next}
      curr_DO_res$do_term_bg_overlap = do_term_bg_overlap[curr_DO_res$ID]
      curr_DO_res$set_mito_overlap_p = mito_overlap_p
      curr_DO_res$set_mito_overlap = mito_overlap
      curr_DO_res$setname = sname
      curr_DO_res$tissue = tissue
      curr_DO_res$ome = ome
      
      # Add enrichment of mito genes in the overlap
      curr_DO_res$mito_overlap_p = 1
      for(ii in 1:nrow(curr_DO_res)){
        curr_overlap_genes = unlist(symb2eg[strsplit(curr_DO_res$geneID[ii],split="\\/")[[1]]])
        curr_overlap_genes = unique(unname(curr_overlap_genes))
        mito_tb = table(
          curr_gene_bg %in% curr_overlap_genes,
          curr_gene_bg %in% mt_entrez_gene
        )
        if(length(mito_tb)>2){
          curr_DO_res$mito_overlap_p[ii] = fisher.test(mito_tb,alt="g")$p.value
        }
      }
      
      do_enrichment_analysis_results = rbind(do_enrichment_analysis_results,curr_DO_res)
      cross_ome_genes = union(cross_ome_genes,curr_genes)
      cross_ome_bg = union(cross_ome_bg,curr_gene_bg)
    }
    
    curr_DO_res = enrichDO(
      cross_ome_genes,
      ont = "DO",
      pvalueCutoff = 1,
      pAdjustMethod = "none",
      universe = cross_ome_bg,
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
    curr_DO_res$set_mito_overlap_p = mito_overlap_p
    curr_DO_res$set_mito_overlap = mito_overlap
    curr_DO_res$setname = sname
    curr_DO_res$tissue = tissue
    curr_DO_res$ome = "cross_ome"
    print(paste(c(sname,curr_DO_res[1,c(2,5)]),collapse="      "))
    
    # Add mito overlap
    curr_DO_res$mito_overlap_p = 1
    for(ii in 1:nrow(curr_DO_res)){
      curr_overlap_genes = unlist(symb2eg[strsplit(curr_DO_res$geneID[ii],split="\\/")[[1]]])
      curr_overlap_genes = unique(unname(curr_overlap_genes))
      mito_tb = table(
        cross_ome_bg %in% curr_overlap_genes,
        cross_ome_bg %in% mt_entrez_gene
      )
      if(length(mito_tb)>2){
        curr_DO_res$mito_overlap_p[ii] = fisher.test(mito_tb,alt="g")$p.value
      }
    }
    do_enrichment_analysis_results = rbind(do_enrichment_analysis_results,curr_DO_res)
  }
}

# is_up = grepl("Up",do_enrichment_analysis_results$setname)
# par(mfrow=c(2,1))
# hist(do_enrichment_analysis_results$pvalue[is_up],
#      main = "Up-regulated sets",col="red")
# hist(do_enrichment_analysis_results$pvalue[!is_up],
#      main = "Down-regulated sets",col="blue")
# qqplot(
#   -log10(do_enrichment_analysis_results$pvalue[is_up]),
#   -log10(do_enrichment_analysis_results$pvalue[!is_up]),
#   xlab = "Up-regulated sets (-log10 p-values)",
#   ylab = "Down-regulated sets (-log10 p-values)",
#   pch=20,col="gray",cex=1.2
# );abline(0,1,col="red",lty=2)

# Adjust p-values
do_enrichment_analysis_results$qvalue = p.adjust(
   do_enrichment_analysis_results$pvalue,method="BH")
do_enrichment_analysis_results = do_enrichment_analysis_results[
  order(do_enrichment_analysis_results$pvalue),]

selected_resuls = do_enrichment_analysis_results[
  do_enrichment_analysis_results$qvalue < 0.05,
]
selected_resuls = selected_resuls[
  selected_resuls$mito_overlap_p < 0.05,
]
dim(selected_resuls)
selected_resuls[,c("Description","tissue","ome","setname","qvalue","mito_overlap_p")]
write.table(
  selected_resuls,
  file=paste(data_dir,"supp_fig_disease_enrichment.tsv"),
  sep="\t",row.names = F,col.names = T,quote=F
)

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
# write.table(do_supp_table,file="supp_table_disease_ontology_enrichment.txt",
#             row.names = F,col.names = T,quote = F,sep="\t")

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


