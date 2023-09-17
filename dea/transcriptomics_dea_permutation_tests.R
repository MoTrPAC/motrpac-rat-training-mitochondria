
library(data.table)
library(DESeq2)
library(MotrpacBicQC)
library(ggplot2)
library(metap)
library(IHW)
library(corrplot)
library(MotrpacRatTraining6moData)
library(parallel)

output_dir = "~/Desktop/MoTrPAC/pass1b_landscape/"
mawg_repo = "~/Desktop/repos/motrpac-mawg/"

# On Stanford's sherlock
# Done in R/4.0.2
output_dir = "/oak/stanford/groups/euan/users/davidama/motrpac/pass1b_landscape/"
mawg_repo = "/oak/stanford/groups/euan/users/davidama/motrpac/pass1b_landscape/motrpac-mawg/"

# source the DEA pipeline
# requires: qvalue, ggcorrplot, reactome.db, fgsea, org.Rn.eg.db, ggplotify
source(paste0(mawg_repo,'pass1b-06/integrative/outliers_covariates/pi1_cook_fx.R'))
source(paste0(mawg_repo,'pass1b-06/tools/get_fx.R'))


##################################################################
# get animal code data
tissue_codes = names(TISSUE_CODE_TO_ABBREV)

# read the RNA-seq metrics to obtain mtRNA read percentages
rnaseq_qa_qc = TRNSCRPT_META
pct_mito_reads = rnaseq_qa_qc$pct_chrM
names(pct_mito_reads) = as.character(rnaseq_qa_qc$vial_label)
rownames(rnaseq_qa_qc) = as.character(rnaseq_qa_qc$vial_label)
colnames(rnaseq_qa_qc) = tolower(colnames(rnaseq_qa_qc))

# Add outliers and phenotypic data
outliers = as.character(OUTLIERS[,1])
pheno = PHENO
# also add vena cava outliers
venacv_outliers = unique(pheno[
  pheno$sex=="female" & pheno$sacrificetime %in% c("1w","2w") & 
    pheno$specimen.processing.sampletypedescription=="Aorta", "viallabel"])
outliers = union(outliers, venacv_outliers)
rownames(pheno) = as.character(pheno$viallabel)

##################################################################
# helper functions

# override the timewise analysis code: add explicit training groups
transcript_timewise_dea_each_sex = function(tissue_code, meta, 
    counts, covariates, curr_outliers, date, training_groups = c('1w','2w','4w','8w'),
    save_rdata=T, write=T){
  
  # fix some inconsistencies 
  if(tissue_code == 't54-hypothalmus'){
    tissue_code = 't54-hypothalamus'
  }
  
  outfile = sprintf('dea/pass1b-06_%s_transcript-rna-seq_timewise-dea_%s.txt',
                    tissue_code,date)
  if(file.exists(outfile)){
    dt = fread(outfile, sep='\t', header=T)
    return(dt)
  }
  
  sex_res = list()
  for(SEX in unique(meta[,sex])){
    
    # subset counts and meta
    curr_samples = meta[sex == SEX, viallabel]
    curr_meta = meta[sex == SEX]
    curr_counts = counts[,curr_samples]
    
    contrasts = list()
    i = 1
    for (tp in training_groups){
      contrasts[[i]] = c('group', tp, 'control')
      i = i+1
    }
    
    # shrunk results 
    # function in pi1_cook_fx.R
    deseq_res_shrunk = run_deseq(curr_counts, # filtered counts
                                 curr_meta, # metadata
                                 covariates, # covariates
                                 'group', # outcome of interest
                                 contrasts, # list of contrasts in format c(outcome_of_interest, numerator_level, denominator_level)
                                 shrink = T)
    
    # non-shrunk results 
    deseq_res = run_deseq(curr_counts, # filtered counts
                          curr_meta, # metadata
                          covariates, # covariates
                          'group', # outcome of interest
                          contrasts, # list of contrasts in format c(outcome_of_interest, numerator_level, denominator_level)
                          shrink = F)
    
    if(save_rdata){
      save(deseq_res, deseq_res_shrunk, file=sprintf('rdata/%s_%s_timewise-dea_%s.RData', tissue_code, SEX, date))
    }
    
    # collect res
    res_shrunk = data.table(deseq_res_shrunk$res)
    res_nonshrunk = data.table(deseq_res$res)
    setnames(res_shrunk, c("log2FoldChange", "lfcSE"), c("shrunk_logFC","shrunk_logFC_se"))
    setnames(res_nonshrunk, c("log2FoldChange", "lfcSE", "stat"), c("logFC","logFC_se", "zscore"))
    res_shrunk = res_shrunk[,.(gene_id, shrunk_logFC, shrunk_logFC_se, numerator, denominator)]
    res = merge(res_nonshrunk, res_shrunk, by=c("gene_id","numerator","denominator"))
    res[,sex := SEX]
    
    setnames(res, c("numerator","pvalue","gene_id"), c("comparison_group","p_value","feature_ID"))
    res[,denominator := NULL]
    
    # add some columns
    res[,tissue := tissue_code]
    res[,assay := 'transcript-rna-seq']
    res[,removed_samples := paste0(curr_outliers, collapse=',')]
    # res[,covariates := paste0(covariates, collapse=',')] added within run_deseq()
    
    # add average intensities 
    norm_counts = as.data.frame(counts(deseq_res$dds, normalized=T))
    ref_sub = norm_counts[,as.character(curr_meta[group == 'control', viallabel])]
    ref_means = rowMeans(ref_sub, na.rm=T)
    ref_se = apply(ref_sub, 1, function(x) sd(x)/sqrt(sum(!is.na(x))) )
    mlist = list()
    i = 1
    for(tp in unique(res[,comparison_group])){
      # get average values
      counts_sub = norm_counts[,as.character(curr_meta[group == tp, viallabel])]
      counts_means = data.table(sex=SEX,
                                comparison_group=tp,
                                comparison_average_intensity=rowMeans(counts_sub, na.rm=T),
                                comparison_average_intensity_se=apply(counts_sub, 1, 
                                                                      function(x) sd(x)/sqrt(sum(!is.na(x))) ),
                                reference_average_intensity=ref_means,
                                reference_average_intensity_se=ref_se,
                                feature_ID=rownames(counts_sub))
      mlist[[i]] = counts_means
      i = i+1
    }
    cmeans = rbindlist(mlist)
    dt = merge(res, cmeans, by=c('feature_ID', 'sex', 'comparison_group'))
    
    sex_res[[SEX]] = dt
  }
  
  dt = rbindlist(sex_res)
  
  dt = dt[,.(
    feature_ID,
    sex,
    comparison_group,
    assay,
    tissue,
    covariates,
    removed_samples,
    logFC,
    logFC_se,
    shrunk_logFC,
    shrunk_logFC_se,
    zscore,
    p_value,
    comparison_average_intensity,
    comparison_average_intensity_se,
    reference_average_intensity,
    reference_average_intensity_se
  )]
  
  # if aorta, remove 1w, 2w F
  if(tissue_code=="t65-aorta"){
    dt = dt[!(sex=='female' & comparison_group %in% c('1w','2w'))]
  }
  
  if(write){
    write.table(dt, file=outfile, sep='\t', col.names=T, row.names=F, quote=F)
  }
  return(dt)
  
}

get_tissue_prepped_data<-function(tissue_name){
  tissue_code = TISSUE_ABBREV_TO_CODE[tissue_name]
  # load tissue data 
  raw_data_name = paste0("TRNSCRPT_",gsub("-","",tissue_name),"_RAW_COUNTS")
  norm_data_name = paste0("TRNSCRPT_",gsub("-","",tissue_name),"_NORM_DATA")
  raw_data = get(raw_data_name)
  norm_data = get(norm_data_name)
  rownames(raw_data) = raw_data$feature_ID
  raw_data = raw_data[,-c(1:4)]
  rownames(norm_data) = norm_data$feature_ID
  norm_data = norm_data[,-c(1:4)]
  # filter unexpressed genes by taking the genes in the norm data
  raw_data = raw_data[rownames(norm_data),]
  curr_vials = colnames(norm_data)
  data = list(
    meta = data.table(cbind(pheno[curr_vials,],rnaseq_qa_qc[curr_vials,])),
    counts = raw_data[,curr_vials],
    tmm = norm_data[,curr_vials]
  )
  # prep data
  prepped = transcript_prep_data(tissue_code, 
    data$meta, data$counts, 
    data$tmm,
    covariates=c('pct_globin', 'rin', 'pct_umi_dup', 'median_5_3_bias'),
    outliers=outliers)
  
  # take training groups that appear in both sexes
  reduced_meta = as.data.frame(prepped$fixed_meta)
  rownames(reduced_meta) = reduced_meta$viallabel
  reduced_meta = reduced_meta[!rownames(reduced_meta) %in% outliers,]
  tb = table(reduced_meta$group,reduced_meta$sex)
  tb = tb[!rownames(tb)=="control",]
  if(!is.null(dim(tb)) && ncol(tb)>1){
    tb = tb[apply(tb>2,1,all),]
    training_groups = rownames(tb)
  }else{
    tb = tb[tb>2]
    training_groups = names(tb)
  }
  
  return(list(
    prepped = prepped, training_groups=training_groups
  ))
}

permute_groups<-function(labels){
  newv = labels[]
  ls = unique(labels)
  for(l in ls){
    inds = labels==l
    l_newv = rep(sample(ls,replace = F),sum(inds))[1:sum(inds)]
    l_newv = sample(l_newv,replace = F)
    newv[inds] = l_newv
  }
  return(newv)
}

perm_group_within_sex_training_dea<-function(prepped,verbose=T, permute = T,
    add_timewise_stats = T,training_groups=NULL,rfile=NULL){
  
  tissue_name = prepped$fixed_meta$tissue[1]
  tissue_code = TISSUE_ABBREV_TO_CODE[tissue_name]
  
  original_group_info = prepped$fixed_meta$group
  if(permute){
    for(sex in unique(prepped$fixed_meta$sex)){
      inds = prepped$fixed_meta$sex == sex
      prepped$fixed_meta$group[inds] = permute_groups(
        prepped$fixed_meta$group[inds])
    }
    new_group_info = prepped$fixed_meta$group
  }else{
    new_group_info = original_group_info
  }
  names(new_group_info) = prepped$fixed_meta$viallabel
  
  perm_training_dea = transcript_training_dea_each_sex(
    tissue_code, 
    prepped$fixed_meta, 
    prepped$fixed_counts, 
    prepped$fixed_covariates, 
    prepped$curr_outliers,
    date = "14122022",write=F)
  
  perm_training_dea = as.data.frame(perm_training_dea)
  is_selected = p.adjust(perm_training_dea$p_value,method="fdr")<0.05
  is_selected[is.na(is_selected)] = F
  selected_features = unique(perm_training_dea[is_selected,"feature_ID"])
  num_selected_features = sum(is_selected,na.rm = T)
  if(verbose){
    print("Number of 0.05 FDR selected features after permutation:")
    print(num_selected_features)
  }
  
  perm_stats = c("num_selected" = num_selected_features)
  is_selected_m = p.adjust(perm_training_dea$p_value_male,method="fdr")<0.05
  is_selected_m[is.na(is_selected_m)] = F
  is_selected_f = p.adjust(perm_training_dea$p_value_female,method="fdr")<0.05
  is_selected_f[is.na(is_selected_f)] = F
  perm_stats["selected_male"] = sum(is_selected_m)
  perm_stats["selected_female"] = sum(is_selected_f)
  perm_stats["selected_both"] = sum(is_selected_f & is_selected_m)
  perm_stats["selected_m_only"] = sum(!is_selected_f & is_selected_m)
  perm_stats["selected_f_only"] = sum(is_selected_f & !is_selected_m)
  
  if(add_timewise_stats){
    if(is.null(training_groups)){
      training_groups = c("1w","2w","4w","8w")
    }
    perm_timewise_dea = transcript_timewise_dea_each_sex(tissue_code,
      prepped$fixed_meta,prepped$fixed_counts,
      prepped$fixed_covariates,prepped$curr_outliers,
      date = "14122022",
      training_groups = training_groups,
      write=F, save_rdata = F)
    perm_timewise_dea = as.data.frame(perm_timewise_dea)
    perm_timewise_dea = perm_timewise_dea[perm_timewise_dea$feature_ID %in% selected_features,]
    selected_logfc_0.25 = unique(perm_timewise_dea[abs(perm_timewise_dea$shrunk_logFC)>0.25,"feature_ID"])
    selected_logfc_0.5 = unique(perm_timewise_dea[abs(perm_timewise_dea$shrunk_logFC)>0.5,"feature_ID"])
    selected_logfc_1 = unique(perm_timewise_dea[abs(perm_timewise_dea$shrunk_logFC)>1,"feature_ID"])
    perm_stats["num_selected_and_logfc_0.25"] = length(selected_logfc_0.25)
    perm_stats["num_selected_and_logfc_0.5"] = length(selected_logfc_0.5)
    perm_stats["num_selected_and_logfc_1"] = length(selected_logfc_1)
  }
  
  # if rfile is not null then save the simulations data
  if (!is.null(rfile)){
    perm_group_labels = new_group_info
    save(
      perm_group_labels,
      perm_stats,
      perm_training_dea,
      file=rfile
    )
  }
  
  return(list(
    "analysis_stats" = perm_stats,
    is_perm = permute,
    labels = list(
      "real" = original_group_info,
      "permuted" = new_group_info
    )
  ))
}

compute_sex_dea_gene_stats<-function(timewise,ids=NULL,
                                     stat_columns = c("zscore","logFC"),add_zero=T){
  timewise = as.data.frame(timewise)
  if(is.null(ids)){
    ids = unique(timewise$feature_ID)
  }
  timewise = timewise[timewise$feature_ID %in% ids,]
  timewise_male = timewise[timewise$sex=="male",]
  timewise_female = timewise[timewise$sex=="female",]
  
  feature_stats = c()
  for(column_name in stat_columns){
    x_male = timewise_male[,c("feature_ID","comparison_group",column_name)]
    x_female = timewise_female[,c("feature_ID","comparison_group",column_name)]
    x_male = reshape(x_male,idvar="feature_ID",timevar = "comparison_group",direction = "wide")
    x_female = reshape(x_female,idvar="feature_ID",timevar = "comparison_group",direction = "wide")
    rownames(x_male) = x_male$feature_ID
    rownames(x_female) = x_female$feature_ID
    x_male = x_male[,-1]
    x_female = x_female[,-1]
    x_female = x_female[rownames(x_male),colnames(x_female)]
    if(add_zero){
      x_male = cbind(0,x_male)
      x_female = cbind(0,x_female)
    }
    sex_diffs = x_male - x_female
    rhos = sapply(1:nrow(x_male),
      function(i,m1,m2)cor(as.numeric(m1[i,]),as.numeric(m2[i,])),m1= x_male, m2= x_female)
    # sanity check
    # all(rhos == diag(cor(t(x_male),t(x_female))),na.rm=T)
    curr_stats = data.frame(
      feature_ID = rownames(sex_diffs),
      stat_name = column_name,
      pearson_rho = rhos,
      max_abs_diff = apply(abs(sex_diffs),1,max,na.rm=T),
      mean_abs_diff = apply(abs(sex_diffs),1,mean,na.rm=T)
    )
    rownames(curr_stats) = NULL
    feature_stats = rbind(feature_stats,curr_stats)
  }
  return(feature_stats)
}

get_timewise_stat_tables<-function(timewise,ids=NULL,column_name = "zscore",add_zero=T){
    timewise = as.data.frame(timewise)
    if(is.null(ids)){
      ids = unique(timewise$feature_ID)
    }
    timewise = timewise[timewise$feature_ID %in% ids,]
    timewise_male = timewise[timewise$sex=="male",]
    timewise_female = timewise[timewise$sex=="female",]
    
    x_male = timewise_male[,c("feature_ID","comparison_group",column_name)]
    x_female = timewise_female[,c("feature_ID","comparison_group",column_name)]
    x_male = reshape(x_male,idvar="feature_ID",timevar = "comparison_group",direction = "wide")
    x_female = reshape(x_female,idvar="feature_ID",timevar = "comparison_group",direction = "wide")
    rownames(x_male) = x_male$feature_ID
    rownames(x_female) = x_female$feature_ID
    x_male = x_male[,-1]
    x_female = x_female[,-1]
    x_female = x_female[rownames(x_male),colnames(x_female)]
    return(list(
      "male" = x_male,
      "female" = x_female
    ))
}

permute_sex_within_group<-function(sex,labels){
  news = sex[]
  for(l in unique(labels)){
    # go over each sex and swap the sex labels of half of the animals
    for(s in c("male","female")){
      inds = which(labels==l & sex==s)
      inds = sample(inds,replace = F)[1:(length(inds)/2)]
      if(s=="male"){
        news[inds] = "female"
      }else{
        news[inds] = "male"
      }
    }
  }
  return(news)
}

perm_sex_get_dea_sex_stats<-function(prepped,training_groups,verbose=T,
                                     permute_sex = T,rfile=NULL){
  
  tissue_name = prepped$fixed_meta$tissue[1]
  tissue_code = TISSUE_ABBREV_TO_CODE[tissue_name]
  
  if (!is.null(rfile) && file.exists(rfile)){
    load(rfile)
    return(sex_perm_stats)
  }
  
  if(permute_sex){
    prepped$fixed_meta$sex = permute_sex_within_group(
      prepped$fixed_meta$sex, prepped$fixed_meta$group)
    if(verbose){print("Sex labels reshuffled")}
  }
  
  if(verbose){print("Computing differential abundance stats")}
  sex_perm_timewise_dea = transcript_timewise_dea_each_sex(tissue_code,
    prepped$fixed_meta,prepped$fixed_counts,
    prepped$fixed_covariates,prepped$curr_outliers,
    date = "14122022",
    training_groups = training_groups,
    write=F, save_rdata = F)
  sex_perm_training_dea = transcript_training_dea_each_sex(
    tissue_code, 
    prepped$fixed_meta, 
    prepped$fixed_counts, 
    prepped$fixed_covariates, 
    prepped$curr_outliers,
    date = "14122022",write=F)
  
  sex_perm_selected_training_genes = sex_perm_training_dea$feature_ID[
    p.adjust(sex_perm_training_dea$p_value,method="fdr")<0.05]
  
  if(verbose){print(paste("found",length(sex_perm_selected_training_genes),"training differential genes"))}
  sex_perm_stats = compute_sex_dea_gene_stats(sex_perm_timewise_dea)
  sex_perm_stats$is_sex_perm = permute_sex
  sex_perm_stats$selected_feature = sex_perm_stats$feature_ID %in% sex_perm_selected_training_genes
  
  if (!is.null(rfile)){
    perm_sex_labels = prepped$fixed_meta$sex
    names(perm_sex_labels) = prepped$fixed_meta$viallabel
    save(
      perm_sex_labels,
      sex_perm_timewise_dea,
      sex_perm_training_dea,
      sex_perm_stats,
      file=rfile
    )
  }
  return(sex_perm_stats)
}

##################################################################
# permutation test: permute groups within sex
perm_dea_results = c()
perm_dea_simulation_labels = c()
for(tissue_name in sort(unique(names(TISSUE_ABBREV_TO_CODE)))){
  tissue_code = TISSUE_ABBREV_TO_CODE[tissue_name]
  if(tissue_name=="" || tissue_code == ""){next}
  if(tissue_name %in% c("VENACV","TESTES","OVARY","PLASMA")){next}
  print(paste("Analyzing tissue:",tissue_name))

  tissue_data_obj = get_tissue_prepped_data(tissue_name)
  prepped = tissue_data_obj$prepped
  training_groups = tissue_data_obj$training_groups
  print("analyzing the following timewise contrasts:")
  print(training_groups)

  unperm_results = perm_group_within_sex_training_dea(
    prepped,verbose=T, permute = F,add_timewise_stats=F,rfile=NULL)
  curr_res = data.frame(
    tissue=tissue_name,
    is_permuted = F
  )
  for(nn in names(unperm_results$analysis_stats)){
    curr_res[[nn]] = unperm_results$analysis_stats[nn]
  }
  print("Number of DEA at 5% FDR:")
  print(curr_res$num_selected[1])
  perm_dea_results = rbind(perm_dea_results,curr_res)

  permuted_counts = mclapply(1:100,
    function(x,y)perm_group_within_sex_training_dea(y,
        permute=T,add_timewise_stats=F,rfile=NULL),
    y=prepped,mc.cores = 25)
  permuted_stats = sapply(permuted_counts,function(x)x[[1]])
  perm_dea_simulation_labels[[tissue_name]] = lapply(permuted_counts,function(x)x[[3]])

  curr_perm_res = data.frame(
    tissue=rep(tissue_name,length(permuted_counts)),
    is_permuted = T
  )
  for(nn in rownames(permuted_stats)){
    curr_perm_res[[nn]] = permuted_stats[nn,]
  }
  perm_dea_results = rbind(perm_dea_results,curr_perm_res)
}
save(perm_dea_results,perm_dea_simulation_labels,
     file=paste0(output_dir,"perm_dea_results_nonnaive_perm.RData"))

##################################################################
# permutation of sex labels: parallelize by tissue
analyze_tissue_sex_perm<-function(tissue_name,output_dir,reps=100){
  tissue_code = TISSUE_ABBREV_TO_CODE[tissue_name]
  if(tissue_name=="" || tissue_code == ""){return(0)}
  if(tissue_name %in% c("VENACV","TESTES","OVARY","PLASMA")){return(0)}
  print(paste("Analyzing tissue:",tissue_name))

  tissue_data_obj = get_tissue_prepped_data(tissue_name)
  prepped = tissue_data_obj$prepped
  training_groups = tissue_data_obj$training_groups
  print("analyzing the following timewise contrasts:")
  print(training_groups)

  # run the unpermuted analysis if not loaded
  print(paste("Running DEA, tissue:",tissue_name))
  currfname = paste0(output_dir,tissue_name,"_real_data_dea_results.RData")
  real_sex_dea_stats = perm_sex_get_dea_sex_stats(
    prepped,training_groups = training_groups,permute_sex = F,
    rfile = currfname)
  print(paste("Done, tissue:",tissue_name))

  # Run the permutations
  for(ii in 1:reps){
    print(paste("Running sex perm dea, rep:",ii,", tissue:",tissue_name))
    try({
      currfname = paste0(output_dir,tissue_name,"_sexperm_dea_results_",ii,".RData")
      tmp = perm_sex_get_dea_sex_stats(prepped,training_groups,rfile=currfname)
    })
  }
  return(1)
}

tissues_for_analysis = sort(unique(names(TISSUE_ABBREV_TO_CODE)))
mclapply(
  tissues_for_analysis,
  analyze_tissue_sex_perm,
  output_dir = output_dir,
  mc.cores = length(tissues_for_analysis)
)

#########################################################################
# Analyze the results
# this code assumes that the working dir has all the RData files
# saved above
setwd("/Users/davidama/Desktop/MoTrPAC/pass1b_landscape")

# group permutation
load("./perm_dea_results_nonnaive_perm.RData")

# names(perm_dea_results)
# p<-ggplot(perm_dea_results,
#           aes(x=tissue,y=log10(num_selected),fill=is_permuted)) + geom_boxplot()
# p

# Supp figure: counts before and after group permutation
b = boxplot(log10(1+num_selected)~tissue,
            data=perm_dea_results[perm_dea_results$is_permuted,],
            las=2,xlab="",col=TISSUE_COLORS[unique(perm_dea_results$tissue)],
            ylim = c(0,0.1+max(log10(1+perm_dea_results$num_selected))),
            ylab = "log10 number of genes at 5% FDR")
real_counts = log10(perm_dea_results[!perm_dea_results$is_permuted,"num_selected"])
emp_pvals = c()
for(tissue in unique(perm_dea_results$tissue)){
  real_score = perm_dea_results[perm_dea_results$tissue==tissue & !perm_dea_results$is_permuted,"num_selected"]
  perm_scores = perm_dea_results[perm_dea_results$tissue==tissue & perm_dea_results$is_permuted,"num_selected"]
  emp_p = (1+sum(perm_scores >= real_score)) / (1+length(perm_scores))
  emp_pvals[tissue] = emp_p
}
pchs = rep(8,length(real_counts))
pchs[emp_pvals > 0.05] = 20
pchs[emp_pvals < 0.05 & emp_pvals > 0.01] = 4
points(x=1:length(real_counts),y = real_counts,lwd=3,col="black",pch=pchs)
legend(x="top",legend=c("p>0.05","p<0.05","p<=0.01"),pch=c(20,4,8),ncol = 1)

# Supp figure: counts before and after group permutation
b = boxplot(log10(1+selected_male)~tissue,
            data=perm_dea_results[perm_dea_results$is_permuted,],
            las=2,xlab="",col=TISSUE_COLORS[unique(perm_dea_results$tissue)],
            ylim = c(0,0.1+max(log10(1+perm_dea_results$selected_male))),
            ylab = "log10 number of genes at 5% FDR")
real_counts = log10(perm_dea_results[!perm_dea_results$is_permuted,"selected_male"])
emp_pvals = c()
for(tissue in unique(perm_dea_results$tissue)){
  real_score = perm_dea_results[perm_dea_results$tissue==tissue & !perm_dea_results$is_permuted,"selected_male"]
  perm_scores = perm_dea_results[perm_dea_results$tissue==tissue & perm_dea_results$is_permuted,"selected_male"]
  emp_p = (1+sum(perm_scores >= real_score)) / (1+length(perm_scores))
  emp_pvals[tissue] = emp_p
}
pchs = rep(8,length(real_counts))
pchs[emp_pvals > 0.05] = 20
pchs[emp_pvals < 0.05 & emp_pvals > 0.01] = 4
points(x=1:length(real_counts),y = real_counts,lwd=3,col="black",pch=pchs)
legend(x="top",legend=c("p>0.05","p<0.05","p<=0.01"),pch=c(20,4,8),ncol = 1)

# Supp figure: counts before and after group permutation
b = boxplot(log10(1+selected_female)~tissue,
            data=perm_dea_results[perm_dea_results$is_permuted,],
            las=2,xlab="",col=TISSUE_COLORS[unique(perm_dea_results$tissue)],
            ylim = c(0,0.1+max(log10(1+perm_dea_results$selected_female))),
            ylab = "log10 number of genes at 5% FDR")
real_counts = log10(perm_dea_results[!perm_dea_results$is_permuted,"selected_female"])
emp_pvals = c()
for(tissue in unique(perm_dea_results$tissue)){
  real_score = perm_dea_results[perm_dea_results$tissue==tissue & !perm_dea_results$is_permuted,"selected_female"]
  perm_scores = perm_dea_results[perm_dea_results$tissue==tissue & perm_dea_results$is_permuted,"selected_female"]
  emp_p = (1+sum(perm_scores >= real_score)) / (1+length(perm_scores))
  emp_pvals[tissue] = emp_p
}
pchs = rep(8,length(real_counts))
pchs[emp_pvals > 0.05] = 20
pchs[emp_pvals < 0.05 & emp_pvals > 0.01] = 4
points(x=1:length(real_counts),y = real_counts,lwd=3,col="black",pch=pchs)
legend(x="top",legend=c("p>0.05","p<0.05","p<=0.01"),pch=c(20,4,8),ncol = 1)


# # Examine outlier permutation results
# wat_data = get_tissue_prepped_data("WAT-SC")
# wat_perm_res = perm_dea_results[perm_dea_results$tissue=="WAT-SC",]
# wat_perm_res = wat_perm_res[wat_perm_res$is_permuted,]
# wat_selected_perm = which(wat_perm_res$num_selected == max(wat_perm_res$num_selected))
# wat_selected_perm_labels = perm_dea_simulation_labels[["WAT-SC"]][[
#   wat_selected_perm
# ]]
# wat_sex_labels = wat_data$prepped$fixed_meta$sex
# table(wat_selected_perm_labels[[1]][wat_sex_labels=="male"],
#       wat_selected_perm_labels[[2]][wat_sex_labels=="male"])

# # helper functions for processing sex perm repeat results
# discretize_scores<-function(x,thr){
#   v = rep(0,length(x))
#   v[x>thr]=1
#   v[x < (-1*thr)] = -1
#   return(v)
# }
# process_sex_perm_rdata<-function(rfile,zthr=3){
#   timewise = NULL
#   try({
#     load(rfile)
#     timewise = sex_perm_timewise_dea
#     selected_features = unique(sex_perm_stats$feature_ID[sex_perm_stats$selected_feature])
#     ztables = get_timewise_stat_tables(timewise)
#     mm = apply(ztables$male,2,discretize_scores,thr=zthr)
#     mf = apply(ztables$female,2,discretize_scores,thr=zthr)
#     diffs = mm-mf
#     diffs = as.data.frame(diffs)
#     diffs$selected = rownames(ztables$male) %in% selected_features
#     rownames(diffs) = rownames(ztables$male)
#     return(diffs)
#   })
#   if(is.null(timewise)){return(NULL)}
# }
# 
# 
# # sex permutation: load and parse the results
# tissues = sort(unique(perm_dea_results$tissue))
# sex_perm_counts = c()
# paired_test_res = c()
# zthr = 2.5
# for(tissue in tissues){
#   print(tissue)
#   real_res_file = paste0(tissue,"_real_data_dea_results.RData")
#   real_res = process_sex_perm_rdata(real_res_file,zthr = zthr)
#   real_counts = sum(real_res$selected & rowSums(abs(real_res[,1:4]))>0)
#   real_gene_counts = rowSums(abs(real_res[,1:4]))
#   curr_v = data.frame(
#     tissue = tissue,
#     is_perm = F,
#     z_diff_count_our_genes = real_counts,
#     z_diff_count_selected_genes = real_counts,
#     num_training_DEA = sum(real_res$selected)
#   )
#   sex_perm_counts = rbind(sex_perm_counts,curr_v)
#   
#   for(rep in 1:100){
#     print(rep)
#     tryCatch({
#       perm_res_file = paste0(tissue,"_sexperm_dea_results_",rep,".RData")
#       perm_res = process_sex_perm_rdata(perm_res_file,zthr = zthr)
#       perm_res = perm_res[rownames(real_res),]
#       perm_counts = sum(perm_res$selected & rowSums(abs(perm_res[,1:4]))>0)
#       perm_counts2 = sum(real_res$selected & rowSums(abs(perm_res[,1:4]))>0)
#       perm_gene_counts = rowSums(abs(perm_res[,1:4]))
#       wilcox_res = wilcox.test(
#         real_gene_counts,perm_gene_counts[names(real_gene_counts)],
#         paired = T,
#         alternative = "greater"
#       )$p.value
#       curr_v = data.frame(
#         tissue = tissue,
#         is_perm = T,
#         z_diff_count_our_genes = perm_counts2,
#         z_diff_count_selected_genes = perm_counts,
#         num_training_DEA = sum(perm_res$selected)
#       )
#       sex_perm_counts = rbind(sex_perm_counts,curr_v)
#       test_v = data.frame(
#         tissue = tissue,
#         paired_test_p = wilcox_res
#       )
#       paired_test_res = rbind(paired_test_res,test_v)
#     },error=function(e){return(NA)})
#   }
# }
# 
# # Supp figure: z-score diff counts before and after group permutation
# b = boxplot(log10(1+z_diff_count_our_genes)~tissue,
#             data=sex_perm_counts[sex_perm_counts$is_perm,],
#             las=2,xlab="",col=TISSUE_COLORS[unique(sex_perm_counts$tissue)],
#             ylab = paste0("log10 num genes (thr=",zthr,")"),
#             ylim = c(0,0.1+log10(max(sex_perm_counts$z_diff_count_our_genes))))
# real_counts = log10(sex_perm_counts[!sex_perm_counts$is_perm,"z_diff_count_our_genes"])
# emp_pvals = c()
# for(tissue in unique(perm_dea_results$tissue)){
#   real_score = sex_perm_counts[
#     sex_perm_counts$tissue==tissue & !sex_perm_counts$is_perm,"z_diff_count_our_genes"]
#   perm_scores = sex_perm_counts[
#     sex_perm_counts$tissue==tissue & sex_perm_counts$is_perm,"z_diff_count_our_genes"]
#   emp_p = (1+sum(perm_scores >= real_score)) / (1+length(perm_scores))
#   emp_pvals[tissue] = emp_p
# }
# pchs = rep(8,length(real_counts))
# pchs[emp_pvals > 0.05] = 20
# pchs[emp_pvals < 0.05 & emp_pvals > 0.01] = 4
# points(x=1:length(real_counts),y = real_counts,lwd=3,col="black",pch=pchs)
# legend(x="top",legend=c("p>0.05","p<0.05","p<0.01"),pch=c(20,4,8),ncol = 1)
# 
# barplot(
#   100*sort(tapply(paired_test_res$paired_test_p,paired_test_res$tissue,
#        function(x)sum(x<0.01)/length(x))),
#   las = 2,
#   ylab = "percent significant vs. permutations"
# )
# 
# # # check outlier simulations
# # tissue = "ADRNL"
# # load(paste0(tissue,"_sex_perm_results.RData"))
# # tissue_counts = sex_perm_counts[sex_perm_counts$tissue==tissue,]
# # tissue_counts = tissue_counts[-1,]
# # outlier_ind = which(tissue_counts$z_diff_count_our_genes ==
# #                               max(tissue_counts$z_diff_count_our_genes))
# # sex_perm_new_comparison_results[[outlier_ind]]
# 
# b = boxplot(log10(1+z_diff_count_selected_genes)~tissue,
#             data=sex_perm_counts[sex_perm_counts$is_perm,],
#             las=2,xlab="",col=TISSUE_COLORS[unique(sex_perm_counts$tissue)],
#             ylab = paste0("log10 num genes (thr=",zthr,")"),
#             ylim = c(0,0.1+log10(max(sex_perm_counts$z_diff_count_selected_genes))))
# real_counts = log10(sex_perm_counts[!sex_perm_counts$is_perm,"z_diff_count_selected_genes"])
# emp_pvals = c()
# for(tissue in unique(perm_dea_results$tissue)){
#   real_score = sex_perm_counts[
#     sex_perm_counts$tissue==tissue & !sex_perm_counts$is_perm,"z_diff_count_selected_genes"]
#   perm_scores = sex_perm_counts[
#     sex_perm_counts$tissue==tissue & sex_perm_counts$is_perm,"z_diff_count_selected_genes"]
#   emp_p = (1+sum(perm_scores >= real_score)) / (1+length(perm_scores))
#   emp_pvals[tissue] = emp_p
# }
# pchs = rep(8,length(real_counts))
# pchs[emp_pvals > 0.05] = 20
# pchs[emp_pvals < 0.05 & emp_pvals > 0.01] = 4
# points(x=1:length(real_counts),y = real_counts,lwd=3,col="black",pch=pchs)
# legend(x="top",legend=c("p>0.05","p<0.05","p<0.01"),pch=c(20,4,8),ncol = 1)

# b = boxplot(log10(1+logfc_diff_count_our_genes)~tissue,
#             data=sex_perm_counts[sex_perm_counts$is_perm,],
#             las=2,xlab="",col=TISSUE_COLORS[unique(sex_perm_counts$tissue)],
#             ylab = paste("Num genes with logFC diff >",logfc_thr),
#             ylim = c(0,0.1+log10(max(sex_perm_counts$logfc_diff_count_our_genes))))
# real_counts = log10(sex_perm_counts[!sex_perm_counts$is_perm,"logfc_diff_count_our_genes"])
# points(x=1:length(real_counts),y = real_counts,lwd=3,col="gray",pch=8)
# 
# 
# b = boxplot(log10(1+logfc_diff_count_selected_genes)~tissue,
#             data=sex_perm_counts[sex_perm_counts$is_perm,],
#             las=2,xlab="",col=TISSUE_COLORS[unique(sex_perm_counts$tissue)],
#             ylab = paste("Num genes with logFC diff >",logfc_thr),
#             ylim = c(0,0.1+log10(max(sex_perm_counts$logfc_diff_count_selected_genes))))
# real_counts = log10(sex_perm_counts[!sex_perm_counts$is_perm,"logfc_diff_count_selected_genes"])
# points(x=1:length(real_counts),y = real_counts,lwd=3,col="gray",pch=8)

# ###################################################################################
# # Simple non-parametric tests
# 
# all_kw_results = c()
# for(tissue_name in sort(unique(names(TISSUE_ABBREV_TO_CODE)))){
#   tissue_code = TISSUE_ABBREV_TO_CODE[tissue_name]
#   if(tissue_name=="" || tissue_code == ""){next}
#   if(tissue_name %in% c("VENACV","TESTES","OVARY","PLASMA")){next}
#   print(paste("Analyzing tissue:",tissue_name))
#   
#   tissue_data_obj = get_tissue_prepped_data(tissue_name)
#   prepped = tissue_data_obj$prepped
#   training_groups = tissue_data_obj$training_groups
#   print("analyzing the following timewise contrasts:")
#   print(training_groups)
#   
#   prepped$fixed_meta = as.data.frame(prepped$fixed_meta)
#   rownames(prepped$fixed_meta) = prepped$fixed_meta$viallabel
#   
#   # test groups within sex
#   for(sex in c("male","female")){
#     curr_samps = prepped$fixed_meta[prepped$fixed_meta$sex==sex,"viallabel"]
#     curr_samps = setdiff(curr_samps,prepped$curr_outliers)
#     m = prepped$fixed_tmm[,curr_samps]
#     g = prepped$fixed_meta[curr_samps,"group"]
#     kw_res = apply(m,1,function(x,g)kruskal.test(x,g=g)$p.value,g=g)
#     kw_res = data.frame(
#       tissue = tissue_name,
#       test_var = "sex",
#       test_group = sex,
#       gene = names(kw_res),
#       p = kw_res
#     )
#     all_kw_results = rbind(all_kw_results,kw_res)
#   }
#   # test sex within groups
#   for(gr in c("control",training_groups)){
#     curr_samps = prepped$fixed_meta[prepped$fixed_meta$group==gr,"viallabel"]
#     curr_samps = setdiff(curr_samps,prepped$curr_outliers)
#     m = prepped$fixed_tmm[,curr_samps]
#     g = prepped$fixed_meta[curr_samps,"sex"]
#     kw_res = apply(m,1,function(x,g)kruskal.test(x,g=g)$p.value,g=g)
#     kw_res = data.frame(
#       tissue = tissue_name,
#       test_var = "group",
#       test_group = gr,
#       gene = names(kw_res),
#       p = kw_res
#     )
#     all_kw_results = rbind(all_kw_results,kw_res)
#   }
# }


