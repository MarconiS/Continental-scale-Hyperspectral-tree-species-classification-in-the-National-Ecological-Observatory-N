#get scores per domain
library(tidyverse)
library(caret)
library(data.table)
accuracy_and_uncertainty = function(vst){
  pairs = readr::read_csv("/Users/sergiomarconi/Downloads/species_results/ALL_final_kld_pairs.csv")
  probabilities = readr::read_csv("/Users/sergiomarconi/Downloads/species_results/ALL_final_kld_probabilities.csv")
  colnames(pairs)=c("id", "individualID", "obs", "pred")
  pairs$domainID = substr(pairs$individualID, 10,12)
  pairs$siteID = substr(pairs$individualID, 14,17)
  
  sp_in_dat = pairs$obs %>% data.frame
  colnames(sp_in_dat) = "taxonID"
  
  pairs = pairs %>% filter(!individualID %in% full_shaded_ids$individualID)
  probabilities = probabilities %>% filter(!individualID %in% full_shaded_ids$individualID)
  full_shaded_ids = vst %>% 
    filter(!str_detect(plantStatus,"Live")) %>%
    select(individualID) %>% unique
  
  vst = vst %>% 
    filter(str_detect(growthForm,"tree")) %>%
    filter(str_detect(plantStatus,"Live")) %>%
    filter(siteID %in% unique(pairs$siteID)) %>%
    group_by(individualID) %>% slice_max(order_by = eventID, n= 1) %>%
    slice(1) %>% ungroup
  
  vst = clean_typos_taxonID(vst)
  
  pairs = pairs %>% filter(!obs %in% c("THUJA", "ILAN", "MORU2", "ABAM")) %>%
    filter(!siteID %in% c("HEAL", "PUUM"))
  #species_per_site = vst %>% filter(taxonID %in% tree_species$taxonID)  
  tot_species = vst %>% select(taxonID, siteID) %>% group_by(siteID) %>% table
  data = data.table::fread("./indir/metadata.csv") %>%
    filter(!taxonID %in% c("ABIES", "BETUL", "FRAXI", "SALIX", "2PLANT", "GLTR", "SASSA", "GYDI",
                           "OXYDE", "HALES", "PINUS", "QUERC", "PICEA", "ULMUS", "MAGNO", "LARIX"))
  #data %>% select(individualID, groupID, siteID, plotID) %>% unique %>% write_csv("./outputs/ids_split.csv")
  #want to check how many species out of the total we have in the dataset for each site
  spst = spdt = spid =summary_freqs = sp_tested_in_site =  needed_missing = list()
  for(ii in 1:ncol(tot_species)){
    sp_in_site =  names(which(tot_species[,ii] >0))
    sort(unique(data$taxonID))
    sp_in_dat = (sp_in_site) %in% unique(data$taxonID)
    sp_tested_in_site[[ii]] = data %>% filter(taxonID %in% sp_in_site, 
                                              siteID == colnames(tot_species)[ii],
                                              groupID == "test") %>% select(taxonID) %>% unique
    needed_missing[[ii]] = sp_in_site[!sp_in_dat]
    sp_frac = sum(sp_in_dat)/length(sp_in_site)
    sp_abundance = tot_species[sp_in_site,ii] 
    itc_frac= sum(sp_abundance[sp_in_dat])/sum(sp_abundance)
    summary_freqs[[ii]] = c(colnames(tot_species)[ii], sp_frac, itc_frac)
    spst[[ii]] = sp_in_site
    spdt[[ii]] = sp_in_dat
    spid[[ii]] = colnames(tot_species)[ii]
  }
  summary_freqs = do.call(rbind.data.frame, summary_freqs) %>% data.frame
  colnames(summary_freqs) = c("siteID", "species_fraction", "abundance_fraction")
  summary_freqs$abundance_fraction = as.numeric(as.character(summary_freqs$abundance_fraction))
  list_missing = do.call(c, needed_missing)
  species_per_site = lapply(1:length(spst), function(x) sum(spdt[[x]])) %>% unlist
  species_per_test = lapply(1:length(sp_tested_in_site), function(x) nrow(sp_tested_in_site[[x]])) %>% unlist
  #reason why wanna go multi site
  overview_species_sites = cbind.data.frame(colnames(tot_species), species_per_site, species_per_test)
  
  
  sp_st = lapply(1:29, function(x) length(spst[[x]]))
  sp_fr = lapply(1:29, function(x) sum(spdt[[x]]))
  sp_st = unlist(sp_st) %>% data.frame
  colnames(sp_st) = "alpha"
  sp_st[["siteID"]] = unlist(spid)
  sp_st["indata"] = unlist(sp_fr)
  
  ordered_dat = sp_st[order(sp_st$alpha, decreasing = T),2]
  
  sp_st = reshape2::melt(sp_st, "siteID")
  sp_st$siteID = factor(sp_st$siteID, levels = ordered_dat)
  ggplot(sp_st, aes(fill=variable, y=value, x=siteID)) + 
    geom_bar(position="dodge", stat="identity") + 
    theme_bw() + theme(axis.text.x = element_text(size=14, angle=45, hjust=1, vjust=1)) 
  
  
  
  
  #overall statistics
  dm_dt = pairs 
  dm_dt$obs = factor(dm_dt$obs, levels = unique(data$taxonID))
  dm_dt$pred = factor(dm_dt$pred, levels = unique(data$taxonID))
  cmdm = confusionMatrix(dm_dt$obs, dm_dt$pred)
  print(cmdm$overall)
  #mcm = as.matrix.data.frame(cmdm$table)
  #rownames(mcm) = colnames(mcm) = colnames(cmdm$table)
  microF1 <- cmdm$overall[1]
  species_per_site = pairs %>% select(obs, siteID) %>% unique
  species_per_site = species_per_site$siteID %>% table %>% data.frame
  
  
  #domain level statistics
  cm = list()
  microF1 = list()
  for(dm in unique(pairs$domainID)){
    dm_dt = pairs %>% filter(domainID == dm)
    
    dm_dt$obs = factor(dm_dt$obs, levels = unique(data$taxonID))
    dm_dt$pred = factor(dm_dt$pred, levels = unique(data$taxonID))
    cmdm = confusionMatrix(dm_dt$obs, dm_dt$pred)
    print(cmdm$overall)
    #mcm = as.matrix.data.frame(cmdm$table)
    #rownames(mcm) = colnames(mcm) = colnames(cmdm$table)
    microF1[[dm]] <- cmdm$overall[1]
    cm[[dm]] = cmdm
  }
  microF1_dom = unlist(microF1)
  cm_dom = cm
  
  #site_level stats
  cm = list()
  microF1 = list()
  for(dm in unique(pairs$siteID)){
    dm_dt = pairs %>% filter(siteID == dm)
    
    dm_dt$obs = factor(dm_dt$obs, levels = unique(data$taxonID))
    dm_dt$pred = factor(dm_dt$pred, levels = unique(data$taxonID))
    cmdm = confusionMatrix(dm_dt$obs, dm_dt$pred)
    print(cmdm$overall)
    #mcm = as.matrix.data.frame(cmdm$table)
    #rownames(mcm) = colnames(mcm) = colnames(cmdm$table)
    microF1[[dm]] <- cmdm$overall[1]
    cm[[dm]] = cmdm
  }
  microF1_site= unlist(microF1)
  microF1_site
  
  site_pairs = pairs %>% group_by(siteID) %>% select(obs) %>% table 
  sp_per_site = apply(site_pairs, 1, function(x)(sum(x>0)))
  entries_per_site = apply(site_pairs, 1, function(x)(sum(x)))
  
  dom_pairs = pairs %>% group_by(domainID) %>% select(obs) %>% table 
  sp_per_domain = apply(dom_pairs, 1, function(x)(sum(x>0)))
  entries_per_domain = apply(dom_pairs, 1, function(x)(sum(x)))
  
  
  #global F1 and accuracy
  dm_dt = pairs
  dm_dt$obs = factor(dm_dt$obs, levels = unique(data$taxonID))
  dm_dt$pred = factor(dm_dt$pred, levels = unique(data$taxonID))
  cmdm = confusionMatrix(dm_dt$obs, dm_dt$pred)
  print(cmdm$overall)
  
  sp_co = cmdm$table
  sp_co = cmdm$table %>% data.frame
  sp_co = sp_co %>% filter(Freq > 0)
  write_csv(sp_co, "./outputs/species_confusion.csv")
  
  # uncertainty quantitative analysis using ranking
  # Can we look at it more continuosly: 0-10 which fraction correctly classified? 
  all_pdist = probabilities %>% select(-one_of("id","X1", "individualID", "effort","max_pid", "confused")) 
  p_max = apply(all_pdist[,-c(1:7)], 1, function(x)(max(x)))
  who_max = lapply(1:nrow(all_pdist), function(x)(colnames(all_pdist[which(all_pdist[x,]==p_max[x])])))
  who_max = unlist(who_max)
  
  final_boxes = cbind.data.frame(all_pdist[["taxonID"]], who_max, p_max)
  bin = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1)
  bin = seq(0,1,by=0.03)
  fraction_ = rep(NA, length(bin))
  for(ii in 2:(length(bin)+1)){
    ith_class = final_boxes %>% filter(final_boxes$p_max < bin[ii] & 
                                         final_boxes$p_max > bin[ii-1])
    fraction_[ii-1] = sum(ith_class[,1] != ith_class[,2])/nrow(ith_class)
    fraction_[is.nan(fraction_)]=0
  }
  
  uncertainty_curve = cbind.data.frame(bin, fraction_)
  colnames(uncertainty_curve) = c("majority_p", "fraction_misclassified")
  ggplot(uncertainty_curve, aes(x = majority_p, y = 1-fraction_misclassified)) + ylim(0,1) + xlim(0,1)+
    geom_point() + theme_bw() + geom_abline(intercept = 0, slope = 1) + stat_smooth(method="lm", se=FALSE)
  ggsave("./Figures/figure_5.png")
  uu = lm(majority_p~fraction_misclassified, data = uncertainty_curve)
  uu = summary(uu)
  uu$adj.r.squared
  
  
  
  
  #family confusion
  pairs = readr::read_csv("/Volumes/Data/speciesClassificationResults/BRDF_April_family_predictions.csv")
  colnames(pairs)=c("id", "obs", "pred")
  dm_dt = pairs
  genuses = c(pairs$obs, pairs$pred) %>% unique
  dm_dt$obs = factor(dm_dt$obs, levels = unique(genuses))
  dm_dt$pred = factor(dm_dt$pred, levels = unique(genuses))
  cmdm = confusionMatrix(dm_dt$obs, dm_dt$pred)
  print(cmdm$overall)
  sp_co = cmdm$table
  sp_co = cmdm$table %>% data.frame
  sp_co = sp_co %>% filter(Freq > 0)
  write_csv(sp_co, "./outputs/genuses_confusion.csv")
}
