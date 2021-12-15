#get scores per domain
library(tidyverse)
library(caret)
library(data.table)
accuracy_and_uncertainty = function(vst_pt = "./indir/vst.csv"
                                    , pairs_sp_pt = "./outdir/site_and_overall_results/ALL_final_kld_pairs.csv"
                                    , pairs_sp_prob = "./outdir/site_and_overall_results/ALL_final_kld_probabilities.csv"
                                    , pairs_fam_pt = "./outdir/site_and_overall_results/ALL_final_family_predictions.csv"
                                    #, pairs_fam_prob = "./outdir/site_and_overall_results/ALL_final_kld_probabilities.csv"
                                    ){
  r_sources = list.files("./src/R/functions/")
  sapply(paste("./src/R/functions/", r_sources, sep=""), source)
  vst = fread(vst_pt) #%>% filter(height > 3)
  
  pairs = readr::read_csv(pairs_sp_pt)
  probabilities = readr::read_csv(pairs_sp_prob)
  colnames(pairs)=c("id", "individualID", "obs", "pred")
  
  pairs$domainID = substr(pairs$individualID, 10,12)
  pairs$siteID = substr(pairs$individualID, 14,17)
  
  sp_in_dat = pairs$obs %>% data.frame
  colnames(sp_in_dat) = "taxonID"
  
  #remove trees that are fully shaded in the target year
  full_shaded_ids = vst %>% 
    filter(!str_detect(plantStatus,"Live")) %>%
    select(individualID) %>% unique
  
  pairs = pairs %>% filter(!individualID %in% full_shaded_ids$individualID)
  probabilities = probabilities %>% filter(!individualID %in% full_shaded_ids$individualID)

  
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
  tot_species = vst %>% filter(!siteID %in% c("HEAL", "PUUM")) %>% select(taxonID, siteID) %>% group_by(siteID) %>% table
  data = data.table::fread("./outdir/metadata.csv") %>%
    filter(!taxonID %in% c("ABIES", "BETUL", "FRAXI", "SALIX", "2PLANT", "GLTR", "SASSA", "GYDI",
                           "OXYDE", "HALES", "PINUS", "QUERC", "PICEA", "ULMUS", "MAGNO", "LARIX")) %>%
    filter(!siteID %in% c("HEAL", "PUUM"))
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
  
  
  sp_st = lapply(1:27, function(x) length(spst[[x]]))
  sp_fr = lapply(1:27, function(x) sum(spdt[[x]]))
  sp_st = unlist(sp_st) %>% data.frame
  colnames(sp_st) = "alpha"
  sp_st[["siteID"]] = unlist(spid)
  sp_st["indata"] = unlist(sp_fr)
  
  ordered_dat = sp_st[order(sp_st$alpha, decreasing = T),2]
  
  sp_st = reshape2::melt(sp_st, "siteID")
  sp_st$siteID = factor(sp_st$siteID, levels = ordered_dat)
  #unused plot: plottin the total number of speices. Now in the bottom of supplement figure as text information
  # ggplot(sp_st, aes(fill=variable, y=value, x=siteID)) + 
  #   geom_bar(position="dodge", stat="identity") + 
  #   theme_bw() + theme(axis.text.x = element_text(size=14, angle=45, hjust=1, vjust=1)) 
  # 
  # 
  
  
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
  microF1 = macroF1 = list()
  for(dm in unique(pairs$domainID)){
    dm_dt = pairs %>% filter(domainID == dm)
    
    dm_dt$obs = factor(dm_dt$obs, levels = unique(data$taxonID))
    dm_dt$pred = factor(dm_dt$pred, levels = unique(data$taxonID))
    cmdm = confusionMatrix(dm_dt$obs, dm_dt$pred)
    print(cmdm$overall)
    #mcm = as.matrix.data.frame(cmdm$table)
    #rownames(mcm) = colnames(mcm) = colnames(cmdm$table)
    microF1[[dm]] <- cmdm$overall[1]
    macroF1[[dm]] = mean(cmdm$byClass[,"F1"], na.rm=T)
    cm[[dm]] = cmdm
  }
  microF1_dom = unlist(microF1)
  macroF1_dom = unlist(macroF1)
  
  cm_dom = cm
  
  #site_level stats
  cm = list()
  microF1 = macroF1 = list()
  for(dm in unique(pairs$siteID)){
    dm_dt = pairs %>% filter(siteID == dm)
    
    dm_dt$obs = factor(dm_dt$obs, levels = unique(data$taxonID))
    dm_dt$pred = factor(dm_dt$pred, levels = unique(data$taxonID))
    cmdm = confusionMatrix(dm_dt$obs, dm_dt$pred)
    print(cmdm$overall)
    #mcm = as.matrix.data.frame(cmdm$table)
    #rownames(mcm) = colnames(mcm) = colnames(cmdm$table)
    microF1[[dm]] <- cmdm$overall[1]
    macroF1[[dm]] = mean(cmdm$byClass[,"F1"], na.rm=T)
    cm[[dm]] = cmdm
  }
  microF1_site= unlist(microF1)
  macroF1_site= unlist(macroF1)
  
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
  write_csv(sp_co, "./outdir/species_confusion.csv")
  
  
  #plot distribution of confusion by health class
  figure_2(pairs, vst)
  
  #plot geographic accuracy
  plot_figure_3(microF1_dom, microF1_site, macroF1_site)

  
  
  # uncertainty quantitative analysis using ranking
  # Can we look at it more continuosly: 0-10 which fraction correctly classified? 
  figure_6(pairs, probabilities)
  #family confusion
  pairs = readr::read_csv(pairs_fam_pt)
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
  write_csv(sp_co, "./outdir/genuses_confusion.csv")
  

    
}
