#get scores per domain
library(tidyverse)
library(caret)
library(data.table)
accuracy_and_uncertainty = function(
  vst_pt = "./indir/vst.csv"
  , pairs_sp_pt = "./outdir/site_and_overall_results/ALL_final_kld_pairs.csv"
  , pairs_sp_prob = "./outdir/site_and_overall_results/ALL_final_kld_probabilities.csv"
  , pairs_fam_pt = "./outdir/site_and_overall_results/ALL_final_family_predictions.csv"
){
  r_sources = list.files("./src/R/functions/")
  sapply(paste("./src/R/functions/", r_sources, sep=""), source)
  vst = fread(vst_pt) #%>% filter(height > 3)

  pairs = pp = readr::read_csv(pairs_sp_pt)
  #pp=pairs
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

  #elle = pairs %>% filter(siteID == "GUAN")
  vst = vst %>%
    #filter(str_detect(growthForm,"tree")) %>%
    filter(str_detect(plantStatus,"Live")) %>%
    #filter(!individualID %in% full_shaded_ids$individualID) %>%
    filter(siteID %in% unique(pairs$siteID)) %>%
    group_by(individualID) %>% #slice_max(order_by = eventID, n= 1) %>%
    slice(1)%>%
    ungroup

  vst = clean_typos_taxonID(vst)
  ee = vst %>% select(taxonID, scientificName) %>% unique
  ee$genus = do.call(rbind.data.frame, strsplit(ee$scientificName, split = " "))[,1]
  ee = unique(ee)
  pairs = pairs %>% filter(!obs %in% c("MORU2")) %>%
    filter(!siteID %in% c("HEAL", "PUUM"))
  #species_per_site = vst %>% filter(taxonID %in% tree_species$taxonID)
  tot_species = vst %>% filter(!siteID %in% c("HEAL", "PUUM")) %>%
    select(taxonID, siteID) %>% group_by(siteID) %>% table
  data = data.table::fread("./outdir/metadata.csv") %>%
    #filter(groupID == "test")
    filter(!taxonID %in% c("ABIES", "BETUL", "FRAXI", "SALIX", "2PLANT", "GLTR", "SASSA", "GYDI",
                           "OXYDE", "HALES", "PINUS", "QUERC", "PICEA", "ULMUS", "MAGNO")) %>%
    filter(!siteID %in% c("HEAL", "PUUM"))
  #data %>% select(individualID, groupID, siteID, plotID) %>% unique %>% write_csv("./outputs/ids_split.csv")
  #want to check how many species out of the total we have in the dataset for each site
  spst = spdt = spid =summary_freqs = sp_tested_in_site = sp_trained_in_site =  needed_missing = list()
  for(ii in 1:ncol(tot_species)){

    sp_in_site =  names(which(tot_species[,ii] >0))
    sort(unique(data$taxonID))
    sp_in_dat = (sp_in_site) %in% unique(data$taxonID)
    sp_tested_in_site[[ii]] = data %>% dplyr::filter(taxonID %in% sp_in_site,
                                                     siteID == colnames(tot_species)[ii],
                                                     groupID == "test") %>% select(taxonID) %>% unique
    sp_trained_in_site[[ii]] = data %>% dplyr::filter(taxonID %in% sp_in_site,
                                                      siteID == colnames(tot_species)[ii],
                                                      groupID == "train") %>% select(taxonID) %>% unique

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
  species_per_train = lapply(1:length(sp_trained_in_site), function(x) nrow(sp_trained_in_site[[x]])) %>% unlist
  species_per_test = lapply(1:length(sp_tested_in_site), function(x) nrow(sp_tested_in_site[[x]])) %>% unlist

  names(species_per_test) = names(species_per_site) = colnames(tot_species)
  #reason why wanna go multi site
  overview_species_sites = cbind.data.frame(colnames(tot_species), species_per_site, species_per_train, species_per_test)


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
  lvls = unique(c(unique(dm_dt$obs), unique(dm_dt$pred)))

  dm_dt$obs = factor(dm_dt$obs, levels = lvls)
  dm_dt$pred = factor(dm_dt$pred, levels = lvls)
  cmdm = confusionMatrix(dm_dt$pred , dm_dt$obs)
  cmdm$table %>% sum
  conf  = as.matrix(cmdm)
  sum(conf)

  #missing_ids = apply(rbind(apply(conf, 1, sum), apply(conf, 2, sum)), 2, sum)==0
  conf = conf %>% t
  #conf = cbind.data.frame(rownames(conf), conf)
  ee = vst %>% select(taxonID, scientificName) %>% unique
  tp = do.call(rbind.data.frame, strsplit(ee$scientificName, split = " "))[,1:2]
  ee$scientificName = paste(tp[,1], tp[,2])
  ee = unique(ee)
  ee = ee %>% filter(scientificName != "Pinus glabra") %>%
    filter(scientificName != "Fraxinus americana/pennsylvanica")
  ee = left_join(data.frame(taxonID = rownames(conf)), ee)
  colnames(conf) = rownames(conf) =  ee$scientificName
  conf = data.frame(taxon = rownames(conf), conf)


  #conf =  conf[!missing_ids, !missing_ids]
  prf1 = data.frame(taxonID = rownames(cmdm$byClass), cmdm$byClass)
  prf1$taxonID = substr(prf1$taxonID, 8, nchar(prf1$taxonID))
  sp_names = vst %>% select(taxonID, scientificName) %>% unique
  prf1_n = left_join(prf1, sp_names) %>% unique
  write_csv(prf1_n, "./outdir/overview_precision_recall_all.csv")
  write_csv((conf), "./outdir/overall_confusion_matrix_all.csv")
  as.matrix(cmdm)
  print(cmdm$overall)
  #mcm = as.matrix.data.frame(cmdm$table)
  #rownames(mcm) = colnames(mcm) = colnames(cmdm$table)
  macroF1= F1_Score_macro(unlist(dm_dt$obs), unlist(dm_dt$pred))
  microF1 = F1_Score_micro(unlist(dm_dt$obs), unlist(dm_dt$pred))
  species_per_site = pairs %>% select(obs, siteID) %>% unique
  species_per_site = species_per_site$siteID %>% table %>% data.frame


  #domain level statistics
  cm = list()
  dir.create("./outdir/cm")

  microF1 = macroF1 = list()
  for(dm in unique(pairs$domainID)){
    dm_dt = pairs %>% filter(domainID == dm)
    lvls = unique(c(unique(dm_dt$obs), unique(dm_dt$pred)))
    dm_dt$obs = factor(dm_dt$obs, levels = lvls)
    dm_dt$pred = factor(dm_dt$pred, levels = lvls)
    cmdm = confusionMatrix(dm_dt$pred, dm_dt$obs)
    print(cmdm$overall)

    conf  = as.matrix(cmdm)
    conf = conf %>% t

    ee = vst %>% select(taxonID, scientificName) %>% unique
    tp = do.call(rbind.data.frame, strsplit(ee$scientificName, split = " "))[,1:2]
    ee$scientificName = paste(tp[,1], tp[,2])
    ee = unique(ee)
    ee = ee %>% filter(scientificName != "Pinus glabra") %>%
      filter(scientificName != "Fraxinus americana/pennsylvanica")
    ee = left_join(data.frame(taxonID = rownames(conf)), ee)

    if(length(conf)==1) {
      colnames(conf) = rownames(conf) =  ee$scientificName
      conf = data.frame(taxon = rownames(conf), conf)
    }else{
      colnames(conf) = rownames(conf) =  ee$scientificName
      conf = data.frame(taxon = rownames(conf), conf)
    }
    write_csv((conf), paste("./outdir/cm/", dm, "confusion_matrix_a_all.csv", sep=""))
    if(ncol(data.frame(cmdm$byClass)) == 1){
      macroF1[[dm]] = mean(cmdm$byClass["F1"], na.rm=T)
    }else{
      macroF1[[dm]] = mean(cmdm$byClass[,"F1"], na.rm=T)
    }
    #mcm = as.matrix.data.frame(cmdm$table)
    #rownames(mcm) = colnames(mcm) = colnames(cmdm$table)
    microF1[[dm]] <- cmdm$overall[1]
    #macroF1[[dm]] = mean(cmdm$byClass[,"F1"], na.rm=T)
    cm[[dm]] = cmdm
  }
  microF1_dom = unlist(microF1)
  macroF1_dom = unlist(macroF1)

  cm_dom = cm

  #site_level stats
  cm = list()
  microF1 = macroF1 = list()
  for(dm in unique(pairs$siteID)[-23]){
    dm_dt = pairs %>% filter(siteID == dm)
    lvls = unique(c(unique(dm_dt$obs), unique(dm_dt$pred)))
    dm_dt$obs = factor(dm_dt$obs, levels = lvls)
    dm_dt$pred = factor(dm_dt$pred, levels = lvls)

    cmdm = confusionMatrix(dm_dt$pred, dm_dt$obs)
    print(cmdm$overall)
    conf  = as.matrix(cmdm)
    conf = conf %>% t
    # missing_ids = apply(rbind(apply(conf, 1, sum), apply(conf, 2, sum)), 2, sum)==0
    # conf =  data.frame(conf)[!missing_ids, !missing_ids]
    ee = vst %>% select(taxonID, scientificName) %>% unique
    tp = do.call(rbind.data.frame, strsplit(ee$scientificName, split = " "))[,1:2]
    ee$scientificName = paste(tp[,1], tp[,2])
    ee = unique(ee)
    ee = ee %>% filter(scientificName != "Pinus glabra") %>%
      filter(scientificName != "Fraxinus americana/pennsylvanica")

    ee = left_join(data.frame(taxonID = rownames(conf)), ee)
    eeb = ee
    colnames(eeb)[1] = "obs"
    ids_and_genus = pairs %>% left_join(eeb)
    colnames(ids_and_genus)[7]='obs_scientific'
    colnames(eeb)[1] = "pred"
    ids_and_genus = ids_and_genus %>% left_join(eeb)
    colnames(ids_and_genus)[c(3,4,8)]=c('obs_sp','pred_sp', 'pred_scientific')

    colnames(conf) = rownames(conf) =  ee$scientificName
    conf = data.frame(taxon = rownames(conf), conf)

    #conf = cbind.data.frame(rownames(conf), conf)
    write_csv(data.frame(conf), paste("./outdir/cm/", dm, "confusion_matrix_a_all.csv", sep=""))

    #mcm = as.matrix.data.frame(cmdm$table)
    #rownames(mcm) = colnames(mcm) = colnames(cmdm$table)
    microF1[[dm]] <- cmdm$overall[1]
    if(ncol(data.frame(cmdm$byClass)) == 1){
      macroF1[[dm]] = mean(cmdm$byClass["F1"], na.rm=T)
      macroF1[[dm]] = F1_Score_macro((dm_dt$obs), (dm_dt$pred))
      microF1[[dm]] = F1_Score_micro((dm_dt$obs), (dm_dt$pred))
    }else{
      macroF1[[dm]] = F1_Score_macro(unlist(dm_dt$obs), unlist(dm_dt$pred))
      microF1[[dm]] = F1_Score_micro(unlist(dm_dt$obs), unlist(dm_dt$pred))
    }
    cm[[dm]] = cmdm
  }
  microF1_site= (unlist(microF1))
  #  microF1_site[["ABBY"]] = 1
  macroF1_site= unlist(macroF1)
  macroF1_site[["ABBY"]] = 1
  microF1_site[["ABBY"]] = 1
  cbind.data.frame(macroF1_site, microF1_site)
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
  #write_csv(sp_co, "./outdir/species_confusion_june23.csv")

  figure_1(pairs, vst)

  #plot distribution of confusion by health class
  figure_2(pairs, vst)

  #plot geographic accuracy
  plot_figure_3(microF1_dom, microF1_site, macroF1_site, species_per_site)

  figure_4()


  # uncertainty quantitative analysis using ranking
  # Can we look at it more continuosly: 0-10 which fraction correctly classified?
  figure_6(pairs, probabilities)
  #family confusion
  pairs = readr::read_csv(pairs_fam_pt)
  pairs = cbind.data.frame(pairs, pp)
  colnames(pairs)=c("id", "obs", "pred", "i", "individualID", "obs_sp", "pred_sp")

  pairs = pairs %>% filter(individualID %in% ids_and_genus$individualID)
  ids_and_genus$obs_scientific = do.call(rbind.data.frame, strsplit(ids_and_genus$obs_scientific, split = " "))[,1]
  ids_and_genus$pred_scientific = do.call(rbind.data.frame, strsplit(ids_and_genus$pred_scientific, split = " "))[,1]
  dm_dt = ids_and_genus
  genuses = c(dm_dt$obs_scientific, dm_dt$pred_scientific) %>% unique
  #lvls = unique(c(unique(genuses$obs), unique(genuses$pred)))
  dm_dt$obs = factor(dm_dt$obs_scientific, levels = unique(genuses))
  dm_dt$pred = factor(dm_dt$pred_scientific, levels = unique(genuses))
  cmdm = confusionMatrix(dm_dt$obs, dm_dt$pred)
  print(cmdm$overall)

  conf  = as.matrix(cmdm)
  conf = conf %>% t
  conf = cbind.data.frame(rownames(conf), conf)
  write_csv(data.frame(conf), paste("./outdir/cm/Family_level_confusion_matrix.csv", sep=""))


  sp_co = cmdm$table
  sp_co = cmdm$table %>% data.frame
  sp_co = sp_co %>% filter(Freq > 0)
  write_csv(sp_co, "./outdir/genuses_confusion.csv")



}
