supplemetary_get_frequencies_species = function(){
  #get scores per domain
  library(tidyverse)
  library(caret)
  library(data.table)
  
  #list_of_species = read_csv("/Volumes//Data/Field_surveys/VST//list_of_species_lifeform.csv")
  pairs = readr::read_csv("/Volumes/Data/speciesClassificationResults/_BRDF_April_kld_pairs.csv")
  probabilities = readr::read_csv("/Volumes/Data/speciesClassificationResults/BRDF_April_kld_probabilities.csv")
  # 
  colnames(probabilities)[which(probabilities[1,-(1:3)] == max(probabilities[1,-(1:3)]))]
  colnames(pairs)=c("id", "individualID", "obs", "pred")
  pairs$domainID = substr(pairs$individualID, 10,12)
  pairs$siteID = substr(pairs$individualID, 14,17)
  
  sp_in_dat = pairs$obs %>% data.frame
  colnames(sp_in_dat) = "taxonID"
  
  vst = fread("/Volumes/Athena/BackUps/Nov2021/actualBKUP/GitHub/neonVegWrangleR/outdir/field_data/neon_vst_data_022021.csv") #%>% filter(height > 3)
  
  vst = vst %>% 
    filter(str_detect(growthForm,"tree")) %>%
    filter(str_detect(plantStatus,"Live")) %>%
    #filter(!canopyPosition %in% c("Full shade", "Mostly shaded")) %>%
    filter(siteID %in% unique(pairs$siteID)) %>%
    #filter(stemDiameter > 5) %>%
    #filter(height > 2) %>%
    group_by(individualID) %>% slice_max(order_by = eventID, n= 1) %>%
    slice(1) %>% ungroup
  
  vst = clean_typos_taxonID(vst)
  
  pairs = pairs %>% filter(!obs %in% c("THUJA", "ILAN", "MORU2", "ABAM")) %>%
    filter(!siteID %in% c("HEAL", "PUUM"))
  #species_per_site = vst %>% filter(taxonID %in% tree_species$taxonID)  
  tot_species = vst %>% select(taxonID, siteID) %>% group_by(siteID) %>% table
  data = data.table::fread("/Volumes/Data/Spectral_library/metadata_1819.csv") %>%
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
  
  tot_number_of_species = vst %>% 
    filter(str_detect(plantStatus,"Live")) %>%
    filter(str_detect(growthForm,"tree")) %>%
    select(taxonID, siteID)
  test_number_of_species = pairs %>% select(obs, siteID) %>% unique
  train_number_of_species = fread("/Volumes/Data/Spectral_library/metadata_0411.csv", select=c("taxonID", "siteID")) %>% unique
  
  train = train_number_of_species$siteID %>% table %>% data.frame
  test = test_number_of_species$siteID %>% table %>% data.frame
  n_sp = tot_number_of_species %>% unique %>% select(siteID) %>% table %>% data.frame
  
  colnames(train) = c("siteID", "Train")
  colnames(test) = c("siteID", "Test")
  colnames(n_sp) = c("siteID", "VST")
  dat = left_join(n_sp, train) %>% left_join(test) %>% left_join(summary_freqs)
  dat$full_train = dat$VST * as.numeric(dat$species_fraction)
  dat
  
  sp_st = reshape2::melt(dat, "siteID")
  sp_st$value = as.numeric(sp_st$value)
  sp_st$value = round(sp_st$value, 2)
  
  sp_st_ = sp_st %>% filter(variable %in% c("abundance_fraction","species_fraction")) %>%
    filter(!siteID == "MOAB")
  
  ordered_sites_by_alpha_div = dat[order(dat$species_fraction),"siteID"][-1]
  num_of_sp = dat %>% select(siteID, VST, full_train)
  num_of_sp$sp_n = paste(as.integer(num_of_sp$VST), as.integer(num_of_sp$full_train), sep="/")
  sp_st_ = left_join(sp_st_, num_of_sp)
  sp_st_$siteID = factor(sp_st_$siteID, levels = ordered_sites_by_alpha_div)
  ggplot(sp_st_, aes(fill=variable, y=value, x=siteID)) + 
    geom_bar(position="dodge", stat="identity") + 
    theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
    geom_text(aes(x = siteID, y = -.05, label = sp_n), size=4) +
    theme(legend.position = "bottom")
  ggsave("./figures/abundance_per_site.png")
  
  table_train = train_number_of_species %>% table %>% data.frame
  table_train$siteID = factor(table_train$siteID, levels = unique(site_coords$siteID[order(site_coords$latitude)]))
  
  ggplot(table_train, aes(x = taxonID, y = siteID))+
    geom_raster(aes(fill=Freq)) + 
    scale_fill_gradient2()+
    #scale_fill_gradient(low="grey90", high="red") +
    labs(x="Taxon", y="Site") + 
    #theme_dark() + 
    theme(panel.grid = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1),
          panel.grid.major.y = element_blank(), legend.position = "bottom")
  ggsave("./figures/speciese_per_site.png")
}
