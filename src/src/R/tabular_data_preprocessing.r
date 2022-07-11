run_data_filtering = function(
  vst_pt = "./indir/vst.csv"
  , spectral_library_pt = "./indir/spectral_library.csv"
  , elevation_data_pt = "./indir/elevation.csv"
){
  library(tidyverse)
  library(data.table)
  r_sources = list.files("./src/R/functions/")
  sapply(paste("./src/R/functions/", r_sources, sep=""), source)
  csv_vst = fread(vst_pt) #%>% filter(height > 3)
  
  #select live trees
  csv_vst = filter_trees(csv_vst)
  #load reflectance
  brick = data.table::fread(spectral_library_pt)
  brick = brick[complete.cases(brick),]
  brick = brick %>% filter(individualID %in% csv_vst$individualID)
  
  #clear typos and remove trees too close to be assigned to a clip
  metadata = clean_typos_taxonID(csv_vst)
  metadata = remove_too_close_stems(metadata, stem_dist = 3)
  
  #filter trees by height, keeping crowns that reach the top of the canopy
  metadata = dplyr::inner_join(metadata, brick)
  metadata$delta_chm = abs(metadata$height - metadata$CHM)
  metadata = metadata %>% group_by(individualID) %>%
    mutate(mean_dh = mean(delta_chm), id_dh = min(delta_chm))
  metadata = metadata %>% filter(delta_chm < 5)
  metadata$individualID %>% unique %>% length
  metadata = metadata %>% filter(!is.na(band_50))
  ids = metadata %>% select(individualID, siteID) %>% unique
  
  
  # metadata = metadata %>% filter(individualID %in% ids)
  sites = ids$siteID %>% table %>% data.frame
  
  #remove noisy pixels
  set.seed(0)
  noise = apply(metadata[,79:424], 1, function(x)any(x > 1 | x < 0))
  metadata[noise,]=NA
  metadata = metadata[complete.cases(metadata[,79:424]),]
  metadata = metadata %>% filter(band_192 <0.45) %>%
    filter(band_193 <0.45) %>%
    filter(band_150 < 0.9) %>%
    filter(band_199 < 0.4) %>%
    filter(band_262 < 0.27) %>%
    filter(band_271 < 0.38) %>%
    filter(band_272 < 0.38) %>%
    filter(band_277 < 0.45) %>%
    filter(band_325 < 0.3) %>%
    filter(band_358 < 0.25) %>%
    filter(band_359 < 0.3) 
  
  ids = metadata %>% select(individualID, siteID, taxonID) %>% unique
  
  #get only one year per site. When 2 available, get 2019
  species_per_site = ids %>% ungroup %>% select(siteID, taxonID) %>%  unique %>%
    ungroup %>% group_by(siteID) %>% mutate(species = n()) %>% select(-one_of("taxonID")) %>%
    unique
  sites = ids$siteID %>% table %>% data.frame
  metadata %>% ungroup %>% select(year, site) %>% table
  sites_2018 = metadata %>% filter(site %in% c("RMNP", "MLBS","NIWO", "GUAN"))%>%
    filter(year == 2018)
  sites_2019 = metadata %>% #filter(!site %in% c("NIWO"))%>%
    filter(year == 2019)
  
  metadata = rbind.data.frame(sites_2018, sites_2019)
  
  # # train test split: optimize for species
  new_df_rfl = remove_outliers(metadata)
  #remove species with 5 or less trees
  new_df_rfl = filter_too_rare_out(new_df_rfl, unique(new_df_rfl$individualID),min_id = 6)
  species_per_site = new_df_rfl %>% select(siteID, taxonID) %>%  unique %>%
    ungroup %>% group_by(siteID) %>% mutate(species = n()) %>% select(-one_of("taxonID")) %>%
    unique
  
  #remove taxa at genus level
  new_df_rfl = new_df_rfl %>% filter(!taxonID %in% c("THUJA", "ULMUS","BETUL","SASSA",
                                                     "PINUS",  "PICEA", "FRAXI"))
  
  #remove pixels 
  new_df_rfl = new_df_rfl %>% filter(delta_chm < 3)
  remove_sites = new_df_rfl %>% select(siteID, taxonID, individualID) %>% unique
  remove_sites$siteID %>% table
  
  set.seed(1987)
  unique_entries = new_df_rfl %>% group_by(individualID)%>% slice(1)
  
  #train test split optimizing number of species in both train and test for different taxa
  new_sp = split_train_test(unique_entries)
  
  new_sp = left_join(new_df_rfl, new_sp)
  #rebalance overabundant pixels by individual tree
  new_sp %>% filter(taxonID == "ACRU") %>%
    select(individualID) %>% unique %>% nrow()
  ACRU  = new_sp %>% filter(taxonID == "ACRU") %>%
    group_by(individualID) %>% 
    slice_min(n=3, order_by = delta_chm) 
  ACSA3  = new_sp %>% filter(taxonID == "ACSA3") %>%
    group_by(individualID) %>% 
    slice_min(n=5, order_by = delta_chm) 
  PSME  = new_sp %>% filter(taxonID == "PSME") %>%
    group_by(individualID) %>% 
    slice_min(n=4, order_by = delta_chm) 
  PIMA  = new_sp %>% filter(taxonID == "PIMA") %>%
    group_by(individualID) %>% 
    slice_min(n=8, order_by = delta_chm)
  POTR5  = new_sp %>% filter(taxonID == "POTR5") %>%
    group_by(individualID) %>% 
    slice_min(n=3, order_by = delta_chm) 
  PIPA2  = new_sp %>% filter(taxonID == "PIPA2") %>%
    group_by(individualID) %>% 
    slice_min(n=3, order_by = delta_chm)
  PICO  = new_sp %>% filter(taxonID == "PICO") %>%
    group_by(individualID) %>% 
    slice_min(n=8, order_by = delta_chm) 
  
  new_sp_ = new_sp %>% filter(!taxonID %in% c("ACSA3", "ACRU", "PSME",
                                              "POTR5", "PIPA2"))
  # 
  new_sp_ = rbind.data.frame(new_sp_, ACSA3, ACRU, PSME, POTR5, PIPA2) 
  colnames(new_sp_)[colnames(new_sp_) == "elevation"] = "field_elevation"
  #new_sp_ = new_sp_ %>% filter(siteID != "BLAN")
  
  #add elevation
  elevation = fread(elevation_data_pt) %>% group_by(individualID, site) %>% 
    summarize_all(mean, na.rm=T) %>% ungroup%>% select(individualID, elevation)
  elevation = elevation %>% group_by(individualID) %>% 
    summarize_all(mean, na.rm=T) %>% ungroup%>% select(individualID, elevation)
  
  new_sp_ = left_join(new_sp_, elevation)
  plot_level_elevation = new_sp_ %>% select(plotID, elevation) %>% group_by(plotID)%>% summarize_all(mean, na.rm=T)
  new_sp_missings =  new_sp_ %>% filter(is.na(elevation)) %>% select(-one_of("elevation"))
  new_sp_missings = left_join(new_sp_missings, plot_level_elevation)
  new_sp_ = new_sp_ %>% filter(!individualID %in% new_sp_missings$individualID)
  new_sp_ = rbind.data.frame(new_sp_, new_sp_missings)
  new_sp_missings =  new_sp_ %>% filter(is.na(elevation)) %>% select(-one_of("elevation"))
  
  check = new_sp_ %>% group_by(siteID) %>% mutate(ave_delta_ele = mean(field_elevation - elevation, na.rm=T), 
                                                  sd_delta_ele = mean(field_elevation - elevation, na.rm = T)) %>%
    select(siteID, ave_delta_ele, sd_delta_ele) %>% unique
  
  #retrieve elevation for those plots there is no information for. Use the field elevation corrected by the average delta with sensor based
  for(id in unique(new_sp_missings$siteID)){
    tmp = new_sp_missings %>% filter(siteID == id)
    new_sp_[new_sp_$individualID %in% unique(tmp$individualID),"elevation"] = 
      new_sp_[new_sp_$individualID %in% unique(tmp$individualID),"field_elevation"] - 
      check[check$siteID == id, "ave_delta_ele"]
  }
  colnames(new_sp_[c(14, 85:431)])
  #new_sp_ = new_sp_ %>% filter(!taxonID %in%c("QUMI", "QUFA"))
  write_csv(new_sp_[c(5:74, 442:447)], "./outdir/metadata.csv")
  write_csv(new_sp_[c(14, 85:431)], "./outdir/features.csv")
  
}
