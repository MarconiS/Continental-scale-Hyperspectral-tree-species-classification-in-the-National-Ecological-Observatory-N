from_txt_to_epsg = function(chr_zn){
  chr_zn = 32600+as.integer(substr(chr_zn, 1, nchar(chr_zn)-1))
  return(chr_zn)
}

l2_norm =  function(brick, ndvi = 0.5, nir = 0.2){
  # filter for no data
  brick = brick %>% ungroup %>% dplyr::select(contains("band"))
  mask1 = apply(brick, 1, function(x)all(x>=0))
  mask2 = apply(brick, 1, function(x)all(x<1))
  brick[!as.logical(mask2), ] = NA
  brick[!as.logical(mask1), ] = NA
  brick[brick<0]=0
  brick[,30:40][brick[,30:40]==0] = NA
  
  #filter for greennes and shadows
  ndvi <- (brick[,"band_90"]- brick[,"band_58"])/(brick[,"band_58"] + brick[,"band_90"]) <ndvi
  nir860 <- (brick[,"band_96"] + brick[,"band_97"])/2 < nir
  mask = as.logical(ndvi * nir860)
  mask[is.na(mask)] = T
  brick[mask,] = NA
  rm(mask, ndvi, nir860)
  brick = brick %>% dplyr::select(one_of(paste("band", 15:360, sep="_")))
  normMat <- sqrt(apply(brick^2,FUN=sum,MAR=1, na.rm=TRUE))
  normMat <- matrix(data=rep(normMat,ncol(brick)),ncol=ncol(brick))
  brick=brick/normMat
  
  
  # rm(cnd,idx)
  
  return(brick)
}

clean_typos_taxonID <- function(species_classification_data){
  #make it a function but for testing ok replicating upper loop
  species_classification_data$taxonID[species_classification_data$taxonID == "2PLANT"] = "NA"
  species_classification_data$taxonID[species_classification_data$taxonID == "FRAM2/FRPE"] = "FRAM2"
  species_classification_data$taxonID[species_classification_data$taxonID == "JUVIV"] = "JUVI"
  species_classification_data$taxonID[species_classification_data$taxonID == "PIGL"] = "PIGL2"
  species_classification_data$taxonID[species_classification_data$taxonID == "PRSES"] = "PRSE2"
  species_classification_data$taxonID[species_classification_data$taxonID == "BEPAP"] = "BEPA"
  species_classification_data$taxonID[species_classification_data$taxonID == "CECAC"] = "CECA4"
  species_classification_data$taxonID[species_classification_data$taxonID == "ACSAS"] = "ACSA3"
  species_classification_data$taxonID[species_classification_data$taxonID == "ACNEN"] = "ACNE2"
  species_classification_data$taxonID[species_classification_data$taxonID == "ACRUR"] = "ACRU"
  species_classification_data$taxonID[species_classification_data$taxonID == "PICOL"] = "PICO"
  species_classification_data$taxonID[species_classification_data$taxonID == "PSMEM"] = "PSME"
  return(species_classification_data)
}

remove_too_close_stems <- function(species_classification_data, stem_dist = 3){
  field_vst = list()
  for (jj in unique(species_classification_data$siteID)){
    new_data = species_classification_data %>% filter(siteID == jj)
    ids = new_data %>% select(individualID, taxonID) %>% data.frame
    #identify tree stems that are closer than 3m apart
    xdist = dist(new_data$itcEasting) %>% as.matrix
    ydist = dist(new_data$itcNorthing) %>% as.matrix
    rownames(xdist) = rownames(ydist) = colnames(xdist) = colnames(ydist) = new_data$individualID
    final_dist = sqrt(xdist^2 + ydist^2)
    index <- which(final_dist < stem_dist, arr.ind=TRUE)
    is_other = index[,1] != index[,2]
    index = index[is_other,]
    #identify if they are diffeernt species
    is_different_species = cbind(rownames(final_dist)[index[,1]], colnames(final_dist)[index[,2]]) %>% data.frame %>% unique
    is_different_species = left_join(is_different_species, ids, by = c("X1"="individualID"))# by = c("X1",  "individualID"))
    is_different_species = left_join(is_different_species, ids, by = c("X2"="individualID"))
    index = is_different_species %>% filter(taxonID.x != taxonID.y) %>% select(X1) %>% unique
    #get only ids with no different species within 3m
    no_too_close = new_data %>% filter(!individualID %in% index$X1)
    field_vst[[jj]] = no_too_close
    #new_data = sf::st_as_sf(no_too_close, coords=c("itcEasting", "itcNorthing"), crs = from_txt_to_epsg(unique(new_data$utmZone)))
    #sf::write_sf(new_data, paste("~/Documents/Data/shp_vst/", jj, ".shp", sep=""))
  }
  species_classification_data = do.call(rbind, field_vst)
  return(species_classification_data)
}

apply_corrections = function(brick){
  brdf_corrections = readr::read_csv("./data/brdf_corrections.csv")
  refl=list()
  for( ii in 1:367){
    y = brick[,ii+1]/10000
    refl[[ii]] = y * bdrf_cor[ii] * topo_cor[ii]
  }
  refl = do.call(cbind.data.frame, refl)
  return(refl)
}
calculate_lat_lon <- function(df){
  latlon = list()
  for(zones in unique(df$utmZone)){
    tmp = df %>% filter(utmZone == zones)
    epsg = 32600+ as.integer(substr(zones, 1, nchar(zones)-1))
    tmp = sf::st_as_sf(tmp, coords = c("itcEasting", "itcNorthing"), crs = epsg)
    tmp = sf::st_transform(tmp, crs=4326)
    latlon[[zones]] = tmp
  }
  latlon = do.call(rbind.data.frame, latlon)
}
clean_reflectance = function(brdf_data, soil_th = NA){
  #filter shaded and NAs
  ancillary = brdf_data %>% select(!contains("band"))
  ids = brdf_data[["individualID"]]
  brdf_data = brdf_data %>% ungroup %>% select(contains("band"))
  if(is.integer(brdf_data[[1]])){
    brdf_data=brdf_data/10000
  }
  if (!is.na(soil_th)){
  }
  #brdf_data[brdf_data < -9000]=NA
  #brdf_data[brdf_data > 12000]=NA
  
  #brdf_data[["individualID"]] = ids
  brdf_data = brdf_data[complete.cases(brdf_data),]
  ancillary = ancillary[complete.cases(brdf_data),]
  
  
  brick = l2_norm(brdf_data, ndvi = 0.4, nir = 0.2)
  brick = cbind.data.frame(ancillary, brick)
  #remove NA bands using a random band as a target
  brick = brick %>% filter(!is.na(band_58))
  #brick = brick[complete.cases(brick),]
  # #apply corrections
  # brick = apply_corrections(brick)
  return(brick)
}

filter_too_rare_out = function(vst, id_in_library, min_id = 5){
  vst = vst %>% data.frame %>% filter(individualID %in% id_in_library)
  vst_ =vst %>% ungroup %>% group_by(individualID) %>% slice(1)
  goodIDS = vst_$taxonID %>% table() %>% data.frame
  goodIDS = goodIDS %>% filter(Freq > min_id) %>% select(".") 
  vst = vst %>% 
    filter(taxonID %in% as.character(goodIDS[[1]]))
  return(vst)
  
}

#df = readr::read_csv("./data/brdf_df_2021_NEON.csv")

plot_spectra <- function(df){
  df = df %>% select(contains("band"))
  df = reshape2::melt(df)
  ggplot(df, aes(x = variable, y = value)) +
    geom_boxplot(size = 0.1, alpha=0.3)+ theme_bw()
}

# df = df %>% group_by(individualID) %>% summarize_all(mean)
# df = df %>% ungroup %>% select(contains("band"))
# df = data.frame(t(df), colnames(df))
# colnames(df) = c("reflectance", "band")
# ggplot(df, aes(x = band, y = reflectance)) + 
#   geom_path(size = 2, color="black") + theme_classic() + theme(text = element_text(size=20),
#          axis.text.x = element_text(angle=45, hjust=1))
# ggsave("~/Documents/spectra_sample.png")
