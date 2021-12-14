library(sf)
clean_spectra <- function(brick){
  # filter for no data 
  mask = brick > 10000
  brick[mask] <-  NA
  mask = brick == -9999
  brick[mask] <-  NA
  
  #filter for greennes and shadows
  ndvi <- (brick[,"band_90"]- brick[,"band_58"])/(brick[,"band_58"] + brick[,"band_90"]) <0.3
  nir860 <- (brick[,"band_96"] + brick[,"band_97"])/2 < 0.1
  mask = as.logical(ndvi | nir860)
  mask[is.na(mask)] = T
  brick[mask,] = NA
  rm(mask, ndvi, nir860)
  brick = brick[,15:365]
  normMat <- sqrt(apply(brick^2,FUN=sum,MAR=1, na.rm=TRUE))
  normMat <- matrix(data=rep(normMat,ncol(brick)),ncol=ncol(brick))
  brick=brick/normMat
  rm(normMat)
  
  #filter for known artifacts
  cnd = (brick[,"band_312"] > 0.03)
  cnd[is.na(cnd)] = T
  #idx <- which(apply(cnd, 1, any))
  brick[cnd,] = NA
  
  cnd = (brick[,24:45] > 0.03)
  idx <- (apply(cnd, 1, any))
  if(length(idx) !=0){
    idx[is.na(idx)] = T
    brick[idx,] = NA
  }
  cnd = (brick[,195:200] > 0.043)
  idx <- (apply(cnd, 1, any))
  if(length(idx) !=0){
    idx[is.na(idx)] = T
    brick[idx,] = NA
  }
  rm(cnd,idx)
  
  # # save pixel positions
  good_pix = !is.na(brick)
  good_pix = (apply(good_pix, 1, all))
  # 
  brick = brick[complete.cases(brick),]
  return(list(refl = brick, good_pix = good_pix))
}


#dataset = readr::read_csv("../NeonSpeciesClassification/dimensionality_reduction/hsi_appended.csv")
hpca_predict <- function(dataset, kld_obj, method = "hist"){
  hsi = dataset %>% dplyr::select(matches("band"))
  kld_array <- kld_obj$bands_grouping
  #perform PCA for each group and select first component
  grp = cbind.data.frame(kld_array, t(hsi))
  if( method=='hist'){
    #get min, max and auc
    kld_refl = list()
    for(gp in unique(kld_array)){
      pcx = grp %>% dplyr::filter(kld_array == gp) 
      pcx = pcx[-1]%>% t 
      min_sband = apply(pcx, 1, min)
      max_sband = apply(pcx, 1, max)
      # scale and sum (value-min)/(max-min)
      pcxgrp =  kld_obj$compression[,gp] #c(min(pcx), max(pcx))
      auc_sband = apply(pcx, 1, function(x)(x-pcxgrp[1])/(pcxgrp[2]-pcxgrp[1]))
      auc_sband =  apply(pcx, 1, sum)
      auc_sband = cbind(min_sband, max_sband, auc_sband)
      colnames(auc_sband) = paste(c("min_kl_", "max_kl_", "auc_kl_"), gp, sep="")
      kld_refl[[gp]] = auc_sband
    }
    kld_refl = do.call(cbind.data.frame, kld_refl)
  }
  return(kld_refl)
}

#get tile with crown predictions from raster
year = 2019
cc = list.files("/Volumes/Stele/crow_maps/predictions", pattern = "727000_4701000_image.shp", full.names = T)
cc = subset(cc, grepl(as.character(year), cc))
tile_itc = sf::read_sf(cc)
tile_itc = st_centroid(tile_itc)

#extract pixels around centroid
brick = brick("~/Documents/Data/HS/HARV_725000_4696000_.tif")
  
#clean spectra and apply transformation
itcs_pixels = raster::extract(brick, tile_itc, buffer = 1, df = T)
colnames(itcs_pixels) = c("individualID", paste("band", 1:369, sep="_"))
data =clean_spectra(itcs_pixels)
#itcs_pixels = itcs_pixels[data$good_pix,]
data = cbind.data.frame(itcs_pixels$individualID[data$good_pix], data$refl)
#load data reduction 
cfc_reduced_dims = readRDS("/Users/sergiomarconi/Documents/GitHub/traitsMaps/indir/kld_hist_trees_30.rds")
foo = hpca_predict(data, cfc_reduced_dims)
foo = cbind.data.frame(data[1], foo)
colnames(foo)[1]="individualID"
foo["siteID"] = siteID
foo["domainID"] = domainID

write_csv(foo, "./weak_label/pred_indir/HARV_725000_4696000.csv")
#itcs_pixels = read_csv("//Volumes/Stele/example_itcextract.csv")
#run prediction




#append reult
preds = read_csv("./weak_label/pred_out/HARV_725000_4696000.csv")
colnames(preds) = c("id", "individualID", "taxonID")
tile_itc$individualID = 1:nrow(tile_itc)
species_predictions = inner_join(tile_itc, preds)
write_sf(species_predictions, "./weak_label/pred_out/HARV_725000_4696000.shp")
