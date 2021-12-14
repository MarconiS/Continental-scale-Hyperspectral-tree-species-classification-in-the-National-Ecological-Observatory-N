remove_outliers = function(metadata){
  new_df_rfl = list()
  for(ii in unique(metadata$taxonID)){
    taxa_foo = metadata %>% data.frame %>% dplyr::filter(taxonID == ii) #%>% unique
    refl_foo = taxa_foo %>% select(contains("band"))
    if(nrow(taxa_foo)> 4){
      pix_ind = apply(refl_foo, 2, function(x){
        out = boxplot.stats(x)$out
        pix_ind = which(x %in% c(out))
        pix_ind
      })
      if(length(pix_ind)>0){
        for(bnd in ncol(refl_foo)){
          taxa_foo[pix_ind[[bnd]],] = NA
        }
      }
      # #check range of spectra
      # 
      max_ideal =  apply(refl_foo[complete.cases(refl_foo),],
                         MARGIN = 2, function(x)quantile(x, 0.99))
      min_ideal =  apply(refl_foo[complete.cases(refl_foo),],
                         MARGIN = 2, function(x)quantile(x, 0.01))
      
      #filter for outliers: too bright
      cnd = apply(refl_foo, 1,function(x)(x > max_ideal))
      idx <- (apply(data.frame(cnd), 2, any))
      if(length(idx) !=0){
        idx[is.na(idx)] = T
        taxa_foo[idx,] = NA
      }
      #filter for outliers: too dark
      cnd = apply(refl_foo, 1,function(x)(x < min_ideal))
      idx <- (apply(data.frame(cnd), 2, any))
      if(length(idx) !=0){
        idx[is.na(idx)] = T
        taxa_foo[idx,] = NA
      }
    }
    new_df_rfl[[ii]] = taxa_foo
    
    #plot spectra of the species
    taxa_foo = taxa_foo %>% filter(!is.na(band_54))
  }
  new_df_rfl = do.call(rbind.data.frame, new_df_rfl)
  new_df_rfl = new_df_rfl %>% filter(!is.na(band_54))
  return(new_df_rfl)
}