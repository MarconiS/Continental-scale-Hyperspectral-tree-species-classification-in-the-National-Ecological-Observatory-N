figure_1 <- function(pairs, data){
  single_site = list.files("./Predictions/", pattern = "pairs.csv")
  nobs = ovrll = macroF1_idst = microF1_all = sp_tested_in_site = sp_trained_in_site = 
    sp_trained_in_general =  macroF1_all = microF1_idst = full_all = full_site = list()
  for(dm in unique(pairs$siteID)){#GUAN, DSNY, YELL, MOAB,  WREF, BONA, [c(10:19, 22:23, 25,27)]
    if(!dm %in% c("KONZ","ABBY", "HEAL", "PUUM")){
      dm_dt = pairs %>% filter(siteID == dm)
      st = single_site[single_site %like% dm]
      dat = fread(paste("./Predictions/", st, sep="")) %>%
        filter(!is.na(V1)) %>% filter(V2 %in% dm_dt$individualID)
      dat = dat %>% filter(V2 %in% dm_dt$individualID)
      #dm_dt = dm_dt  %>% filter(individualID %in% dat$V2)
      colnames(dat) = colnames(dm_dt)[1:4]
      nobs[[dm]] = print(cbind(nrow(dm_dt), nrow(dat)))
      full_all[[dm]] = dm_dt
      full_site[[dm]] = dat
      lvls = unique(c(unique(dat$obs), unique(dat$pred)))
      lvls_ = unique(c(unique(dm_dt$obs), unique(dm_dt$pred)))
      
      sp_tested_in_site[[st]] = dm_dt %>% dplyr::filter(siteID == dm) %>%#, groupID == "test") %>% 
        select(obs) %>% unique
      sp_trained_in_site[[st]] = data %>% dplyr::filter(siteID == dm, groupID == "train") %>% 
        select(taxonID) %>% unique
      sp_trained_in_general[[st]] = data %>% dplyr::filter(siteID == dm) %>% 
        select(taxonID) %>% unique
      # sp_tested_only_site[[st]] = dat %>% #dplyr::filter(siteID == dm, groupID == "train") %>% 
      #   select(taxonID) %>% unique
      if(length(lvls)>1){
        dat$obs = factor(dat$obs, levels = lvls)
        dat$pred = factor(dat$pred, levels = lvls)
        dm_dt$obs = factor(dm_dt$obs, levels = lvls_)
        dm_dt$pred = factor(dm_dt$pred, levels = lvls_)
        cmdm = confusionMatrix(dat$pred, dat$obs)
        cmdm_ = confusionMatrix(dm_dt$pred, dm_dt$obs)
        cmdm_$table
        cmdm$table
        
        ovrll[[dm]] = print(cmdm$overall)
        conf  = as.matrix(cmdm)
        conf = conf %>% t
        
        #reference  = confusionMatrix(dm_dt$pred, dm_dt$obs)
        #missing_ids = apply(rbind(apply(conf, 1, sum), apply(conf, 2, sum)), 2, sum)==0
        #conf =  data.frame(conf)[!missing_ids, !missing_ids]
        
        ee = vst %>% select(taxonID, scientificName) %>% unique
        tp = do.call(rbind.data.frame, strsplit(ee$scientificName, split = " "))[,1:2]
        ee$scientificName = paste(tp[,1], tp[,2])
        ee = unique(ee) 
        ee = ee %>% filter(scientificName != "Pinus glabra") %>%
          filter(scientificName != "Fraxinus americana/pennsylvanica")
        ee = left_join(data.frame(taxonID = rownames(conf)), ee)
        colnames(conf) = rownames(conf) =  ee$scientificName
        conf = data.frame(taxon = rownames(conf), conf)
        colnames(conf)[1] = "rows:observation\\columns:predictions"
        write_csv((conf), paste("./outdir/cm/siteonly/old_", dm, "confusion_matrix_site_only.csv", sep=""))
        
        #write the confusion matrix for each site from the general model
        conf_  = as.matrix(cmdm_)
        conf_ = conf_ %>% t
        ee = vst %>% select(taxonID, scientificName) %>% unique
        tp = do.call(rbind.data.frame, strsplit(ee$scientificName, split = " "))[,1:2]
        ee$scientificName = paste(tp[,1], tp[,2])
        ee = unique(ee) 
        ee = ee %>% filter(scientificName != "Pinus glabra") %>%
          filter(scientificName != "Fraxinus americana/pennsylvanica")
        ee = left_join(data.frame(taxonID = rownames(conf_)), ee)
        colnames(conf_) = rownames(conf_) =  ee$scientificName
        conf_ = data.frame(taxon = rownames(conf_), conf_)
        colnames(conf_)[1] = "rows:observation\\columns:predictions"
        write_csv((conf_), paste("./outdir/cm/general/", dm, "confusion_matrix_site_general.csv", sep=""))
        
        #calculate F1 site
        
        if(dm %in% c("DSNY", "YELL")){
          #in case of 2 classes if one is never predicted caret only returns 1 class. This is not accurate, 
          # because one has 0 F1 (0 precision and 0 recall, returning a NA and in case of 2 classes being )
          macroF1_idst[[st]] = cmdm$overall["Accuracy"]/2
          macroF1_all[[st]] = cmdm_$overall["Accuracy"]/2
          #calulate F1 all
          #dm_dt = dm_dt %>% filter(individualID %in% unlist(dat[,2]))
        } else if(dm %in% c("MOAB",  "WREF", "BONA")){
          cf = table(dm_dt$obs,dm_dt$pred)
          f1_id = list()
          for(tx in colnames(cf)){
            pr = cf[rownames(cf)==tx,colnames(cf)==tx]/sum(cf[,colnames(cf)==tx])
            rc = cf[rownames(cf)==tx,colnames(cf)==tx]/sum(cf[colnames(cf)==tx,])
            f1_id[[tx]] = 2*pr*rc/(pr+rc)
          }
          macroF1_all[[st]] = mean(unlist(f1_id))
          
          cf = table(dat$obs,dat$pred)
          f1_id = list()
          for(tx in colnames(cf)){
            pr = cf[rownames(cf)==tx,colnames(cf)==tx]/sum(cf[,colnames(cf)==tx])
            rc = cf[rownames(cf)==tx,colnames(cf)==tx]/sum(cf[colnames(cf)==tx,])
            f1_id[[tx]] = 2*pr*rc/(pr+rc)
          } 
          macroF1_idst[[st]] = mean(unlist(f1_id))
          # these are perfectly fine, but 2 class cases. Need to calculate the macro F1 manually
        } else {
          cmdm$byClass[is.na(cmdm$byClass[,"F1"]),"F1"] = 0
          cmdm_$byClass[is.na(cmdm_$byClass[,"F1"]),"F1"] = 0
          
          macroF1_all[[st]] = mean(cmdm_$byClass[,"F1"])#F1_Score_macro(unlist(dm_dt$obs), unlist(dm_dt$pred))
          macroF1_idst[[st]] = mean(cmdm$byClass[,"F1"])#F1_Score_macro(unlist(dat[,3]), unlist(dat[,4]))
        } 
          
        microF1_idst[[st]] = (cmdm$overall["Accuracy"])#F1_Score_micro(unlist(dat[,3]), unlist(dat[,4]))
        microF1_all[[st]] = cmdm_$overall["Accuracy"]#F1_Score_micro(unlist(dm_dt$obs), unlist(dm_dt$pred))
      }
      else{
        if(length(lvls_)>1){
          dm_dt$obs = factor(dm_dt$obs, levels = lvls_)
          dm_dt$pred = factor(dm_dt$pred, levels = lvls_)
          cmdm_ = confusionMatrix(dm_dt$pred, dm_dt$obs)
          cmdm_$table
          cf = table(dm_dt$obs,dm_dt$pred)
          f1_id = list()
          for(tx in colnames(cf)){
            pr = cf[rownames(cf)==tx,colnames(cf)==tx]/sum(cf[,colnames(cf)==tx])
            rc = cf[rownames(cf)==tx,colnames(cf)==tx]/sum(cf[colnames(cf)==tx,])
            f1_id[[tx]] = 2*pr*rc/(pr+rc)
            if(is.na(f1_id[[tx]])) f1_id[[tx]] = 0
          }
          macroF1_all[[st]] = mean(unlist(f1_id))
          microF1_all[[st]] = cmdm_$overall["Accuracy"]#F1_Score_macro(unlist(dat[,3]), unlist(dat[,4]))
          
          
          # Score 1 because it was only one class
          microF1_idst[[st]] = 1#F1_Score_micro(unlist(dat[,3]), unlist(dat[,4]))
          macroF1_idst[[st]] = 1#F1_Score_micro(unlist(dm_dt$obs), unlist(dm_dt$pred))
          
          #ovrll[[dm]] = print(cmdm$overall)
          #conf  = as.matrix(cmdm)
          #conf = conf %>% t
          
          #reference  = confusionMatrix(dm_dt$pred, dm_dt$obs)
          #missing_ids = apply(rbind(apply(conf, 1, sum), apply(conf, 2, sum)), 2, sum)==0
          #conf =  data.frame(conf)[!missing_ids, !missing_ids]
          
          # ee = vst %>% select(taxonID, scientificName) %>% unique
          # tp = do.call(rbind.data.frame, strsplit(ee$scientificName, split = " "))[,1:2]
          # ee$scientificName = paste(tp[,1], tp[,2])
          # ee = unique(ee) 
          # ee = ee %>% filter(scientificName != "Pinus glabra") %>%
          #   filter(scientificName != "Fraxinus americana/pennsylvanica")
          # ee = left_join(data.frame(taxonID = rownames(conf)), ee)
          # colnames(conf) = rownames(conf) =  ee$scientificName
          # conf = data.frame(taxon = rownames(conf), conf)
          #write_csv(data.frame(conf), paste("./outdir/cm/", dm, "confusion_matrix_site_june23.csv", sep=""))
          
        }
        
      }
    }
  }
  ids_dat = do.call(rbind.data.frame, full_site)
  #do.call(rbind.data.frame, full_all)
  abby = pairs %>% filter(siteID == "ABBY")
  ids_dat = rbind.data.frame(ids_dat, abby[,1:4])
  species = unique(c(ids_dat$obs, ids_dat$pred))
  ids_dat$obs = factor(ids_dat$obs, levels = species)
  ids_dat$pred = factor(ids_dat$pred, levels = species)
  cmdm_ids = confusionMatrix(ids_dat$obs, ids_dat$pred)
  cmdm_ids$byClass %>% summary
  cmdm_ids$overall
  ccmm = data.frame(cmdm_ids$byClass)
  ccmm = data.frame(taxa = rownames(ccmm), data.frame(cmdm_ids$byClass))
  write_csv( ccmm, "~/Documents/site_only.csv")
  #1 class classification
  microF1_idst[["recheck__ABBY_kld_pairs"]] = macroF1_idst[["recheck__ABBY_kld_pairs"]] =microF1_all[["recheck__ABBY_kld_pairs"]] = macroF1_all[["recheck__ABBY_kld_pairs"]] = 1
  #2 class classification. For some reason Caret trows an error, but we have perfect match 
  trends_comparison = rbind.data.frame(data.frame(Delta = unlist(macroF1_all)-unlist(macroF1_idst), 
                                                  Site = substr(names(macroF1_all),10,13), Score = "MacroF1"),
                                       data.frame(Delta = unlist(microF1_all)-unlist(microF1_idst), 
                                                  Site = substr(names(microF1_all),10,13), Score = "MicroF1")
  )
  
  
  ntaxa = data.frame(Site = substr(names(sp_trained_in_site),10,13), 
                     trained_on_so = unlist(lapply(sp_trained_in_site, nrow)),
                     trained_on_a = unlist(lapply(sp_trained_in_general, nrow)),
                     tested_on = unlist(lapply(sp_tested_in_site, nrow)))
  #add species number
  ntaxa$tr_tst = paste(ntaxa$trained_on_a,  ":", 
                       ntaxa$trained_on_so, "/",
                       ntaxa$tested_on
                       , sep="")
  # trends_comparison$tr_tst[3:4] = paste(trends_comparison$tr_tst[3:4], "/",
  #                                       trends_comparison$Species[3:4]
  #                                       , sep="")
  trends_comparison = trends_comparison %>% left_join(ntaxa) 
  
  sites_ordered_by_micro = trends_comparison %>% filter(Score == "MicroF1") %>% data.table
  sites_ordered_by_micro = sites_ordered_by_micro[order(Delta, decreasing = T)]
  
  trends_comparison[["Site"]] = factor(trends_comparison[["Site"]], 
                                       levels = sites_ordered_by_micro$Site)
  trends_comparison$tr_tst[trends_comparison$Site=="ABBY"] = "1:1/1"
  trends_comparison %>%
    ggplot( aes(x = Site, y = Delta, fill = Score))+
    geom_bar(position="dodge", stat = "identity") + theme_minimal() + 
    theme(panel.grid.major.y = element_blank()) + 
    geom_hline(yintercept = 0)+ ylim(-.77,0.75)+
    ggsci::scale_color_jco()+ggsci::scale_fill_jco() +
    geom_text(aes(x = Site, y = -0.75, label = tr_tst), angle = 45, size=3)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1),
          panel.grid.major.y = element_blank(), legend.position = "bottom")
    ggsave("./Figures/Figure_1_id.png")
    
    unlist(microF1_idst) %>% summary
    micr_site_bysite = mean(unlist(microF1_idst))
    macr_site_bysite = mean(unlist(macroF1_idst))
    micr_site_byall = mean(unlist(microF1_all))
    macr_site_byall = mean(unlist(macroF1_all))
    ggplot(trends_comparison, aes(x = Delta, color = Score)) + 
      geom_density(aes(x = Delta, color = Score, fill=Score, alpha=0.5))+
      geom_rug()+
      #geom_smooth(fullrange=FALSE,  aes(fill = Score), alpha = 0.15)+
      #  geom_point(size=2, aes( aes(x = Delta, y = Tested_om, fill = Score)), colour="black",pch=21) +
      theme_bw() + geom_vline(xintercept = 0) +
      geom_vline(xintercept = micr_site_byall- micr_site_bysite, linetype = "dashed", color = "yellow") +
      geom_vline(xintercept = macr_site_byall -macr_site_bysite, linetype = "dashed", color = "lightblue") +
      theme(panel.grid = element_blank()) + 
      ggsci::scale_color_jco()+ggsci::scale_fill_jco()
    
    ggsave("./Figures/Figure_S4.png")
    
    #plot the delta in a matrix form (supplement s.5)
    full_all = do.call(rbind.data.frame,full_all)
    colnames(full_all)[4] = "pred_full"
    
    full_site = do.call(rbind.data.frame,full_site)
    colnames(full_site)[4] = "pred_site"

    gg = left_join(full_all, full_site, by=c("individualID", "obs"))
    gg$correct = gg$obs == gg$pred_site
    gg$general = gg$obs == gg$pred_full
    
    gg = gg %>% group_by(obs, siteID) %>% mutate(frac_site = sum(correct)/n(), 
                                                 frac_general = sum(general)/n(), 
                                                 tot = n())
    gg$delta = gg$frac_site - gg$frac_general
    double_check = gg %>% select(siteID, obs, delta, tot) %>% unique
    # double_check = left_join(double_check, site_coords)
    # double_check$siteID = factor(double_check$siteID, levels = unique(double_check$siteID[order(double_check$latitude)]))
    # 
    ggplot(double_check, aes(y = siteID, x = obs))+
      geom_raster(aes(fill=delta)) + 
      scale_fill_gradient2()+ coord_flip() + 
      #scale_fill_gradient(low="grey90", high="red") +
      labs(x="Taxon", y="Site") + 
      #theme_dark() + 
      theme(panel.grid = element_blank()) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1),
            panel.grid.major.y = element_blank(), legend.position = "bottom")
    
    
    
    
  accuracy_summary = rbind.data.frame(data.frame(F1_all = unlist(macroF1_all), F1_site = unlist(macroF1_idst), 
                                                 Site = substr(names(macroF1_all),10,13), Score = "MacroF1"),
                                      data.frame(F1_all = unlist(microF1_all), F1_site = unlist(microF1_idst), 
                                                 Site = substr(names(microF1_all),10,13), Score = "MicroF1"))
  return(accuracy_summary)
}
