single_site = list.files("./outdir/site_and_overall_results/", pattern = "pairs.csv")
macroF1_idst = microF1_idst = compiled_dat = list()

all = fread("/Users/sergiomarconi/Downloads/_ALL_kld_pairs.csv") %>% #paste("./outdir/site_and_overall_results/", single_site[29], sep="")) %>%
  filter(!is.na(V1))
sitesp = fread(paste("./outdir/site_and_overall_results/", single_site[30], sep="")) %>%
  filter(!is.na(V1))
colnames(sitesp) = colnames(all) = c("id", "individualID", "obs", "pred")

all$siteID = substr(all$individualID, 14,17)
sitesp$siteID = substr(sitesp$individualID, 14,17)

all = all %>% filter(!siteID %in% c("HEAL", "PUUM"))
sitesp = sitesp %>% filter(!siteID %in% c("HEAL", "PUUM"))

macroF1_site = F1_Score_macro((sitesp$obs), (sitesp$pred))
microF1_site = F1_Score_micro(sitesp$obs, sitesp$pred)

macroF1_all = F1_Score_macro((all$obs), (all$pred))
microF1_all = F1_Score_micro(all$obs, all$pred)

dtt = list()
for(st in unique(all$siteID)[-24]){
  ptt = list.files("./outdir/site_and_overall_results/", pattern = st, full.names = T)[1]
  dtt[[st]] = fread(ptt) %>%
    filter(!is.na(V1))
  
}

dtt2 = do.call(rbind.data.frame, dtt)
colnames(dtt2) = c("id", "individualID", "obs", "pred")
dtt2$siteID = substr(dtt2$individualID, 14,17)

ddt2_all = all %>% filter(individualID %in% dtt2$individualID)

missing =  all %>% filter(siteID =="ABBY")

ddt2_all2 = rbind.data.frame(ddt2_all,
                             missing )
ddt2_sp2 = rbind.data.frame(dtt2,
                             missing )

F1_Score_macro((ddt2_sp2$obs), (ddt2_sp2$pred))
F1_Score_micro((ddt2_sp2$obs), (ddt2_sp2$pred))

F1_Score_macro((ddt2_all2$obs), (ddt2_all2$pred))
F1_Score_micro((ddt2_all2$obs), (ddt2_all2$pred))

missing$siteID %>% table()
smacro = smicro =cmacro = cmicro =list()
for(st in unique(all$siteID)[-24]){
  
  site_only = ddt2_sp2 %>% filter(siteID == st)
  all_d = pairs %>% filter(siteID==st)
  site_only = site_only %>% filter(individualID %in% all_d$individualID)
  smacro[[st]]=F1_Score_macro((site_only$obs), (site_only$pred))
  smicro[[st]]=F1_Score_micro((site_only$obs), (site_only$pred))
  
  combined = ddt2_all2  %>% filter(siteID == st)
  cmacro[[st]]=F1_Score_macro((combined$obs), (combined$pred))
  cmicro[[st]]=F1_Score_micro((combined$obs), (combined$pred))
}
comparisons = cbind.data.frame(unlist(smacro),unlist(smicro),unlist(cmacro),unlist(cmicro))
macros = cbind.data.frame(unlist(smacro),unlist(cmacro))
micros = cbind.data.frame(unlist(smicro),unlist(cmicro))
macros$Delta = macros$`unlist(cmacro)`-macros$`unlist(smacro)`
micros$Delta = micros$`unlist(cmicro)`-micros$`unlist(smicro)`

micros$Site = rownames(micros)
micros$Score = "MicroF1"
macros$Site = rownames(macros)
macros$Score = "MacroF1"
reshuffled_sites = rbind.data.frame(micros[,3:5], macros[,3:5])

#get_list_sp = num_of_sp %>% select(siteID, full_train)
reshuffled_sites = left_join(reshuffled_sites, num_of_sp)
f1_order = reshuffled_sites[["Site"]]
reshuffled_sites = data.table(reshuffled_sites)
reshuffled_sites[order(Delta, decreasing = F)]
sites_ordered_by_micro = reshuffled_sites %>% filter(Score == "MicroF1") 
sites_ordered_by_micro = sites_ordered_by_micro[order(Delta, decreasing = T)]

reshuffled_sites[["Site"]] = factor(reshuffled_sites[["Site"]], 
                                    levels = sites_ordered_by_micro$Site)

reshuffled_sites$tr_tst = paste(#reshuffled_sites$Freq,  ":", 
  reshuffled_sites$Freq, ":",
  reshuffled_sites$Tested_om
  , sep="")
reshuffled_sites$tr_tst[3:4] = paste(reshuffled_sites$tr_tst[3:4], "/",
                                     reshuffled_sites$Species[3:4]
                                     , sep="")   
ggplot(reshuffled_sites, aes(x = Site, y = Delta, fill = Score))+
  geom_bar(position="dodge", stat = "identity") + theme_minimal() + 
  theme(panel.grid.major.y = element_blank()) + 
  geom_hline(yintercept = 0)+ ylim(-.77,0.75)+
  ggsci::scale_color_jco()+ggsci::scale_fill_jco() +
  #geom_text(aes(x = Site, y = -0.75, label = tr_tst), size=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1),
        panel.grid.major.y = element_blank(), legend.position = "bottom")
