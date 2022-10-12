#Figure 1: comparison of individual site - cross site application
library(data.table)
data = fread("/Volumes/Data//NeonSpeciesClassification/Outputs/comparison_single_site_to_all_updated.csv")
data = fread("/Users/sergiomarconi/Desktop/CM/General_toSite_comparison.csv")
vst = fread("./indir/vst.csv")
site_coords =  vst %>% select(siteID, latitude) %>% group_by(siteID) %>% summarize_all(mean)
box_ = data %>% select(-one_of("Genus level", "Species"))
box_ = reshape2::melt(box_)

pairs = readr::read_csv(pairs_sp_pt)
probabilities = readr::read_csv(pairs_sp_prob)
colnames(pairs)=c("id", "individualID", "obs", "pred")

pairs$domainID = substr(pairs$individualID, 10,12)
pairs$siteID = substr(pairs$individualID, 14,17)

# ggplot(box_, aes(y=value, x  =Site, fill=variable))+
#   geom_bar(stat="identity",position=position_dodge()) + 
#   facet_grid(Score~.)+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         panel.grid.major.y = element_blank()) + 
# ggsci::scale_fill_jco()
train_set = fread("./indir/metadata.csv")
train_set = train_set %>% filter(groupID == "train") %>%
  select(taxonID, siteID) %>% unique
dt = train_set$siteID %>% table %>% data.frame
train_set %>% table 
#convert wide to long format
plotDat <- reshape::melt(table(train_set))

#and plot
ggplot(plotDat, aes(taxonID, siteID, fill = c("white", "blue")[ value + 1 ])) +
  geom_tile() +
  scale_fill_identity() +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1))

ggsave(filename = "./Figures/Figure_s10.png", width = 9, height = 4 , dpi = 300)
 colnames(dt)[1] = "Site"
#data$Delta = data$Generalized - data$Site_level
trends_comparison = data #%>% select(No_elevation, d_no_elevation, Score, Species, Tested_om)

ave_delta_macro = mean(unlist(trends_comparison[trends_comparison$Score == "MacroF1", "Delta"]))
ave_delta_micro =  mean(unlist(trends_comparison[trends_comparison$Score == "MicroF1", "Delta"]))

single_site = list.files("./outdir/site_and_overall_results/", pattern = "pairs.csv")[1:28]
macroF1_idst = microF1_idst = compiled_dat = list()
for(st in single_site[-c(10,12)]){
  #dat = read_csv(paste("./Outputs/", st, sep=""))
  dm_dt = pairs %>% dplyr::filter(siteID == substr(st,2,5))
  dat = fread(paste("./outdir/site_and_overall_results/", st, sep="")) %>%
    filter(!is.na(V1))
  dat = dat %>% filter(V2 %in% dm_dt$individualID)
  lvv = unique(c(dat$V3,dat$V4))
  dat$V3 = factor(dat$V3, levels = lvv)
  dat$V4 = factor(dat$V4, levels = lvv)
  
  cmdm = confusionMatrix(dat$V3, dat$V4)
  compiled_dat[[st]] = cbind.data.frame(substr(st, 2,5), dat)
  macroF1_idst[[st]] = F1_Score_macro(unlist(dat[,3]), unlist(dat[,4]))
  microF1_idst[[st]] = F1_Score_micro(unlist(dat[,3]), unlist(dat[,4]))
}
compiled_dat =  do.call(rbind.data.frame, compiled_dat)
macr_site_bysite = F1_Score_macro(unlist(compiled_dat[,3]), unlist(compiled_dat[,4]))
micr_site_bysite = F1_Score_micro(unlist(compiled_dat[,3]), unlist(compiled_dat[,4]))


colnames(compiled_dat) = c("siteID", "n", "ID", "obs", "pred")
compiled_dat = compiled_dat %>% filter(siteID != "BRDF")
compiled_dat$correct  = T
compiled_dat$correct[compiled_dat$obs != compiled_dat$pred]=F

general = compiled_dat %>% filter(siteID == "LL_f")
compiled_dat = compiled_dat %>% filter(siteID !=  "LL_f")

general = general %>% select(ID, correct)
colnames(general)[2] = c("general")
gg = left_join(compiled_dat,  general)



gg = gg %>% group_by(obs, siteID) %>% mutate(frac_site = sum(correct)/n(), 
                                                                 frac_general = sum(general)/n(), 
                                                                 tot = n())
gg$delta = gg$frac_site - gg$frac_general
double_check = gg %>% select(siteID, obs, delta, tot) %>% unique
double_check = left_join(double_check, site_coords)
double_check$siteID = factor(double_check$siteID, levels = unique(double_check$siteID[order(double_check$latitude)]))

ggplot(double_check, aes(y = siteID, x = obs))+
  geom_raster(aes(fill=delta)) + 
  scale_fill_gradient2()+
#scale_fill_gradient(low="grey90", high="red") +
  labs(x="Taxon", y="Site") + 
  #theme_dark() + 
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1),
        panel.grid.major.y = element_blank(), legend.position = "bottom")

  


ggplot(trends_comparison, aes(x = Delta, color = Score)) + 
  geom_density(aes(x = Delta, color = Score, fill=Score, alpha=0.5))+
  geom_rug()+
  #geom_smooth(fullrange=FALSE,  aes(fill = Score), alpha = 0.15)+
  #  geom_point(size=2, aes( aes(x = Delta, y = Tested_om, fill = Score)), colour="black",pch=21) +
  theme_bw() + geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.77- micr_site_bysite, linetype = "dashed", color = "yellow") +
  geom_vline(xintercept = 0.55 -macr_site_bysite, linetype = "dashed", color = "lightblue") +
  theme(panel.grid = element_blank()) + 
  ggsci::scale_color_jco()+ggsci::scale_fill_jco()

reshuffled_sites = trends_comparison
reshuffled_sites = fread("~/Desktop/CM/General_to_site_comparison_updated.csv")
colnames(reshuffled_sites)[1] = "Site"
#get_list_sp = num_of_sp %>% select(siteID, full_train)
reshuffled_sites = left_join(reshuffled_sites, num_of_sp)
f1_order = reshuffled_sites[["Site"]]
reshuffled_sites[order(Delta, decreasing = F)]
sites_ordered_by_micro = reshuffled_sites %>% filter(Score == "MicroF1") 
sites_ordered_by_micro = sites_ordered_by_micro[order(Delta, decreasing = T)]

reshuffled_sites = left_join(reshuffled_sites, dt)
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
  geom_text(aes(x = Site, y = -0.75, label = tr_tst), size=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1),
        panel.grid.major.y = element_blank(), legend.position = "bottom")



#check accuracyby taxa
