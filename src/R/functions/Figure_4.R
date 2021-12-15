figure_4 = function(
  pt_features_importance = "./outdir/features_importance_table.csv"
  , original_feat_pt = "./outdir/features.csv"
  , kld_feat_pt = "./indir/kld_transformed_example.csv"
  , wlg_pt = "./indir/wavelengths.csv"
  , feat_imp_per_feat = "./outdir/features_importance_per_feat.csv"
  , kld_grp_pt = "./indir/kld_groups.csv"
){
  
  library(ggsci)
  variables_importance = readr::read_csv(pt_features_importance)
  variables_importance = variables_importance %>% group_by(Variable) %>%
    mutate(average = mean(Importance), std = sd(Importance))
  variables_importance = variables_importance %>% select(Variable, average, std) %>%unique
  
  

  ggplot(data=variables_importance, aes(x=factor(Variable, levels = c("Site", "Elevation", "Geolocation", "Reflectance")),
                                        y=average, fill = Variable)) + 
    geom_bar(stat="identity")+ylim(0,1)+ scale_fill_jco()+
    geom_errorbar(aes(ymin=average-std, ymax=average+std), width=.2,
                  position=position_dodge(.9))+  theme_bw()+
    theme(legend.position = "none", panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave("./Figures/figure_4_panel_a.png")
  # get spectral transformed
  x_tst = fread(kld_feat_pt, nrows = 3) %>% data.frame
  features_importance = readr::read_csv(feat_imp_per_feat)[1:45,]
  features_importance$run5 = as.numeric(features_importance$run5)
  features_importance[-1]= apply(features_importance[-1], 2, function(x)x/sum(x))
  features_importance$ave_imp= apply(features_importance[-1], 1, function(x)mean(x))
  kld_grp= fread(kld_grp_pt)
  kld_grp$band = paste("band", 11:357, sep="_")
  bands_nm = fread(wlg_pt)
  noise_intervals = bands_nm %>% filter(noise ==1)
  bands_nm = bands_nm %>% filter(noise == 0)
  bands_nm$BandName = paste("band", 1:367, sep="_")
  bands_nm = bands_nm[11:357,]
  bands_nm = cbind.data.frame(bands_nm, kld_grp)
  
  x_mean = apply(x_tst[,c(3:47)], 2, mean)
  x_sd =  apply(x_tst[,c(3:47)], 2, sd)
  V1 = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,
         7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,
         12,12,12,13,13,13,14,14,14,15,15,15)
  kl_reflectance = cbind.data.frame(V1, x_mean, x_sd)
  features_importance = cbind.data.frame(V1, features_importance, x_mean, x_sd)
  #df = inner_join(bands_nm[-1], kl_reflectance)
  df = inner_join(bands_nm[-1], features_importance)
  df = df %>% group_by(V1) %>% 
    mutate(ave_feat = mean(ave_imp)*3)
  
  roi = df %>%ungroup %>%  select(nanometer, V1) %>% unique
  # = reshape2::melt(df)
  #df = df %>% filter(!band %in% c("band_315", "band_316"))
  
  tmp = fread(original_feat_pt, nrows = 3)#[,-c(1:5,353:355)] %>% t
  tmp = tmp[,-c(1)] %>% t
  one_pix = cbind.data.frame(bands_nm[-1], tmp)
  colnames(one_pix)[6] = "pix_reflectance"
  one_pix$pix_reflectance = one_pix$pix_reflectance/sqrt(sum(one_pix$pix_reflectance^2))
  
  df = df %>% filter(!(nanometer> 1.78 & V1==7))
  ggplot(df, aes(x = nanometer, y = x_mean, group=factor(V1))) +
    theme_bw() + geom_hline(yintercept = 0)+
    geom_boxplot(data = df, color = "darkblue", fill = "dodgerblue4", aes(alpha = ave_feat), position = "identity")+
    geom_line(data = one_pix, aes(x = nanometer, y = pix_reflectance,linetype = "dotted"))+
    #scale_x_continuous(breaks = c(0.4,1.33,1.44,1.78,1.95,2.48))+
    theme(legend.position = "bottom", panel.grid = element_blank())
  ggsave("./Figures/figure_4_panel_b.png")
  
}
