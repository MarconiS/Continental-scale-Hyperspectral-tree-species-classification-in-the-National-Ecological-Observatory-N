
figure_6 = function(pairs, probabilities){
  #pairs = readr::read_csv(pairs_sp_pt)
  #probabilities = readr::read_csv(pairs_sp_prob)
  wrong_ids = pairs$individualID[!(pairs$obs == pairs$pred)]
  probabilities = inner_join(pairs, probabilities)
  pdist_wrong = probabilities %>% filter(individualID %in% wrong_ids)
  probabilities$max_pid = apply(probabilities[,-c(1:8)],1, max) 
  probabilities$confused = probabilities$obs!=probabilities$pred

  #check if confused high confidence are field misclassified
  get_wrong_high_confidence = probabilities %>% filter(max_pid >0.8, confused==T)
  get_wrong_high_confidence = left_join(get_wrong_high_confidence, vst)
  get_wrong_high_confidence$utmZone
  
  #calculate lat lon coordinates for vst points
  #get_wrong_high_confidence = calculate_lat_lon(get_wrong_high_confidence)
  
  # median distribution for misclassified species by species
  dist_sp = pdist_wrong %>% select(-one_of("id", "individualID", "obs",
                                           "pred",  "domainID", "siteID", "effort", "X1")) %>% group_by(taxonID) %>%
    summarize_each(median)
  dist_sp = reshape2::melt(dist_sp)
  #ggplot(dist_sp, aes(x = variable, y = value)) + geom_point() + facet_wrap(.~taxonID)
  
  highish = dist_sp %>% filter(value > 0.01)
  highish = highish[!highish$taxonID == highish$variable,]
  
  individual_taxa = sort(unique(highish$taxonID))
  highish$taxonID = factor(highish$taxonID, levels = sort(unique(highish$taxonID)))
  
  ggplot(highish[highish$taxonID %in% individual_taxa[1:35],], aes(x = as.character(variable), y = value)) + geom_point() + ylim(0,1)+
    facet_wrap(.~taxonID, drop = T, scales = "free_x") + theme_bw() + geom_hline(yintercept = 0.1) + 
    theme(axis.text.x=element_text(angle =- 45, vjust = 0.5))+
    xlab("Species confused with") + ylab("P(x)")
  ggsave("./Figures/Figure_S9.png")
  
  # make a quantitative analysis using ranking
  # Can we look at it more continuosly: 0-10 which fraction correctly classified? 
  all_pdist = probabilities %>% select(-one_of("id","X1", "individualID", "effort","max_pid", "confused")) 
  p_max = apply(all_pdist[,-c(1:7)], 1, function(x)(max(x)))
  who_max = lapply(1:nrow(all_pdist), function(x)(colnames(all_pdist[which(all_pdist[x,]==p_max[x])])))
  who_max = unlist(who_max)
  
  final_boxes = cbind.data.frame(all_pdist[["taxonID"]], who_max, p_max)
  bin = seq(0,1,by=0.03)
  fraction_ = rep(NA, length(bin))
  for(ii in 2:(length(bin)+1)){
    ith_class = final_boxes %>% filter(final_boxes$p_max < bin[ii] & 
                                         final_boxes$p_max > bin[ii-1])
    fraction_[ii-1] = sum(ith_class[,1] != ith_class[,2])/nrow(ith_class)
    fraction_[is.nan(fraction_)]=NA
  }
  
  uncertainty_curve = cbind.data.frame(bin, fraction_)
  colnames(uncertainty_curve) = c("majority_p", "fraction_misclassified")
  ggplot(uncertainty_curve, aes(x = majority_p, y = 1-fraction_misclassified)) + ylim(0,1) + xlim(0,1)+
    geom_point() + theme_bw() + geom_abline(intercept = 0, slope = 1) + stat_smooth(method="lm", se=FALSE)
  ggsave("./Figures/Figure_6.png")
  
  uu = lm(majority_p~fraction_misclassified, data = uncertainty_curve)
  uu = summary(uu)
  uu$adj.r.squared
  
}