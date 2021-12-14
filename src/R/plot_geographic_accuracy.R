#Figure results
# US map, colored by domains by macro F1, site size?
plot_figure_3 = function(microF1_dom, microF1_site, macroF1_site, ){
  library(sf)
  library(usmap)
  library("rnaturalearth")
  library("rnaturalearthdata")
  
  us <- ne_countries(scale = "medium", returnclass = "sp") %>% filter( name == "United States")
  us = st_crop(us, c( xmin= -120.1945, ymin= 18.96392, xmax= 179.78, ymax= 71.40767))
  neon_domain = read_sf("./indir/NEON_Domains.shp")
  neon_sites = vst %>% select(siteID, latitude, longitude)%>%
    group_by(siteID)%>% summarize_all(mean)
  
  # append accuracy to site & domain
  
  names(microF1_dom) = substr(names(microF1_dom), 1,3)
  names(microF1_site) = substr(names(microF1_site), 1,4)
  microF1_dom = data.frame(microF1_dom)
  microF1_dom$DomainID = as.integer(substr(rownames(microF1_dom), 2,3))
  microF1_site = data.frame(microF1_site)
  microF1_site$siteID = rownames(microF1_site)
  macroF1_site = data.frame(macroF1_site)
  macroF1_site$siteID = rownames(macroF1_site)
  neon_sites = inner_join(neon_sites, microF1_site)
  neon_sites = inner_join(neon_sites, macroF1_site)
  
  neon_domain_ = left_join(neon_domain, microF1_dom)
  neon_domain_ = neon_domain_ %>% filter(!is.na(microF1_dom))
  neon_domain_ = neon_domain_ %>% filter(OBJECTID !=50)
  
  neon_sites_ = sf::st_as_sf(neon_sites,coords = c("longitude", "latitude"), crs = 4326)
  colnames(species_per_site) = c("siteID", "count")
  neon_sites_ = left_join(neon_sites_, species_per_site)
  neon_domain_$geom2 = st_centroid(st_geometry(neon_domain_))
  neon_domain_ <- neon_domain_ %>%
    mutate(lat = unlist(map(neon_domain_$geom2,1)),
           long = unlist(map(neon_domain_$geom2,2)))
  ggplot() +
    geom_sf(data = neon_domain_, aes(fill = (microF1_dom))) + 
    geom_sf(data = us, alpha = 0) +
    #xlim(-160,-50)+
    geom_sf(data = neon_sites_, size = 2,  alpha=0.8, colour="white") +
    geom_sf(data = neon_sites_, size = 0.5,  colour="black") +
    #geom_sf(data = neon_sites_, aes(size = 0.1)) +
    scale_fill_viridis_c() + 
    # geom_text(aes(x = lat, y = long,label = round(microF1_dom,2)), 
    #                  data = neon_domain_,  size = 5, vjust = 1.1)+
    # geom_point(data = neon_sites, aes(x = longitude, y = latitude, size = microF1_site, color="white", alpha = 0.3)) + 
    #geom_text(data= neon_sites,aes(x=longitude, y=latitude, label=siteID))+
    #color = "darkblue", fontface = "bold", check_overlap = T, ) +
    theme_bw() + theme(legend.position = "bottom")+
    coord_sf(crs = "+proj=laea +lat_0=30 +lon_0=-20 +ellps=GRS80 +units=m +no_defs ")
  
  #figure 3
  ggplot(neon_sites, aes(x = longitude, y = microF1_site)) + geom_point()+
    geom_text(data= neon_sites,aes(x=longitude, y=microF1_site+0.03, label=siteID),
              color = "darkblue", fontface = "bold", check_overlap = T) +theme_bw() + 
    geom_vline(xintercept = -140) + geom_vline(xintercept = -115) + geom_vline(xintercept = -95)
  #geom_vline(xintercept = -140) + geom_vline(xintercept = -115) + geom_vline(xintercept = -95)
  ggsave("./figures/geographic_f1.png")
  
  neon_sites = left_join(neon_sites, species_per_site)
  
  ggplot(neon_sites, aes(x = count, y = microF1_site)) +
    geom_smooth(method = "gam",)+ ylim(0,1.2) + geom_point() +
    geom_text(data= neon_sites,aes(x=count, y=microF1_site+0.03, label=siteID),
              color = "darkblue", fontface = "bold", check_overlap = T) +theme_bw() 
  ggsave("./figures/micro_f1_per_species.png")
  
  ggplot(neon_sites, aes(y = microF1_site))+geom_density(fill = "darkblue", 
                                                         alpha = 0.2) + theme_bw() + geom_hline(yintercept = 0.77)+ ylim(0,1.2) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())
  ggsave("./figures/micro_f1_site.png")
  ggplot(neon_sites, aes(x = count, y = macroF1_site)) + geom_point()+
    geom_smooth(method = "gam", color = "darkorange2")+ ylim(0,1.2) +
    geom_text(data= neon_sites,aes(x=count, y=macroF1_site+0.03, label=siteID),
              color = "darkorange4", fontface = "bold", check_overlap = T) +theme_bw() 
  ggsave("./figures/macro_f1_per_species.png")
  
  ggplot(neon_sites, aes(y = macroF1_site))+geom_density(fill = "darkorange4", 
                                                         alpha = 0.2) + theme_bw() + geom_hline(yintercept = .547)+ ylim(0,1.2) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())
  ggsave("./figures/macro_f1_site.png")
  
}

#which and why are more confused
#1 n species >> confusion
# tree health >> confusion
figure_2 = function(pairs, vst){
  library(data.table)
  vst_on_test_info = left_join(pairs, vst)
  vst_on_test_info[["classified"]] = "Correctly"
  vst_on_test_info[vst_on_test_info$obs != vst_on_test_info$pred, "classified"] = "Incorrectly"
  tmp_for_plot = vst_on_test_info %>% select(classified, nlcdClass, growthForm, plantStatus, canopyPosition)
  
  tmp_for_plot = reshape2::melt(tmp_for_plot, id.vars = "classified")
  tmp_for_plot = tmp_for_plot %>% group_by(value)%>% mutate(total = n()) 
  tmp_for_plot = tmp_for_plot %>% group_by(value, classified)%>% mutate(misclassified = n()) 
  tmp_for_plot = unique(tmp_for_plot)
  
  tmp_for_plot["fraction"] = tmp_for_plot$misclassified / tmp_for_plot$total
  fofo = tmp_for_plot 
  fofo = fofo %>% filter(classified == "Incorrectly", !is.na(value))
  fofo$value = factor(fofo$value, levels = fofo$value[order(fofo$fraction)])
  fofo = fofo %>% filter(!value %like% "Lost") %>%
    filter(!value %like% "Downed") %>% filter(!value %like% "No longer") %>%
    filter(!value %like% "Scrub") %>% filter(!value %like% "Dead") %>%
    filter(!value %like% "Full") %>% filter(!value %like% "dead") 
  ggplot(fofo, aes(y = fraction, x = value, fill = variable)) +
    geom_bar(position="dodge", stat="identity")  + ylim(-0.03,1)+ #coord_flip()+
    theme_minimal()+
    geom_text(aes(x = value, y = -0.03, label = total), size=3)+
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 45,hjust=1)) + ggsci::scale_fill_jco()
  
}
