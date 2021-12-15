filter_trees = function(csv_vst){
  csv_vst =  csv_vst %>% data.frame %>%
    filter(str_detect(plantStatus , "Live")) %>%
    filter(!str_detect(growthForm,"shrub|sapling")) %>%
    group_by(individualID) %>% slice_max(order_by = eventID)
  
  #get the max height for each  tree when available
  csv_no_height = csv_vst %>% filter(is.na(height)) %>% slice_max(stemDiameter)
  csv_no_height = csv_no_height %>% filter(!is.na(canopyPosition)) %>% 
    filter(plantStatus == "Live") %>% filter(stemDiameter > 10)
  csv_vst =  csv_vst %>% slice_max(order_by = height) %>%slice(1)
  csv_vst = rbind.data.frame(csv_vst, csv_no_height) %>%
    filter(canopyPosition %in% c("Full sun", "Open growth", "Partially shaded", NA))
  return(csv_vst)
}
