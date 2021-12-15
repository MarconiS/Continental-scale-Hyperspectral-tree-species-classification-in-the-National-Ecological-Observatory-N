#plot precision/recall/accuracy
precision_recall_accuracy = function(){
library(ggsci)
classes = read_csv("/Users/sergiomarconi/Downloads/Full summary species classification - Sheet1.csv")
classes = classes %>% filter(!taxon %in% c( "accuracy", "macro avg",  "weighted avg"))
classes_ = melt(classes)
classes_ = classes_  %>% filter(variable %in% c("precision", "recall", "f1-score"))
order_f1 = classes$taxon[order(classes$`f1-score`)]
classes_$taxon = factor(classes_$taxon, levels = order_f1)
classes_$genus = substr(classes$taxon, 1,2)
ggplot(classes_, aes(x = taxon, y = value, color = variable, alpha = 0.5)) + 
  geom_point( size = 3) + theme_light() + scale_color_jco() + 
  theme(axis.text.x = element_text(angle = 45, vjust =.5),
        panel.grid.major.y = element_blank(), axis.title.x=element_blank())

}

