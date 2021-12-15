split_train_test = function(unique_entries){
  tr_set = tst_set = list()
  for(st in unique(unique_entries$siteID)){
    token =1
    previous_best = 10
    while(token !=0){
      tmp = unique_entries %>% filter(siteID == st)
      train_plt = tmp %>% ungroup %>% select(plotID) %>% unique  %>% sample_frac(0.7)
      train = tmp %>% filter(plotID %in% unique(train_plt$plotID))
      test = tmp %>% filter(!plotID %in% unique(train_plt$plotID))
      test = test %>% filter(taxonID %in% unique(train$taxonID))
      taxa_in_train = train %>% ungroup %>% select(individualID, taxonID) %>% 
        unique %>% select(taxonID) %>% unique %>% unlist  %>% length
      taxa_in_test = test %>% ungroup %>% select(individualID, taxonID) %>% 
        unique %>% select(taxonID) %>% unique %>% unlist %>% length
      token = token +1
      if(taxa_in_train - taxa_in_test ==0){
        token = 0
        tst_set[[st]] = test
        tr_set[[st]] = train
      }
      if(token >100){
        token = 0
        tst_set[[st]] = test
        tr_set[[st]] = train
      }
      if(taxa_in_train - taxa_in_test < previous_best){
        tst_set[[st]] = test
        tr_set[[st]] = train
        previous_best = taxa_in_train - taxa_in_test
      }
    }
  }
  test = do.call(rbind.data.frame, tst_set)
  train = do.call(rbind.data.frame, tr_set)
  train$groupID = "train"
  test$groupID = "test"
  return(new_sp = rbind(train, test))
}
