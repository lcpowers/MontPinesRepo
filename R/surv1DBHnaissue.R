## Survival==1, is.na(DBH) issue

montpines <- read_csv("../data/fullannual data C.csv") %>% 
  select(-1) %>% ## Get rid of first column that seems like old row numbers from elsewhere
  # filter(demoyr < 2016) %>% 
  arrange(sitetag,yr)  %>% 
  select(sitetag,Site,yr,TagNo,DBH,surv,lags,resurveyarea,newplt)

surv1.DBHna.rows = filter(montpines,is.na(DBH)&surv==1)
write_csv(surv1.DBHna.rows,"surv1_DBHna_rows.csv")

surv1.DBHna.full = filter(montpines,sitetag%in%surv1.DBHna.rows$sitetag)
write_csv(surv1.DBHna.full, "surv1_DBHna_sitetags.csv")
