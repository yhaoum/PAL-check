# nextround_Job.sbat
cores <- 36
mcopts <- list(set.seed=TRUE)
mif_check_paldata_stocksmodel_02 <- readRDS("mif_check_paldata_stocksmodel_02.rds")
RUN_LEVEL = 1
Nmif = switch(RUN_LEVEL, 1, 100)
bake(file = "mif_check_paldata_stocksmodel.rds",{ 
  mifs_global <- foreach(i=1:cores,.packages='pomp', .options.multicore=mcopts) %dopar% {
    continue(mif_check_paldata_stocksmodel_02[[i]],Nmif=1)
  }
  mifs_global
})


