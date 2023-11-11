rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows
# close graphics windows
library(pomp)
library(magrittr)
library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
library(foreach)
library(doParallel)

mif2_panel_Yize_round_01 <- readRDS("mif2_panel_Yize_round_01.rds")

require(doParallel)
run_level <- 1
cores <- switch(run_level, 2, 36)
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)

bake(file = "output/stocks_04(panel_continue)/mif2_panel_Yize_loglik_dummy.rds",{ 
  el <- foreach(i=1:cores,.packages='pomp', .options.multicore=mcopts) %dopar% {
    measlespkg::eval_logLik(mif2_panel_Yize_round_01, ncores =cores, np_pf = 10000,nreps=cores)
  }
  el
})