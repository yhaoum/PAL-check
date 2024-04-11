library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(ggplot2)
sourceCpp('rotavirus_sim_EqEq.cpp')

length = 416

x <- c(beta1 = 12.15,
       beta2 = 0.22,
       beta3 = 0.34)
betas <- t(replicate(length, x))
q <- matrix(0.07, nrow=length, ncol=3)
init_dist <- c(3876549, 30351, 1315221, 57139612, 871, 302852, 19573727, 2550, 131092)


avg_1 <- c()
med_1 <- c()
var_1 <- c()
avg_2 <- c()
med_2 <- c()
var_2 <- c()
avg_3 <- c()
med_3 <- c()
var_3 <- c()
for(i in 1:1000){
  set.seed(123+765*i)
  rota_sim <- rotavirus_sim(length, 82372825, init_dist, betas, q, 0.017, 0.022)
  colnames(rota_sim) <- c("cases1", "cases2", "cases3")
  rota_sim <- data.frame(time = c(1:415), rota_sim)
  
  #mean
  avg_1[i] <- rota_sim$cases1 |> mean()
  avg_2[i] <- rota_sim$cases2 |> mean()
  avg_3[i] <- rota_sim$cases3 |> mean()
  #median
  med_1[i] <- rota_sim$cases1 |> median()
  med_2[i] <- rota_sim$cases2 |> median()
  med_3[i] <- rota_sim$cases3 |> median()
  #variance
  var_1[i] <- rota_sim$cases1 |> var()
  var_2[i] <- rota_sim$cases2 |> var()
  var_3[i] <- rota_sim$cases3 |> var()
}

sim_compare_EqEq <- data.frame(avg_1,avg_2,avg_3,med_1,med_2,med_3,var_1,var_2,var_3)
write.csv(sim_compare_EqEq, file = "sim_compare_WWR_EqEq_old.csv")

#######################################################################################
set.seed(133828)
rota_sim <- rotavirus_sim(length, 82372825, init_dist, betas, q, 0.017, 0.022)
colnames(rota_sim) <- c("cases1", "cases2", "cases3")
rota_sim <- data.frame(time = c(1:415), rota_sim)

p_EqEq_WWR <- ggplot(data=rota_sim) +
  geom_line(aes(x=time, y=cases1),color='green') +
  geom_line(aes(x=time, y=cases2), color='red') +
  geom_line(aes(x=time, y=cases3), color='blue') +
  ylab('Y_sim_WWR')+xlab('date')

# load("~/Study/STATS 489/PAL_check/rotavirus/data/real_rotavirus_metadata.Rdata")
# realdat <- data.frame(time = c(1:416), realdat)
# ggplot(data=realdat) + 
#   geom_line(aes(x=time, y=cases1),color='green') +
#   geom_line(aes(x=time, y=cases2), color='red') +
#   geom_line(aes(x=time, y=cases3), color='blue') +
#   ylab('Y_real')+xlab('date')