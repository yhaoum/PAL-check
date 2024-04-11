##### Pfilter + EqEq
rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows
library(Rcpp)
library(RcppArmadillo)
# install.packages('RcppDist')
sourceCpp('rotavirus_normq.cpp')
load("sim_keep_nonzero_125.Rdata")
prop <- 420

init_dist <- c(3876549, 30351, 1315221, 57139612, 871, 302852, 19573727, 2550, 131092)
# initial_guess <- c(runif(1,5,20),runif(1,0.05,0.5),runif(1,0.2,1),runif(1,0.05,0.4),runif(1,0.05,0.4),runif(1,0.05,0.2),runif(1,10,100))
initial_guess <- c(11.48,0.25,0.35,0.14,0.16,0.021,66.89)

### Just a filter
Sys.time()
y <- c()
for(i in 1:100){
  realdat <- sim_keep_nonzero_125[[i]] |> as.matrix()
  lik_list <- rotavirus_SMC_qropxi(init_dist,realdat , 9, 
                                   c(initial_guess[1],initial_guess[2],initial_guess[3],initial_guess[4],initial_guess[5]), 
                                   gamma_par = c(initial_guess[7],  1/initial_guess[7]), 
                                   norm_par = c(0.07,initial_guess[6]), 
                                   t(prop), 20000,1)
  y[i] <- lik_list$log_lik
}
Sys.time()
save(y, file="palfilter_ovov_100.Rdata")
lik_list$ll_storage |> sum()


# maximization
astep <- c(0.01,0.001,0.001,0.001,0.001,0.001,0.25)
cstep <- 5*c(0.01,0.001,0.001,0.001,0.001,0.001,0.5)

Coordinate_ascent_algorithm2 <- function(a,c,init_params, n_steps){
  tstart <- Sys.time()
  par <- init_params
  traj <- matrix(nrow = 7, ncol= n_steps+1)
  lik <- c()
  lik[1] <- rotavirus_SMC_qropxi(init_dist,realdat , 9, c(par[1],par[2],par[3],par[4],par[5]),norm_par = c(0.07,par[6]), gamma_par = c(par[7],  1/par[7]),t(prop), 10000,1)$log_lik
  for (i in 1:n_steps) {
    t1 <- Sys.time()
    print(paste('iteration =',i))
    for (j in 1:7) {
      print(j)
      #Delta <- sample(c(-1,1),1)
      Delta <- 1
      par_pos <- par
      par_neg <- par
      par_pos[j] <- par[j]+cstep[j]
      par_neg[j] <- par[j]-cstep[j]
      
      grad_pos <- rotavirus_SMC_qropxi(init_dist,realdat , 9, c(par_pos[1],par_pos[2],par_pos[3],par_pos[4],par_pos[5]),norm_par = c(0.07,par_pos[6]), gamma_par = c(par_pos[7],  1/par_pos[7]),t(prop), 1000,1)$log_lik
      grad_neg <- rotavirus_SMC_qropxi(init_dist,realdat , 9, c(par_neg[1],par_neg[2],par_neg[3],par_neg[4],par_neg[5]),norm_par = c(0.07,par_neg[6]), gamma_par = c(par_neg[7],  1/par_neg[7]),t(prop), 1000,1)$log_lik
      print(grad_pos)
      print(grad_neg)
      grad_est <- sign(grad_pos - grad_neg)
      print(grad_est)
      par[j] <- par[j] + astep[j]*grad_est + runif(1,-1,1)*0.1*astep[j]
      if(par[j]==0){par[j] <- 0.01}
    }
    
    traj[,i+1] = par
    lik[i+1] <- rotavirus_SMC_qropxi(init_dist,realdat , 9, c(par[1],par[2],par[3],par[4],par[5]),norm_par = c(0.07,par[6]), gamma_par = c(par[7],  1/par[7]),t(prop), 1000,1)$log_lik
    t2 <- Sys.time()
    print(par)
    print(lik[i+1])
    print(t2-t1)
  }
  tend <- Sys.time()
  time <- tstart - tend
  output <- list(traj = traj, lik =lik, time=time)
  return(output)
}

coordinate_ascent <- Coordinate_ascent_algorithm2(a,c,initial_guess , 6000)

save(coordinate_ascent, file = "output/WWR_filter.Rdata")






