rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
sourceCpp('equidispersed_model.cpp')
load('sim_list_100_eqeq.Rdata')
# load("real_rotavirus_metadata.Rdata")

init_dist <- c(3876549, 30351, 1315221, 57139612, 871, 302852, 19573727, 2550, 131092)

par <- c(11.48,0.25,0.35,0.14,0.16,0.021,66.89)

sim_list_100[[1]][3, 1] <- 0
sim_list_100[[1]][2, 1] <- 0
sim_list_100[[1]][1, 1] <- 0
realdat <- sim_list_100[[1]] |> as.matrix()
rotavirus_equidispersed(init_dist, realdat, 9, 
                        c(par[1],par[2],par[3],par[4],par[5]),
                        0.07)$log_lik

y <- c()
for(i in 1:100){
  realdat <- sim_list_100[[i]] |> as.matrix()
  y[i] <- rotavirus_equidispersed(init_dist, realdat, 9, 
                          c(par[1],par[2],par[3],par[4],par[5]),
                          0.07)$log_lik
}


sim_data <- data.frame(time = c(1:416), sim_list_100[[1]])

ggplot(data=sim_data) + 
  geom_line(aes(x=time, y=cases1),color='green') +
  geom_line(aes(x=time, y=cases2), color='red') +
  geom_line(aes(x=time, y=cases3), color='blue') +
  ylab('Y_sim_pomp')+xlab('date')

initial_guess <- c(runif(1,9,18),runif(1,0.05,0.5),runif(1,0.2,1),runif(1,0.05,0.4),runif(1,0.05,0.4))

astep <- 0.5*c(0.01,0.001,0.001,0.001,0.001,0.001)
cstep <- c(0.01,0.001,0.001,0.001,0.001,0.001)

Coordinate_ascent_algorithm2 <- function(a,c,init_params, n_steps){
  tstart <- Sys.time()
  par <- init_params
  traj <- matrix(nrow = 5, ncol= n_steps+1)
  lik <- c()
  lik[1] <- rotavirus_equidispersed(init_dist, realdat, 9, c(par[1],par[2],par[3],par[4],par[5]),0.07)$log_lik
  # lik[1] <- rotavirus_SMC_obs(init_dist,realdat , 9, c(par[1],par[2],par[3],par[4],par[5]),norm_par = c(0.07,par[6]), gamma_par = c(420,  1/420),t(prop), 1000,1)$log_lik
  for (i in 1:n_steps) {
    if(i%%100==0){print(paste('iteration =',i))}
    t1 <- Sys.time()
    
    for (j in 1:5) {
      # print(j)
      #Delta <- sample(c(-1,1),1)
      Delta <- 1
      par_pos <- par
      par_neg <- par
      par_pos[j] <- par[j]+cstep[j]
      par_neg[j] <- par[j]-cstep[j]
      
      grad_pos <- rotavirus_equidispersed(init_dist, realdat, 9,  c(par_pos[1],par_pos[2],par_pos[3],par_pos[4],par_pos[5]),0.07)$log_lik 
      grad_neg <- rotavirus_equidispersed(init_dist, realdat, 9,  c(par_neg[1],par_neg[2],par_neg[3],par_neg[4],par_neg[5]),0.07)$log_lik 
      # print(grad_pos)
      # print(grad_neg)
      grad_est <- sign(grad_pos - grad_neg)
      # print(grad_est)
      if(i>7500){par[j] <- par[j] + (1/(i-7500)^1.01)*astep[j]*grad_est + (1/(i-7500)^1.01)*runif(1,-1,1)*0.1*astep[j]}
      else{par[j] <- par[j] + astep[j]*grad_est + runif(1,-1,1)*0.1*astep[j]}
      if(par[j]==0){par[j] <- 0.01}
    }
    
    traj[,i+1] = par
    lik[i+1] <- rotavirus_equidispersed(init_dist, realdat, 9, c(par[1],par[2],par[3],par[4],par[5]),0.07)$log_lik
    t2 <- Sys.time()
  }
  tend <- Sys.time()
  time <- tstart - tend
  output <- list(traj = traj, lik =lik, time=time)
  return(output)
}

sign_ascent <- Coordinate_ascent_algorithm2(a,c,initial_guess, 10000)


filename <- paste('sign_ascent_equidispersed',job_id, '.Rdata', sep = '')

save(sign_ascent, file = filename)





