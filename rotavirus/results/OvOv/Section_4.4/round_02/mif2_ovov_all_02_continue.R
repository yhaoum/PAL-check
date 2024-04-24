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
registerDoParallel()
stopifnot(packageVersion("pomp") >= "1.7")

RUN_LEVEL = 2
NP_MIF       = switch(RUN_LEVEL, 4, 50000)
NMIF         = switch(RUN_LEVEL, 4,  100)
ncores       = switch(RUN_LEVEL, 4,  36)
sd           = 2/3
top_n_fits   = 12
# measurement model 
dmeas <- Csnippet("
                  if (ISNA(cases1)) {
                  lik = (give_log) ? 0 : 1;
                  } else {
                        lik =  dbinom(cases1, H1, q1, 1) +
                        dbinom(cases2, H2, q2, 1) +
                        dbinom(cases3, H3, q3, 1);
                      
                    lik = (give_log) ? lik : exp(lik);
                        
                    }")
rmeas <-  Csnippet("
                    cases1 = rbinom(H1, q1);
                    cases2 = rbinom(H2, q2);
                    cases3 = rbinom(H3, q3);
                  ")



rproc <- Csnippet("
    int I = I1 + I2 + I3;
    int trans_S1[3], trans_S2[3], trans_S3[2], trans_I1[3], trans_I2[3], trans_I3[2], trans_R1[3], trans_R2[3], trans_R3[2];
    
    double prob_S1[3],prob_I1[3],prob_R1[3],prob_S2[3],prob_I2[3],prob_R2[3],prob_S3[2],prob_I3[2],prob_R3[2];
    
    double xi = rgamma(sigma_xi, 1/sigma_xi);
    
    double kappa = (1 + beta11*cos(2*3.141593*t/52 + phi)) * xi;
    
    // Define rate
    prob_S1[0] = 1-exp(-dt*beta1*kappa*I/N); // 0->1
    prob_S1[1] = 1-exp(-delta1*dt);
    prob_S1[2] = exp(-delta1*dt) + exp(-dt*beta1*kappa*I/N) - 1;
    
    prob_I1[0] = 1-exp(-gamma*dt);
    prob_I1[1] = 1-exp(-delta1*dt);
    prob_I1[2] = exp(-gamma*dt)+exp(-delta1*dt) - 1;
    
    prob_R1[0] = 1 - exp(-omega*dt);  // E_1,t this goes back to S_1,(t+1)
    prob_R1[1] = 1 - exp(-delta1*dt);
    prob_R1[2] = exp(-omega*dt) + exp(-delta1*dt) - 1;
    
    prob_S2[0] = 1-exp(-dt*beta2*kappa*I/N);
    prob_S2[1] = 1-exp(-delta2*dt);
    prob_S2[2] = exp(-delta2*dt) + exp(-dt*beta2*kappa*I/N) - 1;
    
    prob_I2[0] = 1-exp(-dt*gamma);
    prob_I2[1] = 1-exp(-dt*delta2);
    prob_I2[2] = exp(-dt*gamma)+exp(-dt*delta2) - 1;
    
    prob_R2[0] = 1 - exp(-dt*omega);  // E_1,t this goes back to S_1,(t+1)
    prob_R2[1] = 1 - exp(-dt*delta2);
    prob_R2[2] = exp(-dt*omega) + exp(-dt*delta2) - 1;
    
    // For Age Group (3): Die first before transition;
    
    int S3mD, I3mD, R3mD;
    
    S3mD = rbinom(S3, 1-dt*mu); // S3 minus Death: mu is the death rate, so it's 1-mu here
    I3mD = rbinom(I3, 1-dt*mu);
    R3mD = rbinom(R3, 1-dt*mu);
    
    prob_S3[0] = 1-exp(-dt*beta3*kappa*I/N);
    prob_S3[1] = exp(-dt*beta3*kappa*I/N);
    
    prob_I3[0] = 1 - exp(-dt*gamma);
    prob_I3[1] = exp(-dt*gamma);
    
    prob_R3[0] = 1 - exp(-dt*omega);
    prob_R3[1] = exp(-dt*omega);
    
    // Transition
    // B: S->I
    // C: I->R
    // F: Aging: (1)->(2)->(3)
    // E: R->S
    // D: Death
    //// Note: Here S_1, S_2... are all old value from (t-1)
    rmultinom(S1, &prob_S1, 3, &trans_S1); // B, F, S-B-F
    rmultinom(I1, &prob_I1, 3, &trans_I1); // C, F, I-C-F
    rmultinom(R1, &prob_R1, 3, &trans_R1); // E, F, R-E-F
    
    rmultinom(S2, &prob_S2, 3, &trans_S2); // B, F, S-B-F
    rmultinom(I2, &prob_I2, 3, &trans_I2); // C, F, I-C-F
    rmultinom(R2, &prob_R2, 3, &trans_R2); // E, F, R-E-F
    
    rmultinom(S3mD, &prob_S3, 2, &trans_S3); // B, (S-D)-B
    rmultinom(I3mD, &prob_I3, 2, &trans_I3); // C, (I-D)-C
    rmultinom(R3mD, &prob_R3, 2, &trans_R3); // E, (R-D)-E
    
    S1 = trans_S1[2] + trans_R1[0] + rpois(4*1025.7); // Include Birth
    I1 = trans_I1[2] + trans_S1[0];
    R1 = trans_R1[2] + trans_I1[0];
    
    S2 = trans_S2[2] + trans_R2[0] + trans_S1[1]; // Include Aging
    I2 = trans_I2[2] + trans_S2[0] + trans_I1[1];
    R2 = trans_R2[2] + trans_I2[0] + trans_R1[1];
    
    S3 = trans_S3[1] + trans_R3[0] + trans_S2[1]; // Include Aging
    I3 = trans_I3[1] + trans_S3[0] + trans_I2[1];
    R3 = trans_R3[1] + trans_I3[0] + trans_R2[1];
    
    //Accumvar
    H1 += trans_S1[0];
    H2 += trans_S2[0];
    H3 += trans_S3[0];
    
    q1 = -1; 
    while(q1 < 0 || q1 > 1){
      q1 = rnorm(0.07, sigma_q);
    }
    
    q2 = -1; 
    while(q2 < 0 || q2 > 1){
      q2 = rnorm(0.07, sigma_q);
    }
    
    q3 = -1; 
    while(q3 < 0 || q3 > 1){
      q3 = rnorm(0.07, sigma_q);
    }
")


# define parameters (without betas)
params_fixed <- c(gamma=1, delta1=1/(5*52),delta2=1/(55*52), alpha=1/(78.86912*52), 
                  mu=0, N=82372825, omega=1/(1*52))
# WWR's rinit
rinit <- Csnippet("
    double m = N/(S10+I10+R10+S20+I20+R20+S30+I30+R30);
    I1=nearbyint(m*I10);
    I2=nearbyint(m*I20);
    I3=nearbyint(m*I30);
    S1=nearbyint(m*S10);
    S2=nearbyint(m*S20);
    S3=nearbyint(m*S30);
    R1=nearbyint(m*R10);
    R2=nearbyint(m*R20);
    R3=nearbyint(m*R30);
    H1 = 0;
    H2 = 0;
    H3 = 0;
")

# Set to MLE
params_stocks_stst_mle <- params_fixed

params_stocks_stst_mle["beta1"] <- 11.48
params_stocks_stst_mle["beta2"] <- 0.25
params_stocks_stst_mle["beta3"] <- 0.35
params_stocks_stst_mle["phi"] <- 0.14
params_stocks_stst_mle["beta11"] <- 0.16
params_stocks_stst_mle["sigma_q"] <- 0.021
params_stocks_stst_mle["sigma_xi"] <- 66.89
params_stocks_stst_mle["S10"] <- 0.047061
params_stocks_stst_mle["I10"] <- 0.000368
params_stocks_stst_mle["R10"] <- 0.015967
params_stocks_stst_mle["S20"] <- 0.015967
params_stocks_stst_mle["I20"] <- 0.000011
params_stocks_stst_mle["R20"] <- 0.003677
params_stocks_stst_mle["S30"] <- 0.237624
params_stocks_stst_mle["I30"] <- 0.000031
params_stocks_stst_mle["R30"] <- 0.001591



pt <- pomp::parameter_trans(
  log = c("beta1","beta2","beta3","sigma_q","sigma_xi"),
  logit=c("beta11"),
  barycentric=c("S10","I10","R10",
                "S20","I20","R20",
                "S30","I30","R30"),
  toEst= pomp::Csnippet("T_phi = logit(phi/(M_2PI));"),
  fromEst= pomp::Csnippet("phi = M_2PI*expit(T_phi);")
)

read.table("real_rotavirus_metadata.txt") %>%
  rbind(data.frame(time=0,cases1=NA,cases2=NA,cases3=NA)) %>%
  arrange(time) -> dat


pomp(data = dat,
     times="time",
     t0=0,
     dmeasure = dmeas,
     rmeasure = rmeas,
     rprocess = discrete_time(step.fun = rproc, delta.t = 1/4),
     statenames = c("S1", "I1", "R1", "H1", 
                    "S2", "I2", "R2", "H2",
                    "S3", "I3", "R3", "H3", 
                    "q1", "q2", "q3"),
     paramnames = names(params_stocks_stst_mle),
     accumvars=c("H1", "H2", "H3"),
     rinit=rinit,
     partrans = pt,
     params = params_stocks_stst_mle
) -> sir

sir_panel <- panelPomp::panelPomp(list(unit1=sir),
                                  shared=NULL,
                                  specific=params_stocks_stst_mle |> 
                                    as.matrix() |>
                                    `colnames<-`("unit1")
)

require(doParallel)
cores <- ncores
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)

### Next-round code
el <- readRDS("output_02/round_01/ovov_01_el.rds")
x <- na.omit(el$fits)
score_total = x$logLik
ranking_total = order(score_total, decreasing = TRUE)[1:top_n_fits]

best_fits = dplyr::select(
  x[ranking_total,], -"logLik", -"se"
)

recycle_vec = sort(rep_len(1:top_n_fits, cores))
full_best_fit <- best_fits[recycle_vec, ] 

coef_names <- colnames(full_best_fit)
colnames(full_best_fit) <- gsub(".{7}$","",coef_names)

starting_values <- vector(cores, mode="list")

for(i in 1:cores){
  t(full_best_fit[i, ])  |> 
    as.matrix() |>
    `colnames<-`("unit1") -> starting_values[[i]] 
}

mifs_global <- foreach(i=1:cores,.packages='pomp', .options.multicore=mcopts) %dopar% {
  panelPomp::mif2(
    sir_panel,
    Np = NP_MIF,
    cooling.fraction.50 = 0.5,
    rw.sd = rw_sd(beta1=0.01*sd,beta2=0.01*sd,beta3=0.01*sd,
                  beta11=0.01*sd,phi=0.01*sd,sigma_q=0.01*sd,sigma_xi=0.01*sd,
                  S10=ivp(0.24)*sd,I10=ivp(0.24)*sd,R10=ivp(0.24)*sd,
                  S20=ivp(0.24)*sd,I20=ivp(0.24)*sd,R20=ivp(0.24)*sd,
                  S30=ivp(0.24)*sd,I30=ivp(0.24)*sd,R30=ivp(0.24)*sd),
    cooling.type = "geometric",
    Nmif = NMIF,
    shared.start = numeric(0),
    specific.start = starting_values[[i]],
    block = F
  ) 
}

bake(file = "output_02/round_02/ovov_mif_02.rds",{ 
  mifs_global
})

bake(file = "output_02/round_02/ovov_02_el.rds",{
  el <- measlespkg::eval_logLik(mifs_global, ncores = cores, np_pf = NP_MIF, nreps=cores)
  el
})
