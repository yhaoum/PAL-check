##### Pfilter + OvOv
rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows
# close graphics windows
library(pomp)
library(magrittr)
library(plyr)
library(reshape2)
library(scales)
library(foreach)
library(doParallel)
registerDoParallel()
stopifnot(packageVersion("pomp") >= "1.7")

run_level = 2
np = switch(run_level, 10, 50000)
ncores = switch(run_level, 4, 36)
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
    
    S3mD = rbinom(S3, 1-dt*mu); 
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
params_fixed <- c(gamma=1, delta1=0.003636364,delta2=1/(55*52), alpha=1/(78.86912*52), mu=0, N=82372825, omega=1/(1*52))
rinit <- Csnippet("
    S1=3876549;
    I1=30351;
    R1=1315221;
    S2=57139612;
    I2=871;
    R2=302852;
    S3=19573727;
    I3=2550;
    R3=131092;
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


load("sim_keep_nonzero_125.Rdata")
x <- c()

for(i in 1:50){
  set.seed(438+274*i)
  
  dat <- data.frame(time=c(1:416), sim_keep_nonzero_125[[i]])
  dat |>
    pomp(t0=0,
         time="time",
         rprocess=discrete_time(rproc,delta.t=1/4), ## Obs time every week, but delta.t transition is actually everyday
         rinit=rinit,
         dmeasure=dmeas,
         rmeasure=rmeas,
         accumvars=c("H1","H2","H3"),
         statenames=c("S1","S2","S3",
                      "I1","I2","I3",
                      "R1","R2","R3",
                      "H1","H2","H3",
                      "q1","q2","q3"),
         paramnames = names(params_stocks_stst_mle)
    ) -> m1
  
  theta <- params_stocks_stst_mle
  
  require(doParallel)
  registerDoParallel(ncores)
  # mcopts <- list(set.seed=TRUE)
  doRNG::registerDoRNG()

  foreach(i=1:ncores, .combine=c,
          .options.future=list(seed=998468235L)
  ) %dopar% {
    pfilter(m1,Np=np,params=theta)
  } -> pfs
  
  x[i] <-logmeanexp(logLik(pfs))  
}

save(x, file = "output/pfilter_nonzero_50.Rdata")
