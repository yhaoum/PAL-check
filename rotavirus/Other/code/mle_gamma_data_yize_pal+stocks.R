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

RUN_LEVEL = 1
NP_MIF       = switch(RUN_LEVEL, 4, 5000)
NMIF         = switch(RUN_LEVEL, 4,  100)

# measurement model 
# measurement model 
dmeas <- Csnippet("
                  if (ISNA(cases1)) {
                  lik = (give_log) ? 0 : 1;
                  } else {
                      if(t < 209){
                        lik =  dnbinom_mu(cases1, 1/od, 0.06205271*H1, 1) +
                        dnbinom_mu(cases2, 1/od, 0.06762432*H2, 1) +
                        dnbinom_mu(cases3, 1/od, 0.06988662*H3, 1);
                      }
                      else{
                        lik =  dnbinom_mu(cases1, 1/od, 0.08791755*H1, 1) +
                        dnbinom_mu(cases2, 1/od, 0.09408622*H2, 1) +
                        dnbinom_mu(cases3, 1/od, 0.09649926*H3, 1);
                      }
                    
                    lik = (give_log) ? lik : exp(lik);
                        
                    }")
rmeas <-  Csnippet("
                  if(t < 209){
                   cases1 = rnbinom_mu(1/od,0.06205271*H1);
                   cases2 = rnbinom_mu(1/od,0.06762432*H2);
                   cases3 = rnbinom_mu(1/od,0.06988662*H3);
                  }
                  else{
                    cases1 = rnbinom_mu(1/od,0.08791755*H1);
                    cases2 = rnbinom_mu(1/od,0.09408622*H2);
                    cases3 = rnbinom_mu(1/od,0.09649926*H3);
                  }
                  ")




# process model is Markovian SIRS with 3 age classes
sir.step <- Csnippet("
                     double rate[19];
                     double dN[19];
                     double Beta1;
                     double Beta2;
                     double Beta3;
                     double I;
                     double dW;
                     
                     // compute the environmental stochasticity
                     dW = rgammawn(sigma,dt);
                     I= I1+I2+I3;
                     Beta1 = beta1*(1 + beta11 * cos(M_2PI/52*t + phi));
                     Beta2 = beta2*(1 + beta11 * cos(M_2PI/52*t + phi));
                     Beta3 = beta3*(1 + beta11 * cos(M_2PI/52*t + phi));
                     rate[0] = alpha*N;         // Birth
                     rate[1] = Beta1*I/N*dW/dt; // S1 Trans: dW is gamma white noise, i.e. st in X.
                     rate[2] = delta1;          // S1 Aging
                     rate[3] = Beta2*I/N*dW/dt; // S2 Trans
                     rate[4] = delta2;          // S2 Aging
                     rate[5] = Beta3*I/N*dW/dt; // S3 Trans
                     rate[6] = mu;              // S3 Death
                     rate[7] = gamma;           // I1 Trans
                     rate[8] = delta1;          // I1 Aging
                     rate[9] = gamma;           // I2 Trnas
                     rate[10] = delta2;         // I2 Aging
                     rate[11] = gamma;          // I3 Trans
                     rate[12] = mu;             // I3 Death
                     rate[13] = delta1;         // R1 Aging
                     rate[14] = omega;          // R1 Trans
                     rate[15] = delta2;         // R2 Aging
                     rate[16] = omega;          // R2 Trans
                     rate[17] = mu;             // R3 Death
                     rate[18] = omega;          // R3 Trans
                     dN[0] = rpois(rate[0]*dt); // alpha*N*dt
                     reulermultinom(2, S1, &rate[1], dt, &dN[1]);
                     reulermultinom(2, S2, &rate[3], dt, &dN[3]);
                     reulermultinom(2, S3, &rate[5], dt, &dN[5]);
                     reulermultinom(2, I1, &rate[7], dt, &dN[7]);
                     reulermultinom(2, I2, &rate[9], dt, &dN[9]);
                     reulermultinom(2, I3, &rate[11], dt, &dN[11]);
                     reulermultinom(2, R1, &rate[13], dt, &dN[13]);
                     reulermultinom(2, R2, &rate[15], dt, &dN[15]);
                     reulermultinom(2, R3, &rate[17], dt, &dN[17]);
                     S1 += dN[0] - dN[1] - dN[2] + dN[14];
                     S2 += dN[2] - dN[3] - dN[4]  + dN[16];
                     S3 += dN[4] - dN[5] - dN[6] + dN[18];
                     I1 += dN[1]          - dN[7] - dN[8];
                     I2 += dN[3] + dN[8]  - dN[9] - dN[10];
                     I3 += dN[5] + dN[10] - dN[11] - dN[12];
                     R1 += dN[7]           - dN[13] - dN[14];
                     R2 += dN[9]  + dN[13] - dN[15] - dN[16];
                     R3 += dN[11] + dN[15] - dN[17] - dN[18];
                     H1 += dN[1];
                     H2 += dN[3];
                     H3 += dN[5];
                     ")



# ------------ deterministic skeleton-----------------------------
sir.skel <- "
double rate[19];
double term[19];
double Beta1;
double Beta2;
double Beta3;

Beta1 = beta1*(1 + beta11 * cos(M_2PI/52*t + phi)); //seasonal forcing
Beta2 = beta2*(1 + beta11 * cos(M_2PI/52*t + phi)); //seasonal forcing
Beta3 = beta3*(1 + beta11 * cos(M_2PI/52*t + phi)); //seasonal forcing

rate[0] = alpha*N;
rate[1] = Beta1*(I1+I2+I3)/N;
rate[2] = delta1;

rate[3] = Beta2*(I1+I2+I3)/N;
rate[4] = delta2;

rate[5] = Beta3*(I1+I2+I3)/N;
rate[6] = mu;

rate[7] = gamma;
rate[8] = delta1;

rate[9] = gamma;
rate[10] = delta2;

rate[11] = gamma;
rate[12] = mu;

rate[13] = delta1;
rate[14] = omega;  

rate[15] = delta2;  
rate[16] = omega;  

rate[17] = mu;  
rate[18] = omega;  


// compute the several terms
term[0] = rate[0];

term[1] = rate[1] * S1;
term[2] = rate[2] * S1;

term[3] = rate[3] * S2;
term[4] = rate[4] * S2;

term[5] = rate[5] * S3;
term[6] = rate[6] * S3;

term[7] = rate[7] * I1;
term[8] = rate[8] * I1;

term[9] = rate[9] * I2;
term[10] = rate[10] * I2;

term[11] = rate[11] * I3;
term[12] = rate[12] * I3;

term[13] = rate[13] * R1;
term[14] = rate[14] * R1;

term[15] = rate[15] * R2;
term[16] = rate[16] * R2;

term[17] = rate[17] * R3;
term[18] = rate[18] * R3;


DS1 = term[0] - term[1] - term[2] + term[14];
DI1 = term[1]          - term[7] - term[8];
DR1 = term[7]          - term[13] - term[14];
DH1 = term[1];

DS2 = term[2] - term[3] - term[4]  + term[16];
DI2 = term[3] + term[8]  - term[9] - term[10];
DR2 = term[9]  + term[13] - term[15] - term[16];
DH2 = term[3];

DS3 = term[4] - term[5] - term[6] + term[18];
DI3 = term[5] + term[10] - term[11] - term[12];
DR3 = term[11] + term[15] - term[17] - term[18];
DH3 = term[5];

" 

## ----rinit-------------------------------------------------
rinit_raw <- Csnippet("
  I1 = nearbyint(y1/(gamma+delta1));
  I2 = nearbyint((y2+delta1*I1)/((delta2+gamma)));
  I3 = nearbyint((y3+delta2*I2)/((mu+gamma)));
  S1 = nearbyint((alpha*N-(gamma+delta1)*I1+omega*(N*alpha/delta1-I1))/(delta1+omega));
  S2 = nearbyint((delta1*S1-(delta2+gamma)*I2+delta1*I1+omega*(N*alpha/delta2-I2))/(delta2+omega));
  S3 = nearbyint((delta2*S2-(gamma+mu)*I3+delta2*I2+omega*(N*alpha/mu-I3))/(omega+mu));
  R1 = nearbyint((N*alpha/delta1-S1-I1));
  R2 = nearbyint((N*alpha/delta2-S2-I2));
  R3 = nearbyint((N*alpha/mu-S3-I3));
")

rinit <- Csnippet("
    I1=2871;
    I2=639;
    I3=174;
    S1=5.09484e+06;
    S2=5.73856e+07;
    S3=1.96976e+07;
    R1=124410;
    R2=57072;
    R3=9578;
    H1 = 0;
    H2 = 0;
    H3 = 0;
")

read.table("real_rotavirus_metadata.txt") %>%
  rbind(data.frame(time=0,cases1=NA,cases2=NA,cases3=NA)) %>%
  arrange(time) -> dat

# define parameters (without betas)
params_fixed <- c(gamma=1, delta1=1/(5*52),delta2=1/(55*52), alpha=1/(78.86912*52), 
                  mu=1/(18.86912*52), N=82372825, omega=1/(1*52))
params <- c(params_fixed)
params["beta1"] <- 11.298
params["beta2"] <- 0.267
params["beta3"] <- 0.433
params["beta11"] <- 0.148
params["phi"] <- 0.085
params["od"] <- 0.111
params["sigma"] <- 0.091
params["y1"] <- 2882
params["y2"] <- 628
params["y3"] <- 174

sir_fixed_params <- c(params_fixed, params["y1"], params["y2"], params["y3"])

pt <- pomp::parameter_trans(
  log = c("beta1","beta2","beta3","sigma","od"),
  logit=c("beta11"),
  toEst= pomp::Csnippet("T_phi = logit(phi/(M_2PI));"),
  fromEst= pomp::Csnippet("phi = M_2PI*expit(T_phi);")
)
                            

pomp(data = dat,
     times="time",
     t0=1-6*52,
     dmeasure = dmeas,
     rmeasure = rmeas,
     rprocess = euler(step.fun = sir.step, delta.t = 1/10),
     statenames = c("S1", "I1", "R1", "H1", "S2", "I2", "R2", "H2","S3","I3", "R3", "H3"),
     paramnames = names(params),
     accumvars=c("H1", "H2", "H3"),
     skeleton=vectorfield(Csnippet(sir.skel)),
     rinit=rinit_raw,
     partrans = pt,
     params = params
) -> sir

sir_panel <- panelPomp::panelPomp(list(unit1=sir),
                     shared=NULL,
                     specific=params |> 
                       as.matrix() |>
                       `colnames<-`("unit1")
                     )

require(doParallel)
cores <- 36
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)


sir_box <- rbind(
  beta1=c(10,15),
  beta2=c(0.2,0.40),
  beta3=c(0.3,0.5),
  beta11=c(0.11,0.16),
  phi=c(0.01,0.3),
  od=c(0.001,0.3),
  sigma=c(0.001,0.2)
)

sobol_design(
  lower=c(beta1=10,beta2=0.2,beta3=0.3,beta11=0.11,phi=0.01,od=0.001,sigma=0.001),
  upper=c(beta1=15,beta2=0.4,beta3=0.5,beta11=0.16,phi=0.3,od=0.3,sigma=0.2),
  nseq=100
) -> guesses

## Make it panelPomp
c(apply(sir_box,1,function(x)runif(1,min=x[1],max=x[2])),sir_fixed_params) |> 
  as.matrix() |>
  `colnames<-`("unit1") -> starting

bake(file = "mif_Yize_panel_01.rds",{ 
     mifs_global <- foreach(i=1:cores,.packages='pomp', .options.multicore=mcopts) %dopar% {
       panelPomp::mif2(
         sir_panel,
         Np = 2,
         cooling.fraction.50 = 0.5,
         rw.sd = rw_sd(beta1=0.002,beta2=0.002,beta3=0.002,
                       beta11=0.001,phi=0.01,od=0.01, sigma=0.01),
         cooling.type = "geometric",
         Nmif = 2,
         shared.start = numeric(0),
         specific.start = starting,
         block = F
       )
     }
    mifs_global
})





