### This .R file can generate log-lik=-7000 by using Pal's data and Stocks' model
##### Simulate
#### WWR-Potential Model

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
                     rate[8] = delta1;          // I2 Aging
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

# read in the data
# add at t=0 a row of NAs to not have problems with the accumulator variable since
# t0 is much less than t[1]
# rinit <- Csnippet("
#     I1=2871;
#     I2=639;
#     I3=174;
#     S1=5.09484e+06;
#     S2=5.73856e+07;
#     S3=1.96976e+07;
#     R1=124410;
#     R2=57072;
#     R3=9578;
#     H1 = 0;
#     H2 = 0;
#     H3 = 0;
# ")

read.table("real_rotavirus_metadata.txt") %>%
  rbind(data.frame(time=0,cases1=NA,cases2=NA,cases3=NA)) %>%
  arrange(time) -> dat

params_fixed <- c(gamma=1, delta1=1/(5*52),delta2=1/(55*52), alpha=1/(78.86912*52), 
                  mu=1/(18.86912*52), N=82372825, omega=1/(1*52))
params <- params_fixed
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

params_stocks_stst_mle <- params

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
     paramnames = names(params_stocks_stst_mle),
     accumvars=c("H1", "H2", "H3"),
     partrans = pt,
     rinit=rinit_raw,
     params = params_stocks_stst_mle
) -> sir

## Pfilter
theta <- params_stocks_stst_mle

plan(multicore)
foreach(i=1:4, .combine=c,
        .options.future=list(seed=998468235L)
) %dopar% {
  pfilter(sir,Np=10000,params=theta)
} -> pfs

logmeanexp(logLik(pfs),se=TRUE)


### -7021.9882607


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

initial_parameters_tbl = dplyr::bind_cols(
  pomp::runif_design(lower=c(beta1=10,beta2=0.2,beta3=0.3,beta11=0.11,phi=0.01,od=0.001,sigma=0.001),
                     upper=c(beta1=15,beta2=0.4,beta3=0.5,beta11=0.16,phi=0.3,od=0.3,sigma=0.2),
                     nseq=100)
  )
lapply(1:nrow(initial_parameters_tbl), function(z)
  coef_to_pparams(initial_parameters_tbl[z,])
)

list(shared = shared_params, specific = specific_params)


bake(file = "mif_check_paldata_stocksmodel.rds",{ 
  mifs_global <- foreach(i=1:cores,.packages='pomp', .options.multicore=mcopts) %dopar% {
    mif2(
      sir,
      start=c(apply(sir_box,1,function(x)runif(1,min=x[1],max=x[2])),sir_fixed_params),
      Np=2,
      Nmif=2,
      cooling.type="geometric",
      cooling.fraction.50=0.5,
      transform=TRUE,
      rw.sd=rw_sd(beta1=0.002,beta2=0.002,beta3=0.002,beta11=0.001,phi=0.01,od=0.01, sigma=0.01)
    )
  }
  mifs_global
})
