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

# measurement model 
dmeas <- Csnippet("
                  if (ISNA(cases1)) {
                  lik = (give_log) ? 0 : 1;
                  } else {
                  lik =  dnbinom_mu(cases1, 1/od, H1, 1) +
                  dnbinom_mu(cases2, 1/od, H2, 1) +
                  dnbinom_mu(cases3, 1/od, H3, 1);
                  lik = (give_log) ? lik : exp(lik);}") # mean = r(1-p)/p
rmeas <-  Csnippet("
                   cases1 = rnbinom_mu(1/od,H1);
                   cases2 = rnbinom_mu(1/od,H2);
                   cases3 = rnbinom_mu(1/od,H3);")


# process model is Markovian SIRS with 3 age classes
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
                     rate[0] = alpha*N;
                     rate[1] = Beta1*I/N*dW/dt;
                     rate[2] = delta1;
                     rate[3] = Beta2*I/N*dW/dt;
                     rate[4] = delta2;
                     rate[5] = Beta3*I/N*dW/dt;
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
                     dN[0] = rpois(rate[0]*dt); // births are Poisson
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

rate[0] = alpha*N; // Birth
rate[1] = Beta1*(I1+I2+I3)/N; // S -> I
rate[2] = delta1;

rate[3] = Beta2*(I1+I2+I3)/N;
rate[4] = delta2;

rate[5] = Beta3*(I1+I2+I3)/N;
rate[6] = mu;   // Death = Delta

rate[7] = gamma;  // I -> R
rate[8] = delta1; // (1) -> (2)

rate[9] = gamma;
rate[10] = delta2;

rate[11] = gamma;
rate[12] = mu; // I's Death

rate[13] = delta1; // R_1 -> R_2
rate[14] = omega; // R_1 -> S_1

rate[15] = delta2;  // R_2 -> R_3
rate[16] = omega;  // R_2 -> S_2

rate[17] = mu;   // R'S Death
rate[18] = omega; // R_3 -> S_3 


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

# read in the data
# add at t=0 a row of NAs to not have problems with the accumulator variable since
# t0 is much less than t[1]
read.table("rota_study.txt") %>%
  rbind(data.frame(time=0,cases1=NA,cases2=NA,cases3=NA)) %>%
  arrange(time) -> dat


# define parameters (without betas)
params_fixed <- c(gamma=1, delta1=1/(5*52),delta2=1/(55*52), alpha=1/(78.86912*52), 
                  mu=1/(18.86912*52), N=82372825, omega=1/(1*52))
first_data <- c(y1=dat$cases1[2], y2=dat$cases2[2], y3=dat$cases3[2])


# initializer
init <- function(params, t0, ...) {
  x0 <- c(S1=0,I1=0,R1=0,H1=0,S2=0,I2=0,R2=0,H2=0,S3=0,I3=0,R3=0,H3=0)
  y <- params[c("y1","y2","y3")]
  x0["I1"] <- y[1]/((params["gamma"]+params["delta1"]))
  x0["I2"] <- (y[2]+params["delta1"]*x0["I1"])/((params["delta2"]+params["gamma"]))
  x0["I3"] <- (y[3]+params["delta2"]*x0["I2"])/((params["mu"]+params["gamma"]))
  x0["S1"] <- (params["alpha"]*params["N"]-(params["gamma"]+params["delta1"])*x0["I1"]+
                 params["omega"]*(params["N"]*params["alpha"]/params["delta1"]-x0["I1"]))/(params["delta1"]+params["omega"])
  I_bar    <-  x0["I1"]+x0["I2"]+x0["I3"]
  x0["S2"] <- (params["delta1"]*x0["S1"]-(params["delta2"]+params["gamma"])*x0["I2"]+params["delta1"]*x0["I1"]+
                 params["omega"]*(params["N"]*params["alpha"]/params["delta2"]-x0["I2"]))/(params["delta2"]+params["omega"])
  x0["S3"] <- (params["delta2"]*x0["S2"]-(params["gamma"]+params["mu"])*x0["I3"]+params["delta2"]*x0["I2"]+
                 params["omega"]*(params["N"]*params["alpha"]/params["mu"]-x0["I3"]))/(params["omega"]+params["mu"])
  x0["R1"] <- (params["N"]*params["alpha"]/params["delta1"]-x0["S1"]-x0["I1"])
  x0["R2"] <- (params["N"]*params["alpha"]/params["delta2"]-x0["S2"]-x0["I2"])
  x0["R3"] <- (params["N"]*params["alpha"]/params["mu"]-x0["S3"]-x0["I3"])
  round(x0) 
}

#help parameters with different data ie the mean data
mean_data <- c(y1=mean(dat$cases1[-1]), y2=mean(dat$cases2[-1]), y3=mean(dat$cases3[-1])) 
help_param <- c(params_fixed,mean_data)
# analytic guess for the betas
beta_ana <-  function(params){
  beta_ana <- c(beta1=0, beta2=0, beta3=0)
  I_bar <- init(params)["I1"]+init(params)["I2"]+init(params)["I3"]
  beta_ana["beta1"] <- ((params["gamma"]+params["delta1"])*init(params)["I1"]*params["N"])/(init(params)["S1"]*I_bar)
  beta_ana["beta2"] <- ((params["delta2"]+params["gamma"])*init(params)["I2"]-params["delta1"]*init(params)["I1"])*params["N"]/(init(params)["S2"]*I_bar)
  beta_ana["beta3"] <- ((params["mu"]+params["gamma"])*init(params)["I3"]-params["delta2"]*init(params)["I2"])*params["N"]/(init(params)["S3"]*I_bar)
  return(beta_ana)
}
dat <- realdat |> as.data.frame()
params_fixed <- c(gamma=1, delta1=1/(5*52),delta2=1/(55*52), alpha=1/(78.86912*52),
                  mu=1/(18.86912*52), N=82372825, omega=1/(1*52))
first_data <- c(y1=dat$cases1[2], y2=dat$cases2[2], y3=dat$cases3[2])
# paramtervector with betas and inital data 
params <- c(beta_ana(help_param), beta11=0.15, phi=0.1, params_fixed,first_data,od=0.3, sigma=0.05)

toEst<- Csnippet("
                 Tbeta1  = log(beta1);
                 Tbeta2  = log(beta2);
                 Tbeta3  = log(beta3);
                 Tbeta11 = logit(beta11);
                 Tsigma = log(sigma);
                 Tphi    = logit(phi/(M_2PI));
                 Tod = log(od);")

fromEst <-Csnippet("
                   Tbeta1  = exp(beta1);
                   Tbeta2  = exp(beta2);
                   Tbeta3  = exp(beta3);
                   Tsigma = exp(sigma);
                   Tbeta11 = expit(beta11);
                   Tphi    = M_2PI*expit(phi);
                   Tod = exp(od);")

pomp(data = dat,
     times="time",
     t0=1-6*52,
     dmeasure = dmeas,
     rmeasure = rmeas,
     rprocess = euler.sim(step.fun = sir.step, delta.t = 1/10),
     statenames = c("S1", "I1", "R1", "H1", "S2", "I2", "R2", "H2","S3","I3", "R3", "H3"),
     paramnames = names(params),
     zeronames=c("H1", "H2", "H3"),
     skeleton=vectorfield(Csnippet(sir.skel)),
     initializer=init,
     toEstimationScale=toEst,
     fromEstimationScale=fromEst,
     params = params
) -> sir

# # # ##how to create study data!! with rotavirus.txt instead
# a<- simulate(sir,seed=12345, obs=TRUE, as.data.frame=TRUE)
# rota_study<- data.frame(time=a$time[-1],cases1=a$cases1[-1],cases2=a$cases2[-1],cases3=a$cases3[-1])
# plot(rota_study$cases1, type="l")
# write.table(rota_study, file = "/afs/pdc.kth.se/home/t/tstoc/Private/paper/gamma_study/mle_gamma_study/rota_study.txt" )

################### GLOBAL SEARCH ############################
require(doParallel)
cores <- 20
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)



sir_box <- rbind(
  beta1=c(12,16),
  beta2=c(0.2,0.40),
  beta3=c(0.3,0.5),
  beta11=c(0.12,0.16),
  phi=c(0,0.2),
  od=c(0.1,0.5),
  sigma=c(0.001,0.2)
)

sir_fixed_params <- c(params_fixed, y1=dat$cases1[2], y2=dat$cases2[2], y3=dat$cases3[2])

stew(file="mle_gamma_study-%d.rda",{
  
  t_global <- system.time({
    mifs_global <- foreach(i=1:cores,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      mif2(
        sir,
        start=c(apply(sir_box,1,function(x)runif(1,min=x[1],max=x[2])),sir_fixed_params),
        Np=5000,
        Nmif=300,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(beta1=0.002,beta2=0.002,beta3=0.002,beta11=0.001,phi=0.01,od=0.01, sigma=0.01)
      )
    }
  })
},seed=1270401374,kind="L'Ecuyer")



mifs_global %>%
  conv.rec(c("loglik","nfail","beta1","beta2","beta3","beta11","phi","od","sigma")) %>%
  melt() %>%subset(iteration>175)%>%
  ggplot(aes(x=iteration,y=value,color=variable,group=L1))+
  geom_line()+
  guides(color=FALSE)+
  labs(x="MIF2 Iteration",y="")+
  facet_wrap(~variable,scales="free_y",ncol=2)+
  theme_bw()


stew(file="mle_gamma_study_lik-%d.rda",{
  
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")


liks_global
best <- which.max(liks_global[,1])
round(liks_global[best,],3)
round(coef(mifs_global[[best]]),3)


