# Section 7.2 & Supplement C.2

rproc <- Csnippet("
    double rate[3], trans[7];
    double N = S+E+I+R;         // population size
    double xi = rgammawn(1,dt); // Is this mean 1 gamma noise??
    // True death is not by expon rate, but binomial probability.
    rate[0] = -xi*beeta*I/N;    // stochastic force of infection
    rate[1] = -rho;         // rho is our original sigma
    rate[2] = -gamma;      // recovery
    
    // die first before transit
    trans[0] = rbinom(S, delta);
    trans[2] = rbinom(E, delta);
    trans[4] = rbinom(I, delta);
    trans[6] = rbinom(R, delta);
    
    trans[1] = rbinom(trans[0], 1-exp(rate[0]));
    trans[3] = rbinom(trans[2], 1-exp(rate[1]));
    trans[5] = rbinom(trans[4], 1-exp(rate[2]));
    
    // with birth case
    S = trans[0] - trans[1] + rpois(0.05*pi_0S*N);
    E = trans[2] + trans[1] - trans[3];
    I = trans[4] + trans[3] - trans[5] + rpois(0.05*pi_0I*N);
    R = trans[6] + trans[5]; 
    
    // Truncated Normal reporting rate
    q_I = -1;
  
    while(q_I<0 || q_I>1){
    q_I = rnorm(mu_q, sigma_q);
    }
")

## ----rinit-------------------------------------------------
rinit <- Csnippet("
  double pi_0[4]={pi_0S,0,pi_0I,0}; // pi_0 should be specified
  int trans_0[4];

  rmultinom(n, &pi_0, 4, &trans_0);
  S = trans_0[0];
  E = trans_0[1];
  I = trans_0[2];
  R = trans_0[3];
")


## ----rmeasure-------------------------------------------------
rmeas <- Csnippet("
  Y_I = rbinom(I, q_I); // this is observed cases
")

dmeas <- Csnippet("
  lik = dbinom(Y_I, I, q_I, 1);
")

sim1 <- simulate(t0=0, 
                 times = 1:100, 
                 paramnames = c("beeta", "rho", "gamma", "delta",
                                "pi_0S", "pi_0I", "n",
                                "mu_q","sigma_q"),
                 params = c(beeta=0.8, rho=0.1, gamma=0.2, delta=0.95, 
                            pi_0S=0.99, pi_0I=0.01, 
                            n=10000,
                            mu_q=0.5, sigma_q=0.1),
                 rinit = rinit,
                 rprocess = discrete_time(rproc),
                 rmeasure = rmeas,
                 statenames=c("S","E","I","R","q_I"),
                 obsnames=c("Y_I"))

sim_data <- as.data.frame(sim1)

dat <- subset(sim_data, select = c("time", "Y_I"))

dat |>
  pomp(t0=0,
       time="time",
       rprocess=discrete_time(rproc),
       rinit=rinit,
       dmeasure=dmeas,
       rmeasure=rmeas,
       statenames=c("S","E","I","R","q_I"),
       paramnames=c("beeta", "rho", "gamma", "delta",
                    "pi_0S", "pi_0I", "n",
                    "mu_q","sigma_q")
  ) -> m1

theta <- c(beeta=0.8, rho=0.1, gamma=0.2, delta=0.95, 
           pi_0S=0.99, pi_0I=0.01, 
           n=10000,
           mu_q=0.5, sigma_q=0.1)

library(doFuture)
plan(multicore)
foreach(i=1:4, .combine=c,
        .options.future=list(seed=998468235L)
) %dopar% {
  pfilter(m1,Np=10000,params=theta)
} -> pfs
## ----pfilter1b-----------------------------------------------
logmeanexp(logLik(pfs),se=TRUE)