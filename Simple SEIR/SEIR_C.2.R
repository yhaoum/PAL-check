#This model can generate the plot in C.2

#simulate
library(pomp)

rproc <- Csnippet("
    double rate[3], trans[6];
    double alpha[4]={a_0,a_1,a_2,a_3};
  
    double xi = rgammawn(1,dt);
    // True death is not by expon rate, but binomial probability.
    rate[0] = -xi*beeta*I/n;    // stochastic force of infection
    rate[1] = -rho;         // rho is our original sigma
    rate[2] = -gamma;      // recovery
    
    // die first before transit
    trans[1] = rbinom(S, delta);
    trans[3] = rbinom(E, delta);
    trans[5] = rbinom(I, delta);
    
    trans[0] = rbinom(trans[1], 1-exp(rate[0]));
    trans[2] = rbinom(trans[3], 1-exp(rate[1]));
    trans[4] = rbinom(trans[5], 1-exp(rate[2]));
    
    // without birth case
    S = trans[1] - trans[0] + rpois(a_0);
    E = trans[3] + trans[0] - trans[2] + rpois(a_1);
    I = trans[5] + trans[2] - trans[4] + rpois(a_2);
    R = n - S - E - I + rpois(a_3); 
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
  int obs[4], arr_S[4], arr_E[4], arr_I[4], arr_R[4];
  int Y[4];
  double kappa[4]={k_0,k_1,k_2,k_3};
  
  //Detect then misreport;
  obs[0] = rbinom(S, q_S);
  obs[1] = rbinom(E, q_E);
  obs[2] = rbinom(I, q_I); // this is observed cases
  obs[3] = rbinom(R, q_R);
  
  //misreporting matrix:
  double prob_S[4] = {g00,g01,g02,g03};
  double prob_E[4] = {g10,g11,g12,g13};
  double prob_I[4] = {g20,g21,g22,g23};
  double prob_R[4] = {g30,g31,g32,g33};
  
  //Observed Y:
  rmultinom(obs[0], &prob_S[0], 4, &arr_S);
  rmultinom(obs[1], &prob_E[0], 4, &arr_E);
  rmultinom(obs[2], &prob_I[0], 4, &arr_I);
  rmultinom(obs[3], &prob_R[0], 4, &arr_R);
  
  
  for(int i=0; i<4; i++){
    Y[i] = arr_S[i]+arr_E[i]+arr_I[i]+arr_R[i];
  }
  
  Y_S = Y[0]+rpois(k_0);
  Y_E = Y[1]+rpois(k_1);
  Y_I = Y[2]+rpois(k_2);
  Y_R = Y[3]+rpois(k_3);
")

sim1 <- simulate(t0=0, 
                 times = 1:200, 
                 paramnames = c("beeta", "rho", "gamma", "delta",
                                "pi_0S", "pi_0I", "n",
                                "a_0","a_1","a_2","a_3",
                                "k_0","k_1","k_2","k_3",
                                "g00","g01","g02","g03",
                                "g10","g11","g12","g13",
                                "g20","g21","g22","g23",
                                "g30","g31","g32","g33",
                                "q_S","q_E","q_I","q_R"),
                 params = c(beeta=0.5, rho=0.05, gamma=0.1, delta=1, 
                            pi_0S=0.99, pi_0I=0.01, n=100000,
                            a_0=0,a_1=0,a_2=0,a_3=0,
                            k_0=0,k_1=0,k_2=0,k_3=0,
                            g00=0.95,g01=0.0,g02=0.05,g03=0.0,
                            g10=0.3,g11=0.0,g12=0.7,g13=0.0,
                            g20=0.15,g21=0.0,g22=0.85,g23=0.0,
                            g30=0.0,g31=0.0,g32=0.0,g33=1.0,
                            q_S=0.1,q_E=0.1,q_I=0.3,q_R=0.2),
                 rinit = rinit,
                 rprocess = discrete_time(rproc),
                 rmeasure = rmeas,
                 statenames=c("S","E","I","R"),
                 obsnames=c("Y_S","Y_E","Y_I","Y_R"))

sim_data <- as.data.frame(sim1)
