---
  title: "Test Here"
author: "Yize Hao"
date: "2023-09-23"
output: html_document
---
  
  ```{r}
library(pomp)
library(ggplot2)
```

```{r}
rproc <- Csnippet("
    double rate[3], trans[7];
    double alpha[4]={a_0,a_1,a_2,a_3};
  
    double xi = rgammawn(1,dt);
    // True death is not by expon rate, but binomial probability.
    rate[0] = -beeta*I/(S+E+I+R);    // stochastic force of infection
    rate[1] = -rho;         // rho is our original sigma
    rate[2] = -gamma;      // recovery
    
    // die first before transit
    trans[1] = rbinom(S, delta);
    trans[3] = rbinom(E, delta);
    trans[5] = rbinom(I, delta);
    trans[6] = rbinom(R, delta);
    
    trans[0] = rbinom(trans[1], 1-exp(rate[0]));
    trans[2] = rbinom(trans[3], 1-exp(rate[1]));
    trans[4] = rbinom(trans[5], 1-exp(rate[2]));
    
    // with birth case
    S = trans[1] - trans[0] + rpois(a_0);
    E = trans[3] + trans[0] - trans[2] + rpois(a_1);
    I = trans[5] + trans[2] - trans[4] + rpois(a_2);
    R = trans[6] + trans[4] + rpois(a_3); 
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
  rmultinom(obs[0], &prob_S, 4, &arr_S);
  rmultinom(obs[1], &prob_E, 4, &arr_E);
  rmultinom(obs[2], &prob_I, 4, &arr_I);
  rmultinom(obs[3], &prob_R, 4, &arr_R);
  
  
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
                 params = c(beeta=0.5, rho=0.05, gamma=0.1, delta=0.98, 
                            pi_0S=0.99, pi_0I=0.01, 
                            n=100000,
                            a_0=4000,   a_1=0,    a_2=0,      a_3=0,
                            k_0=1000,   k_1=1000, k_2=1000,   k_3=1000,
                            g00=0.95,   g01=0.0,  g02=0.05,   g03=0.0,
                            g10=0.3,    g11=0,    g12=0.7,    g13=0.0,
                            g20=0.15,   g21=0.0,  g22=0.85,   g23=0.0,
                            g30=0.0,    g31=0.0,  g32=0.0,    g33=1.0,
                            q_S=0.1,    q_E=0.1,  q_I=0.3,    q_R=0.2),
                 rinit = rinit,
                 rprocess = discrete_time(rproc),
                 rmeasure = rmeas,
                 statenames=c("S","E","I","R"),
                 obsnames=c("Y_S","Y_E","Y_I","Y_R"))

(sim_data <- as.data.frame(sim1))

ggplot(data=sim_data) + 
  geom_line(aes(x=time, y=S),color='green') +
  geom_line(aes(x=time, y=E), color='orange') +
  geom_line(aes(x=time, y=I), color='red') +
  geom_line(aes(x=time, y=R), color='blue') +
  ylab('Values')+xlab('date')

ggplot(data=sim_data) + 
  geom_line(aes(x=time, y=Y_S),color='green') +
  geom_line(aes(x=time, y=Y_E), color='orange') +
  geom_line(aes(x=time, y=Y_I), color='red') +
  geom_line(aes(x=time, y=Y_R), color='blue') +
  ylab('Values')+xlab('date')

S_pal <- c(100525, 102086, 103762, 105456, 107013, 108600, 110174, 111667,
           113014, 114496, 115876, 117254, 118500, 119436, 120605, 121626,
           122701, 123779, 124648, 125466, 126258, 126873, 127729, 128447,
           129143, 129725, 130341, 130859, 131321, 131679, 131874, 132228,
           132399, 132544, 132486, 132411, 132118, 131808, 131546, 131166,
           130852, 130242, 129805, 129071, 128351, 127664, 126769, 125827,
           124683, 123670, 122634, 121426, 120445, 119216, 117762, 116361,
           114932, 113570, 112343, 110832, 109162, 107364, 105721, 104133,
           102343, 100643,  98877,  97023,  95260,  93501,  91736,  90048,
           88327,  86659,  85151,  83570,  81896,  80331,  78959,  77568,
           76201,  74857,  73540,  72159,  71034,  69973,  68810,  67663,
           66723,  66007,  65031,  64181,  63397,  62608,  61930,  61150,
           60548,  60065,  59489,  58790,  58254,  57863,  57504,  57244,
           57065,  56790,  56596,  56394,  56163,  56050,  55899,  55686,
           55514,  55308,  55212,  55303,  55243,  55422,  55433,  55467,
           55546,  55799,  56116,  56329,  56308,  56260,  56409,  56524,
           56709,  56991,  57198,  57423,  57533,  57784,  57911,  58149,
           58123,  58443,  58587,  58794,  58950,  59232,  59459,  59713,
           59912,  60210,  60359,  60615,  60850,  60981,  61100,  61181,
           61430,  61531,  61850,  61950,  62134,  62259,  62271,  62396,
           62499,  62716,  62938,  62994,  63206,  63285,  63478,  63535,
           63529,  63706,  63913,  64087,  64186,  64321,  64528,  64612,
           64603,  64729,  64836,  65046,  65040,  65176,  65250,  65345,
           65421,  65471,  65662,  65814,  65943,  66038,  66129,  66201,
           66294,  66458,  66574,  66661,  66682,  66722,  66750,  66790)
S_pomp <- sim_data$S
ggplot() + 
  geom_line(aes(x=time, y=S_pal), color='purple') +
  geom_line(aes(x=time, y=S_pomp),color='green') +labs(x="time", y="S")

tot_pomp <- sim_data$S+sim_data$E+sim_data$I+sim_data$R
tot_pal <- c(101967, 103872, 105928, 107930, 109811, 111713, 113557, 115331,
             116980, 118769, 120447, 122184, 123758, 125072, 126554, 127963,
             129387, 130877, 132190, 133514, 134786, 135988, 137359, 138676,
             139923, 141146, 142360, 143519, 144711, 145853, 146939, 148106,
             149128, 150176, 151058, 152093, 152966, 153782, 154680, 155603,
             156572, 157398, 158309, 159054, 159781, 160630, 161451, 162208,
             162724, 163445, 164246, 164937, 165753, 166489, 167035, 167648,
             168214, 168887, 169676, 170277, 170790, 171300, 171870, 172458,
             172981, 173371, 173905, 174496, 174862, 175317, 175717, 176111,
             176615, 176954, 177513, 178052, 178559, 178970, 179455, 179946,
             180499, 180991, 181381, 181688, 182159, 182617, 183000, 183343,
             183707, 184055, 184333, 184722, 185145, 185443, 185747, 186013,
             186349, 186649, 186943, 187202, 187452, 187791, 188006, 188355,
             188694, 189021, 189307, 189550, 189759, 189957, 190142, 190329,
             190506, 190726, 190867, 191165, 191278, 191568, 191618, 191778,
             191865, 192265, 192663, 192962, 193094, 193213, 193309, 193481,
             193715, 193883, 193976, 194065, 193992, 194033, 194035, 194200,
             194191, 194277, 194315, 194505, 194547, 194721, 194981, 195053,
             195083, 195303, 195439, 195635, 195820, 195882, 195958, 195905,
             196000, 196162, 196157, 196181, 196222, 196247, 196215, 196252,
             196324, 196412, 196437, 196486, 196633, 196685, 196779, 196842,
             196781, 196827, 197020, 197035, 197066, 197189, 197247, 197326,
             197265, 197352, 197302, 197306, 197337, 197369, 197405, 197390,
             197355, 197400, 197522, 197704, 197892, 198019, 198036, 198109,
             198199, 198312, 198358, 198357, 198386, 198364, 198332, 198453)

ggplot() + 
  geom_line(aes(x=time, y=tot_pomp), color='purple') +
  geom_line(aes(x=time, y=tot_pal),color='green')+labs(x="time",y="Total #")
````