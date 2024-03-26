// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;

IntegerVector rmultinom_1(unsigned int size, NumericVector probs, unsigned int N) {
  IntegerVector outcome(N);
  rmultinom(size, probs.begin(), N, outcome.begin());
  return outcome;
}


mat rmultinom_rcpp(unsigned int n, unsigned int size, NumericVector probs) {
  unsigned int N = probs.length();
  IntegerMatrix sim(N, n);
  for (unsigned int i = 0; i < n; i++) {
    sim(_,i) = rmultinom_1(size, probs, N);
  }
  mat out = as<arma::mat>(sim);
  return out;
}

void update_K(mat& K,
              vec state,
              vec betas,
              int N,
              int t,
              double phase,
              double rho
){
  // seasonal component
  
  mat B(3,3, fill::zeros);
  for(int i = 0; i<3;i++){
    for(int j = 0; j<3;j++){
      B(i,j) = betas(i); 
    }
  }
  
  double kappa = (1 + rho*cos(2*3.141593*t/52 + phase));
  
  vec I;
  I = {state(1),state(4),state(7)};
  vec lambda = B*I*kappa/N;
  K(0,1) = 1-exp(-0.25*lambda(0));
  K(3,4) = 1-exp(-0.25*lambda(1));
  K(6,7) = 1-exp(-0.25*lambda(2));
  
  
  K(0,0)= 1 - K(0,1)-  K(0,2) - K(0,3) - K(0,4) - K(0,5) - K(0,6) - K(0,7) - K(0,8);
  K(1,1)= 1 - K(1,0) - K(1,2) - K(1,3) - K(1,4) - K(1,5) - K(1,6) - K(1,7) - K(1,8);
  K(2,2)= 1 - K(2,0) - K(2,1) - K(2,3) - K(2,4) - K(2,5) - K(2,6) - K(2,7) - K(2,8);
  K(3,3)= 1 - K(3,0) - K(3,1) - K(3,2) - K(3,4) - K(3,5) - K(3,6) - K(3,7) - K(3,8);
  K(4,4)= 1 - K(4,0) - K(4,1) - K(4,2) - K(4,3) - K(4,5) - K(4,6) - K(4,7) - K(4,8);
  K(5,5)= 1 - K(5,0) - K(5,1) - K(5,2) - K(5,3) - K(5,4) - K(5,6) - K(5,7) - K(5,8);
  K(6,6)= 1 - K(6,0) - K(6,1) - K(6,2) - K(6,3) - K(6,4) - K(6,5) - K(6,7) - K(6,8);
  K(7,7)= 1 - K(7,0) - K(7,1) - K(7,2) - K(7,3) - K(7,4) - K(7,5) - K(7,6) - K(7,8);
  K(8,8)= 1 - K(8,0) - K(8,1) - K(8,2) - K(8,3) - K(8,4) - K(8,5) - K(8,6) - K(8,7);
  // if(t==1){cout << K(0,1) << endl;}
}


// [[Rcpp::export]]
mat rotavirus_sim(int length,
                  int N,
                  vec init_dist,
                  mat betas,
                  mat q,
                  double phase,
                  double amp
){
  
  mat state(length*4, 9, fill::zeros);
  state.row(0) = init_dist.t();
  
  mat transitions(9,9,fill::zeros);
  mat weekly_incidence(length,3, fill::zeros);
  mat K(9,9,fill::value(0));
  double h = 0.25;
  // Recovery rates;
  double gamma = 1;
  K(1,2) = 1 - exp(-h*gamma);
  K(4,5) = 1 - exp(-h*gamma);
  K(7,8) = 1 - exp(-h*gamma);
  // aging rates
  double delta_1 = 0.003636364;
  double delta_2 = 0.0003496503;
  K(0,3) = 1 - exp(-h*delta_1);
  K(1,4) = 1 - exp(-h*delta_1);
  K(2,5) = 1 - exp(-h*delta_1);
  K(3,6) = 1 - exp(-h*delta_2);
  K(4,7) = 1 - exp(-h*delta_2);
  K(5,8) = 1 - exp(-h*delta_2);
  // immunity waning rates;
  double omega =  0.01923077;
  K(2,0) = 1 - exp(-h*omega);
  K(5,3) = 1 - exp(-h*omega);
  K(8,6) = 1 - exp(-h*omega);
  
  for(int t = 1; t<length; t++){
    // if(t==52){cout << state.row(t*4-5)<<endl;}
    for(int i = 0; i<4; i++){
      
      update_K(K,state.row(4*(t)-4+i).t(), betas.row(t-1).t(),N, t, phase, amp);
      //cout << K << endl;
      for(int m = 0; m < 9; m++){
        transitions.row(m) = rmultinom_rcpp(1,state(4*(t)-4+i,m), as<Rcpp::NumericVector>(wrap(K.row(m)))).t();
      }
      
      state.row(4*(t)-4+i+1) = sum(transitions,0);
      weekly_incidence(t-1,0) += transitions(0,1);
      weekly_incidence(t-1,1) += transitions(3,4);
      weekly_incidence(t-1,2) += transitions(6,7);
      state(4*(t)-4+i+1,0) += R::rpois(4*1025.7);
    }
    
  }
  //under-reporting
  mat observations(length-1,3,fill::zeros);
  
  for(int j = 0; j<3;j++){
    for(int t = 0; t<length-1; t++){
      observations(t,j) = R::rbinom(weekly_incidence(t,j),q(t,j)); 
    }
  }
  
  return(observations);
}