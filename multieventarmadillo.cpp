#include <RcppArmadillo.h>

using namespace Rcpp;

 // [[Rcpp::depends(RcppArmadillo)]]
 // [[Rcpp::export]]
 
/* multievent likelihood */
arma::mat multieventrcpp(NumericVector b, IntegerMatrix ch, IntegerVector fc, IntegerVector fs) {
  /* b = parameters */
  /* ch = capture-recapture histories (individual format) */
  /* fc = date of first capture */
  /* fs = state at first capture */
  
// OBSERVATIONS
// 0 = non-detected
// 1 = seen and ascertained as non-breeder
// 2 = seen and ascertained as breeder
// 3 = not ascertained
//   
// STATES
// 1 = alive non-breeder
// 2 = alive breeder
// 3 = dead
//   
// PARAMETERS
// phiNB  survival prob. of non-breeders
// phiB  survival prob. of breeders
// pNB  detection prob. of non-breeders
// pB  detection prob. of breeders
// psiNBB transition prob. from non-breeder to breeder
// psiBNB transition prob. from breeder to non-breeder
// piNB prob. of being in initial state non-breeder
// deltaNB prob to ascertain the breeding status of an individual encountered as non-breeder
// deltaB prob to ascertain the breeding status of an individual encountered as breeder
//   
// logit link for all parameters
// note: below, we decompose the state and obs process in two steps composed of binomial events, 
// which makes the use of the logit link appealing; 
// if not, a multinomial (aka generalised) logit link should be used

  int km = ch.nrow();
  int nh = ch.ncol();  
  int npar = b.size();
  NumericVector par(npar);
  for (int i = 0; i < npar; i++) {
    par(i) = 1.0 / (1.0 + exp(-b(i)));
  }
  double piNB = par(0); // careful, indexing starts at 0 in Rcpp!
  double phiNB = par(1);
  double phiB = par(2);
  double psiNBB = par(3);
  double psiBNB = par(4);
  double pNB = par(5);
  double pB = par(6);
  double deltaNB = par(7);
  double deltaB = par(8);
  
  /* prob of obs (rows) cond on states (col) */
  arma::mat B1(3, 3);
  B1(0,0) = 1.0-pNB;
  B1(0,1) = pNB;
  B1(0,2) = 0.0;
  B1(1,0) = 1.0-pB;
  B1(1,1) = 0.0;
  B1(1,2) = pB;
  B1(2,0) = 1.0;
  B1(2,1) = 0.0;
  B1(2,2) = 0.0;
  
  arma::mat B2(3, 4);
  B2(0,0) = 1.0;
  B2(0,1) = 0.0;
  B2(0,2) = 0.0;
  B2(0,3) = 0.0;
  B2(1,0) = 0.0;
  B2(1,1) = deltaNB;
  B2(1,2) = 0.0;
  B2(1,3) = 1-deltaNB;
  B2(2,0) = 0.0;
  B2(2,1) = 0.0;
  B2(2,2) = deltaB;
  B2(2,3) = 1-deltaB;
  
  arma::mat B(3, 4);
  B = trans(B1 * B2);
  /* B = t(matrix(c(1-pNB,pNB*deltaNB,0,pNB*(1-deltaNB),1-pB,0,pB*deltaB,pB*(1-deltaB),1,0,0,0),nrow=3,ncol=4,byrow=T)) */
  
  /* first encounter */
  arma::mat BE1(3, 3);
  BE1(0,0) = 0.0;
  BE1(0,1) = 1.0;
  BE1(0,2) = 0.0;
  BE1(1,0) = 0.0;
  BE1(1,1) = 0.0;
  BE1(1,2) = 1.0;
  BE1(2,0) = 1.0;
  BE1(2,1) = 0.0;
  BE1(2,2) = 0.0;
  
  arma::mat BE2(3, 4);
  BE2(0,0) = 1.0;
  BE2(0,1) = 0.0;
  BE2(0,2) = 0.0;
  BE2(0,3) = 0.0;
  BE2(1,0) = 0.0;
  BE2(1,1) = deltaNB;
  BE2(1,2) = 0.0;
  BE2(1,3) = 1.0-deltaNB;
  BE2(2,0) = 0.0;
  BE2(2,1) = 0.0;
  BE2(2,2) = deltaB;
  BE2(2,3) = 1.0-deltaB;
  
  arma::mat BE(4,3);
  BE = trans(BE1 * BE2);

  /* prob of states at t+1 given states at t */
  arma::mat A1(3, 3);
  A1(0,0) = phiNB;
  A1(0,1) = 0.0;
  A1(0,2) = 1.0-phiNB;
  A1(1,0) = 0.0;
  A1(1,1) = phiB;
  A1(1,2) = 1.0-phiB;
  A1(2,0) = 0.0;
  A1(2,1) = 0.0;
  A1(2,2) = 1.0;
  
  arma::mat A2(3, 3);
  A2(0,0) = 1.0-psiNBB;
  A2(0,1) = psiNBB;
  A2(0,2) = 0.0;
  A2(1,0) = psiBNB;
  A2(1,1) = 1.0-psiBNB;
  A2(1,2) = 0.0;
  A2(2,0) = 0.0;
  A2(2,1) = 0.0;
  A2(2,2) = 1.0;

  arma::mat A(3, 3);
  A = A1 * A2;
  /* A <- matrix(c(phiNB*(1-psiNBB),phiNB*psiNBB,1-phiNB,phiB*psiBNB,phiB*(1-psiBNB),1-phiB,0,0,1),nrow=3,ncol=3,byrow=T) */
    
  /* init states */
  arma::rowvec PROP(3);
  PROP(0) = piNB;
  PROP(1) = 1.0-piNB;
  PROP(2) = 0.0;
  
  /* likelihood */
  arma::mat lik(1,3);
  arma::mat ALPHA(1,3);
    for (int i = 0; i < nh; i++) {
    double ei = fc(i)-1;
    int oe = fs(i);
    IntegerVector evennt = ch(_,i);
    ALPHA = PROP % BE.row(oe);
        for (int j = ei+1; j < km; j++) {
        ALPHA = (ALPHA * A) % B.row(evennt(j));
        }
    lik += log(sum(ALPHA.row(0)));
        }
  lik = -lik(0,0);
  return lik;
}

/*** R

set.seed(1)

# read in data
data = read.table('titis2.txt')
#data = rbind(data,data,data,data,data) # increase sample size artificially

# define various quantities
nh <- dim(data)[1]
k <- dim(data)[2]
km1 <- k-1

# counts
eff <- rep(1,nh)
  
# compute the date of first capture fc, and state at initial capture init.state
fc <- NULL
init.state <- NULL
for (i in 1:nh){
  temp <- 1:k
  fc <- c(fc,min(which(data[i,]!=0)))
  init.state <- c(init.state,data[i,fc[i]])
}

# init values
binit <- runif(9)
  
# transpose data
data <- t(data)

#---------- standard implementation of likelihood in R

devMULTIEVENT <- function(b,data,eff,e,garb,nh,km1){
    
# data encounter histories, eff counts
# e vector of dates of first captures
# garb vector of initial states 
# km1 nb of recapture occasions (nb of capture occ - 1)
# nh nb ind
    
# OBSERVATIONS (+1)
# 0 = non-detected
# 1 = seen and ascertained as non-breeder
# 2 = seen and ascertained as breeder
# 3 = not ascertained
    
# STATES
# 1 = alive non-breeder
# 2 = alive breeder
# 3 = dead
    
# PARAMETERS
# phiNB  survival prob. of non-breeders
# phiB  survival prob. of breeders
# pNB  detection prob. of non-breeders
# pB  detection prob. of breeders
# psiNBB transition prob. from non-breeder to breeder
# psiBNB transition prob. from breeder to non-breeder
# piNB prob. of being in initial state non-breeder
# deltaNB prob to ascertain the breeding status of an individual encountered as non-breeder
# deltaB prob to ascertain the breeding status of an individual encountered as breeder
    
# logit link for all parameters
# note: below, we decompose the state and obs process in two steps composed of binomial events, 
# which makes the use of the logit link appealing; 
# if not, a multinomial (aka generalised) logit link should be used
    par = plogis(b)
    piNB <- par[1]
    phiNB <- par[2]
    phiB <- par[3]
    psiNBB <- par[4]
    psiBNB <- par[5]
    pNB <- par[6]
    pB <- par[7]
    deltaNB <- par[8]
    deltaB <- par[9]
    
# prob of obs (rows) cond on states (col)
    B1 = matrix(c(1-pNB,pNB,0,1-pB,0,pB,1,0,0),nrow=3,ncol=3,byrow=T)
    B2 = matrix(c(1,0,0,0,0,deltaNB,0,1-deltaNB,0,0,deltaB,1-deltaB),nrow=3,ncol=4,byrow=T)
    B = t(B1 %*% B2)
#B = t(matrix(c(1-pNB,pNB*deltaNB,0,pNB*(1-deltaNB),1-pB,0,pB*deltaB,pB*(1-deltaB),1,0,0,0),nrow=3,ncol=4,byrow=T))
    
# first encounter
    BE1 = matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3,byrow=T)
    BE2 = matrix(c(1,0,0,0,0,deltaNB,0,1-deltaNB,0,0,deltaB,1-deltaB),nrow=3,ncol=4,byrow=T)
    BE = t(BE1 %*% BE2) 
#BE = t(matrix(c(0,deltaNB,0,(1-deltaNB),0,0,deltaB,(1-deltaB),1,0,0,0),nrow=3,ncol=4,byrow=T))
    
# prob of states at t+1 given states at t
    A1 <- matrix(c(phiNB,0,1-phiNB,0,phiB,1-phiB,0,0,1),nrow=3,ncol=3,byrow=T)
    A2 <- matrix(c(1-psiNBB,psiNBB,0,psiBNB,1-psiBNB,0,0,0,1),nrow=3,ncol=3,byrow=T)
    A <- A1 %*% A2
#A <- matrix(c(phiNB*(1-psiNBB),phiNB*psiNBB,1-phiNB,phiB*psiBNB,phiB*(1-psiBNB),1-phiB,0,0,1),nrow=3,ncol=3,byrow=T)
    
# init states
    PI <- c(piNB,1-piNB,0)
    
# likelihood
    l <- 0
    for (i in 1:nh) # loop on ind
   {
      ei <- e[i] # date of first det
      oe <- garb[i] + 1 # init obs
      evennt <- data[,i] + 1 # add 1 to obs to avoid 0s in indexing
      ALPHA <- PI*BE[oe,]
     for (j in (ei+1):(km1+1)) # cond on first capture
     {
        if ((ei+1)>(km1+1)) {break} # sous MATLAB la commande >> 8:7 rend >> null, alors que sous R, ça rend le vecteur c(8,7)!
        ALPHA <- (ALPHA %*% A)*B[evennt[j],]
      }
      l <- l + log(sum(ALPHA))#*eff[i]
    }
    l <- -l
    l
  }

# evaluate the likelihood
devMULTIEVENT(binit,data,eff,fc,init.state,nh,km1) # standard R
multieventrcpp(binit,data,fc,init.state) # Rcpp

 # fit model with standard lik
 deb=Sys.time()
 tmpmin1 <- optim(binit,devMULTIEVENT,NULL,hessian=T,data,eff,fc,init.state,nh,km1,method="BFGS")
 fin=Sys.time()
 res1 = fin-deb
 
 # fit model with Rcpp lik
 deb=Sys.time()
 tmpmin2 <- optim(binit,multieventrcpp,NULL,hessian=T,data,fc,init.state,method="BFGS")
 fin=Sys.time()
 res2 = fin-deb
 
 # standard implementation
 res1
 tmpmin1$par
 
 # Rcpp implementation
 res2
 tmpmin2$par

 # The optimization is not stochastic, but depending on what else I'm doing, comptation times may vary, hence a benchmark
 library(microbenchmark)
 res = microbenchmark(optim(binit,devMULTIEVENT,NULL,hessian=T,data,eff,fc,init.state,nh,km1,method="BFGS"),
   optim(binit,multieventrcpp,NULL,hessian=T,data,fc,init.state,method="BFGS"),
   times=5
 )
 
 res2 = summary(res)
 
*/