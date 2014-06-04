## cqrp.R

rm(list=ls());
source("fxns.R");

env <- as.numeric(Sys.getenv("LSB_JOBINDEX"));
set.seed(env);

## **************************************************************
## generate data
## **************************************************************
cqr <- function(n,b,nz,r,c,tau)
{
  Z <- t(replicate(n,arima.sim(list(order=c(1,0,0),ar=r),n=length(b))))*sqrt(1-r^2);
  Z[Z>=2] <- 2;
  Z[Z<=-2] <- -2;
  T <- Z%*%b+rnorm(n)*(9+(Z[,nz[11]]-Z[,nz[1]])*b[nz[1]]/qnorm(tau));
  cen <- rexp(n,c);
  d <- as.numeric(T<=cen);
  X <- drop(pmin(T,cen));
  return(list(Z=Z,T=T,cen=cen,d=d,X=X));
}

## parameters for generating data
n <- 400;
p <- 10000;
nz <- seq(5,55,by=5);
b <- rep(0,p);
## at median, nz[1:10] are important; at tau, nz[2:11] are important
## 6th cov is marginally unimportant
b[nz] <- c(1.5,1.5,1.5,1.5,1.5,-1.5,1.5,1.5,1.5,1.5,0);
r <- 0.8;
c <- 0.15; # ~30% censoring
tau <- 0.25;

## generate data
int <- cqr(n,b,nz,r,c,tau); ## for calculating intercept
train <- cqr(n,b,nz,r,c,tau);

## parameters for iterative screening
R <- 20;
a <- NULL;
tol <- 1e-4;
maxit <- 250;
standardize <- TRUE;
silent <- FALSE;
U <- U.cqr;
h <- function(x){x};

## initial betas for score test screening
init.5 <- c(KMquant(int$X,int$d,0.5),rep(0,ncol(train$Z)));
init.tau <- c(KMquant(int$X,int$d,tau),rep(0,ncol(train$Z)));

toremove <- setdiff(ls(),"env");

## **************************************************************
## wald
## **************************************************************
library(quantreg);
time.wald <- system.time(wald <- apply(train$Z,2,function(z)
    {
        fit <- crq(Surv(train$X,train$d)~z,
                   method="PengHuang");
        return(coef(fit,taus=c(tau,0.5))[-1,]);
    }))[3];
wald.5 <- order(-abs(wald[2,]));
wald.tau <- order(-abs(wald[1,]));

## **************************************************************
## score
## **************************************************************
## median
time.ss <- system.time(
    ss.5 <- ss(cbind(train$X,train$d),cbind(1,train$Z),init.5,
               U.cqr,h,0.5))[3];
rs.5 <- rs(cbind(train$X,train$d),cbind(1,train$Z),step=10,B=100,
           init=init.5,U.cqr,h,0.5);

## tau
time.ss <- time.ss+system.time(
    ss.tau <- ss(cbind(train$X,train$d),cbind(1,train$Z),init.tau,
               U.cqr,h,tau))[3];
rs.tau <- rs(cbind(train$X,train$d),cbind(1,train$Z),step=10,B=100,
             init=init.tau,U.cqr,h,tau);

## **************************************************************
## iterative
## **************************************************************
## median
psm.5 <- psm(cbind(train$X,train$d),cbind(1,train$Z),
             R,a,tol,maxit,init=init.5,
             standardize=TRUE,silent,U,h,tau=0.5);

## tau
psm.tau <- psm(cbind(train$X,train$d),cbind(1,train$Z),
               R,a,tol,maxit,init=init.tau,
               standardize,silent,U,h,tau=tau);

rm(list=c(toremove,"toremove"));
save.image(paste("results/cqrp",env,".RData",sep=""));
