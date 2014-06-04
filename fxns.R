## fxns.R

## **************************************************************
## score test screening
## if the first column of x all have the same value, this
## assumes there's an intercept and removes it from screening
## y = outcome
## x = covariates, n x p, include intercept as column of 1's
## init = initial beta, set to all 0 by default
## U = score equation, assume marginal score equations are given
## by U(0)
## ... = additional arguments
## **************************************************************
ss <- function(y,x,init=NULL,U,...)
{
    sds <- apply(x,2,sd,na.rm=TRUE);
    sds[which(sds==0)] <- 1;
    x <- scale(x,center=FALSE,scale=sds);
    b <- init; if(is.null(init)){ b <- rep(0,ncol(x)); }
    ret <- order(-abs(U(y,x,b,...)));
    if(length(unique(x[,1]))==1){ ret <- setdiff(ret,1)-1; }
    return(ret);
}

## **************************************************************
## reproducibility selection
## if the first column of x all have the same value, this
## assumes there's an intercept and removes it from screening
## y = outcome
## x = covariates, n x p, include intercept as column of 1's
## screen = function used for screening
## step = step size for number to retain, to reduce computation
## B = number of bootstrap samples
## init = initial beta, set to all 0 by default
## U = score equation, assume marginal score equations are given
## by U(0)
## ... = additional arguments
## **************************************************************
## return ranks after score screening
ssr <- function(y,x,init=NULL,U,...)
{
    sds <- apply(x,2,sd,na.rm=TRUE);
    sds[which(sds==0)] <- 1;
    x <- scale(x,center=FALSE,scale=sds);
    b <- init; if(is.null(init)){ b <- rep(0,ncol(x)); }
    r <- rank(-abs(U(y,x,b,...))); ## ranks of original data
    if(length(unique(x[,1]))==1){ r <- r[-1]; } ## remove intercept
    return(r);
}

rs <- function(y,x,step=1,B=100,init=NULL,U,...)
{
    n <- nrow(x);
    b <- init; if(is.null(init)){ b <- rep(0,ncol(x)); }
    U.m <- ssr(y,x,b,U,...);

    ## bootstrap
    cat("bootstrapping and screening...");
    U.B <- matrix(NA,nrow=B,ncol=length(U.m));
    for(i in 1:B)
    {
        ind <- sample(1:n,n,replace=TRUE);
        U.B[i,] <- ssr(as.matrix(y)[ind,],x[ind,],b,U,...);
    }

    ## find standardized difference between obs and expected overlaps
    ## difference must be positive
    cat("calculating overlaps...\n");
    ## reorder covariates to make calculating easier
    ord <- order(U.m);
    U.B <- U.B[,ord];
    p <- length(U.m);
    j <- seq(step,p,by=step);
    oj <- rep(0,length(j));
    for(k in 1:length(j))
    {
        if(j[k]%%500==0){ print(j[k]); }
        ## calculate overlap
        oj[k] <- mean(apply(as.matrix(U.B[,1:j[k]])<=j[k],1,sum));
    }

    ## standardize
    v <- j^2/p*(p-j)/p*(p-j)/(p-1); ## variance
    diff <- (oj-j^2/p)/sqrt(v/B); ## diff must be positive
    diff[is.na(diff)] <- 0; ## var = 0 if keeping all covs
    
    return(list(order=ord,oj=oj,diff=diff,cutoff=j[which.max(diff)]));
}

## **************************************************************
## iterative score test screening using projected subgradient
## y = outcome
## x = covariates, n x p, include intercept as column of 1's
## R = l1 norm of true beta
## a = square summable step size scaling, default is "optimal"
## tol = convergence of U(b) to 0
## maxit = max number of iterations
## init = initial value of parameter vector
## standardize = standardize covariates to variance 1
## silent = report progress
## U = estimating equation
## ... = additional arguments to be passed to U
## **************************************************************
## soft thresholding function
Sl <- function(x,l){ return(pmax(0,abs(x)-l)*sign(x)); }
## euclidean project of x onto l1 ball of radius R
Pr <- function(x,R)
{
    if(sum(abs(x))<=R){ return(x); }
    xstar <- sort(abs(x),decreasing=TRUE);
    Sks <- sapply(xstar,function(z){ sum(abs(Sl(xstar,z))); });
    k <- max(which(Sks<=R));
    mu <- xstar[k]-(R-Sks[k])/k;
    return(Sl(x,mu));
}

## project subgradient method
psm <- function(y,x,R,a=NULL,tol=1e-4,maxit=100,init=NULL,
                standardize=TRUE,silent=TRUE,U,...)
{
    ## get rid of missing data
    cc <- complete.cases(y,x);
    y <- as.matrix(y)[cc,];
    x <- as.matrix(as.matrix(x)[cc,]);

    ## standardize covariates
    if(standardize)
    {
        if(!silent){ cat("standardizing..."); }
        sds <- apply(x,2,sd);
        sdnz <- which(sds!=0); ## sd of intercept covariate will be 0
        init[sdnz] <- init[sdnz]*sds[sdnz];
        sds[-sdnz] <- 1;
        x <- scale(x,center=FALSE,scale=sds);
    }

    ## iteration
    if(!silent){ cat("iterating...\n"); }
    b <- init; if(is.null(b)){ b <- rep(0,ncol(x)); }
    oldb <- b+1;
    t <- 1;
    ## return best solution
    g2best <- Inf; bbest <- b;
    ## "optimal" step size, treat R as l2 norm of true beta
    if(is.null(a))
    { a <- R/pi*sqrt(6/sum(U(y,x,b,...)^2)); }
    while(t<=maxit)
    {
        g <- U(y,x,b,...);
        g2 <- sum(g^2);
        if(g2<tol){ cat("gradient=0\n"); break; }
        if(g2<=g2best){ g2best <- g2; bbest <- b; }
        b <- Pr(b+a/t*g,R);
        if(!silent)
        {
            cat(t,"step",round(a/t,4),"old gradient:",round(g2,4),
                "current variables:",sum(b!=0),"\n");
        }
        if(max(abs(b-oldb))<tol){ cat("beta converged\n"); break; }
        oldb <- b;
        if(t==maxit){ cat("max reached\n"); }
        t <- t+1;
    }

    ## report the beta that best solves the estimating equation
    ## unstandardize
    if(standardize){ bbest <- bbest/sds; }

    ## last threshold on the small values
    bbest[abs(bbest)<=tol] <- 0;
    return(b=bbest);
}

## **************************************************************
## example estimating equations U
## **************************************************************
ols <- function(y,x,b)
{
    y <- as.matrix(y); x <- as.matrix(x);
    return(t(x)%*%(y-x%*%b)/length(y));
}

## ==============================================================
## censored quantile regression
## h = transformation function
## tau = quantile
## ==============================================================
U.cqr <- function(outcomes,covariates,b,h,tau)
{
  hY <- h(outcomes[,1]);
  eta <- covariates%*%b;
  SchY <- KM(hY,1-outcomes[,2],hY);
  Sceta <- KM(hY,1-outcomes[,2],eta);

  term1 <- tau*as.numeric(hY>eta);
  term2 <- (1-tau)*Sceta/SchY;
  term2[hY>eta|outcomes[,2]==0] <- 0;

  return(t(covariates)%*%(term1-term2)/length(hY));
}

## ==============================================================
## cox w normal measurement error
## for normal measurement error
## Db.f = function that returns D(beta)%*%beta, takes args (b,...)
## writing this fxn will be faster than the matrix multiplication
## when S is huge.
## ==============================================================
U.coxme <- function(outcomes,Z,b,Db.f,...)
{
  X <- outcomes[,1];
  d <- outcomes[,2];
  Z <- as.matrix(Z);

  ## rows = death times, cols = subjs, 1 if at risk, 0 else
  atrisk <- drop(outer(drop(X[d==1]),drop(X),"<="));
  score <- as.matrix(exp(Z%*%b));
  ## value of Sk at each deathtime (row)
  S0 <- as.matrix(atrisk%*%score);
  S1 <- as.matrix(atrisk%*%(drop(score)*Z));

  ## correction terms
  cZ <- as.matrix(Z[d==1,])+
    matrix(drop(Db.f(b,...)),nrow=sum(d),ncol=ncol(Z),byrow=TRUE);
  cS0 <- -exp(as.matrix(Z[d==1,])%*%b)+exp(cZ%*%b);
  cS1 <- -as.matrix(Z[d==1,])*drop(exp(as.matrix(Z[d==1,])%*%b))+
    cZ*drop(exp(cZ%*%b));
  S0 <- S0+cS0; S1 <- S1+cS1;

  U <- as.matrix(apply(cZ,2,sum))-
    as.matrix(apply(S1/drop(S0),2,sum));
  
  #### ordinary cox model
  ##U <- as.matrix(apply(as.matrix(Z[d==1,]),2,sum))-
  ##    as.matrix(apply(S1/drop(S0),2,sum));

  return(U/nrow(Z));
}

## ===============================================================
## aft rank estimating equation
## survival time should be log-transformed
## code from Lu Tian
## ===============================================================
U.aft <- function(obs,covs,beta)
{
  # get rid of missing data
  if(length(covs)==0){ cc <- complete.cases(obs); } else
  { cc <- complete.cases(cbind(obs,covs)); }
  X <- obs[cc,1];
  d <- obs[cc,2];
  Z <- matrix(as.matrix(covs)[cc,],nrow=sum(cc));
  n <- length(X);
  
  # calculate residuals
  e <- drop(X-Z%*%beta);
  index <- order(e);
  e <- e[index];
  d <- d[index];
  Z <- as.matrix(Z[index,]);
  
  w1 <- n:1
  w2 <- cumsum(d);
  
  score=Z*d*w1-Z*w2;
  score=apply(score,2,sum);
  return(-score/n^2);
}

## **************************************************************
## helper functions for estimating equation functions
## **************************************************************
library(survival);
KM <- function(X,d,x)
{
  S <- survfit(Surv(X,d)~1);
  # for numerical precision
  Stime <- round(S$time,12);
  x <- round(x,12);
  ret <- sapply(x,function(xx)
    {
      if(xx<min(Stime)){ return(1); } else
      { ti <- max(which(Stime<=xx)); return(S$surv[ti]); }
    });
  return(ret);
}

## quantile function for survival dist using KM estimate
## q = P(T<=t), this function returns t
KMquant <- function(X,d,q)
{
  S <- survfit(Surv(X,d)~1);
  # for numerical precision
  Ssurv <- round(S$surv,12);
  q <- round(q,12);

  ret <- sapply(q,function(qq)
    { return(min(S$time[Ssurv<=(1-qq)])); });
  return(ret);
}
