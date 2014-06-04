**************************************************************
README file for the R code for "Score test variable
screening", by Sihai Dave Zhao and Yi Li.
**************************************************************
==============================================================
fxns.R: functions for implementing score screening
==============================================================
-- Screening functions; instructions are included as comments
   in the R code
   -- ss(y,x,init,U,...): marginal score screening function
   -- rs(y,x,step,B,init,U,...): calculates
      reproducibile screening threshold
   -- psm(y,x,R,a,tol,maxit,init,standardize,silient,U,...):
      iterative score screening using projected subgradient
      method
-- Example estimating equations
   -- U.ols: linear model
   -- U.cqr: censored linear regression
   -- U.coxme: Cox model with measurement error
   -- U.aft: accelerated failure time model

==============================================================
cqrp.R: code from Example 1 of the simulations
==============================================================
-- Generates data from a linear censored quantile regression
   model
-- Implements Wald and score screening
