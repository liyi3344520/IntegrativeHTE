
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IntegrativeHTE

The goal of *IntegrativeHET* is to implement integrative analyses for
the heterogenous treatment effect combining a randomized trial and a
real-world evidence study.

Two datasets

  - The randomized trial contains observations on (A,X,Y), where the
    treatment assignment A is randomized.

  - The real-world evidence study contains observations on (A,X,Y),
    where the treatment assignment A may be confounded.

## Installation with `devtools`:

``` r
devtools::install_github("shuyang1987/IntegrativeHTE")
```

### Main Paper: coming soon

The reference paper will come soon.

### Usage

IntHTE(A, X, Y, delta, family.oc, gamma1=0.90, gamma2=2, nptb=50,
nboots=50,
alpha.CI=0.05)

### Arguments

| Argument  |                                                                                                                                                                  |
| --------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| A         | is a vector of treatment ( (n+m) x 1), where n is the RT sample size and m is the RW sample size                                                                 |
| X         | is a matrix of covariates without intercept ( (n+m) x p)                                                                                                         |
| Y         | is a vector of outcome ( (n+m) x 1)                                                                                                                              |
| delta     | is a vector of the binary indicator of belonging to the randomized trial (RT) sample; i.e., 1 if the unit belongs to the RT sample, and 0 otherwise ( (n+m) x 1) |
| family.oc | specifies the family for the outcome model                                                                                                                       |
| –         | “gaussian”: a linear regression model for the continuous outcome                                                                                                 |
| –         | “binomial”: a logistic regression model for the binary outcome                                                                                                   |
| nptb      | is the number for replication-based method for estimating Sigma\_SS                                                                                              |
| gamma1    | is the parameter in the violation test such that T \> chisq\_{1-gamma1} indicates the bias of the RW analysis is significant                                     |
| gamma2    | is a smooth parameter in using the normal cummulative distribution with mean 0 and variance gamma2^2 function to approximate the indicator function              |
| nboots    | is the number of bootstrap samples                                                                                                                               |
| alpha.CI  | is the level of percentile bootstrap confidence interval                                                                                                         |

### Value

|         |                                                                                                          |
| ------- | -------------------------------------------------------------------------------------------------------- |
| est     | the HTE estimators including                                                                             |
| –       | i) the covariate adjustment estimator;                                                                   |
| –       | ii) the RT estimating equation estimator;                                                                |
| –       | iii) the RT\&RW estimating equation estimator;                                                           |
| –       | iv) the elastic integrative estimator.                                                                   |
| ve      | the bootstrap variance estimate for est. Note: ve is conservative for the elastic integrative estimator. |
| boot.CI | 1-alpha.CI percentiel bootstrap confidence interval for the elastic integrative estimator                |

## Example

This is an example for illustration.

``` r

set.seed(1234)

library(MASS)
library(rootSolve)
library(stats)
library(utils)
library(IntegrativeHTE)

n <- 1e5
beta0<- c(0,1,1)
psi0 <- c(1,1,1)

family.oc <- "binomial" #family.oc <- "gaussian"

X.pop <- cbind(rnorm(n,0,1),rnorm(n,0,1))
if( family.oc=="gaussian" ){
  Y0   <-  cbind(1,X.pop)%*%beta0+rnorm(n,0,1)
  Y1   <-  cbind(1,X.pop)%*%beta0+cbind(1,X)%*%psi0+rnorm(n,0,1)
}
if( family.oc=="binomial" ){
  lY0   <- -cbind(1,X.pop)%*%psi0
  lY1   <-  cbind(1,X.pop)%*%psi0
  pY0 <- exp(lY0)/(1+exp(lY0) )
  pY1 <- exp(lY1)/(1+exp(lY1) )
  Y0    <- rbinom(n,1, pY0 )
  Y1    <- rbinom(n,1, pY1 )
}

## RCT data
eS <- expoit(-cbind(1,X.pop)%*%c(5.5,1,1))
S <- sapply(eS,rbinom,n = 1, size = 1)
S.ind <- which(S==1)
n.t <- length(S.ind)
X.t <- X.pop[S.ind,]
Y0.t <- Y0[S.ind]
Y1.t <- Y1[S.ind]
A.t  <- rbinom(n.t,1,0.5) # randomly assign treatment to trial participant
Y.t  <-  Y1.t*A.t+Y0.t*(1-A.t)

## OS data
m <- 5000
P.ind <- sample(1:n,size = m)
X.os <- X.pop[P.ind,]
Y0.os <- Y0[P.ind]
Y1.os <- Y1[P.ind]

eA <- expoit( cbind(1,X.os)%*%c(1,-1,-1))

A.os <- sapply(eA,rbinom,n=1,size=1)
Y.os <- Y1.os*A.os+Y0.os*(1-A.os)

A<-c(A.t,A.os)
X<-rbind(X.t,X.os)
Y<-c(Y.t,Y.os)
delta<-c(rep(1,n.t),rep(0,m))

IntegrativeHTE::IntHTE(A, X, Y, delta, family.oc, gamma1=0.90, gamma2=2, nptb=50, nboots=50, alpha.CI=0.05)
#> $est
#>     covj.t.1     covj.t.2     covj.t.3  ee.rt(ml).1  ee.rt(ml).2 
#>    0.9828114    1.0873996    0.9158159    1.0019063    1.0694406 
#>  ee.rt(ml).3 opt.ee(ml).1 opt.ee(ml).2 opt.ee(ml).3       elas.1 
#>    0.9999815    0.9767846    0.9699336    1.0375235    0.9940706 
#>       elas.2       elas.3 
#>    1.0384034    1.0116912 
#> 
#> $ve
#>  ee.rt(ml).1  ee.rt(ml).2  ee.rt(ml).3 opt.ee(ml).1 opt.ee(ml).2 
#>  0.021681578  0.013069224  0.012383619  0.001755750  0.002564049 
#> opt.ee(ml).3       elas.1       elas.2       elas.3 
#>  0.002583618  0.017412566  0.013943164  0.021976628 
#> 
#> $boot.CI
#>           bootq1   bootq2
#> elas.1 0.8381128 1.323909
#> elas.2 0.8697130 1.366554
#> elas.3 0.7771233 1.541844
```
