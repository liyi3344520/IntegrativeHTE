#' Integrative analysis of the HTE (IntHTE)
#'
#' Implements integrative HTE analysis combing randomized trial and real-world data.
#' @param A is a vector of treatment  ( (n+m) x 1).
#' @param X is a matrix of covariates without intercept  ( (n+m) x p).
#' @param Y is a vector of outcome  ( (n+m) x 1).
#' @param delta is a vector of the binary indicator of belonging to the randomized trial (RT) sample;
#' i.e., 1 if the unit belongs to the RT sample, and 0 otherwise  ( (n+m) x 1).
#' @param family.oc specifies the family for the outcome model.
#' \code{"gaussian"}: a linear regression model for the continuous outcome.
#'
#' \code{"binomial"}: a logistic regression model for the binary outcome.
#' @param nptb is the number for replication-based method for estimating Sigma_SS.
#' @param gamma1 is the parameter in the violation test such that T > chisq_{1-gamma1} indicates the bias of the RW analysis is significant.
#' @param gamma2 is a smooth parameter in using the normal cummulative distribution with mean 0 and variance gamma2^2 function to approximate the indicator function.
#' @param nboots is the number of bootstrap samples.
#' @param alpha.CI is the level of percentile bootstrap confidence interval.
#'
#' @return
#'
#' \itemize{
#'
#' \item \code{est}: the HTE estimators including
#' i) the covariate adjustment estimator;
#' ii) the RT estimating equation estimator;
#' iii) the RT&RW estimating equation estimator;
#' iv) the elastic integrative estimator.
#'
#' \item \code{ve}:  the bootstrap variance estimate for \code{est}. Note: ve is conservative for the elastic integrative estimator.
#'
#' \item \code{boot.CI}: 1-alpha.CI percentiel bootstrap confidence interval for the elastic integrative estimator.
#'
#' }
#'
#' @details
#' Details will be provided in the reference paper.
#'
#' @import utils MASS rootSolve stats
#'
#' @references
#'
#'A reference paper will come up soon.
#'
#' @examples
#'
#'library(MASS)
#'library(rootSolve)
#'library(stats)
#'library(utils)
#'
#'set.seed(1234)
#'
#'## population size
#'n <- 1e5
#'beta0<- c(0,1,1)
#'psi0 <- c(1,1,1)
#'family.oc <- "binomial"
#'
#'X.pop <- cbind(rnorm(n,0,1),rnorm(n,0,1))
#'if( family.oc=="Gaussian" ){
#'  Y0   <-  cbind(1,X.pop)%*%beta0+rnorm(n,0,1)
#'  Y1   <-  cbind(1,X.pop)%*%beta0+cbind(1,X)%*%psi0+rnorm(n,0,1)
#'}
#'if( family.oc=="binomial" ){
#'  lY0   <- -cbind(1,X.pop)%*%psi0
#'  lY1   <-  cbind(1,X.pop)%*%psi0
#'  pY0 <- exp(lY0)/(1+exp(lY0) )
#'  pY1 <- exp(lY1)/(1+exp(lY1) )
#'  Y0    <- rbinom(n,1, pY0 )
#'  Y1    <- rbinom(n,1, pY1 )
#'}
#'
#'## RCT data
#'eS <- expoit(-cbind(1,X.pop)%*%c(5.5,1,1))
#'mean(eS)*n
#'S <- sapply(eS,rbinom,n = 1, size = 1)
#'S.ind <- which(S==1)
#'n.t <- length(S.ind)
#'n.t
#'X.t <- X.pop[S.ind,]
#'Y0.t <- Y0[S.ind]
#'Y1.t <- Y1[S.ind]
#'A.t  <- rbinom(n.t,1,0.5) # randomly assign treatment to trial participant
#'Y.t  <-  Y1.t*A.t+Y0.t*(1-A.t)
#'
#'## OS data
#'m <- 5000
#'P.ind <- sample(1:n,size = m)
#'X.os <- X.pop[P.ind,]
#'Y0.os <- Y0[P.ind]
#'Y1.os <- Y1[P.ind]
#'
#'eA <- expoit( cbind(1,X.os)%*%c(1,-1,-1))
#'mean(eA)
#'
#'A.os <- sapply(eA,rbinom,n=1,size=1)
#'Y.os <- Y1.os*A.os+Y0.os*(1-A.os)
#'
#'A<-c(A.t,A.os)
#'X<-rbind(X.t,X.os)
#'Y<-c(Y.t,Y.os)
#'delta<-c(rep(1,n.t),rep(0,m))
#'
#'gamma1=0.9
#'gamma2=2
#'nptb=5
#'nboots=2
#'alpha.CI=0.05
#'
#'IntHTE(A, X, Y, delta, family.oc,
#'gamma1, gamma2, nptb, nboots, alpha.CI)
#'
#' @export

IntHTE <- function(A, X, Y, delta, family.oc="gaussian", gamma1=0.90, gamma2=2, nptb=50, nboots=50, alpha.CI=0.05){
  options(warn=-1)

  est<-NULL
  ve <-NULL

  X<-as.matrix(X)
  loc.t<-which(delta==1)
  loc.os<-which(delta==0)
  n.t<-length(loc.t)
  m<-length(loc.os)

  ## RT data
  A.t<-A[loc.t]
  X.t<-X[loc.t,]
  Y.t<-Y[loc.t]
  dat.t  <- list(Y = Y.t, A = A.t, X = X.t, Y.t = Y.t, A.t = A.t, X.t = X.t, q = rep(1,n.t),  ml.ps = rep(0.5,n.t))

  ## RW data
  A.os<-A[loc.os]
  X.os<-X[loc.os,]
  Y.os<-Y[loc.os]
  dat.os  <- list(Y = Y.os, A = A.os, X = X.os, Y.os = Y.os, A.os = A.os, X.os = X.os, q = rep(1,m))

  ml.X.t<-cbind(X.t,X.t^2,apply(X.t, 1, utils::combn, 2, prod))
  ml.X.os<-cbind(X.os,X.os^2,apply(X.os, 1, utils::combn, 2, prod))
  dat.t$ml.X.t<-ml.X.t
  dat.os$ml.X.os<-ml.X.os

  glm.out <- glm(A.t~ml.X.t,family=quasibinomial,weights=q,data=dat.t)
  ml.ps.t <- glm.out$fitted.values
  dat.t$ml.ps <- ml.ps.t

  glm.out <- glm(A.os~ml.X.os,family=quasibinomial,weights=q,data=dat.os)
  ml.ps.os <- glm.out$fitted.values
  dat.os$ml.ps <- ml.ps.os

  if(family.oc=="gaussian"){
    mu0.out <- glm(Y.t[which(A.t==0)]~ml.X.t[which(A.t==0),],weights=q[which(A.t==0)],data=dat.t)
    ml.mu0.t <- cbind(1,ml.X.t)%*%mu0.out$coeff
    dat.t$ml.mu0 <- ml.mu0.t

    mu0.out <- glm(Y.os[which(A.os==0)]~ml.X.os[which(A.os==0),],weights=q[which(A.os==0)],data=dat.os)
    ml.mu0.os <- cbind(1,ml.X.os)%*%mu0.out$coeff
    dat.os$ml.mu0 <- ml.mu0.os

    # covariate adjustment (trial)
    X.tadj<- cbind(1,X.t)*(2*A.t-1)/2
    jj<-glm(Y.t~X.tadj-1)
    est[paste("covj.t.",1:3,sep="" )]<-reg.t<-jj$coeff/2
    reg.t0<-reg.t
  }

  if(family.oc=="binomial"){
    mu0.out <- glm(Y.t[which(A.t==0)]~ml.X.t[which(A.t==0),],weights=q[which(A.t==0)],data=dat.t,family="quasibinomial")
    lml.mu0.t <- cbind(1,ml.X.t)%*%mu0.out$coeff
    ml.mu0.t <- exp(lml.mu0.t)/( 1+exp(lml.mu0.t) )
    dat.t$ml.mu0 <- ml.mu0.t

    mu0.out <- glm(Y.os[which(A.os==0)]~ml.X.os[which(A.os==0),],weights=q[which(A.os==0)],data=dat.os,family="quasibinomial")
    lml.mu0.os <- cbind(1,ml.X.os)%*%mu0.out$coeff
    ml.mu0.os <- exp(lml.mu0.os)/(1+exp(lml.mu0.os))
    dat.os$ml.mu0 <- ml.mu0.os

    # covariate adjustment (trial)
    X.tadj<- cbind(1,X.t)*(2*A.t-1)/2
    jj<-glm(Y.t~X.tadj-1,family="quasibinomial")
    est[paste("covj.t.",1:3,sep="" )]<-reg.t<-jj$coeff/2
    reg.t0<-reg.t
  }

  dat.integ<-list( A=c(A.os,A.t),
                   X=rbind(X.os,X.t),
                   Y = c(Y.os,Y.t),
                   q =c(rep(1,m),rep(1,n.t)),
                   ml.mu0=c(ml.mu0.os,ml.mu0.t),
                   ml.ps=c(ml.ps.os,ml.ps.t),
                   ml.X=rbind(ml.X.os,ml.X.t))

  psi0<-rep(1,dim(X)[2]+1)
  n.par<-length(psi0)

  if(family.oc=="gaussian"){
    ## ee (trial) ml
    jj<-rootSolve::multiroot(f = ee1.ml, start = psi0, dat = dat.t,par0=reg.t0)
    est[paste("ee.rt(ml).",1:n.par,sep="" )]<-opt.t<-jj$root

    ## ee (integration) ml
    jj<-rootSolve::multiroot(f = ee1.ml, start = psi0, dat = dat.integ,par0=reg.t0)
    est[paste("opt.ee(ml).",1:n.par,sep="" )]<-jj$root

    ## trying a test before combining RCT and OS
    ## perturbation based variance estimation

    S.os1<-ee1.ml(par=opt.t,dat=dat.os, par0=reg.t0)
    #S.os1
  }

  if(family.oc=="binomial"){
    ## ee (trial) ml
    jj<-rootSolve::multiroot(f = ee2.ml, start = psi0, dat = dat.t,par0=reg.t0)
    est[paste("ee.rt(ml).",1:n.par,sep="" )]<-opt.t<-jj$root

    ## ee (integration) ml
    jj<-rootSolve::multiroot(f = ee2.ml, start = psi0, dat = dat.integ,par0=reg.t0)
    est[paste("opt.ee(ml).",1:n.par,sep="" )]<-jj$root


    ## trying a test before combining RCT and OS
    ## perturbation based variance estimation

    S.os1<-ee2.ml(par=opt.t,dat=dat.os, par0=reg.t0)
    #S.os1
  }

  ptb.S.os1<-matrix(0,nptb,length(S.os1))
  cnames<-c(paste("ee.rt(ml).",1:n.par,sep="" ),paste("opt.ee(ml).",1:n.par,sep="" ),paste("elas.",1:n.par,sep="" ))
  ptb<-matrix(0,nptb,length(cnames))
  colnames(ptb)<-cnames

  for(kkk in 1:nptb){
    dat.t$q <-rexp(n.t,1)
    dat.os$q<-rexp(m,1)

    glm.out <- glm(A.os~ml.X.os,family=quasibinomial,weights=q,data=dat.os)
    ml.ps.os <- glm.out$fitted.values
    dat.os$ml.ps <- ml.ps.os

    if(family.oc=="gaussian"){
      mu0.out <- glm(Y.t[which(A.t==0)]~ml.X.t[which(A.t==0),],weights=q[which(A.t==0)],data=dat.t)
      ml.mu0.t <- cbind(1,ml.X.t)%*%mu0.out$coeff
      dat.t$ml.mu0 <- ml.mu0.t

      mu0.out <- glm(Y.os[which(A.os==0)]~ml.X.os[which(A.os==0),],weights=q[which(A.os==0)],data=dat.os)
      ml.mu0.os <- cbind(1,ml.X.os)%*%mu0.out$coeff
      dat.os$ml.mu0 <- ml.mu0.os
    }

    if(family.oc=="binomial"){
      mu0.out <- glm(Y.t[which(A.t==0)]~ml.X.t[which(A.t==0),],weights=q[which(A.t==0)],data=dat.t,family="quasibinomial")
      lml.mu0.t <- cbind(1,ml.X.t)%*%mu0.out$coeff
      ml.mu0.t <- exp(lml.mu0.t)/( 1+exp(lml.mu0.t) )
      dat.t$ml.mu0 <- ml.mu0.t

      mu0.out <- glm(Y.os[which(A.os==0)]~ml.X.os[which(A.os==0),],weights=q[which(A.os==0)],data=dat.os,family="quasibinomial")
      lml.mu0.os <- cbind(1,ml.X.os)%*%mu0.out$coeff
      ml.mu0.os <- exp(lml.mu0.os)/(1+exp(lml.mu0.os))
      dat.os$ml.mu0 <- ml.mu0.os
    }

    dat.integ<-list( Y = c(Y.os,Y.t),A=c(A.os,A.t),X=rbind(X.os,X.t),
                     q =c(dat.os$q,dat.t$q),
                     ml.mu0=c(ml.mu0.os,ml.mu0.t),ml.ps=c(ml.ps.os,ml.ps.t),
                     ml.X=rbind(ml.X.os,ml.X.t))
    if(family.oc=="gaussian"){
      jj<-rootSolve::multiroot(f = ee1.ml, start = psi0, dat = dat.t, par0=reg.t0)
      ptb[kkk,paste("ee.rt(ml).",1:n.par,sep="" )]<-opt.t<-jj$root

      jj<-rootSolve::multiroot(f = ee1.ml, start = psi0, dat = dat.integ, par0=reg.t0)
      ptb[kkk,paste("opt.ee(ml).",1:n.par,sep="" )]<-jj$root

      ptb.S.os1[kkk,]<-ee1.ml(par=opt.t,dat=dat.os, par0=reg.t0)
    }
    if(family.oc=="binomial"){
      jj<-rootSolve::multiroot(f = ee2.ml, start = psi0, dat = dat.t, par0=reg.t0)
      ptb[kkk,paste("ee.rt(ml).",1:n.par,sep="" )]<-opt.t<-jj$root

      jj<-rootSolve::multiroot(f = ee2.ml, start = psi0, dat = dat.integ, par0=reg.t0)
      ptb[kkk,paste("opt.ee(ml).",1:n.par,sep="" )]<-jj$root

      ptb.S.os1[kkk,]<-ee2.ml(par=opt.t,dat=dat.os, par0=reg.t0)
    }


  }

  ve<-apply(ptb,2,var)
  mu.S.os1<-apply(ptb.S.os1,2,mean,na.rm=TRUE)
  dimS<-length(S.os1)
  Sigma.S1<-matrix(0,dimS,dimS)
  for(jj in 1:dimS){
    for(kk in jj:dimS){
      Sigma.S1[jj,kk]<-Sigma.S1[kk,jj]<-
        mean( (ptb.S.os1[,kk]-mu.S.os1[kk])*(ptb.S.os1[,jj]-mu.S.os1[jj]),na.rm=TRUE )
    }
  }
  Tstat1<- t(S.os1)%*%MASS::ginv(Sigma.S1)%*%(S.os1)
  Tstat1

  Icomb1 =  (1-pnorm(Tstat1-qchisq(1-gamma1,df=n.par),0,gamma2))
  est[paste("elas.",1:n.par,sep="" )]<-  (1-Icomb1)*est[paste("ee.rt(ml).",1:n.par,sep="" )]+Icomb1*est[paste("opt.ee(ml).",1:n.par,sep="" )]


  est0<-est
  A.t0 <-A.t
  X.t0 <-X.t
  Y.t0 <-Y.t
  A.os0 <-A.os
  X.os0 <-X.os
  Y.os0 <-Y.os

  Tstat1
  estboot<-matrix(NA,nrow=nboots,ncol=n.par )
  colnames(estboot)<-paste("elas.",1:n.par,sep="" )
  if(nboots>1){
    for(kk in 1:nboots){
      indices.os<-sample(m,m,replace = TRUE)
      indices.t<-sample(n.t,n.t,replace = TRUE)
      A.t <-A.t0[indices.t]
      X.t <-X.t0[indices.t,]
      Y.t <-Y.t0[indices.t]
      A.os <-A.os0[indices.os]
      X.os <-X.os0[indices.os,]
      Y.os <-Y.os0[indices.os]
      kktemp <-estforboot(A.t,X.t,Y.t,A.os,X.os,Y.os,n.t,m,cnames,family.oc,nptb,gamma1,gamma2)
      estboot[kk,]<-kktemp$est
    }

    bootq1<-apply(estboot,2,quantile,probs=alpha.CI/2,  type=5,na.rm = TRUE)
    bootq2<-apply(estboot,2,quantile,probs=1-alpha.CI/2,type=5,na.rm = TRUE)
  }

  if(nboots<2){
    bootq1<-apply(estboot,2,quantile,probs=alpha.CI/2,  type=5,na.rm = TRUE)
    bootq2<-apply(estboot,2,quantile,probs=1-alpha.CI/2,type=5,na.rm = TRUE)
  }

  ve[c(paste("elas.",1:n.par,sep="" )) ]<- apply(estboot,2,var,na.rm = TRUE)



  boot.CI<-cbind(bootq1,bootq2)
  return(list(est=est,ve=ve,boot.CI=boot.CI))


}


####################################################################################
## function list
####################################################################################

#' expoit
#'
#' A function calculate e(x)/{1+e(x)}.
#' @param x a numeric
#' @return ex/(1+ex)
#' @export

expoit <- function(x) {return (exp(x)/(1+exp(x)))}

#' logit
#'
#' A function calculate log{x/(1-x)}.
#' @param x a numeric
#' @return log{x/(1-x)}
#' @export

logit <- function(x) {return (log(x/(1-x)))}

#' estforboot
#'
#' A function computes the bootstrap statistics.
#' @param A.t is a vector of treatment  ( n x 1).
#' @param X.t is a matrix of covariates without intercept  ( n x p).
#' @param Y.t is a vector of outcome  ( n x 1).
#' @param A.os is a vector of treatment  ( m x 1).
#' @param X.os is a matrix of covariates without intercept  ( m x p).
#' @param Y.os is a vector of outcome  ( m x 1).
#' @param n.t is the RT sample size.
#' @param m is the RW sample size.
#' @param cnames is the names for the estimators.
#' @param family.oc specifies the family for the outcome model.
#' \code{"gaussian"}: a linear regression model for the continuous outcome.
#' \code{"binomial"}: a logistic regression model for the binary outcome.
#' @param nptb is the number for replication-based method for estimating Sigma_SS.
#' @param gamma1 is the parameter in the violation test such that T > chisq_{1-gamma1} indicates the bias of the RW analysis is significant.
#' @param gamma2 is a smooth parameter in using the normal cummulative distribution with mean 0 and variance gamma2^2 function to approximate the indicator function.
#' @return the boostrap statistics
#' @import utils MASS rootSolve stats
#' @export


estforboot<-function(A.t,X.t,Y.t,A.os,X.os,Y.os,n.t,m,cnames,family.oc,nptb,gamma1,gamma2){
  est <-NULL

  psi0<-rep(1,dim(X.t)[2]+1)
  n.par<-length(psi0)

  dat.t  <- list(Y.t = Y.t, A.t = A.t, X.t = X.t, Y = Y.t,
                 A = A.t, X = X.t, q = rep(1,n.t), ml.ps = rep(0.5,n.t))
  dat.os <- list(Y.os = Y.os, A.os = A.os, X.os = X.os , Y = Y.os,
                 A = A.os, X = X.os, q = rep(1,m))

  ml.X.t<-cbind(X.t,X.t^2,apply(X.t, 1, utils::combn, 2, prod))
  ml.X.os<-cbind(X.os,X.os^2,apply(X.os, 1, utils::combn, 2, prod))
  dat.t$ml.X.t<-ml.X.t
  dat.os$ml.X.os<-ml.X.os

  glm.out <- glm(A.t~ml.X.t,family=quasibinomial,weights=q,data=dat.t)
  ml.ps.t <- glm.out$fitted.values
  dat.t$ml.ps <- ml.ps.t

  glm.out <- glm(A.os~ml.X.os,family=quasibinomial,weights=q,data=dat.os)
  ml.ps.os <- glm.out$fitted.values
  dat.os$ml.ps <- ml.ps.os

  if(family.oc=="gaussian"){
    mu0.out <- glm(Y.t[which(A.t==0)]~ml.X.t[which(A.t==0),],weights=q[which(A.t==0)],data=dat.t)
    ml.mu0.t <- cbind(1,ml.X.t)%*%mu0.out$coeff
    dat.t$ml.mu0 <- ml.mu0.t

    mu0.out <- glm(Y.os[which(A.os==0)]~ml.X.os[which(A.os==0),],weights=q[which(A.os==0)],data=dat.os)
    ml.mu0.os <- cbind(1,ml.X.os)%*%mu0.out$coeff
    dat.os$ml.mu0 <- ml.mu0.os

    # covariate adjustment (trial)
    X.tadj<- cbind(1,X.t)*(2*A.t-1)/2
    jj<-glm(Y.t~X.tadj-1)
    reg.t<-jj$coeff/2
    reg.t0<-reg.t
  }

  if(family.oc=="binomial"){
    mu0.out <- glm(Y.t[which(A.t==0)]~ml.X.t[which(A.t==0),],weights=q[which(A.t==0)],data=dat.t,family="quasibinomial")
    lml.mu0.t <- cbind(1,ml.X.t)%*%mu0.out$coeff
    ml.mu0.t <- exp(lml.mu0.t)/( 1+exp(lml.mu0.t) )
    dat.t$ml.mu0 <- ml.mu0.t

    mu0.out <- glm(Y.os[which(A.os==0)]~ml.X.os[which(A.os==0),],weights=q[which(A.os==0)],data=dat.os,family="quasibinomial")
    lml.mu0.os <- cbind(1,ml.X.os)%*%mu0.out$coeff
    ml.mu0.os <- exp(lml.mu0.os)/(1+exp(lml.mu0.os))
    dat.os$ml.mu0 <- ml.mu0.os

    # covariate adjustment (trial)
    X.tadj<- cbind(1,X.t)*(2*A.t-1)/2
    jj<-glm(Y.t~X.tadj-1,family="quasibinomial")
    reg.t<-jj$coeff/2
    reg.t0<-reg.t
  }

  dat.integ<-list( Y = c(Y.os,Y.t),A=c(A.os,A.t),X=rbind(X.os,X.t),
                   q =c(rep(1,m),rep(1,n.t)),
                   ml.mu0=c(ml.mu0.os,ml.mu0.t),ml.ps=c(ml.ps.os,ml.ps.t),
                   ml.X=rbind(ml.X.os,ml.X.t))

  if(family.oc=="gaussian"){
    # ee (trial) ml
    jj<-rootSolve::multiroot(f = ee1.ml, start = psi0, dat = dat.t , par0=reg.t0)
    est.ee.rt<-opt.t<-jj$root

    # ee (integration) ml
    jj<-rootSolve::multiroot(f = ee1.ml, start = psi0, dat = dat.integ , par0=reg.t0)
    est.opt.ee<-jj$root

    S.os1<-ee1.ml(par=opt.t,dat=dat.os, par0=reg.t0)
    S.os1
  }
  if(family.oc=="binomial"){
    # ee (trial) ml
    jj<-rootSolve::multiroot(f = ee2.ml, start = psi0, dat = dat.t , par0=reg.t0)
    est.ee.rt<-opt.t<-jj$root

    # ee (integration) ml
    jj<-rootSolve::multiroot(f = ee2.ml, start = psi0, dat = dat.integ , par0=reg.t0)
    est.opt.ee<-jj$root

    S.os1<-ee2.ml(par=opt.t,dat=dat.os, par0=reg.t0)
    S.os1
  }

  ## trying a test before combining RCT and OS
  ## perturbation based variance estimation

  ptb.S.os1<-matrix(0,nptb,length(S.os1))
  cnames<-c(paste("elas.",1:n.par,sep="" ))
  ptb<-matrix(0,nptb,length(cnames))
  colnames(ptb)<-cnames

  for(kkk in 1:nptb){
    dat.t$q <-rexp(n.t,1)
    dat.os$q<-rexp(m,1)

    glm.out <- glm(A.os~ml.X.os,family=quasibinomial,weights=q,data=dat.os)
    ml.ps.os <- glm.out$fitted.values
    dat.os$ml.ps <- ml.ps.os

    if(family.oc=="gaussian"){
      mu0.out <- glm(Y.t[which(A.t==0)]~ml.X.t[which(A.t==0),],weights=q[which(A.t==0)],data=dat.t)
      ml.mu0.t <- cbind(1,ml.X.t)%*%mu0.out$coeff
      dat.t$ml.mu0 <- ml.mu0.t

      mu0.out <- glm(Y.os[which(A.os==0)]~ml.X.os[which(A.os==0),],weights=q[which(A.os==0)],data=dat.os)
      ml.mu0.os <- cbind(1,ml.X.os)%*%mu0.out$coeff
      dat.os$ml.mu0 <- ml.mu0.os
    }

    if(family.oc=="binomial"){
      mu0.out <- glm(Y.t[which(A.t==0)]~ml.X.t[which(A.t==0),],weights=q[which(A.t==0)],data=dat.t,family="quasibinomial")
      lml.mu0.t <- cbind(1,ml.X.t)%*%mu0.out$coeff
      ml.mu0.t <- exp(lml.mu0.t)/( 1+exp(lml.mu0.t) )
      dat.t$ml.mu0 <- ml.mu0.t

      mu0.out <- glm(Y.os[which(A.os==0)]~ml.X.os[which(A.os==0),],weights=q[which(A.os==0)],data=dat.os,family="quasibinomial")
      lml.mu0.os <- cbind(1,ml.X.os)%*%mu0.out$coeff
      ml.mu0.os <- exp(lml.mu0.os)/(1+exp(lml.mu0.os))
      dat.os$ml.mu0 <- ml.mu0.os
    }

    dat.integ<-list( Y = c(Y.os,Y.t),A=c(A.os,A.t),X=rbind(X.os,X.t),
                     q =c(dat.os$q,dat.t$q),
                     ml.mu0=c(ml.mu0.os,ml.mu0.t),ml.ps=c(ml.ps.os,ml.ps.t),
                     ml.X=rbind(ml.X.os,ml.X.t))

    if(family.oc=="gaussian"){
      jj<-rootSolve::multiroot(f = ee1.ml, start = psi0, dat = dat.t , par0=reg.t0)
      opt.t<-jj$root

      ptb.S.os1[kkk,]<-ee1.ml(par=opt.t,dat=dat.os, par0=reg.t0)
    }
    if(family.oc=="binomial"){
      jj<-rootSolve::multiroot(f = ee2.ml, start = psi0, dat = dat.t , par0=reg.t0)
      opt.t<-jj$root

      ptb.S.os1[kkk,]<-ee2.ml(par=opt.t,dat=dat.os, par0=reg.t0)
    }

  }

  ve<-apply(ptb,2,var)
  mu.S.os1<-apply(ptb.S.os1,2,mean,na.rm=TRUE)
  dimS<-length(S.os1)
  Sigma.S1<-matrix(0,dimS,dimS)
  for(jj in 1:dimS){
    for(kk in jj:dimS){
      Sigma.S1[jj,kk]<-Sigma.S1[kk,jj]<-
        mean( (ptb.S.os1[,kk]-mu.S.os1[kk])*(ptb.S.os1[,jj]-mu.S.os1[jj]),na.rm=TRUE )
    }
  }
  Tstat1<- t(S.os1)%*%MASS::ginv(Sigma.S1)%*%(S.os1)
  Tstat1

  Icomb1<- (1-pnorm(Tstat1-qchisq(1-gamma1,df=n.par),0,gamma2))
  est[paste("elas.",1:n.par,sep="" )]<-  (1-Icomb1)*est.ee.rt+Icomb1*est.opt.ee

  return( list(est=est) )
}

#' ee1.ml
#'
#' An estimating function for the HTE with a continuous outcome.
#' @param par parameter.
#' @param dat data list.
#' @param par0 required parameters.
#' @return estimating function values.
#' @export

ee1.ml<-function(par,dat,par0){
  ## opt ee for the continuous outcome
  psi<-par
  Y<-dat$Y
  A<-dat$A
  X<-dat$X
  q<-dat$q
  ps<-dat$ml.ps
  mu0<-dat$ml.mu0
  H<-Y-A*cbind(1,X)%*%psi
  apply(cbind(1,X)* matrix( (H-mu0)*(A-ps)*q,length(Y),1+dim(X)[2], byrow=FALSE),2,mean)
}

#' ee2.ml
#'
#' An estimating function for the HTE with a binary outcome.
#' @param par parameter
#' @param dat data list
#' @param par0 required parameters
#' @return estimating function values
#' @export

ee2.ml<-function(par,dat,par0){
  psi<-par
  Y<-dat$Y
  A<-dat$A
  X<-dat$X
  q<-dat$q
  ps<-dat$ml.ps
  mu0<-dat$ml.mu0
  le<-cbind(1,X)%*%psi
  le0<-cbind(1,X)%*%par0
  tau<- (exp(le)-1 )/(exp(le)+1 )
  tau0<- (exp(le0)-1 )/(exp(le0)+1 )
  H<-Y-A*tau
  partialtau<-2*exp(le0)/(( exp(le0)+1 )^2)
  apply(cbind(1,X)* matrix( partialtau*( mu0*(1-mu0) )^(-1)*(H-mu0)*(A-ps)*q,length(Y),1+dim(X)[2], byrow=FALSE),2,sum)
}
