library(mvtnorm)
library(truncnorm)
library(fields)

## Likelihood for top level
lgcp_likelihood<-function(data,params){
  lambda<-exp(data$X%*%params$beta+params$y)
  ll<-sum(data$ni*log(lambda)-(1/data$nx)^2*lambda)
  return(ll)
}

## Elliptical slice sampling
update_y<-function(data,params){
  Sig<-exp(params$log_sig2)*Matern(data$dist,range=exp(params$log_phi),nu=params$nu)
  v<-rmvnorm(n=1,sigma=Sig)[1,]
  v<-v-mean(v)
  
  thresh<-lgcp_likelihood(data,params)+log(runif(1))
  
  theta<-runif(1,0,2*pi)
  theta.min<-min(theta,theta-2*pi)
  theta.max<-max(theta,theta-2*pi)
  
  params_star<-params
  params_star$y<-params$y*cos(theta)+v*sin(theta)
  
  if(lgcp_likelihood(data,params_star)>thresh){
    params<-params_star
  } else{
    accept<-FALSE
    while(!accept){
      if(theta<0){
        theta.min<-theta
      } else{theta.max<-theta}
      theta<-runif(1,theta.min,theta.max)
      
      params_star$y<-params$y*cos(theta)+v*sin(theta)
      
      if(lgcp_likelihood(data,params_star)>thresh){
        params<-params_star
        accept<-TRUE
      }
    }
  }
  
  return(params)
}

update_log_sig2_full<-function(data,params,priors,likelihood){
  params_star<-params
  params_star$log_sig2 <- rnorm(n=1,mean=params$log_sig2,sd=sqrt(priors$log_sig2.pvar))
  
  R <- likelihood(data,params_star) - likelihood(data,params) +
    dnorm(params_star$log_sig2,mean=priors$log_sig2.mean,sd=sqrt(priors$log_sig2.var),log=TRUE) -
    dnorm(params$log_sig2,mean=priors$log_sig2.mean,sd=sqrt(priors$log_sig2.var),log=TRUE)
  
  
  if(runif(1)<exp(R)){params<-params_star}
  return(params)
}

log_GP_full<-function(data,params){
  Sig<-exp(params$log_sig2)*Matern(data$dist,range=exp(params$log_phi),nu=params$nu)
  return(dmvnorm(params$y,sigma = Sig,log = TRUE))
}

update_log_phi_full<-function(data,params,priors,likelihood){
  params_star<-params
  params_star$log_phi <- rtruncnorm(n=1,mean=params$log_phi,sd=sqrt(priors$log_phi.pvar),a=-3,b=-1.6)
  
  R <- likelihood(data,params_star) - likelihood(data,params) +
    log(dtruncnorm(params_star$log_phi,mean=priors$log_phi.mean,sd=sqrt(priors$log_phi.var),a=-3,b=-1.6))+
    log(dtruncnorm(params$log_phi,mean=params_star$log_phi,sd=sqrt(priors$log_phi.pvar),a=-3,b=-1.6))-
    log(dtruncnorm(params_star$log_phi,mean=priors$log_phi.mean,sd=sqrt(priors$log_phi.var),a=-3,b=-1.6))-
    log(dtruncnorm(params_star$log_phi,mean=params$log_phi,sd=sqrt(priors$log_phi.pvar),a=-3,b=-1.6))
  
  
  if(runif(1)<exp(R)){params<-params_star}
  
  return(params)
}


update_beta_full<-function(data,params,priors,likelihood){
  params_star<-params
  params_star$beta<-rnorm(n=length(params$beta),mean=params$beta,sd=sqrt(priors$beta.pvar))
  
  R <- likelihood(data,params_star) + sum(dnorm(params_star$beta,mean=priors$beta.beta0,sd=sqrt(priors$beta.sig20),log=TRUE))-
    (likelihood(data,params) + sum(dnorm(params$beta,mean=priors$beta.beta0,sd=sqrt(priors$beta.sig20),log=TRUE)))
  
  if(runif(1)<exp(R)){params<-params_star}
  
  return(params)
}


mcmc_lgcp<-function(data,params,priors,iters,burn,thin){
  out<-list()
  out$y_full<-matrix(nrow=length(params$y),ncol=iters)
  out$log_sig2_full<-NULL
  out$log_phi_full<-NULL
  out$beta_full<-matrix(ncol=length(params$beta),nrow=iters)
  out$time<-Sys.time()
  
  for(k in 1:iters){
    params<-update_y(data,params)
    params<-update_log_sig2_full(data,params,priors,log_GP_full)
    params<-update_log_phi_full(data,params,priors,log_GP_full)
    params<-update_beta_full(data,params,priors,lgcp_likelihood)
    
    
    out$y_full[,k]<-params$y
    out$log_sig2_full[k]<-params$log_sig2
    out$log_phi_full[k]<-params$log_phi
    out$beta_full[k,]<-params$beta
    
    ## adaptive proposals
    if(k%%100==0){
      beta_acc<-length(unique(out$beta_full[(k-99):k,]))/100
      sig2_acc<-length(unique(out$log_sig2_full[(k-99):k]))/100
      phi_acc<-length(unique(out$log_phi_full[(k-99):k]))/100
      
      if(beta_acc>0.4){priors$beta.pvar<-priors$beta.pvar*1.2}
      if(beta_acc<0.2){priors$beta.pvar<-priors$beta.pvar*0.8}
      
      if(sig2_acc>0.4){priors$log_sig2.pvar<-priors$log_sig2.pvar*1.2}
      if(sig2_acc<0.2){priors$log_sig2.pvar<-priors$log_sig2.pvar*0.8}
      
      if(phi_acc>0.4){priors$log_phi.pvar<-priors$log_phi.pvar*1.2}
      if(phi_acc<0.2){priors$log_phi.pvar<-priors$log_phi.pvar*0.8}
    }
    
  }
  
  out$y<-out$y_full[,seq(burn+1,iters,thin)]
  out$log_sig2<-out$log_sig2_full[seq(burn+1,iters,thin)]
  out$log_phi<-out$log_phi_full[seq(burn+1,iters,thin)]
  out$beta<-out$beta_full[seq(burn+1,iters,thin),]
  
  
  out$time<-Sys.time()-out$time
  return(out)
}


