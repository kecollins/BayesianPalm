palm_intensity<-function(u,beta,sig2,phi,X){
  return(exp(X%*%beta+sig2/2+sig2*exp(-u/phi)))
}

palm_likelihood<-function(data,params){
  beta<-params$beta
  sig2<-exp(params$log_sig2)
  phi<-exp(params$log_phi)
  
  R<-data$rmax/2
  
  return(sum(log(palm_intensity(data$dist_pairs,beta,sig2,phi,data$X1))) - 
           sum(palm_intensity(data$dist_cross,beta,sig2,phi,data$X2))*25)
}

update_beta0<-function(data,params,priors){
  params_star<-params
  params_star$beta[1]<-rnorm(n=1,mean=params$beta[1],sd=sqrt(priors$beta0.pvar))
  
  R<-params$eta*palm_likelihood(data,params_star)-
    params$eta*palm_likelihood(data,params)+
    dnorm(params_star$beta[1],sd=1000,log=TRUE)-
    dnorm(params$beta[1],sd=1000,log=TRUE)
  
  if(runif(1)<exp(R)){params<-params_star}
  return(params)
}

update_beta1<-function(data,params,priors){
  params_star<-params
  params_star$beta[2]<-rnorm(n=1,mean=params$beta[2],sd=sqrt(priors$beta1.pvar))
  
  R<-params$eta*palm_likelihood(data,params_star)-
    params$eta*palm_likelihood(data,params)+
    dnorm(params_star$beta[2],sd=1000,log=TRUE)-
    dnorm(params$beta[2],sd=1000,log=TRUE)
  
  if(runif(1)<exp(R)){params<-params_star}
  return(params)
}

update_beta2<-function(data,params,priors){
  params_star<-params
  params_star$beta[3]<-rnorm(n=1,mean=params$beta[3],sd=sqrt(priors$beta2.pvar))
  
  R<-params$eta*palm_likelihood(data,params_star)-
    params$eta*palm_likelihood(data,params)+
    dnorm(params_star$beta[3],sd=1000,log=TRUE)-
    dnorm(params$beta[3],sd=1000,log=TRUE)
  
  if(runif(1)<exp(R)){params<-params_star}
  return(params)
}

update_log_sig2<-function(data,params,priors){
  params_star<-params
  params_star$log_sig2<-rnorm(n=1,mean=params$log_sig2,sd=sqrt(priors$sig2.pvar))
  
  R<-params$eta*palm_likelihood(data,params_star)-
    params$eta*palm_likelihood(data,params)+
    dnorm(params_star$log_sig2,sd=sqrt(priors$sig2.var),log=TRUE)-
    dnorm(params$log_sig2,sd=sqrt(priors$sig2.var),log=TRUE)
  
  if(runif(1)<exp(R)){params<-params_star}
  
  return(params)
}

update_log_phi<-function(data,params,priors){
  params_star<-params
  params_star$log_phi<-rtruncnorm(n=1,mean=params$log_phi,sd=sqrt(priors$phi.pvar),
                                  a=priors$phi.lower,b=priors$phi.upper)
  
  R<-params$eta*palm_likelihood(data,params_star)-
    params$eta*palm_likelihood(data,params)+
    #log(dtruncnorm(params_star$log_phi,sd=sqrt(priors$phi.var),a=priors$phi.lower,b=priors$phi.upper))-
    #log(dtruncnorm(params$log_phi,sd=sqrt(priors$phi.var),a=priors$phi.lower,b=priors$phi.upper))+
    log(dtruncnorm(params$log_phi,mean=params_star$log_phi,sd=sqrt(priors$phi.pvar),a=priors$phi.lower,b=priors$phi.upper))-
    log(dtruncnorm(params_star$log_phi,mean=params$log_phi,sd=sqrt(priors$phi.pvar),a=priors$phi.lower,b=priors$phi.upper))
  
  if(runif(1)<exp(R)){params<-params_star}
  
  return(params)
}

mcmc_palm<-function(data,params,priors,burn,iters,thin){
  out<-list()
  out$beta0_full<-NULL
  out$beta1_full<-NULL
  out$beta2_full<-NULL
  out$log_sig2_full<-NULL
  out$log_phi_full<-NULL
  
  j<-1
  t<-Sys.time()
  for(k in 1:iters){
    params<-update_beta0(data,params,priors)
    params<-update_beta1(data,params,priors)
    params<-update_beta2(data,params,priors)
    params<-update_log_sig2(data,params,priors)
    params<-update_log_phi(data,params,priors)
    
    out$beta0_full[k]<-params$beta[1]
    out$beta1_full[k]<-params$beta[2]
    out$beta2_full[k]<-params$beta[3]
    out$log_sig2_full[k]<-params$log_sig2
    out$log_phi_full[k]<-params$log_phi
    
    
    if(k%%100==0){
      beta0_acc<-length(unique(out$beta0_full[(k-99):k]))/100
      beta1_acc<-length(unique(out$beta1_full[(k-99):k]))/100
      beta2_acc<-length(unique(out$beta2_full[(k-99):k]))/100
      sig2_acc<-length(unique(out$log_sig2_full[(k-99):k]))/100
      phi_acc<-length(unique(out$log_phi_full[(k-99):k]))/100
      
      if(beta0_acc>0.4){priors$beta0.pvar<-priors$beta0.pvar*1.2}
      if(beta0_acc<0.2){priors$beta0.pvar<-priors$beta0.pvar*0.8}
      
      if(beta1_acc>0.4){priors$beta1.pvar<-priors$beta1.pvar*1.2}
      if(beta1_acc<0.2){priors$beta1.pvar<-priors$beta1.pvar*0.8}
      
      if(beta2_acc>0.4){priors$beta2.pvar<-priors$beta2.pvar*1.2}
      if(beta2_acc<0.2){priors$beta2.pvar<-priors$beta2.pvar*0.8}
      
      if(sig2_acc>0.4){priors$sig2.pvar<-priors$sig2.pvar*1.2}
      if(sig2_acc<0.2){priors$sig2.pvar<-priors$sig2.pvar*0.8}
      
      if(phi_acc>0.4){priors$phi.pvar<-priors$phi.pvar*1.2}
      if(phi_acc<0.2){priors$phi.pvar<-priors$phi.pvar*0.8}
    }
  }
  
  out$beta0<-out$beta0_full[seq(burn+1,iters,thin)]
  out$beta1<-out$beta1_full[seq(burn+1,iters,thin)]
  out$beta2<-out$beta2_full[seq(burn+1,iters,thin)]
  out$log_sig2<-out$log_sig2_full[seq(burn+1,iters,thin)]
  out$log_phi<-out$log_phi_full[seq(burn+1,iters,thin)]
  
  out$time<-Sys.time()-t
  return(out)
}


