library(truncnorm)

palm_intensity<-function(u,lambda,log_sig2,log_phi){
  return(lambda*exp(exp(log_sig2)*exp(-u/exp(log_phi))))
}

palm_likelihood<-function(data,params){
  lambda<-params$lambda
  log_sig2<-params$log_sig2
  log_phi<-params$log_phi
  
  R<-data$rmax/2
  
  return(2*sum(log(palm_intensity(data$dist_pairs,lambda,log_sig2,log_phi))) - 
           sum(palm_intensity(data$dist_cross,lambda,log_sig2,log_phi))/data$nx^2)
}

update_lambda<-function(data,params,priors){
  params_star<-params
  
  params_star$lambda<-rtruncnorm(n=1,mean=params$lambda,sd=sqrt(priors$lambda.pvar),a=0)
  
  R<-params$eta*palm_likelihood(data,params_star) -params$eta*palm_likelihood(data,params) +
    log(dtruncnorm(params_star$lambda,mean=nrow(data$points),sd=sqrt(priors$lambda.var),a=0)) +
    log(dtruncnorm(params$lambda,mean=params_star$lambda,sd=sqrt(priors$lambda.pvar),a=0)) -
    log(dtruncnorm(params$lambda,mean=nrow(data$points),sd=sqrt(priors$lambda.var),a=0)) -
    log(dtruncnorm(params_star$lambda,mean=params$lambda,sd=sqrt(priors$lambda.pvar),a=0))
  
  if(!is.na(R)){if(runif(1)<exp(R)){params<-params_star}}
  
  return(params)
}

update_log_sig2<-function(data,params,priors){
  params_star<-params
  
  params_star$log_sig2<-rnorm(n=1,mean=params$log_sig2,sd=sqrt(priors$log_sig2.pvar))
  
  R<-params$eta*palm_likelihood(data,params_star) - params$eta*palm_likelihood(data,params) +
    dnorm(params_star$log_sig2,mean=0,sd=sqrt(priors$log_sig2.var),log=TRUE) -
    dnorm(params$log_sig2,mean=0,sd=sqrt(priors$log_sig2.var),log=TRUE)
  
  if(!is.na(R)){if(runif(1)<exp(R)){params<-params_star}}
  return(params)
}

update_log_phi<-function(data,params,priors){
  params_star<-params
  
  params_star$log_phi<-rtruncnorm(n=1,mean=params$log_phi,sd=sqrt(priors$log_phi.pvar),a=-3,b=-1.6)
  
  R<-params$eta*palm_likelihood(data,params_star) - params$eta*palm_likelihood(data,params) +
    log(dtruncnorm(params$log_phi,mean=params_star$log_phi,sd=sqrt(priors$log_phi.pvar),a=-3,b=-1.6)) -
    log(dtruncnorm(params_star$log_phi,mean=params$log_phi,sd=sqrt(priors$log_phi.pvar),a=-3,b=-1.6))
  
  if(!is.na(R)){if(runif(1)<exp(R)){params<-params_star}}
  return(params)
}

mcmc_palm<-function(data,params,priors,burn,iters,thin){
  out<-list()
  out$lambda<-NULL
  out$log_sig2<-NULL
  out$log_phi<-NULL
  out$lambda_full<-NULL
  out$log_sig2_full<-NULL
  out$log_phi_full<-NULL

  t<-Sys.time()
  for(k in 1:iters){
    params<-update_lambda(data,params,priors)
    params<-update_log_sig2(data,params,priors)
    params<-update_log_phi(data,params,priors)
    
    out$lambda_full[k]<-params$lambda
    out$log_sig2_full[k]<-params$log_sig2
    out$log_phi_full[k]<-params$log_phi
    
    if(k%%100==0){
      lambda_acc<-length(unique(out$lambda_full[(k-99):k]))/100
      sig2_acc<-length(unique(out$log_sig2_full[(k-99):k]))/100
      phi_acc<-length(unique(out$log_phi_full[(k-99):k]))/100
      
      if(lambda_acc>0.4){priors$lambda.pvar<-priors$lambda.pvar*1.2}
      if(lambda_acc<0.2){priors$lambda.pvar<-priors$lambda.pvar*0.8}
      
      if(sig2_acc>0.4){priors$log_sig2.pvar<-priors$log_sig2.pvar*1.2}
      if(sig2_acc<0.2){priors$log_sig2.pvar<-priors$log_sig2.pvar*0.8}
      
      if(phi_acc>0.4){priors$log_phi.pvar<-priors$log_phi.pvar*1.2}
      if(phi_acc<0.2){priors$log_phi.pvar<-priors$log_phi.pvar*0.8}
    }
  }

  out$lambda<-out$lambda_full[seq(burn+1,iters,thin)]
  out$log_sig2<-out$log_sig2_full[seq(burn+1,iters,thin)]
  out$log_phi<-out$log_phi_full[seq(burn+1,iters,thin)]
  out$time<-Sys.time()-t
  return(out)
}


d_log_palm<-function(u,lambda,log_sig2,log_phi){
  return(cbind(1/lambda,
               exp(log_sig2)*exp(-u/exp(log_phi)),
               exp(log_sig2)*u/exp(log_phi)*exp(-u/exp(log_phi))))
}

d_palm<-function(u,lambda,log_sig2,log_phi){
  d1<-exp(exp(log_sig2)*exp(-u/exp(log_phi)))
  d2<-exp(log_sig2)*exp(-u/exp(log_phi))*lambda*exp(exp(log_sig2)*exp(-u/exp(log_phi)))
  d3<-exp(log_sig2)*u/exp(log_phi)*exp(-u/exp(log_phi))*lambda*exp(exp(log_sig2)*exp(-u/exp(log_phi)))
  
  return(cbind(d1,d2,d3))
}


palm_score<-function(data,params){
  return(2*apply(d_log_palm(data$dist_pairs,params$lambda,params$log_sig2,params$log_phi),2,sum) - 
           apply(d_palm(data$dist_cross,params$lambda,params$log_sig2,params$log_phi),2,sum)/data$nx^2)
}



