
palm_intensity<-function(u,lambda,mu,sig2){
  return(lambda*(1+1/(4*pi*sig2*mu)*exp(-u^2/(4*sig2))))
}

palm_likelihood<-function(data,params){
  lambda<-params$lambda
  log_mu<-params$log_mu
  log_sig2<-params$log_sig2
  
  R<-data$rmax/2
  
  return(2*sum(log(palm_intensity(data$dist_pairs,lambda,exp(log_mu),exp(log_sig2)))) - 
           sum(palm_intensity(data$dist_cross,lambda,exp(log_mu),exp(log_sig2)))/data$nx^2)
}


update_lambda<-function(data,params,priors){
  params_star<-params
  
  params_star$lambda<-rtruncnorm(n=1,mean=params$lambda,sd=sqrt(priors$lambda.pvar),a=0)
  
  R <- params$eta*palm_likelihood(data,params_star) - params$eta*palm_likelihood(data,params) +
    log(dtruncnorm(params_star$lambda,mean=nrow(data$points),sd=sqrt(priors$lambda.var),a=0)) +
    log(dtruncnorm(params$lambda,mean=params_star$lambda,sd=sqrt(priors$lambda.pvar),a=0)) -
    log(dtruncnorm(params$lambda,mean=nrow(data$points),sd=sqrt(priors$lambda.var),a=0)) -
    log(dtruncnorm(params_star$lambda,mean=params$lambda,sd=sqrt(priors$lambda.pvar),a=0))
  
  if(runif(1)<exp(R)){params<-params_star}
  
  return(params)
}

update_log_mu<-function(data,params,priors){
  params_star<-params
  
  params_star$log_mu<-rnorm(n=1,mean=params$log_mu,sd=sqrt(priors$log_mu.pvar))
  
  R <- params$eta*palm_likelihood(data,params_star) - params$eta*palm_likelihood(data,params) +
    dnorm(params_star$log_mu,mean=0,sd=sqrt(priors$log_mu.var),log=TRUE) -
    dnorm(params$log_mu,mean=0,sd=sqrt(priors$log_mu.var),log=TRUE)
  
  if(runif(1)<exp(R)){params<-params_star}
  
  return(params)
}

update_log_sig2<-function(data,params,priors){
  params_star<-params
  
  params_star$log_sig2<-rnorm(n=1,mean=params$log_sig2,sd=sqrt(priors$log_sig2.pvar))
  
  R <- params$eta*palm_likelihood(data,params_star) - params$eta*palm_likelihood(data,params) +
    dnorm(params_star$log_sig2,mean=-5,sd=sqrt(priors$log_sig2.var),log=TRUE) -
    dnorm(params$log_sig2,mean=-5,sd=sqrt(priors$log_sig2.var),log=TRUE)
  
  if(runif(1)<exp(R)){params<-params_star}
  
  return(params)
}

mcmc_palm<-function(data,params,priors,iters,burn,thin){
  out<-list()
  out$lambda<-NULL
  out$log_mu<-NULL
  out$log_sig2<-NULL
  out$lambda_full<-NULL
  out$log_mu_full<-NULL
  out$log_sig2_full<-NULL
  
  t<-Sys.time()
  for(k in 1:iters){
    params<-update_lambda(data,params,priors)
    params<-update_log_mu(data,params,priors)
    params<-update_log_sig2(data,params,priors)
    
    out$lambda_full[k]<-params$lambda
    out$log_mu_full[k]<-params$log_mu
    out$log_sig2_full[k]<-params$log_sig2
    
    if(k%%100==0){
      lambda_acc<-length(unique(out$lambda_full[(k-99):k]))/100
      mu_acc<-length(unique(out$log_mu_full[(k-99):k]))/100
      sig2_acc<-length(unique(out$log_sig2_full[(k-99):k]))/100
      
      if(lambda_acc>0.4){priors$lambda.pvar<-priors$lambda.pvar*1.2}
      if(lambda_acc<0.2){priors$lambda.pvar<-priors$lambda.pvar*0.8}
      
      if(mu_acc>0.4){priors$log_mu.pvar<-priors$log_mu.pvar*1.2}
      if(mu_acc<0.2){priors$log_mu.pvar<-priors$log_mu.pvar*0.8}
      
      if(mu_acc>0.4){priors$log_sig2.pvar<-priors$log_sig2.pvar*1.2}
      if(mu_acc<0.2){priors$log_sig2.pvar<-priors$log_sig2.pvar*0.8}
    }
  }

  out$lambda<-out$lambda_full[seq(burn+1,iters,thin)]
  out$log_mu<-out$log_mu_full[seq(burn+1,iters,thin)]
  out$log_sig2<-out$log_sig2_full[seq(burn+1,iters,thin)]
  out$time<-Sys.time()-t
  return(out)
}



d_log_palm<-function(u,lambda,log_mu,log_sig2){
  return(cbind(1/lambda,
               (-1/(4*pi*exp(log_sig2)*exp(log_mu))*exp(-u^2/(4*exp(log_sig2))))/
                 (1+1/(4*pi*exp(log_sig2)*exp(log_mu))*exp(-u^2/(4*exp(log_sig2)))),
               ((u^2/(4*exp(log_sig2))-1)/(4*pi*exp(log_sig2)*exp(log_mu))*exp(-u^2/(4*exp(log_sig2))))/
                 ((1+1/(4*pi*exp(log_sig2)*exp(log_mu))*exp(-u^2/(4*exp(log_sig2)))))))
}

d_palm<-function(u,lambda,log_mu,log_sig2){
  d1<-1+1/(4*pi*exp(log_sig2))
  d2<-lambda*(-1/(4*pi*exp(log_sig2)*exp(log_mu))*exp(-u^2/(4*exp(log_sig2))))
  d3<-lambda*((u^2/(4*exp(log_sig2))-1)/(4*pi*exp(log_sig2)*exp(log_mu))*exp(-u^2/(4*exp(log_sig2))))
  
  return(cbind(d1,d2,d3))
}


palm_score<-function(data,params){
  return(2*apply(d_log_palm(data$dist_pairs,params$lambda,params$log_mu,params$log_sig2),2,sum) - 
           apply(d_palm(data$dist_cross,params$lambda,params$log_mu,params$log_sig2),2,sum)/data$nx^2)
}


