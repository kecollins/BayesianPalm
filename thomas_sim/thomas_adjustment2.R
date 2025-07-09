source('thomas_functions.R')
library(spatstat)
library(doParallel)
cl<-makeCluster(2)
registerDoParallel(cl)
registerDoParallel(cores=10)


pp_gen<-readRDS('pp_gen_thomas.rds')

palm_thomas0.4<-readRDS('palm_thomas0.4.rds')
palm_empirical_thomas0.4<-readRDS('palm_empirical_thomas0.4.rds')

nx<-20
grid<-expand.grid(seq((1/nx)/2,1,1/nx),seq((1/nx)/2,1,1/nx))


mu<-30
nu<-10
sig2<-0.05^2
rmax<-0.4

n_sims<-100

adj_out<-foreach(k=c(1:n_sims),.errorhandling = 'pass')%do%{
  est<-list()
  est$lambda<-mean(palm_thomas0.4[[k]]$lambda)
  est$log_mu<-mean(palm_thomas0.4[[k]]$log_mu)
  est$log_sig2<-mean(palm_thomas0.4[[k]]$log_sig2)
  est$nu<-mean(palm_thomas0.4[[k]]$lambda/exp(palm_thomas0.4[[k]]$log_mu))
  
  test_post<-foreach(k=c(1:100),.errorhandling = 'remove')%dopar%{
    pp<-rThomas(kappa=exp(est$log_mu),scale=sqrt(exp(est$log_sig2)),mu=est$nu)
    
    data<-list(points=cbind(pp$x,pp$y),nx=nx,rmax=rmax)
    dist_mat<-as.matrix(dist(data$points))
    data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
    data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
    data$dist_cross<-as.matrix(rdist(data$points,grid))
    data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
    
    params<-list(lambda=nrow(data$points),log_mu=log(mu),log_sig2=log(sig2),eta=1)
    
    priors<-list(lambda.var=1000,log_mu.var=10,log_sig2.var=10,
                 lambda.pvar=10,log_mu.pvar=0.1,log_sig2.pvar=0.1)
    
    burn<-2000
    iters<-20000
    thin<-18
    
    post_palm<-mcmc_palm(data,params,priors,iters,burn,thin)
    
    return(post_palm)
  }
  
  nu_unadj<-lapply(test_post,function(x) log(x$lambda/exp(x$log_mu)))
  nu_adj<-nu_unadj
  stop<-FALSE
  eta<-0
  while(!stop){
    nu_qs<-lapply(nu_adj,function(x) quantile(x,c(0.025,0.975)))
    cov<-mean(unlist(lapply(nu_qs,function(x) (x[1]<log(est$nu))&(x[2]>log(est$nu)))))
    if(cov>=0.95){stop<-TRUE}
    
    eta<-eta+0.1
    nu_adj<-lapply(nu_unadj,function(x) mean(x)+eta*(x-mean(x)))
  }
  
  
  nu_unadj<-log(palm_thomas0.4[[k]]$lambda/exp(palm_thomas0.4[[k]]$log_mu))
  nu_adj<-mean(nu_unadj)+eta*(nu_unadj-mean(nu_unadj))
  
  
  sig2_adj<-lapply(test_post,function(x) x$log_sig2)
  sig2_unadj<-sig2_adj
  stop<-FALSE
  eta<-0
  while(!stop){
    sig2_qs<-lapply(sig2_adj,function(x) quantile(x,c(0.025,0.975)))
    cov<-mean(unlist(lapply(sig2_qs,function(x) (x[1]<est$log_sig2)&(x[2]>est$log_sig2))))
    if(cov>=0.95){stop<-TRUE}
    
    eta<-eta+0.1
    sig2_adj<-lapply(sig2_unadj,function(x) mean(x)+eta*(x-mean(x)))
  }
  
  
  sig2_unadj<-palm_thomas0.4[[k]]$log_sig2
  sig2_adj<-mean(sig2_unadj)+eta*(sig2_unadj-mean(sig2_unadj))
  
  
  mu_adj<-lapply(test_post,function(x) x$log_mu)
  mu_unadj<-mu_adj
  stop<-FALSE
  eta<-0
  while(!stop){
    mu_qs<-lapply(mu_adj,function(x) quantile(x,c(0.025,0.975)))
    cov<-mean(unlist(lapply(mu_qs,function(x) (x[1]<est$log_mu)&(x[2]>est$log_mu))))
    if(cov>=0.95){stop<-TRUE}
    
    eta<-eta+0.1
    mu_adj<-lapply(mu_unadj,function(x) mean(x)+eta*(x-mean(x)))
  }
  
  
  mu_unadj<-palm_thomas0.4[[k]]$log_mu
  mu_adj<-mean(mu_unadj)+eta*(mu_unadj-mean(mu_unadj))
  
  
  out<-list()
  out$log_nu<-nu_adj
  out$log_sig2<-sig2_adj
  out$log_mu<-mu_adj
  
  return(out)
}
saveRDS(adj_out,paste0('adj_thomas2.rds'))




adj_out<-foreach(k=c(1:n_sims),.errorhandling = 'pass')%do%{
  est<-list()
  est$lambda<-mean(palm_empirical_thomas0.4[[k]]$lambda)
  est$log_mu<-mean(palm_empirical_thomas0.4[[k]]$log_mu)
  est$log_sig2<-mean(palm_empirical_thomas0.4[[k]]$log_sig2)
  est$nu<-mean(palm_empirical_thomas0.4[[k]]$lambda/exp(palm_empirical_thomas0.4[[k]]$log_mu))
  
  test_post<-foreach(k=c(1:100),.errorhandling = 'remove')%dopar%{
    pp<-rThomas(kappa=exp(est$log_mu),scale=sqrt(exp(est$log_sig2)),mu=est$nu)
    
    data<-list(points=cbind(pp$x,pp$y),nx=nx,rmax=rmax)
    dist_mat<-as.matrix(dist(data$points))
    data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
    data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
    data$dist_cross<-as.matrix(rdist(data$points,grid))
    data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
    
    params<-list(lambda=nrow(data$points),log_mu=log(mu),log_sig2=log(sig2),eta=1)
    
    priors<-list(lambda.var=10,log_mu.var=10,log_sig2.var=10,
                 lambda.pvar=10,log_mu.pvar=0.1,log_sig2.pvar=0.1)
    
    burn<-2000
    iters<-20000
    thin<-18
    
    post_palm<-mcmc_palm(data,params,priors,iters,burn,thin)
    
    return(post_palm)
  }
  
  nu_unadj<-lapply(test_post,function(x) log(x$lambda/exp(x$log_mu)))
  nu_adj<-nu_unadj
  stop<-FALSE
  eta<-0
  while(!stop){
    nu_qs<-lapply(nu_adj,function(x) quantile(x,c(0.025,0.975)))
    cov<-mean(unlist(lapply(nu_qs,function(x) (x[1]<log(est$nu))&(x[2]>log(est$nu)))))
    if(cov>=0.95){stop<-TRUE}
    
    eta<-eta+0.1
    nu_adj<-lapply(nu_unadj,function(x) mean(x)+eta*(x-mean(x)))
  }
  
  
  nu_unadj<-log(palm_empirical_thomas0.4[[k]]$lambda/exp(palm_empirical_thomas0.4[[k]]$log_mu))
  nu_adj<-mean(nu_unadj)+eta*(nu_unadj-mean(nu_unadj))
  
  
  sig2_adj<-lapply(test_post,function(x) x$log_sig2)
  sig2_unadj<-sig2_adj
  stop<-FALSE
  eta<-0
  while(!stop){
    sig2_qs<-lapply(sig2_adj,function(x) quantile(x,c(0.025,0.975)))
    cov<-mean(unlist(lapply(sig2_qs,function(x) (x[1]<est$log_sig2)&(x[2]>est$log_sig2))))
    if(cov>=0.95){stop<-TRUE}
    
    eta<-eta+0.1
    sig2_adj<-lapply(sig2_unadj,function(x) mean(x)+eta*(x-mean(x)))
  }
  
  
  sig2_unadj<-palm_empirical_thomas0.4[[k]]$log_sig2
  sig2_adj<-mean(sig2_unadj)+eta*(sig2_unadj-mean(sig2_unadj))
  
  
  mu_adj<-lapply(test_post,function(x) x$log_mu)
  mu_unadj<-mu_adj
  stop<-FALSE
  eta<-0
  while(!stop){
    mu_qs<-lapply(mu_adj,function(x) quantile(x,c(0.025,0.975)))
    cov<-mean(unlist(lapply(mu_qs,function(x) (x[1]<est$log_mu)&(x[2]>est$log_mu))))
    if(cov>=0.95){stop<-TRUE}
    
    eta<-eta+0.1
    mu_adj<-lapply(mu_unadj,function(x) mean(x)+eta*(x-mean(x)))
  }
  
  
  mu_unadj<-palm_empirical_thomas0.4[[k]]$log_mu
  mu_adj<-mean(mu_unadj)+eta*(mu_unadj-mean(mu_unadj))
  
  
  out<-list()
  out$log_nu<-nu_adj
  out$log_sig2<-sig2_adj
  out$log_mu<-mu_adj
  
  return(out)
}
saveRDS(adj_out,paste0('adj_empirical_thomas2.rds'))




