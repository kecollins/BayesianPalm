library(doParallel)
cl<-makeCluster(2)
registerDoParallel(cl)
registerDoParallel(cores=10)


pp_gen<-readRDS('pp_gen_homogeneous.rds')

palm_hom_sims0.4<-readRDS('palm_hom_sims0.4.rds')
palm_hom_sims_empirical0.4<-readRDS('palm_hom_sims_empirical0.4.rds')

nx<-20
grid<-expand.grid(seq((1/nx)/2,1,1/nx),seq((1/nx)/2,1,1/nx))

mu<-log(300)-1/2
sig2<-1
phi<-0.1
nu<-0.5
nx<-20
rmax<-0.4

n_sims<-100

adj_out<-foreach(k=c(1:n_sims),.errorhandling = 'pass')%do%{
  est<-list()
  est$lambda<-mean(palm_hom_sims0.4[[k]]$lambda)
  est$log_sig2<-mean(palm_hom_sims0.4[[k]]$log_sig2)
  est$log_phi<-mean(palm_hom_sims0.4[[k]]$log_phi)
  est$mu<-mean(log(palm_hom_sims0.4[[k]]$lambda)-exp(palm_hom_sims0.4[[k]]$log_sig2)/2)
  
  X<-1
  
  test_post<-foreach(k=c(1:100),.errorhandling = 'remove')%dopar%{
    pp_gen_temp<-generate_LGCP(X,beta=est$mu,sig2=exp(est$log_sig2),phi=exp(est$log_phi),grid=grid,plot=FALSE)
    
    pp<-list(points=pp_gen_temp$pts,grid=grid)
    
    data<-list(points=pp$points,nx=nx,rmax=rmax)
    dist_mat<-as.matrix(dist(pp$points))
    data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
    data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
    data$dist_cross<-as.matrix(rdist(pp$points,pp$grid))
    data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
    
    params<-list(lambda=nrow(data$points),log_sig2=log(sig2),log_phi=log(phi),eta=1)
    
    priors<-list(log_sig2.var=1000,log_phi.var=0.1,log_phi.mean=log(phi),
                 lambda.var=100000,lambda.pvar=1,log_sig2.pvar=0.1/4,log_phi.pvar=0.005)
    
    burn<-2000
    iters<-20000
    thin<-18
    
    post_palm<-mcmc_palm(data,params,priors,burn,iters,thin)
    
    return(post_palm)
  }
  
  mu_unadj<-lapply(test_post,function(x) log(x$lambda)-exp(x$log_sig2)/2)
  mu_adj<-mu_unadj
  stop<-FALSE
  eta<-0
  while(!stop){
    mu_qs<-lapply(mu_adj,function(x) quantile(x,c(0.025,0.975)))
    cov<-mean(unlist(lapply(mu_qs,function(x) (x[1]<est$mu)&(x[2]>est$mu))))
    if(cov>=0.95){stop<-TRUE}
    
    eta<-eta+0.1
    mu_adj<-lapply(mu_unadj,function(x) mean(x)+eta*(x-mean(x)))
  }
  
  
  mu_unadj<-log(palm_hom_sims0.4[[k]]$lambda)-exp(palm_hom_sims0.4[[k]]$log_sig2)/2
  mu_adj<-mean(mu_unadj)+eta*(mu_unadj-mean(mu_unadj))
  
  
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
  
  
  sig2_unadj<-palm_hom_sims0.4[[k]]$log_sig2
  sig2_adj<-mean(sig2_unadj)+eta*(sig2_unadj-mean(sig2_unadj))
  
  
  phi_adj<-lapply(test_post,function(x) x$log_phi)
  phi_unadj<-phi_adj
  stop<-FALSE
  eta<-0
  while(!stop){
    phi_qs<-lapply(phi_adj,function(x) quantile(x,c(0.025,0.975)))
    cov<-mean(unlist(lapply(phi_qs,function(x) (x[1]<est$log_phi)&(x[2]>est$log_phi))))
    if(cov>=0.95){stop<-TRUE}
    
    eta<-eta+0.1
    phi_adj<-lapply(phi_unadj,function(x) mean(x)+eta*(x-mean(x)))
  }
  
  
  phi_unadj<-palm_hom_sims0.4[[k]]$log_phi
  phi_adj<-mean(phi_unadj)+eta*(phi_unadj-mean(phi_unadj))
  
  
  out<-list()
  out$mu<-mu_adj
  out$log_sig2<-sig2_adj
  out$log_phi<-phi_adj
  
  return(out)
}
saveRDS(adj_out,paste0('adj_hom2.rds'))




adj_out<-foreach(k=c(1:n_sims),.errorhandling='remove')%do%{
  est<-list()
  est$lambda<-mean(palm_hom_sims_empirical0.4[[k]]$lambda)
  est$log_sig2<-mean(palm_hom_sims_empirical0.4[[k]]$log_sig2)
  est$log_phi<-mean(palm_hom_sims_empirical0.4[[k]]$log_phi)
  est$mu<-mean(log(palm_hom_sims_empirical0.4[[k]]$lambda)-exp(palm_hom_sims_empirical0.4[[k]]$log_sig2)/2)
  
  X<-1
  
  test_post<-foreach(k=c(1:100),.errorhandling = 'remove')%dopar%{
    pp_gen_temp<-generate_LGCP(X,beta=est$mu,sig2=exp(est$log_sig2),phi=exp(est$log_phi),grid=grid,plot=FALSE)
    
    pp<-list(points=pp_gen_temp$pts,grid=grid)
    
    data<-list(points=pp$points,nx=nx,rmax=rmax)
    dist_mat<-as.matrix(dist(pp$points))
    data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
    data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
    data$dist_cross<-as.matrix(rdist(pp$points,pp$grid))
    data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
    
    params<-list(lambda=nrow(data$points),log_sig2=log(sig2),log_phi=log(phi),eta=1)
    
    priors<-list(log_sig2.var=1000,log_phi.var=0.1,log_phi.mean=log(phi),
                 lambda.var=10,lambda.pvar=1,log_sig2.pvar=0.1/4,log_phi.pvar=0.005)
    
    burn<-2000
    iters<-20000
    thin<-18
    
    post_palm<-mcmc_palm(data,params,priors,burn,iters,thin)
    
    return(post_palm)
  }
  
  mu_adj<-lapply(test_post,function(x) log(x$lambda)-exp(x$log_sig2)/2)
  stop<-FALSE
  eta<-0
  while(!stop){
    mu_qs<-lapply(mu_adj,function(x) quantile(x,c(0.025,0.975)))
    cov<-mean(unlist(lapply(mu_qs,function(x) (x[1]<est$mu)&(x[2]>est$mu))))
    if(cov>=0.95){stop<-TRUE}
    
    eta<-eta+0.1
    mu_adj<-lapply(mu_adj,function(x) mean(x)+eta*(x-mean(x)))
  }
  
  
  mu_unadj<-lapply(test_post,function(x) log(x$lambda)-exp(x$log_sig2)/2)
  mu_adj<-mu_unadj
  stop<-FALSE
  eta<-0
  while(!stop){
    mu_qs<-lapply(mu_adj,function(x) quantile(x,c(0.025,0.975)))
    cov<-mean(unlist(lapply(mu_qs,function(x) (x[1]<est$mu)&(x[2]>est$mu))))
    if(cov>=0.95){stop<-TRUE}
    
    eta<-eta+0.1
    mu_adj<-lapply(mu_unadj,function(x) mean(x)+eta*(x-mean(x)))
  }
  
  
  mu_unadj<-log(palm_hom_sims_empirical0.4[[k]]$lambda)-exp(palm_hom_sims_empirical0.4[[k]]$log_sig2)/2
  mu_adj<-mean(mu_unadj)+eta*(mu_unadj-mean(mu_unadj))
  
  
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
  
  
  sig2_unadj<-palm_hom_sims_empirical0.4[[k]]$log_sig2
  sig2_adj<-mean(sig2_unadj)+eta*(sig2_unadj-mean(sig2_unadj))
  
  
  phi_adj<-lapply(test_post,function(x) x$log_phi)
  phi_unadj<-phi_adj
  stop<-FALSE
  eta<-0
  while(!stop){
    phi_qs<-lapply(phi_adj,function(x) quantile(x,c(0.025,0.975)))
    cov<-mean(unlist(lapply(phi_qs,function(x) (x[1]<est$log_phi)&(x[2]>est$log_phi))))
    if(cov>=0.95){stop<-TRUE}
    
    eta<-eta+0.1
    phi_adj<-lapply(phi_unadj,function(x) mean(x)+eta*(x-mean(x)))
  }
  
  
  phi_unadj<-palm_hom_sims_empirical0.4[[k]]$log_phi
  phi_adj<-mean(phi_unadj)+eta*(phi_unadj-mean(phi_unadj))
  
  
  out<-list()
  out$mu<-mu_adj
  out$log_sig2<-sig2_adj
  out$log_phi<-phi_adj
  
  return(out)
}
saveRDS(adj_out,paste0('adj_hom_empirical2.rds'))





