source('palm_lambda_1.21.R')
library(doParallel)
library(MASS)
library(fields)
cl<-makeCluster(2)
registerDoParallel(cl)
registerDoParallel(cores=10)


mu<-log(300)-1/2
sig2<-1
phi<-0.1
nu<-0.5
nx<-20
rmax<-0.4

n_sims<-100

grid<-expand.grid(seq((1/nx)/2,1,1/nx),seq((1/nx)/2,1,1/nx))

pp_gen<-readRDS('pp_gen_homogeneous.rds')

sims_out<-foreach(k=c(1:n_sims))%dopar%{
  pp<-list(points=pp_gen[[k]]$pts,grid=grid)
  
  data<-list(points=pp$points,nx=nx,rmax=rmax)
  dist_mat<-as.matrix(dist(pp$points))
  data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
  data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
  data$dist_cross<-as.matrix(rdist(pp$points,pp$grid))
  data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
  
  
  ########################
  ### Initial MCMC run ###
  ########################
  
  params<-list(lambda=nrow(data$points),log_sig2=log(sig2),log_phi=log(phi),eta=1)
  
  priors<-list(log_sig2.var=1,log_phi.var=0.1,log_phi.mean=log(phi),
               lambda.var=10,lambda.pvar=10,log_sig2.pvar=0.1/4,log_phi.pvar=0.1)
  
  burn<-2000
  iters<-20000
  thin<-18
  
  post_palm<-mcmc_palm(data,params,priors,burn,iters,thin)

  return(post_palm)
}

saveRDS(sims_out,paste0('palm_hom_sims_empirical',rmax,'.rds'))

sims_out<-foreach(k=c(1:n_sims))%dopar%{
  pp<-list(points=pp_gen[[k]]$pts,grid=grid)
  
  data<-list(points=pp$points,nx=nx,rmax=rmax)
  dist_mat<-as.matrix(dist(pp$points))
  data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
  data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
  data$dist_cross<-as.matrix(rdist(pp$points,pp$grid))
  data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
  
  
  ########################
  ### Initial MCMC run ###
  ########################
  
  params<-list(lambda=nrow(data$points),log_sig2=log(sig2),log_phi=log(phi),eta=1)
  
  priors<-list(log_sig2.var=1,log_phi.var=0.1,log_phi.mean=log(phi),
               lambda.var=1000,lambda.pvar=10,log_sig2.pvar=0.1,log_phi.pvar=0.1)
  
  burn<-2000
  iters<-20000
  thin<-18
  
  post_palm<-mcmc_palm(data,params,priors,burn,iters,thin)

  return(post_palm)
}

saveRDS(sims_out,paste0('palm_hom_sims',rmax,'.rds'))
