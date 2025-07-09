source('thomas_functions.R')
library(doParallel)
cl<-makeCluster(2)
registerDoParallel(cl)
registerDoParallel(cores=10)

pp_gen<-readRDS('pp_gen_thomas.rds')

n_sims<-100
rmax<-0.4

nx<-20
grid<-expand.grid(seq((1/nx)/2,1,1/nx),seq((1/nx)/2,1,1/nx))

sims_out<-foreach(k=c(1:n_sims))%dopar%{
  pp<-pp_gen[[k]]
  
  data<-list(points=cbind(pp$x,pp$y),nx=nx,rmax=rmax)
  dist_mat<-as.matrix(dist(data$points))
  data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
  data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
  data$dist_cross<-as.matrix(rdist(data$points,grid))
  data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
  
  params<-list(lambda=nrow(data$points),log_mu=log(mu),log_sig2=log(sig2),eta=1)
  
  priors<-list(lambda.var=1000,log_mu.var=10,log_sig2.var=10,
               lambda.pvar=10,log_mu.pvar=0.1,log_sig2.pvar=0.1)
  
  post_palm<-mcmc_palm(data,params,priors,20000,2000,18)
  
  return(post_palm)
}

saveRDS(sims_out,paste0('palm_thomas',rmax,'.rds'))


sims_out<-foreach(k=c(1:n_sims))%dopar%{
  pp<-pp_gen[[k]]
  
  data<-list(points=cbind(pp$x,pp$y),nx=nx,rmax=rmax)
  dist_mat<-as.matrix(dist(data$points))
  data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
  data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
  data$dist_cross<-as.matrix(rdist(data$points,grid))
  data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
  
  params<-list(lambda=nrow(data$points),log_mu=log(mu),log_sig2=log(sig2),eta=1)
  
  priors<-list(lambda.var=10,log_mu.var=10,log_sig2.var=10,
               lambda.pvar=10,log_mu.pvar=0.1,log_sig2.pvar=0.1)
  
  post_palm<-mcmc_palm(data,params,priors,20000,2000,18)
  
  return(post_palm)
}

saveRDS(sims_out,paste0('palm_empirical_thomas',rmax,'.rds'))

