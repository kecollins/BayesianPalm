source('simulate_LGCP.R')
pp_gen<-readRDS('pp_gen_homogeneous.rds')

palm_hom_sims0.4<-readRDS('palm_hom_sims0.4.rds')

nx<-20
grid<-expand.grid(seq((1/nx)/2,1,1/nx),seq((1/nx)/2,1,1/nx))

rmax<-0.4
n_sims<-100

adj_out<-list()

for(k in 1:n_sims){
  H_inv<-cov(cbind(palm_hom_sims0.4[[k]]$lambda,palm_hom_sims0.4[[k]]$log_sig2,palm_hom_sims0.4[[k]]$log_phi))
  
  est<-list()
  est$lambda<-mean(palm_hom_sims0.4[[k]]$lambda)
  est$log_sig2<-mean(palm_hom_sims0.4[[k]]$log_sig2)
  est$log_phi<-mean(palm_hom_sims0.4[[k]]$log_phi)
  est$mu<-mean(log(palm_hom_sims0.4[[k]]$lambda)-exp(palm_hom_sims0.4[[k]]$log_sig2)/2)
  
  X<-1
  score_out<-foreach(k=c(1:1000),.combine='rbind')%dopar%{
    
    pp_gen_temp<-generate_LGCP(X,beta=est$mu,sig2=exp(est$log_sig2),phi=exp(est$log_phi),grid=grid,plot=FALSE)
    
    pp<-list(points=pp_gen_temp$pts,grid=grid)
    
    data<-list(points=pp$points,nx=nx,rmax=rmax)
    dist_mat<-as.matrix(dist(pp$points))
    data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
    data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
    data$dist_cross<-as.matrix(rdist(pp$points,pp$grid))
    data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
    
    return(palm_score(data,est))
  }
  
  J<-cov(score_out)
  
  eta<-1/sum(eigen(H_inv%*%J)$values)
  
  
  
  pp<-list(points=pp_gen[[k]]$pts,grid=grid)
  
  data<-list(points=pp$points,nx=nx,rmax=rmax)
  dist_mat<-as.matrix(dist(pp$points))
  data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
  data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
  data$dist_cross<-as.matrix(rdist(pp$points,pp$grid))
  data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
  
  
  params<-list(lambda=nrow(data$points),log_sig2=log(sig2),log_phi=log(phi),eta=eta)
  
  priors<-list(log_sig2.var=1,log_phi.var=0.1,log_phi.mean=log(phi),
               lambda.var=1000,lambda.pvar=10,log_sig2.pvar=0.1,log_phi.pvar=0.1)
  
  burn<-2000
  iters<-20000
  thin<-18
  
  post_palm<-mcmc_palm(data,params,priors,burn,iters,thin)
  
  adj_out[[k]]<-post_palm
  print(k)
}

saveRDS(adj_out,'adj_hom1.rds')

palm_hom_sims0.4<-readRDS('palm_hom_sims_empirical0.4.rds')

adj_out<-list()

for(k in 1:n_sims){
  H_inv<-cov(cbind(palm_hom_sims0.4[[k]]$lambda,palm_hom_sims0.4[[k]]$log_sig2,palm_hom_sims0.4[[k]]$log_phi))
  
  est<-list()
  est$lambda<-mean(palm_hom_sims0.4[[k]]$lambda)
  est$log_sig2<-mean(palm_hom_sims0.4[[k]]$log_sig2)
  est$log_phi<-mean(palm_hom_sims0.4[[k]]$log_phi)
  est$mu<-mean(log(palm_hom_sims0.4[[k]]$lambda)-exp(palm_hom_sims0.4[[k]]$log_sig2)/2)
  
  X<-1
  score_out<-foreach(k=c(1:1000),.combine='rbind')%dopar%{
    
    pp_gen_temp<-generate_LGCP(X,beta=est$mu,sig2=exp(est$log_sig2),phi=exp(est$log_phi),grid=grid,plot=FALSE)
    
    pp<-list(points=pp_gen_temp$pts,grid=grid)
    
    data<-list(points=pp$points,nx=nx,rmax=rmax)
    dist_mat<-as.matrix(dist(pp$points))
    data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
    data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
    data$dist_cross<-as.matrix(rdist(pp$points,pp$grid))
    data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
    
    return(palm_score(data,est))
  }
  
  J<-cov(score_out)
  
  eta<-1/sum(eigen(H_inv%*%J)$values)
  
  
  
  pp<-list(points=pp_gen[[k]]$pts,grid=grid)
  
  data<-list(points=pp$points,nx=nx,rmax=rmax)
  dist_mat<-as.matrix(dist(pp$points))
  data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
  data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
  data$dist_cross<-as.matrix(rdist(pp$points,pp$grid))
  data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
  
  
  params<-list(lambda=nrow(data$points),log_sig2=log(sig2),log_phi=log(phi),eta=eta)
  
  priors<-list(log_sig2.var=1,log_phi.var=0.1,log_phi.mean=log(phi),
               lambda.var=10,lambda.pvar=10,log_sig2.pvar=0.1,log_phi.pvar=0.1)
  
  burn<-2000
  iters<-20000
  thin<-18
  
  post_palm<-mcmc_palm(data,params,priors,burn,iters,thin)
  
  adj_out[[k]]<-post_palm
  print(k)
}

saveRDS(adj_out,'adj_hom_empirical1.rds')
