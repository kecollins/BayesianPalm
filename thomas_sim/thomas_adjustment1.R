source('thomas_functions.R')
pp_gen<-readRDS('pp_gen_thomas.rds')
palm_thomas0.4<-readRDS('palm_thomas0.4.rds')
palm_empirical_thomas0.4<-readRDS('palm_empirical_thomas0.4.rds')

nx<-20
grid<-expand.grid(seq((1/nx)/2,1,1/nx),seq((1/nx)/2,1,1/nx))

rmax<-0.4
n_sims<-100

adj_out<-list()

for(k in 1:n_sims){
  H_inv<-cov(cbind(palm_thomas0.4[[k]]$lambda,palm_thomas0.4[[k]]$log_mu,palm_thomas0.4[[k]]$log_sig2))
  
  est<-list()
  est$lambda<-mean(palm_thomas0.4[[k]]$lambda)
  est$log_mu<-mean(palm_thomas0.4[[k]]$log_mu)
  est$log_sig2<-mean(palm_thomas0.4[[k]]$log_sig2)
  est$nu<-mean(palm_thomas0.4[[k]]$lambda/exp(palm_thomas0.4[[k]]$log_mu))
  
  X<-1
  score_out<-foreach(k=c(1:1000),.combine='rbind')%dopar%{
    
    pp<-rThomas(kappa=exp(est$log_mu),scale=sqrt(exp(est$log_sig2)),mu=est$nu)
    
    data<-list(points=cbind(pp$x,pp$y),nx=nx,rmax=rmax)
    dist_mat<-as.matrix(dist(data$points))
    data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
    data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
    data$dist_cross<-as.matrix(rdist(data$points,grid))
    data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
    
    return(palm_score(data,est))
  }
  
  J<-cov(score_out)
  
  eta<-1/sum(eigen(H_inv%*%J)$values)
  
  pp<-pp_gen[[k]]
  
  data<-list(points=cbind(pp$x,pp$y),nx=nx,rmax=rmax)
  dist_mat<-as.matrix(dist(data$points))
  data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
  data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
  data$dist_cross<-as.matrix(rdist(data$points,grid))
  data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
  
  params<-list(lambda=nrow(data$points),log_mu=log(mu),log_sig2=log(sig2),eta=eta)
  
  priors<-list(lambda.var=1000,log_mu.var=10,log_sig2.var=10,
               lambda.pvar=10,log_mu.pvar=0.1,log_sig2.pvar=0.1)
  
  burn<-2000
  iters<-20000
  thin<-18
  
  post_palm<-mcmc_palm(data,params,priors,iters,burn,thin)
  
  adj_out[[k]]<-post_palm
  print(k)
}

saveRDS(adj_out,'thomas_adj_hom1.rds')

palm_thomas0.4<-readRDS('palm_empirical_thomas0.4.rds')

adj_out<-list()

for(k in 1:n_sims){
  H_inv<-cov(cbind(palm_thomas0.4[[k]]$lambda,palm_thomas0.4[[k]]$log_mu,palm_thomas0.4[[k]]$log_sig2))
  
  est<-list()
  est$lambda<-mean(palm_thomas0.4[[k]]$lambda)
  est$log_mu<-mean(palm_thomas0.4[[k]]$log_mu)
  est$log_sig2<-mean(palm_thomas0.4[[k]]$log_sig2)
  est$nu<-mean(palm_thomas0.4[[k]]$lambda/exp(palm_thomas0.4[[k]]$log_mu))
  
  X<-1
  score_out<-foreach(k=c(1:1000),.combine='rbind')%dopar%{
    
    pp<-rThomas(kappa=exp(est$log_mu),scale=sqrt(exp(est$log_sig2)),mu=est$nu)
    
    data<-list(points=cbind(pp$x,pp$y),nx=nx,rmax=rmax)
    dist_mat<-as.matrix(dist(data$points))
    data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
    data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
    data$dist_cross<-as.matrix(rdist(data$points,grid))
    data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
    
    return(palm_score(data,est))
  }
  
  J<-cov(score_out)
  
  eta<-1/sum(eigen(H_inv%*%J)$values)
  
  pp<-pp_gen[[k]]
  
  data<-list(points=cbind(pp$x,pp$y),nx=nx,rmax=rmax)
  dist_mat<-as.matrix(dist(data$points))
  data$dist_pairs<-dist_mat[upper.tri(dist_mat)] 
  data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
  data$dist_cross<-as.matrix(rdist(data$points,grid))
  data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
  
  params<-list(lambda=nrow(data$points),log_mu=log(mu),log_sig2=log(sig2),eta=eta)
  
  priors<-list(lambda.var=10,log_mu.var=10,log_sig2.var=10,
               lambda.pvar=10,log_mu.pvar=0.1,log_sig2.pvar=0.1)
  
  burn<-2000
  iters<-20000
  thin<-18
  
  post_palm<-mcmc_palm(data,params,priors,iters,burn,thin)
  
  adj_out[[k]]<-post_palm
  print(k)
}

saveRDS(adj_out,'thomas_adj_hom_empirical1.rds')
