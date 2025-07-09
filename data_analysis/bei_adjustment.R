library(FastGP)
library(doParallel)
cl<-makeCluster(2)
registerDoParallel(cl)
registerDoParallel(cores=5)

bei_post<-readRDS('bei_post.rds')

est<-list()
est$beta<-c(mean(bei_post$beta0),mean(bei_post$beta1),mean(bei_post$beta2))
est$log_sig2<-mean(bei_post$log_sig2)
est$log_phi<-mean(bei_post$log_phi)


bei.extra_scaled<-bei.extra
bei.extra_scaled$elev$v<-(bei.extra_scaled$elev$v-mean(bei.extra_scaled$elev$v))/sd(bei.extra_scaled$elev$v)
bei.extra_scaled$grad$v<-(bei.extra_scaled$grad$v-mean(bei.extra_scaled$grad$v))/sd(bei.extra_scaled$grad$v)

pp_gen<-list()
for(k in 1:100){
  pp<-generate_LGCP(data$X,beta=est$beta,sig2=exp(est$log_sig2),phi=exp(est$log_phi),grid)
  pp_gen[[k]]<-as.ppp(pp$pts,bei$window)
  print(k)
}
saveRDS(pp_gen,'bei_pp_gen.rds')

pp_gen<-readRDS('bei_pp_gen.rds')




for(j in 1:20){
  test_post<-foreach(k=c((j-1)*5+1:5))%dopar%{
    pp<-pp_gen[[k]]
    
    pp$points<-cbind(pp$x,pp$y)
    pp$grid<-grid
    
    data<-list(X=cbind(1,as.vector(t(bei.extra$elev$v)),as.vector(t(bei.extra$grad$v))),
               rmax=400,points=pp$points)
    
    data$X[,2]<-(data$X[,2]-mean(data$X[,2]))/sd(data$X[,2])
    data$X[,3]<-(data$X[,3]-mean(data$X[,3]))/sd(data$X[,3])
    
    point_cov<-matrix(nrow=nrow(data$points),ncol=3)
    for(i in 1:nrow(data$points)){
      temp_ind<-which.min(sqrt((pp$grid[,1]-pp$points[i,1])^2+(pp$grid[,2]-pp$points[i,2])^2))
      point_cov[i,]<-data$X[temp_ind,]
    }
    
    
    
    dist_mat<-as.matrix(dist(pp$points))
    ind_mat<-matrix(rep(1:nrow(pp$points),ncol(dist_mat)),ncol=ncol(dist_mat))
    data$ind_pairs<-ind_mat[which(dist_mat!=0)]
    data$dist_pairs<-dist_mat[which(dist_mat!=0)] 
    data$ind_pairs<-data$ind_pairs[which(data$dist_pairs<(data$rmax/2))]
    data$dist_pairs<-data$dist_pairs[which(data$dist_pairs<(data$rmax/2))]
    
    data$X1<-point_cov[data$ind_pairs,]
    
    data$dist_cross<-as.matrix(rdist(pp$points,pp$grid))
    ind_mat2<-matrix(rep(1:nrow(pp$grid),nrow(pp$points)),nrow=nrow(pp$points),byrow=TRUE)
    data$ind_pairs2<-ind_mat2[which(data$dist_cross<(data$rmax/2))]
    data$dist_cross<-data$dist_cross[which(data$dist_cross<(data$rmax/2))]
    
    data$X2<-data$X[data$ind_pairs2,]
    
    
    ########################
    ### Initial MCMC run ###
    ########################
    
    params<-list(beta=c(-5.5686,0.15013,0.12973),log_sig2=0.4,log_phi=3.5,eta=1)
    
    priors<-list(phi.lower=log(20),phi.upper=log(200),phi.var=10,sig2.var=10,
                 beta0.pvar=1,beta1.pvar=0.01,beta2.pvar=1,
                 sig2.pvar=0.005,phi.pvar=0.005)
    
    
    burn<-2000
    iters<-10000
    thin<-8
    
    post_palm<-mcmc_palm(data,params,priors,burn,iters,thin)
    return(post_palm)
  }
  saveRDS(test_post,paste0('test_post_2',j,'.rds'))
  print(j)
}


test_post<-NULL
for(j in 1:20){
  temp<-readRDS(paste0('test_post_2',j,'.rds'))
  test_post<-c(test_post,temp)
}



beta0_adj<-lapply(test_post,function(x) x$beta0_full[5001:10000])
beta0_unadj<-beta0_adj
stop<-FALSE
eta<-0
while(!stop){
  beta0_hpd<-lapply(beta0_adj,function(x) HPDinterval(mcmc(x)))
  cov<-mean(unlist(lapply(beta0_hpd,function(x) (x[1]<est$beta[1])&(x[2]>est$beta[1]))))
  
  #beta0_qs<-lapply(beta0_adj,function(x) quantile(x,c(0.025,0.975)))
  #cov<-mean(unlist(lapply(beta0_qs,function(x) (x[1]<est$beta[1])&(x[2]>est$beta[1]))))
  if(cov>=0.95){stop<-TRUE}
  
  eta<-eta+0.1
  beta0_adj<-lapply(beta0_unadj,function(x) mean(x)+eta*(x-mean(x)))
}


beta0_unadj<-bei_post$beta0_full[5001:10000]
beta0_adj<-mean(beta0_unadj)+eta*(beta0_unadj-mean(beta0_unadj))

beta1_adj<-lapply(test_post,function(x) x$beta1_full[5001:10000])
beta1_unadj<-beta1_adj
stop<-FALSE
eta<-0
while(!stop){
  beta1_hpd<-lapply(beta1_adj,function(x) HPDinterval(mcmc(x)))
  cov<-mean(unlist(lapply(beta1_hpd,function(x) (x[1]<est$beta[2])&(x[2]>est$beta[2]))))
  
  #beta1_qs<-lapply(beta1_adj,function(x) quantile(x,c(0.025,0.975)))
  #cov<-mean(unlist(lapply(beta1_qs,function(x) (x[1]<est$beta[2])&(x[2]>est$beta[2]))))
  if(cov>=0.95){stop<-TRUE}
  
  eta<-eta+0.1
  beta1_adj<-lapply(beta1_unadj,function(x) mean(x)+eta*(x-mean(x)))
}


beta1_unadj<-bei_post$beta1_full[5001:10000]
beta1_adj<-mean(beta1_unadj)+eta*(beta1_unadj-mean(beta1_unadj))

beta2_adj<-lapply(test_post,function(x) x$beta2_full[5001:10000])
beta2_unadj<-beta2_adj
stop<-FALSE
eta<-0
while(!stop){
  beta2_hpd<-lapply(beta2_adj,function(x) HPDinterval(mcmc(x)))
  cov<-mean(unlist(lapply(beta2_hpd,function(x) (x[1]<est$beta[3])&(x[2]>est$beta[3]))))
  
  #beta2_qs<-lapply(beta2_adj,function(x) quantile(x,c(0.025,0.975)))
  #cov<-mean(unlist(lapply(beta2_qs,function(x) (x[1]<est$beta[3])&(x[2]>est$beta[3]))))
  if(cov>=0.95){stop<-TRUE}
  
  eta<-eta+0.1
  beta2_adj<-lapply(beta2_unadj,function(x) mean(x)+eta*(x-mean(x)))
}


beta2_unadj<-bei_post$beta2_full[5001:10000]
beta2_adj<-mean(beta2_unadj)+eta*(beta2_unadj-mean(beta2_unadj))



sig2_adj<-lapply(test_post,function(x) x$log_sig2_full[5001:10000])
sig2_unadj<-sig2_adj
stop<-FALSE
eta<-0
while(!stop){
  sig2_hpd<-lapply(sig2_adj,function(x) HPDinterval(mcmc(x)))
  cov<-mean(unlist(lapply(sig2_hpd,function(x) (x[1]<est$log_sig2)&(x[2]>est$log_sig2))))
  
  #sig2_qs<-lapply(sig2_adj,function(x) quantile(x,c(0.025,0.975)))
  #cov<-mean(unlist(lapply(sig2_qs,function(x) (x[1]<est$log_sig2)&(x[2]>est$log_sig2))))
  if(cov>=0.95){stop<-TRUE}
  
  eta<-eta+0.1
  sig2_adj<-lapply(sig2_unadj,function(x) mean(x)+eta*(x-mean(x)))
}


sig2_unadj<-bei_post$log_sig2_full[5001:10000]
sig2_adj<-mean(sig2_unadj)+eta*(sig2_unadj-mean(sig2_unadj))


phi_adj<-lapply(test_post,function(x) x$log_phi_full[5001:10000])
phi_unadj<-phi_adj
stop<-FALSE
eta<-0
while(!stop){
  phi_hpd<-lapply(phi_adj,function(x) HPDinterval(mcmc(x)))
  cov<-mean(unlist(lapply(phi_hpd,function(x) (x[1]<est$log_phi)&(x[2]>est$log_phi))))
  
  #phi_qs<-lapply(phi_adj,function(x) quantile(x,c(0.025,0.975)))
  #cov<-mean(unlist(lapply(phi_qs,function(x) (x[1]<est$log_phi)&(x[2]>est$log_phi))))
  if(cov>=0.95){stop<-TRUE}
  
  eta<-eta+0.1
  phi_adj<-lapply(phi_unadj,function(x) mean(x)+eta*(x-mean(x)))
}


phi_unadj<-bei_post$log_phi_full[5001:10000]
phi_adj<-mean(phi_unadj)+eta*(phi_unadj-mean(phi_unadj))

#### Posterior means
mean(beta0_adj-mean(bei.extra$elev$v)*beta1_adj/sd(bei.extra$elev$v)-mean(bei.extra$grad$v)*beta2_adj/sd(bei.extra$grad$v))
mean(beta1_adj/sd(bei.extra$elev$v))
mean(beta2_adj/sd(bei.extra$grad$v))
mean(exp(sig2_adj))
mean(exp(phi_adj))

#### 95% CI
quantile(beta0_adj-mean(bei.extra$elev$v)*beta1_adj/sd(bei.extra$elev$v)-mean(bei.extra$grad$v)*beta2_adj/sd(bei.extra$grad$v),c(0.025,0.975))
quantile(beta1_adj/sd(bei.extra$elev$v),c(0.025,0.975))
quantile(beta2_adj/sd(bei.extra$grad$v),c(0.025,0.975))
quantile(exp(sig2_adj),c(0.025,0.975))
quantile(exp(phi_adj),c(0.025,0.975))


#### 90% CI
quantile(beta0_adj-mean(bei.extra$elev$v)*beta1_adj/sd(bei.extra$elev$v)-mean(bei.extra$grad$v)*beta2_adj/sd(bei.extra$grad$v),c(0.05,0.95))
quantile(beta1_adj/sd(bei.extra$elev$v),c(0.05,0.95))
quantile(beta2_adj/sd(bei.extra$grad$v),c(0.05,0.95))
quantile(exp(sig2_adj),c(0.05,0.95))
quantile(exp(phi_adj),c(0.05,0.95))


#### 95% HPD
HPDinterval(mcmc(beta0_adj-mean(bei.extra$elev$v)*beta1_adj/sd(bei.extra$elev$v)-mean(bei.extra$grad$v)*beta2_adj/sd(bei.extra$grad$v)),prob=0.95)
HPDinterval(mcmc(beta1_adj/sd(bei.extra$elev$v)),prob=0.95)
HPDinterval(mcmc(beta2_adj/sd(bei.extra$grad$v)),prob=0.95)
HPDinterval(mcmc(exp(sig2_adj)),prob=0.95)
HPDinterval(mcmc(exp(phi_adj)),prob=0.95)


#### 90% HPD
HPDinterval(mcmc(beta0_adj-mean(bei.extra$elev$v)*beta1_adj/sd(bei.extra$elev$v)-mean(bei.extra$grad$v)*beta2_adj/sd(bei.extra$grad$v)),prob=0.9)
HPDinterval(mcmc(beta1_adj/sd(bei.extra$elev$v)),prob=0.9)
HPDinterval(mcmc(beta2_adj/sd(bei.extra$grad$v)),prob=0.9)
HPDinterval(mcmc(exp(sig2_adj)),prob=0.9)
HPDinterval(mcmc(exp(phi_adj)),prob=0.9)





#### Pretty pictures
pdf("bei_points.pdf",width=10,height=5)
plot(bei,pch=19,main="")
dev.off()

pdf("bei_cov.pdf",width=10,height=5)
plot(bei.extra,main="")
dev.off()

