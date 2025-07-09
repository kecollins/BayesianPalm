source('full_lgcp_mcmc_1.2.R')
library(doParallel)
cl<-makeCluster(2)
registerDoParallel(cl)
registerDoParallel(cores=25)

pp_gen<-readRDS('pp_gen_homogeneous.rds')
X<-1
beta<-log(300)-1/2
sig2<-1
phi<-0.1
nu<-0.5
nx<-20

n_sims<-100

grid<-expand.grid(seq((1/nx)/2,1,1/nx),seq((1/nx)/2,1,1/nx))

sims_out<-foreach(k=c(1:n_sims))%dopar%{

  pp<-list()
  pp$points<-pp_gen[[k]]$pts
  pp$grid<-grid
  
  ##### FULL MCMC
  y0<-log(pp_gen[[k]]$lam)-X%*%beta
  
  
  ## Count up points in grid squares
  sq.index=NULL
  for(i in 1:nrow(pp$points)){
    sq.index[i]=which.min(sqrt((pp$points[i,1]-pp$grid[,1])^2+(pp$points[i,2]-pp$grid[,2])^2))
  }
  
  ni<-NULL
  for(i in 1:nrow(pp$grid)){
    ni[i]=sum(sq.index==i)
  }
  
  data<-list(ni=ni,nx=nx,dist=as.matrix(dist(pp$grid)),X=X)
  
  ## Parameters
  params<-list(beta=beta,log_sig2=log(sig2),log_phi=log(phi),nu=nu,y=y0)
  
  ## Priors
  priors<-list(beta.beta0=c(0,0),beta.sig20=1000,beta.pvar=0.1,
               log_sig2.mean=log(sig2),log_sig2.var=1,log_sig2.pvar=0.01/4,
               log_phi.mean=log(phi),log_phi.var=0.1,log_phi.pvar=0.005)
  
  
  post<-mcmc_lgcp(data,params,priors,20000,2000,18)
  
  return(post)
}

saveRDS(sims_out,'full_hom_sims.rds')
