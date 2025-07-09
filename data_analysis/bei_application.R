library(spatstat)
source('bei_functions.R')

pp<-bei
grid<-expand.grid(bei.extra$elev$xcol,bei.extra$elev$yrow)

plot(pp,pch=19,cex=0.5)
m1<-kppm(pp ~.,data=bei.extra_scaled,clusters="LGCP",method="palm",rmax=400)
logLik(m1)

m1

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

params<-list(beta=c(-4.9905810,0.1727164,0.33973),log_sig2=log(1.32184),log_phi=log(42.62346),eta=1)

priors<-list(phi.lower=log(20),phi.upper=log(200),phi.var=10,sig2.var=10,
             beta0.pvar=1,beta1.pvar=0.01,beta2.pvar=1,
             sig2.pvar=0.005,phi.pvar=0.005)

burn<-5000
iters<-10000
thin<-1

post_palm<-mcmc_palm(data,params,priors,burn,iters,thin)
saveRDS(post_palm,'bei_post.rds')

plot(bei_post$beta2_full[seq((burn:iters))/thin],type="l")

par(mfrow=c(1,3))
plot(post_palm$beta0_full,type="l")
plot(post_palm$beta1_full,type="l")
plot(post_palm$beta2_full,type="l")
plot(post_palm$log_sig2_full,type="l")
plot(post_palm$log_phi_full,type="l")

plot(post_palm$beta0,type="l")
plot(post_palm$beta1,type="l")
plot(post_palm$beta2,type="l")
plot(post_palm$log_sig2,type="l")
plot(post_palm$log_phi,type="l")


length(unique(post_palm$beta0_full[9000:10000]))/1000




