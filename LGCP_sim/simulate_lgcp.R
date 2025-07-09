library(mvtnorm)
library(fields)
library(spatstat)


generate_RF<-function(sig2,phi,dist_mat){
  
  cov_mat<-sig2*exp(-dist_mat/phi)
  
  y<-rcpp_rmvnorm(n=1,mu=rep(0,nrow(dist_mat)),S = cov_mat)
  
  return(y[1,])
}


generate_LGCP<-function(X,beta,sig2,phi,grid,area,plot=FALSE){
  dist_mat<-as.matrix(dist(grid))
  
  Y<-generate_RF(sig2,phi,dist_mat)
  Y<-Y-mean(Y) #center
  
  lambda<-exp(X%*%beta+Y)
  N<-rpois(n=1,lambda=max(lambda))
  pts<-cbind(runif(N,0,1),runif(N,0,1))
  pts.keep<-NULL
  for(i in 1:N){
    
    ind<-which.min(sqrt((pts[i,1]-grid[,1])^2+(pts[i,2]-grid[,2])^2))
    
    p<-lambda[ind]/(max(lambda))
    if(runif(1)<p){
      pts.keep<-rbind(pts.keep,pts[i,])
      }
  }
  
  if(plot==TRUE){
    quilt.plot(grid[,1],grid[,2],log(lambda),nx=sqrt(nrow(grid)),ny=sqrt(nrow(grid)))
    points(pts.keep,pch=19)
                           }
  
  return(list(pts=pts.keep,lam=lambda))
}



