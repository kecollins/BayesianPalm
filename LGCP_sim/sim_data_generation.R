source('simulate_lgcp.R')

###################
### HOMOGENEOUS ###
###################

set.seed(80)
mu<-log(300)-1/2
sig2<-1
phi<-0.1
nu<-0.5
nx<-20
rmax=0.5

n_sims<-100

grid<-expand.grid(seq((1/nx)/2,1,1/nx),seq((1/nx)/2,1,1/nx))

X<-1
pp_gen<-list()
for(i in 1:n_sims){
  pp_gen[[i]]<-generate_LGCP(X,mu,sig2,phi,grid,plot=FALSE)
}

saveRDS(pp_gen,'pp_gen_homogeneous.rds')


#################
### COVARIATE ###
#################

set.seed(80)
nx<-20
grid<-expand.grid(seq((1/nx)/2,1,1/nx),seq((1/nx)/2,1,1/nx))

#### Create spatial covariate
X<-cbind(1,grid[,2])

beta<-c(4.5,2)
sig2<-1
phi<-0.1
rmax<-0.5

n_sims<-100

grid<-expand.grid(seq((1/nx)/2,1,1/nx),seq((1/nx)/2,1,1/nx))

pp_gen<-list()
for(i in 1:n_sims){
  pp_gen[[i]]<-generate_LGCP(X,beta,sig2,phi,grid,plot=FALSE)
}

saveRDS(pp_gen,'pp_gen_covariate.rds')

