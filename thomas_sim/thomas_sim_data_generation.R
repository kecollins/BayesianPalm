library(spatstat)

set.seed(80)
mu<-30
nu<-10
scale<-0.05

n_sims<-100

pp_gen<-list()
for(k in 1:n_sims){
  pp_gen[[k]]<-rThomas(kappa=mu,scale,mu=nu)
}

saveRDS(pp_gen,'pp_gen_thomas.rds')
