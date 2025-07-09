library(coda)
library(kableExtra)

pp_gen<-readRDS('pp_gen_homogeneous.rds')

palm_hom_sims0.4<-readRDS('palm_hom_sims0.4.rds')
palm_hom_sims0.6<-readRDS('palm_hom_sims0.6.rds')
palm_hom_sims0.8<-readRDS('palm_hom_sims0.8.rds')
palm_hom_sims1<-readRDS('palm_hom_sims1.rds')

palm_hom_sims_empirical0.4<-readRDS('palm_hom_sims_empirical0.4.rds')
palm_hom_sims_empirical0.6<-readRDS('palm_hom_sims_empirical0.6.rds')
palm_hom_sims_empirical0.8<-readRDS('palm_hom_sims_empirical0.8.rds')
palm_hom_sims_empirical1<-readRDS('palm_hom_sims_empirical1.rds')


full_hom_sims<-readRDS('full_hom_sims.rds')

mu<-log(300)-1/2
sig2<-1
phi<-0.1

### BIASES

# mu
mu_bias_full<-mean(unlist(lapply(full_hom_sims,function(x) mean(x$beta)))-mu)

mu_bias_0.4<-mean(unlist(lapply(palm_hom_sims0.4,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)
mu_bias_0.6<-mean(unlist(lapply(palm_hom_sims0.6,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)
mu_bias_0.8<-mean(unlist(lapply(palm_hom_sims0.8,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)
mu_bias_1<-mean(unlist(lapply(palm_hom_sims1,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)


mu_bias_0.4E<-mean(unlist(lapply(palm_hom_sims_empirical0.4,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)
mu_bias_0.6E<-mean(unlist(lapply(palm_hom_sims_empirical0.6,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)
mu_bias_0.8E<-mean(unlist(lapply(palm_hom_sims_empirical0.8,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)
mu_bias_1E<-mean(unlist(lapply(palm_hom_sims_empirical1,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)


# sig2
sig2_bias_full<-mean(unlist(lapply(full_hom_sims,function(x) mean(exp(x$log_sig2))))-sig2)

sig2_bias_0.4<-mean(unlist(lapply(palm_hom_sims0.4,function(x) mean(exp(x$log_sig2))))-sig2)
sig2_bias_0.6<-mean(unlist(lapply(palm_hom_sims0.6,function(x) mean(exp(x$log_sig2))))-sig2)
sig2_bias_0.8<-mean(unlist(lapply(palm_hom_sims0.8,function(x) mean(exp(x$log_sig2))))-sig2)
sig2_bias_1<-mean(unlist(lapply(palm_hom_sims1,function(x) mean(exp(x$log_sig2))))-sig2)

sig2_bias_0.4E<-mean(unlist(lapply(palm_hom_sims_empirical0.4,function(x) mean(exp(x$log_sig2))))-sig2)
sig2_bias_0.6E<-mean(unlist(lapply(palm_hom_sims_empirical0.6,function(x) mean(exp(x$log_sig2))))-sig2)
sig2_bias_0.8E<-mean(unlist(lapply(palm_hom_sims_empirical0.8,function(x) mean(exp(x$log_sig2))))-sig2)
sig2_bias_1E<-mean(unlist(lapply(palm_hom_sims_empirical1,function(x) mean(exp(x$log_sig2))))-sig2)


# phi
phi_bias_full<-mean(unlist(lapply(full_hom_sims,function(x) mean(exp(x$log_phi))))-phi)

phi_bias_0.4<-mean(unlist(lapply(palm_hom_sims0.4,function(x) mean(exp(x$log_phi))))-phi)
phi_bias_0.6<-mean(unlist(lapply(palm_hom_sims0.6,function(x) mean(exp(x$log_phi))))-phi)
phi_bias_0.8<-mean(unlist(lapply(palm_hom_sims0.8,function(x) mean(exp(x$log_phi))))-phi)
phi_bias_1<-mean(unlist(lapply(palm_hom_sims1,function(x) mean(exp(x$log_phi))))-phi)

phi_bias_0.4E<-mean(unlist(lapply(palm_hom_sims_empirical0.4,function(x) mean(exp(x$log_phi))))-phi)
phi_bias_0.6E<-mean(unlist(lapply(palm_hom_sims_empirical0.6,function(x) mean(exp(x$log_phi))))-phi)
phi_bias_0.8E<-mean(unlist(lapply(palm_hom_sims_empirical0.8,function(x) mean(exp(x$log_phi))))-phi)
phi_bias_1E<-mean(unlist(lapply(palm_hom_sims_empirical1,function(x) mean(exp(x$log_phi))))-phi)


### RMSE

# mu
mu_rmse_full<-sqrt(mean((unlist(lapply(full_hom_sims,function(x) mean(x$beta)))-mu)^2))

mu_rmse_0.4<-sqrt(mean((unlist(lapply(palm_hom_sims0.4,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)^2))
mu_rmse_0.6<-sqrt(mean((unlist(lapply(palm_hom_sims0.6,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)^2))
mu_rmse_0.8<-sqrt(mean((unlist(lapply(palm_hom_sims0.8,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)^2))
mu_rmse_1<-sqrt(mean((unlist(lapply(palm_hom_sims1,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)^2))


mu_rmse_0.4E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical0.4,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)^2))
mu_rmse_0.6E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical0.6,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)^2))
mu_rmse_0.8E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical0.8,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)^2))
mu_rmse_1E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical1,function(x) mean(log(x$lambda)-exp(x$log_sig2)/2)))-mu)^2))


# sig2
sig2_rmse_full<-sqrt(mean((unlist(lapply(full_hom_sims,function(x) mean(exp(x$log_sig2))))-sig2)^2))

sig2_rmse_0.4<-sqrt(mean((unlist(lapply(palm_hom_sims0.4,function(x) mean(exp(x$log_sig2))))-sig2)^2))
sig2_rmse_0.6<-sqrt(mean((unlist(lapply(palm_hom_sims0.6,function(x) mean(exp(x$log_sig2))))-sig2)^2))
sig2_rmse_0.8<-sqrt(mean((unlist(lapply(palm_hom_sims0.8,function(x) mean(exp(x$log_sig2))))-sig2)^2))
sig2_rmse_1<-sqrt(mean((unlist(lapply(palm_hom_sims1,function(x) mean(exp(x$log_sig2))))-sig2)^2))

sig2_rmse_0.4E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical0.4,function(x) mean(exp(x$log_sig2))))-sig2)^2))
sig2_rmse_0.6E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical0.6,function(x) mean(exp(x$log_sig2))))-sig2)^2))
sig2_rmse_0.8E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical0.8,function(x) mean(exp(x$log_sig2))))-sig2)^2))
sig2_rmse_1E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical1,function(x) mean(exp(x$log_sig2))))-sig2)^2))


# phi
phi_rmse_full<-sqrt(mean((unlist(lapply(full_hom_sims,function(x) mean(exp(x$log_phi))))-phi)^2))

phi_rmse_0.4<-sqrt(mean((unlist(lapply(palm_hom_sims0.4,function(x) mean(exp(x$log_phi))))-phi)^2))
phi_rmse_0.6<-sqrt(mean((unlist(lapply(palm_hom_sims0.6,function(x) mean(exp(x$log_phi))))-phi)^2))
phi_rmse_0.8<-sqrt(mean((unlist(lapply(palm_hom_sims0.8,function(x) mean(exp(x$log_phi))))-phi)^2))
phi_rmse_1<-sqrt(mean((unlist(lapply(palm_hom_sims1,function(x) mean(exp(x$log_phi))))-phi)^2))

phi_rmse_0.4E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical0.4,function(x) mean(exp(x$log_phi))))-phi)^2))
phi_rmse_0.6E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical0.6,function(x) mean(exp(x$log_phi))))-phi)^2))
phi_rmse_0.8E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical0.8,function(x) mean(exp(x$log_phi))))-phi)^2))
phi_rmse_1E<-sqrt(mean((unlist(lapply(palm_hom_sims_empirical1,function(x) mean(exp(x$log_phi))))-phi)^2))




posterior_mean<-data.frame(round(
  rbind(c(mu_bias_full,sig2_bias_full,phi_bias_full,mu_rmse_full,sig2_rmse_full,phi_rmse_full),
        c(mu_bias_0.4,sig2_bias_0.4,phi_bias_0.4,mu_rmse_0.4,sig2_rmse_0.4,phi_rmse_0.4),
        c(mu_bias_0.6,sig2_bias_0.6,phi_bias_0.6,mu_rmse_0.6,sig2_rmse_0.6,phi_rmse_0.6),
        c(mu_bias_0.8,sig2_bias_0.8,phi_bias_0.8,mu_rmse_0.8,sig2_rmse_0.8,phi_rmse_0.8),
        c(mu_bias_1,sig2_bias_1,phi_bias_1,mu_rmse_1,sig2_rmse_1,phi_rmse_1),
        c(mu_bias_0.4E,sig2_bias_0.4E,phi_bias_0.4E,mu_rmse_0.4E,sig2_rmse_0.4E,phi_rmse_0.4E),
        c(mu_bias_0.6E,sig2_bias_0.6E,phi_bias_0.6E,mu_rmse_0.6E,sig2_rmse_0.6E,phi_rmse_0.6E),
        c(mu_bias_0.8E,sig2_bias_0.8E,phi_bias_0.8E,mu_rmse_0.8E,sig2_rmse_0.8E,phi_rmse_0.8E),
        c(mu_bias_1E,sig2_bias_1E,phi_bias_1E,mu_rmse_1E,sig2_rmse_1E,phi_rmse_1E)
  )
,2))

rownames(posterior_mean)<-c("FL","PL_0.2","PL_0.3","PL_0.4","PL_0.5","PLE_0.2","PLE_0.3","PLE_0.4","PLE_0.5")

colnames(posterior_mean)<-c("mu","sigma2","phi","mu","sigma2","phi")

kable(posterior_mean,format='latex')


### ESS/Minute
time_full<-unlist(lapply(full_hom_sims,function(x) as.numeric(x$time,units="secs")))
time_0.4<-unlist(lapply(palm_hom_sims0.4,function(x) as.numeric(x$time,units="secs")))
time_0.6<-unlist(lapply(palm_hom_sims0.6,function(x) as.numeric(x$time,units="secs")))
time_0.8<-unlist(lapply(palm_hom_sims0.8,function(x) as.numeric(x$time,units="secs")))
time_1<-unlist(lapply(palm_hom_sims1,function(x) as.numeric(x$time,units="secs")))

mu_ess_full<-unlist(lapply(full_hom_sims,function(x) effectiveSize(x$beta_full)))
mu_ess0.4<-unlist(lapply(palm_hom_sims0.4,function(x) effectiveSize(log(x$lambda_full)-exp(x$log_sig2_full)/2)))
mu_ess0.6<-unlist(lapply(palm_hom_sims0.6,function(x) effectiveSize(log(x$lambda_full)-exp(x$log_sig2_full)/2)))
mu_ess0.8<-unlist(lapply(palm_hom_sims0.8,function(x) effectiveSize(log(x$lambda_full)-exp(x$log_sig2_full)/2)))
mu_ess1<-unlist(lapply(palm_hom_sims1,function(x) effectiveSize(log(x$lambda_full)-exp(x$log_sig2_full)/2)))

sig2_ess_full<-unlist(lapply(full_hom_sims,function(x) effectiveSize(x$log_sig2_full)))
sig2_ess0.4<-unlist(lapply(palm_hom_sims0.4,function(x) effectiveSize(x$log_sig2_full)))
sig2_ess0.6<-unlist(lapply(palm_hom_sims0.6,function(x) effectiveSize(x$log_sig2_full)))
sig2_ess0.8<-unlist(lapply(palm_hom_sims0.8,function(x) effectiveSize(x$log_sig2_full)))
sig2_ess1<-unlist(lapply(palm_hom_sims1,function(x) effectiveSize(x$log_sig2_full)))

phi_ess_full<-unlist(lapply(full_hom_sims,function(x) effectiveSize(x$log_phi_full)))
phi_ess0.4<-unlist(lapply(palm_hom_sims0.4,function(x) effectiveSize(x$log_phi_full)))
phi_ess0.6<-unlist(lapply(palm_hom_sims0.6,function(x) effectiveSize(x$log_phi_full)))
phi_ess0.8<-unlist(lapply(palm_hom_sims0.8,function(x) effectiveSize(x$log_phi_full)))
phi_ess1<-unlist(lapply(palm_hom_sims1,function(x) effectiveSize(x$log_phi_full)))


kable(round(cbind(
  c(mean(mu_ess_full/time_full),mean(mu_ess0.4/time_0.4),mean(mu_ess0.6/time_0.6),
    mean(mu_ess0.8/time_0.8),mean(mu_ess1/time_1)),
  c(mean(sig2_ess_full/time_full),mean(sig2_ess0.4/time_0.4),mean(sig2_ess0.6/time_0.6),
    mean(sig2_ess0.8/time_0.8),mean(sig2_ess1/time_1)),
  c(mean(phi_ess_full/time_full),mean(phi_ess0.4/time_0.4),mean(phi_ess0.6/time_0.6),
    mean(phi_ess0.8/time_0.8),mean(phi_ess1/time_1))
),2),format='latex')

## check acceptance rates


mean(unlist(lapply(palm_hom_sims0.4,function(x) length(unique(x$lambda_full))/20000)))
mean(unlist(lapply(palm_hom_sims0.6,function(x) length(unique(x$lambda_full))/20000)))
mean(unlist(lapply(palm_hom_sims0.8,function(x) length(unique(x$lambda_full))/20000)))
mean(unlist(lapply(palm_hom_sims1,function(x) length(unique(x$lambda_full))/20000)))


##################################
### CREDIBLE INTERVAL COVERAGES ##
##################################

adj_hom1<-readRDS('adj_hom1.rds')
adj_hom_empirical1<-readRDS('adj_hom_empirical1.rds')

adj_hom2<-readRDS('adj_hom2.rds')
adj_hom_empirical2<-readRDS('adj_hom_empirical2.rds')

######### EMPIRICAL COVERAGE/CI LENGTH
###### FULL MCMC

mu_full_qs<-lapply(full_hom_sims,function(x) quantile(x$beta,c(0.025,0.975)))
mu_full_cov<-mean(unlist(lapply(mu_full_qs,function(x) (x[1]<mu)*(x[2]>mu))))
mu_full_length<-mean(unlist(lapply(mu_full_qs,function(x) x[2]-x[1])))
mu_full_length_med<-median(unlist(lapply(mu_full_qs,function(x) x[2]-x[1])))

log_sig2_full_qs<-lapply(full_hom_sims,function(x) quantile(x$log_sig2,c(0.025,0.975)))
log_sig2_full_cov<-mean(unlist(lapply(log_sig2_full_qs,function(x) (x[1]<log(sig2))*(x[2]>log(sig2)))))
log_sig2_full_length<-mean(unlist(lapply(log_sig2_full_qs,function(x) x[2]-x[1])))
log_sig2_full_length_med<-median(unlist(lapply(log_sig2_full_qs,function(x) x[2]-x[1])))

log_phi_full_qs<-lapply(full_hom_sims,function(x) quantile(x$log_phi,c(0.025,0.975)))
log_phi_full_cov<-mean(unlist(lapply(log_phi_full_qs,function(x) (x[1]<log(phi))*(x[2]>log(phi)))))
log_phi_full_length<-mean(unlist(lapply(log_phi_full_qs,function(x) x[2]-x[1])))
log_phi_full_length_med<-median(unlist(lapply(log_phi_full_qs,function(x) x[2]-x[1])))

###### NON EMPIRICAL
### palm_hom_0.4 unadjusted
mu_unadj_qs1<-lapply(palm_hom_sims0.4,function(x) quantile(log(x$lambda)-exp(x$log_sig2)/2,c(0.025,0.975)))
mu_unadj_cov1<-mean(unlist(lapply(mu_unadj_qs1, function(x) (x[1]<mu)*(x[2]>mu))))
mu_unadj_length1<-mean(unlist(lapply(mu_unadj_qs1,function(x) x[2]-x[1])))
mu_unadj_length1_med<-median(unlist(lapply(mu_unadj_qs1,function(x) x[2]-x[1])))

log_sig2_unadj_qs1<-lapply(palm_hom_sims0.4,function(x) quantile(x$log_sig2,c(0.025,0.975)))
log_sig2_unadj_cov1<-mean(unlist(lapply(log_sig2_unadj_qs1, function(x) (x[1]<log(sig2))*(x[2]>log(sig2)))))
log_sig2_unadj_length1<-mean(unlist(lapply(log_sig2_unadj_qs1, function(x) x[2]-x[1])))
log_sig2_unadj_length1_med<-median(unlist(lapply(log_sig2_unadj_qs1, function(x) x[2]-x[1])))

log_phi_unadj_qs1<-lapply(palm_hom_sims0.4,function(x) quantile(x$log_phi,c(0.025,0.975)))
log_phi_unadj_cov1<-mean(unlist(lapply(log_phi_unadj_qs1, function(x) (x[1]<log(phi))*(x[2]>log(phi)))))
log_phi_unadj_length1<-mean(unlist(lapply(log_phi_unadj_qs1, function(x) x[2]-x[1])))
log_phi_unadj_length1_med<-median(unlist(lapply(log_phi_unadj_qs1, function(x) x[2]-x[1])))

### palm_hom_0.4 adjusted (magnitude)
mu_adj_qs1<-lapply(adj_hom1,function(x) quantile(log(x$lambda)-exp(x$log_sig2)/2,c(0.025,0.975)))
mu_adj_cov1<-mean(unlist(lapply(mu_adj_qs1, function(x) (x[1]<mu)*(x[2]>mu))))
mu_adj_length1<-mean(unlist(lapply(mu_adj_qs1,function(x) x[2]-x[1])))
mu_adj_length1_med<-median(unlist(lapply(mu_adj_qs1,function(x) x[2]-x[1])))

log_sig2_adj_qs1<-lapply(adj_hom1,function(x) quantile(x$log_sig2,c(0.025,0.975)))
log_sig2_adj_cov1<-mean(unlist(lapply(log_sig2_adj_qs1, function(x) (x[1]<log(sig2))*(x[2]>log(sig2)))))
log_sig2_adj_length1<-mean(unlist(lapply(log_sig2_adj_qs1, function(x) x[2]-x[1])))
log_sig2_adj_length1_med<-median(unlist(lapply(log_sig2_adj_qs1, function(x) x[2]-x[1])))

log_phi_adj_qs1<-lapply(adj_hom1,function(x) quantile(x$log_phi,c(0.025,0.975)))
log_phi_adj_cov1<-mean(unlist(lapply(log_phi_adj_qs1, function(x) (x[1]<log(phi))*(x[2]>log(phi)))))
log_phi_adj_length1<-mean(unlist(lapply(log_phi_adj_qs1, function(x) x[2]-x[1])))
log_phi_adj_length1_med<-median(unlist(lapply(log_phi_adj_qs1, function(x) x[2]-x[1])))

### palm_hom_0.4 adjusted (GPC)
mu_adj_qs2<-lapply(adj_hom2,function(x) quantile(x$mu,c(0.025,0.975)))
mu_adj_cov2<-mean(unlist(lapply(mu_adj_qs2, function(x) (x[1]<mu)*(x[2]>mu))))
mu_adj_length2<-mean(unlist(lapply(mu_adj_qs2,function(x) x[2]-x[1])))
mu_adj_length2_med<-median(unlist(lapply(mu_adj_qs2,function(x) x[2]-x[1])))

log_sig2_adj_qs2<-lapply(adj_hom2,function(x) quantile(x$log_sig2,c(0.025,0.975)))
log_sig2_adj_cov2<-mean(unlist(lapply(log_sig2_adj_qs2, function(x) (x[1]<log(sig2))*(x[2]>log(sig2)))))
log_sig2_adj_length2<-mean(unlist(lapply(log_sig2_adj_qs2, function(x) x[2]-x[1])))
log_sig2_adj_length2_med<-median(unlist(lapply(log_sig2_adj_qs2, function(x) x[2]-x[1])))

log_phi_adj_qs2<-lapply(adj_hom2,function(x) quantile(x$log_phi,c(0.025,0.975)))
log_phi_adj_cov2<-mean(unlist(lapply(log_phi_adj_qs2, function(x) (x[1]<log(phi))*(x[2]>log(phi)))))
log_phi_adj_length2<-mean(unlist(lapply(log_phi_adj_qs2, function(x) x[2]-x[1])))
log_phi_adj_length2_med<-median(unlist(lapply(log_phi_adj_qs2, function(x) x[2]-x[1])))

###### EMPIRICAL
mu_unadj_qs2<-lapply(palm_hom_sims_empirical0.4,function(x) quantile(log(x$lambda)-exp(x$log_sig2)/2,c(0.025,0.975)))
mu_unadj_cov2<-mean(unlist(lapply(mu_unadj_qs2, function(x) (x[1]<mu)*(x[2]>mu))))
mu_unadj_length2<-mean(unlist(lapply(mu_unadj_qs2,function(x) x[2]-x[1])))
mu_unadj_length2_med<-median(unlist(lapply(mu_unadj_qs2,function(x) x[2]-x[1])))

log_sig2_unadj_qs2<-lapply(palm_hom_sims_empirical0.4,function(x) quantile(x$log_sig2,c(0.025,0.975)))
log_sig2_unadj_cov2<-mean(unlist(lapply(log_sig2_unadj_qs2, function(x) (x[1]<log(sig2))*(x[2]>log(sig2)))))
log_sig2_unadj_length2<-mean(unlist(lapply(log_sig2_unadj_qs2, function(x) x[2]-x[1])))
log_sig2_unadj_length2_med<-median(unlist(lapply(log_sig2_unadj_qs2, function(x) x[2]-x[1])))

log_phi_unadj_qs2<-lapply(palm_hom_sims_empirical0.4,function(x) quantile(x$log_phi,c(0.025,0.975)))
log_phi_unadj_cov2<-mean(unlist(lapply(log_phi_unadj_qs2, function(x) (x[1]<log(phi))*(x[2]>log(phi)))))
log_phi_unadj_length2<-mean(unlist(lapply(log_phi_unadj_qs2, function(x) x[2]-x[1])))
log_phi_unadj_length2_med<-median(unlist(lapply(log_phi_unadj_qs2, function(x) x[2]-x[1])))

### palm_hom_empirical0.4 adjusted (magnitude)
mu_adj_qs3<-lapply(adj_hom_empirical1,function(x) quantile(log(x$lambda)-exp(x$log_sig2)/2,c(0.025,0.975)))
mu_adj_cov3<-mean(unlist(lapply(mu_adj_qs3, function(x) (x[1]<mu)*(x[2]>mu))))
mu_adj_length3<-mean(unlist(lapply(mu_adj_qs3,function(x) x[2]-x[1])))
mu_adj_length3_med<-median(unlist(lapply(mu_adj_qs3,function(x) x[2]-x[1])))

log_sig2_adj_qs3<-lapply(adj_hom_empirical1,function(x) quantile(x$log_sig2,c(0.025,0.975)))
log_sig2_adj_cov3<-mean(unlist(lapply(log_sig2_adj_qs3, function(x) (x[1]<log(sig2))*(x[2]>log(sig2)))))
log_sig2_adj_length3<-mean(unlist(lapply(log_sig2_adj_qs3, function(x) x[2]-x[1])))
log_sig2_adj_length3_med<-median(unlist(lapply(log_sig2_adj_qs3, function(x) x[2]-x[1])))

log_phi_adj_qs3<-lapply(adj_hom_empirical1,function(x) quantile(x$log_phi,c(0.025,0.975)))
log_phi_adj_cov3<-mean(unlist(lapply(log_phi_adj_qs3, function(x) (x[1]<log(phi))*(x[2]>log(phi)))))
log_phi_adj_length3<-mean(unlist(lapply(log_phi_adj_qs3, function(x) x[2]-x[1])))
log_phi_adj_length3_med<-median(unlist(lapply(log_phi_adj_qs3, function(x) x[2]-x[1])))

### palm_hom_0.4 adjusted (GPC)
mu_adj_qs4<-lapply(adj_hom_empirical2,function(x) quantile(x$mu,c(0.025,0.975)))
mu_adj_cov4<-mean(unlist(lapply(mu_adj_qs4, function(x) (x[1]<mu)*(x[2]>mu))))
mu_adj_length4<-mean(unlist(lapply(mu_adj_qs4,function(x) x[2]-x[1])))
mu_adj_length4_med<-median(unlist(lapply(mu_adj_qs4,function(x) x[2]-x[1])))

log_sig2_adj_qs4<-lapply(adj_hom_empirical2,function(x) quantile(x$log_sig2,c(0.025,0.975)))
log_sig2_adj_cov4<-mean(unlist(lapply(log_sig2_adj_qs4, function(x) (x[1]<log(sig2))*(x[2]>log(sig2)))))
log_sig2_adj_length4<-mean(unlist(lapply(log_sig2_adj_qs4, function(x) x[2]-x[1])))
log_sig2_adj_length4_med<-median(unlist(lapply(log_sig2_adj_qs4, function(x) x[2]-x[1])))

log_phi_adj_qs4<-lapply(adj_hom_empirical2,function(x) quantile(x$log_phi,c(0.025,0.975)))
log_phi_adj_cov4<-mean(unlist(lapply(log_phi_adj_qs4, function(x) (x[1]<log(phi))*(x[2]>log(phi)))))
log_phi_adj_length4<-mean(unlist(lapply(log_phi_adj_qs4, function(x) x[2]-x[1])))
log_phi_adj_length4_med<-median(unlist(lapply(log_phi_adj_qs4, function(x) x[2]-x[1])))



kable(round(rbind(
c(mu_full_cov,mu_unadj_cov1,mu_unadj_cov2,mu_adj_cov1,mu_adj_cov2,mu_adj_cov3,mu_adj_cov4),
c(log_sig2_full_cov,log_sig2_unadj_cov1,log_sig2_unadj_cov2,log_sig2_adj_cov1,log_sig2_adj_cov2,log_sig2_adj_cov3,log_sig2_adj_cov4),
c(log_phi_full_cov,log_phi_unadj_cov1,log_phi_unadj_cov2,log_phi_adj_cov1,log_phi_adj_cov2,log_phi_adj_cov3,log_phi_adj_cov4)
),2),format="latex")

kable(round(rbind(
  c(mu_full_length,mu_unadj_length1,mu_unadj_length2,mu_adj_length1,mu_adj_length2,mu_adj_length3,mu_adj_length4),
  c(log_sig2_full_length,log_sig2_unadj_length1,log_sig2_unadj_length2,log_sig2_adj_length1,log_sig2_adj_length2,log_sig2_adj_length3,log_sig2_adj_length4),
  c(log_phi_full_length,log_phi_unadj_length1,log_phi_unadj_length2,log_phi_adj_length1,log_phi_adj_length2,log_phi_adj_length3,log_phi_adj_length4)
),2),format="latex")

kable(round(rbind(
  c(mu_full_length_med,mu_unadj_length1_med,mu_unadj_length2_med,mu_adj_length1_med,mu_adj_length2_med,mu_adj_length3_med,mu_adj_length4_med),
  c(log_sig2_full_length_med,log_sig2_unadj_length1_med,log_sig2_unadj_length2_med,log_sig2_adj_length1_med,log_sig2_adj_length2_med,log_sig2_adj_length3_med,log_sig2_adj_length4_med),
  c(log_phi_full_length_med,log_phi_unadj_length1_med,log_phi_unadj_length2_med,log_phi_adj_length1_med,log_phi_adj_length2_med,log_phi_adj_length3_med,log_phi_adj_length4_med)
),2),format="latex")


######### PLOTS
par(mfrow=c(3,3))
for(k in 1:100){
  plot(density(log(palm_hom_sims_empirical0.4[[k]]$lambda)-exp(palm_hom_sims_empirical0.4[[k]]$log_sig2)/2),xlim=c(3,6),lty=3,col="red",main=k)
  lines(density(adj_hom_empirical2[[k]]$mu),col="red")
  lines(density(adj_hom2[[k]]$mu),col="blue")
  lines(density(log(palm_hom_sims0.4[[k]]$lambda)-exp(palm_hom_sims0.4[[k]]$log_sig2)/2),lty=3,col="blue")
  lines(density(full_hom_sims[[k]]$beta))
  
  plot(density(palm_hom_sims_empirical0.4[[k]]$log_sig2),xlim=c(-1,1),col="red",lty=3)
  lines(density(adj_hom_empirical2[[k]]$log_sig2),col="red")
  lines(density(palm_hom_sims0.4[[k]]$log_sig2),lty=3,col="blue")
  lines(density(adj_hom2[[k]]$log_sig2),col="blue")
  lines(density(full_hom_sims[[k]]$log_sig2))
  
  plot(density(palm_hom_sims_empirical0.4[[k]]$log_phi),xlim=c(-4,-1),col="red",lty=3)
  lines(density(adj_hom_empirical2[[k]]$log_phi),col="red")
  lines(density(palm_hom_sims0.4[[k]]$log_phi),lty=3,col="blue")
  lines(density(adj_hom2[[k]]$log_phi),col="blue")
  lines(density(full_hom_sims[[k]]$log_phi))
}


pdf("example_posterior_problem.pdf",width=15,height=5)
par(mfrow=c(1,3))
k<-83
plot(density(log(palm_hom_sims0.4[[k]]$lambda)-exp(palm_hom_sims0.4[[k]]$log_sig2)/2),xlim=c(3,6),lty=3,col="blue",main=expression(beta[0]),lwd=2,cex.main=2,cex.lab=1.5,xlab="")
lines(density(full_hom_sims[[k]]$beta),lwd=2)

plot(density(palm_hom_sims0.4[[k]]$log_sig2),xlim=c(-1,1),col="blue",lty=3,main=expression(log(sigma^2)),lwd=2,cex.main=2,cex.lab=1.5,xlab="")
lines(density(full_hom_sims[[k]]$log_sig2),lwd=2)

plot(density(palm_hom_sims0.4[[k]]$log_phi),xlim=c(-4,-1),col="blue",lty=3,main=expression(log(phi)),lwd=2,cex.main=2,cex.lab=1.5,xlab="")
lines(density(full_hom_sims[[k]]$log_phi),lwd=2)
dev.off()

pdf("example_adjustment.pdf",width=15,height=5)
par(mfrow=c(1,3))
k<-83
plot(density(full_hom_sims[[k]]$beta),xlim=c(3,6),main=expression(beta[0]),lwd=2,cex.main=2,cex.lab=1.5,xlab="")
lines(density(adj_hom_empirical2[[k]]$mu),col="red",lwd=2)
lines(density(adj_hom2[[k]]$mu),col="blue",lwd=2)

plot(density(full_hom_sims[[k]]$log_sig2),xlim=c(-1,1),main=expression(log(sigma^2)),lwd=2,cex.main=2,cex.lab=1.5,xlab="")
lines(density(adj_hom_empirical2[[k]]$log_sig2),col="red",lwd=2)
lines(density(adj_hom2[[k]]$log_sig2),col="blue",lwd=2)

plot(density(full_hom_sims[[k]]$log_phi),xlim=c(-5,0),main=expression(log(phi)),lwd=2,cex.main=2,cex.lab=1.5,xlab="")
lines(density(adj_hom_empirical2[[k]]$log_phi),col="red",lwd=2)
lines(density(adj_hom2[[k]]$log_phi),col="blue",lwd=2)

dev.off()
