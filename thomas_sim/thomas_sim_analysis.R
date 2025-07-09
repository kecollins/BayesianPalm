thomas_0.4<-readRDS('palm_thomas0.4.rds')
thomas_empirical_0.4<-readRDS('palm_empirical_thomas0.4.rds')
pp_gen<-readRDS('pp_gen_thomas.rds')

mu<-30
nu<-10
sig2<-0.05^2

### BIASES
bias_mu<-mean(unlist(lapply(thomas_0.4,function(x) mean(exp(x$log_mu)-mu))))
bias_muE<-mean(unlist(lapply(thomas_empirical_0.4,function(x) mean(exp(x$log_mu)-mu))))

bias_nu<-mean(unlist(lapply(thomas_0.4,function(x) mean(x$lambda/exp(x$log_mu)-nu))))
bias_nuE<-mean(unlist(lapply(thomas_empirical_0.4,function(x) mean(x$lambda/exp(x$log_mu)-nu))))

bias_sig2<-mean(unlist(lapply(thomas_0.4,function(x) mean(exp(x$log_sig2)-sig2))))
bias_sig2E<-mean(unlist(lapply(thomas_empirical_0.4,function(x) mean(exp(x$log_sig2)-sig2))))

### MSE
rmse_mu<-sqrt(mean(unlist(lapply(thomas_0.4,function(x) mean(exp(x$log_mu)-mu)))^2))
rmse_muE<-sqrt(mean(unlist(lapply(thomas_empirical_0.4,function(x) mean(exp(x$log_mu)-mu)))^2))

rmse_nu<-sqrt(mean(unlist(lapply(thomas_0.4,function(x) mean(x$lambda/exp(x$log_mu)-nu)))^2))
rmse_nuE<-sqrt(mean(unlist(lapply(thomas_empirical_0.4,function(x) mean(x$lambda/exp(x$log_mu)-nu)))^2))

rmse_sig2<-sqrt(mean(unlist(lapply(thomas_0.4,function(x) mean(exp(x$log_sig2)-sig2)))^2))
rmse_sig2E<-sqrt(mean(unlist(lapply(thomas_empirical_0.4,function(x) mean(exp(x$log_sig2)-sig2)))^2))


bias_rmse_data <- data.frame(
  Model = c("PL$_{0.2}$", "PLE$_{0.2}$"),
  Bias_mu = round(c(bias_mu,bias_muE),2),
  Bias_nu = round(c(bias_nu,bias_nuE),2),
  Bias_sig2 = round(c(bias_sig2,bias_sig2E),5),
  RMSE_mu = round(c(rmse_mu,rmse_muE),2),
  RMSE_nu = round(c(rmse_nu,rmse_nuE),2),
  RMSE_sig2 = round(c(rmse_sig2,rmse_sig2E),5)
)

# Generate the LaTeX table
bias_rmse_data %>%
  kable("latex", booktabs = TRUE, escape = FALSE,
        col.names = c("Model", "$\\mu$", "$\\nu2$", "$\\sigma^2$", "$\\mu$", "$\\nu$", "$\\sigma^2$"),
        caption = "Posterior mean bias and RMSE for the Thomas process across models.") %>%
  add_header_above(c(" " = 1, "Bias" = 3, "RMSE" = 3)) %>%
  kable_styling(latex_options = c("hold_position"))


pdf("thomas.pdf",width=15,height=5)
par(mfrow=c(1,3))
boxplot(unlist(lapply(thomas_0.4,function(x) mean(exp(x$log_mu)))),
        unlist(lapply(thomas_empirical_0.4,function(x) mean(exp(x$log_mu)))),
        main=expression(paste("Posterior means of ",mu)),
        names=c("PL","PLE"))
abline(h=mu,col="blue",lwd=3)
boxplot(unlist(lapply(thomas_0.4,function(x) mean(x$lambda/exp(x$log_mu)))),
        unlist(lapply(thomas_empirical_0.4,function(x) mean(x$lambda/exp(x$log_mu)))),
        main=expression(paste("Posterior means of ",nu)),
        names=c("PL","PLE"))
abline(h=nu,col="blue",lwd=3)
boxplot(unlist(lapply(thomas_0.4,function(x) mean(exp(x$log_sig2)))),
        unlist(lapply(thomas_empirical_0.4,function(x) mean(exp(x$log_sig2)))),
        main=expression(paste("Posterior means of ",sigma^2)),
        names=c("PL","PLE"))
abline(h=sig2,col="blue",lwd=3)
dev.off()



#### COVERAGES
thomas_adj1<-readRDS('thomas_adj_hom1.rds')
thomas_empirical_adj1<-readRDS('thomas_adj_hom_empirical1.rds')
thomas_adj2<-readRDS('adj_thomas2.rds')
thomas_empirical_adj2<-readRDS('adj_empirical_thomas2.rds')

nu_qs<-lapply(thomas_adj1,function(x) quantile(x$lambda/exp(x$log_mu),c(0.025,0.975)))
mu_qs<-lapply(thomas_adj1,function(x) quantile(x$log_mu,c(0.025,0.975)))
sig2_qs<-lapply(thomas_adj1,function(x) quantile(x$log_sig2,c(0.025,0.975)))

mean(unlist(lapply(nu_qs,function(x) (x[1]<(nu)) *(x[2]>(nu)))))
mean(unlist(lapply(mu_qs,function(x) (x[1]<log(mu)) *(x[2]>log(mu)))))
mean(unlist(lapply(sig2_qs,function(x) (x[1]<log(sig2)) *(x[2]>log(sig2)))))

nu_qsE<-lapply(thomas_empirical_adj1,function(x) quantile(x$lambda/exp(x$log_mu),c(0.025,0.975)))
mu_qsE<-lapply(thomas_empirical_adj1,function(x) quantile(x$log_mu,c(0.025,0.975)))
sig2_qsE<-lapply(thomas_empirical_adj1,function(x) quantile(x$log_sig2,c(0.025,0.975)))

mean(unlist(lapply(nu_qsE,function(x) (x[1]<log(nu)) *(x[2]>log(nu)))))
mean(unlist(lapply(mu_qsE,function(x) (x[1]<log(mu)) *(x[2]>log(mu)))))
mean(unlist(lapply(sig2_qsE,function(x) (x[1]<log(sig2)) *(x[2]>log(sig2)))))


nu_qs2<-lapply(thomas_adj2,function(x) quantile(x$log_nu,c(0.025,0.975)))
mu_qs2<-lapply(thomas_adj2,function(x) quantile(x$log_mu,c(0.025,0.975)))
sig2_qs2<-lapply(thomas_adj2,function(x) quantile(x$log_sig2,c(0.025,0.975)))

mean(unlist(lapply(nu_qs2,function(x) (x[1]<log(nu)) *(x[2]>log(nu)))))
mean(unlist(lapply(mu_qs2,function(x) (x[1]<log(mu)) *(x[2]>log(mu)))))
mean(unlist(lapply(sig2_qs2,function(x) (x[1]<log(sig2)) *(x[2]>log(sig2)))))

nu_qs2E<-lapply(thomas_empirical_adj2,function(x) quantile(x$log_nu,c(0.025,0.975)))
mu_qs2E<-lapply(thomas_empirical_adj2,function(x) quantile(x$log_mu,c(0.025,0.975)))
sig2_qs2E<-lapply(thomas_empirical_adj2,function(x) quantile(x$log_sig2,c(0.025,0.975)))

mean(unlist(lapply(nu_qs2E,function(x) (x[1]<log(nu)) *(x[2]>log(nu)))))
mean(unlist(lapply(mu_qs2E,function(x) (x[1]<log(mu)) *(x[2]>log(mu)))))
mean(unlist(lapply(sig2_qs2E,function(x) (x[1]<log(sig2)) *(x[2]>log(sig2)))))

#undadjusted coverage
nu_qs_full<-lapply(thomas_0.4, function(x) quantile(x$lambda/exp(x$log_mu),c(0.025,0.975)))
mu_qs_full<-lapply(thomas_0.4, function(x) quantile(x$log_mu,c(0.025,0.975)))
sig2_qs_full<-lapply(thomas_0.4, function(x) quantile(x$log_sig2,c(0.025,0.975)))

mean(unlist(lapply(nu_qs_full,function(x) (x[1]<(nu)) *(x[2]>(nu)))))
mean(unlist(lapply(mu_qs_full,function(x) (x[1]<log(mu)) *(x[2]>log(mu)))))
mean(unlist(lapply(sig2_qs_full,function(x) (x[1]<log(sig2)) *(x[2]>log(sig2)))))

nu_qs_fullE<-lapply(thomas_empirical_0.4, function(x) quantile(x$lambda/exp(x$log_mu),c(0.025,0.975)))
mu_qs_fullE<-lapply(thomas_empirical_0.4, function(x) quantile(x$log_mu,c(0.025,0.975)))
sig2_qs_fullE<-lapply(thomas_empirical_0.4, function(x) quantile(x$log_sig2,c(0.025,0.975)))

mean(unlist(lapply(nu_qs_fullE,function(x) (x[1]<(nu)) *(x[2]>(nu)))))
mean(unlist(lapply(mu_qs_fullE,function(x) (x[1]<log(mu)) *(x[2]>log(mu)))))
mean(unlist(lapply(sig2_qs_fullE,function(x) (x[1]<log(sig2)) *(x[2]>log(sig2)))))
