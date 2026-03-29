
# Collect data in list
N_days <- length(nee_obs)
M <- 30
log_rw_obs <- log(rw_obs+1e-6)
N_rings <- length(log_rw_obs)
data <- mget(c('N_days', 'N_stands', 'N_stand_years', 
               'N_years_perstand', 'N_days_peryear', 
               'stand_start_idxs', 'stand_end_idxs',
               'stand_years', 
               'year_padded_start_idxs', 'year_padded_end_idxs',  
               'year_true_start_idxs', 'year_true_end_idxs',  
               'nee_obs', 'M',
               
               'N_rings', 'N_trees', 'N_all_years', 'stand_idxs',
               "N_years", 'all_years', 'years', 'all_years_idxs',
               'tree_start_idxs', 'tree_end_idxs', 'log_rw_obs'
               
))

# data$N_years_perstand <- array(N_years_perstand, dim = 1)
# data$stand_start_idxs <- array(stand_start_idxs, dim = 1)
# data$stand_end_idxs <- array(stand_end_idxs, dim = 1)


# ---------------
# 1. Fixed kappa
# ---------------

# Posterior quantification
modelstan <- stan_model(file.path(wd, 'stan/growth', 'simplemodel_growth_GSL_paddeddata.stan'))
fit <- sampling(modelstan, data, chains = 4, cores = 4, 
                seed = 12345, iter = 2000, warmup = 1000)
saveRDS(fit, file.path(wd, 'output/fit', 'QUAL_5stands.rds'))
samples <- util$extract_expectand_vals(fit)

# Retrodictive checks
par(mar = c(4,4,1,1), mfrow = c(3,2))
for(t in sample(1:data$N_trees,12)){
  
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")
  
  util$plot_conn_pushforward_quantiles(samples, names, data$years[idxs],
                                       xlab="Year", ylab="Log ring width (per mm)", 
                                       display_ylim=c(-2, 2), display_xlim = range(data$all_years))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.5, col="black")
  
}

# Posterior inference
par(mfrow = c(1,2))
util$plot_expectand_pushforward(samples[['rho_sh']], 30, display_name = bquote('Stand-level GP'~rho), flim = c(0,20))
util$plot_expectand_pushforward(samples[['rho_sp']], 30, display_name = bquote('Tree-level GP'~rho), flim = c(0,20))


par(mfrow = c(2,3), cex.main = 0.8)
for(s in 1:data$N_stands){
  idxs <- stand_start_idxs[s]:stand_end_idxs[s]
  util$plot_disc_pushforward_quantiles(samples, paste0('gsl_est[', idxs, ']'), 
                                       xlab = 'Year', ylab = 'GSL (weeks)', main = bquote(kappa~'= 0.05'),
                                       display_ylim = c(22, 32))
  
}

util$plot_expectand_pushforward(samples[['beta_gsl']], 20, display_name = bquote(beta[GSL]~'(1/week)'))


par(mfrow = c(2,1), cex.main = 0.8)
i <- 6
idxs <- data$year_padded_start_idxs[i]:data$year_padded_end_idxs[i ]
util$plot_conn_pushforward_quantiles(samples, paste0('f[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'f(t)')
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "white", cex = 0.5)
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "black", cex = 0.25)
abline(h = 0, lwd = 1.5, col = util $c_mid_teal, lty = 2)
abline(v = c(data$year_true_start_idxs[i], data$year_true_end_idxs[i]), lwd = 1, col = util $c_mid_teal, lty = 2)
util$plot_conn_pushforward_quantiles(samples, paste0('g_sharp[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'g(t)')
abline(v = c(data$year_true_start_idxs[i], data$year_true_end_idxs[i]), lwd = 1, col = util $c_mid_teal, lty = 2)





# --------------
# 2. Free kappa
# --------------

# Posterior quantification
modelstan <- stan_model(file.path(wd, 'stan/growth', 'simplemodel_growth_GSL_paddeddata_freekappa.stan'))
fit <- sampling(modelstan, data, chains = 4, cores = 4, 
                seed = 12345, iter = 2000, warmup = 1000)
saveRDS(fit, file.path(wd, 'output/fit', 'QUAL_5stands_freekappa_differentgsl0_2.rds'))
samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha_nee', 'rho_nee', 'gamma_nee', 'sigma_nee',
                                         'kappa', 'f_tilde',
                                         'alpha', 'beta_gsl', 'rho_sp', 'gamma_sp',
                                         'f_tilde_sh', 'rho_sh', 'gamma_sh', 'sigma'),
                                       check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples)

util$plot_pairs_by_chain(samples[['alpha']], 'alpha',
                         samples[['beta_gsl']], 'beta_gsl')

util$plot_pairs_by_chain(samples[['gamma_sh']], 'gamma_sh',
                         samples[['beta_gsl']], 'beta_gsl')

par(mar = c(4,4,1,1), mfrow = c(1,2))
util$plot_expectand_pushforward(samples[['alpha']], 40, display_name = bquote(alpha),
                                ylim = c(0,2), flim = c(-2,2))
prior <- rnorm(1e6,  0, 0.69)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples[['beta_gsl']], 30, display_name =  bquote(beta[GSL]~'(1/week)'),
                                ylim = c(0,30), flim = c(-0.4,0.4))
prior <- rnorm(1e6,  0, log(1.8) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)

# Retrodictive checks
par(mar = c(4,4,1,1), mfrow = c(3,2))
for(t in sample(1:data$N_trees,12)){
  
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")
  
  util$plot_conn_pushforward_quantiles(samples, names, data$years[idxs],
                                       xlab="Year", ylab="Log ring width (per mm)", 
                                       display_ylim=c(-2, 2), display_xlim = range(data$all_years))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.5, col="black")
  
}

# Posterior inference
par(mfrow = c(1,3), cex.main = 0.8)
util$plot_expectand_pushforward(samples[['kappa']], 20, display_name = bquote('Sharpness'~kappa),
                                ylim = c(0,10))
prior <- rnorm(1e6, 0, 0.6/2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
par(mfrow = c(2,3), cex.main = 0.8)
for(s in 1:data$N_stands){
  idxs <- stand_start_idxs[s]:stand_end_idxs[s]
  util$plot_disc_pushforward_quantiles(samples, paste0('gsl_est[', idxs, ']'), 
                                       xlab = 'Year', ylab = 'GSL (weeks)',
                                       display_ylim = c(24, 32))
  
}
util$plot_expectand_pushforward(samples[['beta_gsl']], 20, display_name = bquote(beta[GSL]~'(1/week)'), flim = c(-0.3,0.3))
prior <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)

par(mfrow = c(2,1), cex.main = 0.8)
i <- 40
idxs <- data$year_padded_start_idxs[i]:data$year_padded_end_idxs[i ]
util$plot_conn_pushforward_quantiles(samples, paste0('f[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'f(t)')
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "white", cex = 0.5)
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "black", cex = 0.25)
abline(h = 0, lwd = 1.5, col = util $c_mid_teal, lty = 2)
abline(v = c(data$year_true_start_idxs[i], data$year_true_end_idxs[i]), lwd = 1, col = util $c_mid_teal, lty = 2)
util$plot_conn_pushforward_quantiles(samples, paste0('g[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'g(t)')
abline(v = c(data$year_true_start_idxs[i], data$year_true_end_idxs[i]), lwd = 1, col = util $c_mid_teal, lty = 2)
text(x = mean(idxs), y = 0.3, labels = 
       paste0('GSL = [', round(util$ensemble_mcmc_quantile_est(samples[[paste0('gsl_est[',i,']')]], c(0.05)),1),
              ' - ', round(util$ensemble_mcmc_quantile_est(samples[[paste0('gsl_est[',i,']')]], c(0.95)),1), ']'))




# ---------------------------
# #. Free kappa and fixed nu
# ---------------------------

# Posterior quantification
modelstan <- stan_model(file.path(wd, 'stan/growth', 'simplemodel_growth_GSL_paddeddata_freekappa_fixednu.stan'))
fit <- sampling(modelstan, data, chains = 4, cores = 4, 
                seed = 12345, iter = 2000, warmup = 1000)
saveRDS(fit, file.path(wd, 'output/fit', 'ACSH_15stands_freekappa_fixednu.rds'))
samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha_nee', 'rho_nee', 'gamma_nee', 'sigma_nee',
                                         'kappa', 
                                         # 'f_tilde',
                                         'alpha', 'beta_gsl', 'rho_sp', 'gamma_sp',
                                         # 'f_tilde_sh', 
                                         'rho_sh', 'gamma_sh', 'sigma'),
                                       check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples)

util$plot_pairs_by_chain(samples[['alpha']], 'alpha',
                         samples[['beta_gsl']], 'beta_gsl')
util$plot_pairs_by_chain(samples[['gamma_sh']], 'gamma_sh',
                         samples[['rho_sh']], 'rho_sh')
util$plot_pairs_by_chain(samples[['gamma_sp']], 'gamma_sp',
                         samples[['rho_sp']], 'rho_sp')
util$plot_pairs_by_chain(samples[['kappa']], 'kappa',
                         samples[['beta_gsl']], 'beta_gsl')
par(mfrow = c(1,1), cex.main = 0.8)
util$plot_expectand_pushforward(samples[['kappa']], 30, display_name =  bquote(kappa), flim = c(0,0.5))
prior <- rnorm(1e6,  0, 0.6/2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)

par(mfrow = c(2,1), cex.main = 0.8, mar = c(4,4,1,1))
i <- 40
idxs <- data$year_padded_start_idxs[i]:data$year_padded_end_idxs[i ]
util$plot_conn_pushforward_quantiles(samples, paste0('f[',idxs, ']'), idxs,
                                     xlab = ' ', ylab = 'f(t)')
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "white", cex = 0.5)
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "black", cex = 0.25)
abline(h = 0, lwd = 1.5, col = util $c_mid_teal, lty = 2)
abline(v = c(data$year_true_start_idxs[i], data$year_true_end_idxs[i]), lwd = 1, col = util $c_mid_teal, lty = 2)
util$plot_conn_pushforward_quantiles(samples, paste0('g[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'g(t)')
abline(v = c(data$year_true_start_idxs[i], data$year_true_end_idxs[i]), lwd = 1, col = util $c_mid_teal, lty = 2)
text(x = mean(idxs)+25, y = 0.3, labels = 
       paste0('GSL =\n[', round(util$ensemble_mcmc_quantile_est(samples[[paste0('gsl_est[',i,']')]], c(0.05)),1),
              ' - ', round(util$ensemble_mcmc_quantile_est(samples[[paste0('gsl_est[',i,']')]], c(0.95)),1), '] weeks'),
     cex = 0.8)

par(mfrow = c(3,2), cex.main = 0.8)
util$plot_expectand_pushforward(samples[['kappa']], 30, display_name =  bquote(kappa), flim = c(0,2))
prior <- rnorm(1e6,  0, 0.6/2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples[['alpha_nee']], 30, display_name =  bquote(kappa), flim = c(-2,2))
prior <- rnorm(1e6,  0, 1)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples[['rho_nee']], 30, display_name =  bquote(kappa), flim = c(0,100))
prior <- rnorm(1e6,  30, 10)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples[['rho_sh']], 30, display_name =  bquote(rho[short]), flim = c(0,30))
prior <- rlnorm(1e6,  1.7, 0.26)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)



# ---------------------------
# #. Free kappa and free nu
# ---------------------------

# Posterior quantification
modelstan <- stan_model(file.path(wd, 'stan/growth', 'simplemodel_growth_GSL_paddeddata_freekappa_freenu.stan'))
fit <- sampling(modelstan, data, chains = 4, cores = 4, 
                seed = 12345, iter = 2000, warmup = 1000)
saveRDS(fit, file.path(wd, 'output/fit', 'QUAL_11stands_freekappa_freenu_gsl0_21.rds'))
fit <- readRDS(file.path(wd, 'output/fit', 'QUAL_11stands_freekappa_freenu.rds'))
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)
samples <- util$extract_expectand_vals(fit)
gc()
base_samples <- util$filter_expectands(samples,
                                       c('alpha_nee', 'rho_nee', 'gamma_nee', 'sigma_nee',
                                         'kappa', 'nu',
                                         # 'f_tilde',
                                         'alpha', 'beta_gsl', 'rho_sp', 'gamma_sp',
                                         # 'f_tilde_sh', 
                                         'rho_sh', 'gamma_sh', 'sigma'),
                                       check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples)

util$plot_pairs_by_chain(samples[['kappa']], 'kappa',
                         samples[['nu']], 'nu')
util$plot_pairs_by_chain(samples[['alpha']], 'alpha',
                         samples[['beta_gsl']], 'beta_gsl')
util$plot_pairs_by_chain(samples[['gamma_sh']], 'gamma_sh',
                         samples[['rho_sh']], 'rho_sh')


par(mfrow = c(3,2), cex.main = 0.8)
util$plot_expectand_pushforward(samples[['kappa']], 30, display_name =  bquote(kappa), flim = c(0,2))
prior <- rnorm(1e6,  0, 0.2/2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
prior <- rnorm(1e6,  0, 0.6/2.57)
lines(density(prior), col = 'darkgreen', lwd = 1.5, lty = 1)
util$plot_expectand_pushforward(samples[['nu']], 30, display_name =  bquote(nu), flim = c(-2,0))
prior <- rnorm(1e6,  0, 1/2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples[['alpha_nee']], 30, display_name =  bquote(alpha[NEE]), flim = c(-2,2))
prior <- rnorm(1e6,  0, 1)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples[['rho_nee']], 30, display_name =  bquote(rho[NEE]), flim = c(0,60))
prior <- rnorm(1e6,  30, 10)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples[['gamma_nee']], 30, display_name =  bquote(gamma[NEE]), flim = c(0,10))
prior <- rnorm(1e6,  0, 1)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
prior <- rnorm(1e6,  0, 6/2.32)
lines(density(prior), col = 'darkgreen', lwd = 1.5, lty = 1)
util$plot_expectand_pushforward(samples[['sigma_nee']], 30, display_name =  bquote(sigma[NEE]), flim = c(0,5))
prior <- rnorm(1e6,  0, 1)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)

par(mfrow = c(3,2), cex.main = 0.8)
util$plot_expectand_pushforward(samples[['rho_sh']], 30, display_name =  bquote(rho[short]), flim = c(0,30))
prior <- rlnorm(1e6,  2.6, 0.18)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples[['gamma_sh']], 30, display_name =  bquote(gamma[short]), flim = c(0,5))
prior <- rnorm(1e6, 0,  log(3) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples[['sigma']], 30, display_name =  bquote(sigma), flim = c(0,1))
prior <- rnorm(1e6, 0, 0.095 / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)

par(mfrow = c(1,1), cex.main = 0.8)
util$plot_disc_pushforward_quantiles(samples, paste0('gsl_est[', 1:data$N_stand_years, ']'), 
                                     xlab = 'Year', ylab = 'GSL (weeks)',
                                     display_ylim = c(14.5, 29.5))


par(mfrow = c(2,1), cex.main = 0.8)
i <- 50
idxs <- data$year_padded_start_idxs[i]:data$year_padded_end_idxs[i ]
util$plot_conn_pushforward_quantiles(samples, paste0('f[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'f(t)')
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "white", cex = 0.5)
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "black", cex = 0.25)
abline(h = 0, lwd = 1.5, col = util $c_mid_teal, lty = 2)
abline(v = c(data$year_true_start_idxs[i], data$year_true_end_idxs[i]), lwd = 1, col = util $c_mid_teal, lty = 2)
util$plot_conn_pushforward_quantiles(samples, paste0('g[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'g(t)')
abline(v = c(data$year_true_start_idxs[i], data$year_true_end_idxs[i]), lwd = 1, col = util $c_mid_teal, lty = 2)
text(x = mean(idxs), y = 0.3, labels = 
       paste0('GSL = [', round(util$ensemble_mcmc_quantile_est(samples[[paste0('gsl_est[',i,']')]], c(0.05)),1),
              ' - ', round(util$ensemble_mcmc_quantile_est(samples[[paste0('gsl_est[',i,']')]], c(0.95)),1), ']'))

util$plot_pairs_by_chain(samples[['kappa']], 'kappa',
                         samples[[paste0('gsl_est[',i,']')]], paste0('gsl_est[',i,']'))

par(mfrow = c(4,3), cex.main = 0.8, mgp = c(1.5,0.5,0), cex.axis = 0.8)
for(s in 1:data$N_stands){
  idxs <- stand_start_idxs[s]:stand_end_idxs[s]
  util$plot_disc_pushforward_quantiles(samples, paste0('gsl_est[', idxs, ']'), 
                                       xlab = 'Year', ylab = 'GSL (weeks)',
                                       display_ylim = c(14.5, 30.5), main = uniq_stand_ids[s])
  
}


# Retrodictive checks
par(mar = c(4,4,1,1), mfrow = c(3,2))
for(t in sample(1:data$N_trees,12)){
  
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")
  
  util$plot_conn_pushforward_quantiles(samples, names, data$years[idxs],
                                       xlab="Year", ylab="Log ring width (per mm)", 
                                       display_ylim=c(-2, 2), display_xlim = range(data$all_years))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.5, col="black")
  
}
