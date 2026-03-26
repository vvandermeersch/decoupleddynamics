
# Collection data into list
N_days <- length(nee_obs)
M <- 30
data <- mget(c('N_days', 'N_stands', 'N_stand_years', 
               'N_years_perstand', 'N_days_peryear', 
               'stand_start_idxs', 'stand_end_idxs',
               'year_padded_start_idxs', 'year_padded_end_idxs',  
               'nee_obs', 'M'))

data$N_years_perstand <- array(N_years_perstand, dim = 1)
data$stand_start_idxs <- array(stand_start_idxs, dim = 1)
data$stand_end_idxs <- array(stand_end_idxs, dim = 1)

modelstan <- stan_model(file.path(wd, 'stan/growth', 'simplemodel_onlyGSL.stan'))
fit <- sampling(modelstan, data, chains = 4, cores = 4, 
                seed = 12345, iter = 2000, warmup = 1000)
samples <- util$extract_expectand_vals(fit)

par(mfrow = c(2,1), cex.main = 0.8)
idxs <- seq(1, 365, 1)
util$plot_conn_pushforward_quantiles(samples, paste0('f[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'NEE')
idxs <- seq(1, 365+15, 1)
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "white", cex = 0.5)
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "black", cex = 0.25)
abline(h = 0, lwd = 1, col = util $c_mid_teal, lty = 2)
util$plot_conn_pushforward_quantiles(samples, paste0('g[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'NEE')


# Padding +/-15days before and after
N_days <- length(nee_obs)
M <- 30
data <- mget(c('N_days', 'N_stands', 'N_stand_years', 
               'N_years_perstand', 'N_days_peryear', 
               'stand_start_idxs', 'stand_end_idxs',
               'year_padded_start_idxs', 'year_padded_end_idxs',  
               'year_true_start_idxs', 'year_true_end_idxs',  
               'nee_obs', 'M'))
data$N_years_perstand <- array(N_years_perstand, dim = 1)
data$stand_start_idxs <- array(stand_start_idxs, dim = 1)
data$stand_end_idxs <- array(stand_end_idxs, dim = 1)

modelstan <- stan_model(file.path(wd, 'stan/growth', 'simplemodel_onlyGSL_paddeddata.stan'))
fit <- sampling(modelstan, data, chains = 4, cores = 4, 
                seed = 12345, iter = 2000, warmup = 1000)
samples <- util$extract_expectand_vals(fit)
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


N_days <- length(nee_obs)
M <- 30
N_rings <- length(log_rw_obs)
log_rw_obs <- log(rw_obs)
data <- mget(c('N_days', 'N_stands', 'N_stand_years', 
               'N_years_perstand', 'N_days_peryear', 
               'stand_start_idxs', 'stand_end_idxs',
               'year_start_idxs', 'year_end_idxs',  
               'nee_obs', 'M',
               
               'N_rings', 'N_trees', 'N_all_years', 'stand_idxs',
               "N_years", 'all_years', 'years', 'all_years_idxs',
               'tree_start_idxs', 'tree_end_idxs', 'log_rw_obs'
               
               ))

data$N_years_perstand <- array(N_years_perstand, dim = 1)
data$stand_start_idxs <- array(stand_start_idxs, dim = 1)
data$stand_end_idxs <- array(stand_end_idxs, dim = 1)

modelstan <- stan_model(file.path(wd, 'stan/growth', 'simplemodel_growth_GSL.stan'))
fit <- sampling(modelstan, data, chains = 4, cores = 4, 
                seed = 12345, iter = 2000, warmup = 1000)
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha_nee', 'rho_nee', 'gamma_nee', 
                                         'sigma_nee', 'f_tilde',
                                         
                                         'alpha', 'beta_gsl', 'rho_sp', 'gamma_sp',
                                         'f_tilde_sh', 'rho_sh', 'gamma_sh',
                                         'sigma'), check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples)

par(mfrow = c(1,2))
util$plot_expectand_pushforward(samples[['rho_sh']], 30, display_name = bquote('Stand-level GP'~rho), flim = c(0,15))
util$plot_expectand_pushforward(samples[['rho_sp']], 30, display_name = bquote('Tree-level GP'~rho), flim = c(0,15))

par(mar = c(4,4,1,1))
for(t in 1:data$N_trees){
  
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- paste0("log_rw_pred[",idxs,"]")
  
  util$plot_conn_pushforward_quantiles(samples, names, data$years[idxs],
                                       xlab="Year", ylab="Log ring width (per mm)", 
                                       display_ylim=c(-2, 2), display_xlim = range(data$all_years))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.5, col="black")

}



par(mfrow = c(1,2), cex.main = 0.8)
util$plot_disc_pushforward_quantiles(samples, paste0('gsl_est[', 1:data$N_stand_years, ']'), 
                                     xlab = 'Year', ylab = 'GSL (weeks)', main = bquote(kappa~'= 0.05'),
                                     display_ylim = c(16, 24), xticklabs = data$all_years)
util$plot_expectand_pushforward(samples[['beta_gsl']], 20, display_name = bquote(beta[GSL]~'(1/week)'))


par(mfrow = c(2,1), cex.main = 0.8)
idxs <- seq(365*5+1, 365*8, 1)
util$plot_conn_pushforward_quantiles(samples, paste0('f[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'NEE')
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "white", cex = 0.5)
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "black", cex = 0.25)
abline(h = 0, lwd = 1, col = util $c_mid_teal, lty = 2)
util$plot_conn_pushforward_quantiles(samples, paste0('g[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'NEE')



par(mfrow = c(2,1), cex.main = 0.8)
idxs <- seq(1, 365, 1)
util$plot_conn_pushforward_quantiles(samples, paste0('f[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'NEE')
idxs <- seq(1, 365+15, 1)
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "white", cex = 0.5)
points(x = idxs, y = data$nee_obs[idxs],
       pch = 20, col = "black", cex = 0.25)
abline(h = 0, lwd = 1, col = util $c_mid_teal, lty = 2)
util$plot_conn_pushforward_quantiles(samples, paste0('g[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'NEE')






modelstan <- stan_model(file.path(wd, 'stan/growth', 'simplemodel_growth_GSL_yearlyalphas.stan'))
fit <- sampling(modelstan, data, chains = 4, cores = 4, 
                seed = 12345, iter = 2000, warmup = 1000)
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
util$plot_div_pairs(paste0('alpha_nee[',1:data$N_stand_years,']'), 'sigma_alpha_nee', samples, diagnostics,
                    transforms = list('sigma_alpha_nee' = 1))
