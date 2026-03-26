
# Collect data in list
N_days <- length(nee_obs)
M <- 30
log_rw_obs <- log(rw_obs)
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
saveRDS(fit, file.path(wd, 'output/fit', 'QUAL_5stands_freekappa.rds'))
samples <- util$extract_expectand_vals(fit)
base_samples <- 
util$check_all_expectand_diagnostics()


real alpha_nee;
real<lower=0> rho_nee;
real<lower=0> gamma_nee;
real<lower=0> sigma_nee;
real<lower=0> kappa;
array[N_stand_years] vector[M] f_tilde;


// ---------------------
  // 2. Tree growth model
// ---------------------
  
  real alpha; // log ring width baseline

real beta_gsl; // GSL slope (1/week)

real<lower=1> rho_sp; // lifetime proportional growth scale
real<lower=0> gamma_sp; // lifetime proportional growth variation

// Short-term proportional growth functional behavior
array[N_stands] vector[N_all_years] f_tilde_sh; // Non-centered functional behavior
real<lower=1> rho_sh;
real<lower=0> gamma_sh; 

real<lower=0> sigma; // Proportional measurement error


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
prior <- rlnorm(1e6,  -0.57, 0.6)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)
util$plot_disc_pushforward_quantiles(samples, paste0('gsl_est[', 1:data$N_stand_years, ']'), 
                                     xlab = 'Year', ylab = 'GSL (weeks)',
                                     display_ylim = c(25, 32), xticklabs = data$all_years)
util$plot_expectand_pushforward(samples[['beta_gsl']], 20, display_name = bquote(beta[GSL]~'(1/week)'), flim = c(-0.1,0.3))
prior <- rnorm(1e6, 0, log(1.8) / 2.57)
lines(density(prior), col = util$c_mid_teal, lwd = 2, lty = 2)

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
util$plot_conn_pushforward_quantiles(samples, paste0('g[',idxs, ']'), idxs,
                                     xlab = 'Days across 21 years', ylab = 'g(t)')
abline(v = c(data$year_true_start_idxs[i], data$year_true_end_idxs[i]), lwd = 1, col = util $c_mid_teal, lty = 2)


