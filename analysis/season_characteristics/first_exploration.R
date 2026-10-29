library(rstan)

# 1. A simple Gaussian (symetrical growing season)

## 1.1. Plot the curve
gaussian <- function(x, gpp0, gamma, peakday, width) {
  gpp0 + gamma * exp(-0.5 * ((x - peakday) / width)^2)
}
g_sym  <- gaussian(1:365, gpp0 = 0.1, gamma = 7, peakday = 200, width = 55)
r <- (gpp_threshold-0.1)/7        
z <- sqrt(-2*log(r))  
SOS <- 200-55*z
EOS <- 200+55*z

par(mfrow = c(1,1), mar = c(3.5,3.5,1,1), cex.lab = 0.8, cex.axis = 0.8, mgp = c(2.5,1,0))
plot(NULL, type = "l", lwd = 2, xlim = c(1, 365), ylim = c(0, 8),
     xlab = "DOY", ylab = "GPP (gC/m²/day)", yaxs = "i", xaxt = 'n', bty = 'n')
axis(1, at = seq(0,365,60))
abline(h = gpp_threshold, lty = 2, col = "grey40")   # threshold
polygon(c(SOS, EOS, EOS, SOS), c(0, 0, gpp_threshold, gpp_threshold),
        col = paste0(util$c_light_teal, '30'), border = NA)
x <- seq(SOS, EOS, length.out = 1e4)
y <- gaussian(x, gpp0 = 0.1, gamma = 7, peakday = 200, width = 55)
polygon(c(x, rev(x)), c(y, rep(gpp_threshold, 1e4)),
        col = paste0(util$c_light_teal, '50'), border = NA) 
segments(x0 = SOS, y0 = 0, y1 = gpp_threshold, lty = 2, col = util$c_mid_teal)
segments(x0 = EOS, y0 = 0, y1 = gpp_threshold, lty = 2, col = util$c_mid_teal)
lines(doy, g_sym, lwd = 2, col = util$c_mid_teal)

## 1.2. Posteriro quantification
model1 <- stan_model(file.path(wd, 'stan', 'gs_charac', 'model1_simplegaussian.stan'))

init_fun <- function() {

  stand_idxs <- rep(1:gpp_data$N_stands, times = gpp_data$N_years_perstand)
  
  mu_peakday <- rnorm(gpp_data$N_stands, 180, 80 / 2.57)
  mu_log_gamma <- rnorm(gpp_data$N_stands, log(5), 0.75)
  mu_log_width <- rnorm(gpp_data$N_stands, log(40), log(2) / 2.57)
  
  tau_peakday <- abs(rnorm(1, 0, 10))
  tau_log_gamma <- abs(rnorm(1, 0, 0.25))
  tau_log_width <- abs(rnorm(1, 0, 0.25))
  
  list(
    gpp0  = rnorm(gpp_data$N_stands, 0, 1),
    sigma = abs(rnorm(gpp_data$N_stands, 0, 1)),
    mu_peakday = mu_peakday, mu_log_gamma = mu_log_gamma, mu_log_width = mu_log_width,
    tau_peakday = tau_peakday, tau_log_gamma = tau_log_gamma, tau_log_width = tau_log_width,
    peakday   = rnorm(gpp_data$N_years, mu_peakday[stand_idxs], tau_peakday),
    log_gamma = rnorm(gpp_data$N_years, mu_log_gamma[stand_idxs], tau_log_gamma),
    log_width = rnorm(gpp_data$N_years, mu_log_width[stand_idxs], tau_log_width)
  )
}

gpp_data$gpp_thresh <- 1
fit1 <- sampling(model1, gpp_data, chains = 4, cores = 4, init = init_fun)
samples1 <- util$extract_expectand_vals(fit1)

## 1.3. Retrodictive checks
sy <- 58
s <- 1 + 365 * (sy - 1);
e <- s + 364
names <- paste0('gpp_pred[', s:e, ']')
par(mfrow = c(1,1), mar = c(4,4,1,1), cex.lab = 0.8, cex.axis = 0.8)
util$plot_conn_pushforward_quantiles(samples1, names, plot_xs = 1:365, ylab = 'GPP (gC/m²/day)', xlab = 'DOY')
points(x = 1:365, y = gpp_data$gpp_obs[s:e], col = 'white', cex = 1.2, pch = 20)
points(x = 1:365, y = gpp_data$gpp_obs[s:e], col = 'black', cex = 0.5, pch = 20)

## 1.4. Posterior inference
par(mfrow = c(4,1), mar = c(4,4,1,1), cex.lab = 0.8, cex.axis = 0.8)
for(s in 1:gpp_data$N_stands){
  idxs <- gpp_data$stand_start_idxs[s]:gpp_data$stand_end_idxs[s]
  names <- paste0('GSL[', idxs, ']')
  util$plot_disc_pushforward_quantiles(samples1, names, ylab = 'GSL (days)', xlab = 'Year',
                                       display_ylim = c(165, 245))
  text(x = 0, y = 240, labels = paste0('Stand ', s), adj = 0, cex = 0.8)
}

par(mfrow = c(4,1), mar = c(4,4,1,1), cex.lab = 0.8, cex.axis = 0.8)
for(s in 1:gpp_data$N_stands){
  idxs <- gpp_data$stand_start_idxs[s]:gpp_data$stand_end_idxs[s]
  names <- paste0('AUC[', idxs, ']')
  util$plot_disc_pushforward_quantiles(samples1, names, ylab = 'AUC (gC/m²)', xlab = 'Year',
                                       display_ylim = c(400, 1600))
  text(x = 0, y = 1550, labels = paste0('Stand ', s), adj = 0, cex = 0.8)
}


# 2. An asymmetrical Gaussian

## 2.1 Plot the curve!
## Instead of having one unique width, we allow the width of the spring to be different than the with of the autumn
asym_gaussian <- function(t, gpp0, gamma, peakday, width_s, width_a) {
  width <- ifelse(t < peakday, width_s, width_a)
  gpp0 + gamma * exp(-0.5 * ((t - peakday) / width)^2)
}
g_asym <- asym_gaussian(doy, gpp0 = 0.1, gamma = 7, peakday = 180, width_s = 40, width_a = 70)
r <- (gpp_threshold-0.1)/7        
z <- sqrt(-2*log(r))  
SOS <- 180-40*z # spring width
EOS <- 180+70*z # autumn width
par(mfrow = c(1,1), mar = c(3.5,3.5,1,1), cex.lab = 0.8, cex.axis = 0.8, mgp = c(2.5,1,0))
plot(NULL, type = "l", lwd = 2, xlim = c(1, 365), ylim = c(0, 8),
     xlab = "DOY", ylab = "GPP (gC/m²/day)", yaxs = "i", xaxt = 'n', bty = 'n')
axis(1, at = seq(0,365,60))
abline(h = gpp_threshold, lty = 2, col = "grey40")   # threshold
polygon(c(SOS, EOS, EOS, SOS), c(0, 0, gpp_threshold, gpp_threshold),
        col = paste0(util$c_light_teal, '30'), border = NA)
x <- seq(SOS, EOS, length.out = 1e4)
y <- asym_gaussian(x, gpp0 = 0.1, gamma = 7, peakday = 180, width_s = 40, width_a = 70)
polygon(c(x, rev(x)), c(y, rep(gpp_threshold, 1e4)),
        col = paste0(util$c_light_teal, '50'), border = NA) 
segments(x0 = SOS, y0 = 0, y1 = gpp_threshold, lty = 2, col = util$c_mid_teal)
segments(x0 = EOS, y0 = 0, y1 = gpp_threshold, lty = 2, col = util$c_mid_teal)
lines(doy, g_asym, lwd = 2, col = util$c_mid_teal)
segments(x0 = 180, y0 = 0, y1 = 7.1, lty = 1, col = 'white', lwd = 2.5)
segments(x0 = 180, y0 = 0, y1 = 7.1, lty = 2, col = util$c_dark_teal, lwd = 1)
points(x = 180, y = 7.1, pch = 20, col = 'white', cex = 1.5)
points(x = 180, y = 7.1, pch = 20, col = util$c_dark_teal, cex = 0.9)


## 2.2. Posterior quantification
model2 <- stan_model(file.path(wd, 'stan', 'gs_charac', 'model2_assymgaussian.stan'))

init_fun <- function() {
  
  stand_idxs <- rep(1:gpp_data$N_stands, times = gpp_data$N_years_perstand)
  
  mu_peakday <- rnorm(gpp_data$N_stands, 180, 80 / 2.57)
  mu_log_gamma <- rnorm(gpp_data$N_stands, log(5), 0.75)
  mu_log_width <- rnorm(gpp_data$N_stands, log(40), log(2) / 2.57)
  mu_log_width_s <- rnorm(gpp_data$N_stands, log(40), log(2) / 2.57)
  mu_log_width_a <- rnorm(gpp_data$N_stands, log(40), log(2) / 2.57)
  
  tau_peakday <- abs(rnorm(1, 0, 10))
  tau_log_gamma <- abs(rnorm(1, 0, 0.25))
  tau_log_width <- abs(rnorm(1, 0, 0.25))
  tau_log_width_s <- abs(rnorm(1, 0, 0.25))
  tau_log_width_a <- abs(rnorm(1, 0, 0.25))
  
  list(
    gpp0  = rnorm(gpp_data$N_stands, 0, 1),
    sigma = abs(rnorm(gpp_data$N_stands, 0, 1)),
    mu_peakday = mu_peakday, mu_log_gamma = mu_log_gamma, mu_log_width = mu_log_width,
    tau_peakday = tau_peakday, tau_log_gamma = tau_log_gamma, tau_log_width = tau_log_width,
    peakday   = rnorm(gpp_data$N_years, mu_peakday[stand_idxs], tau_peakday),
    log_gamma = rnorm(gpp_data$N_years, mu_log_gamma[stand_idxs], tau_log_gamma),
    log_width_s = rnorm(gpp_data$N_years, mu_log_width_s[stand_idxs], tau_log_width_s),
    log_width_a = rnorm(gpp_data$N_years, mu_log_width_a[stand_idxs], tau_log_width_a)
  )
}

gpp_data$gpp_thresh <- 1
fit2 <- sampling(model2, gpp_data, chains = 4, cores = 4, init = init_fun)
samples2 <- util$extract_expectand_vals(fit2)

## 1.3. Retrodictive checks
sy <- 58
s <- 1 + 365 * (sy - 1);
e <- s + 364
names <- paste0('gpp_pred[', s:e, ']')
par(mfrow = c(2,1), mar = c(4,4,1,1), cex.lab = 0.8, cex.axis = 0.8)
util$plot_conn_pushforward_quantiles(samples2, names, plot_xs = 1:365, ylab = 'GPP (gC/m²/day)', xlab = 'DOY')
points(x = 1:365, y = gpp_data$gpp_obs[s:e], col = 'white', cex = 1.2, pch = 20)
points(x = 1:365, y = gpp_data$gpp_obs[s:e], col = 'black', cex = 0.5, pch = 20)
util$plot_conn_pushforward_quantiles(samples2, names, plot_xs = 1:365, ylab = 'GPP (gC/m²/day)', xlab = 'DOY')
y_sym <- sapply(names, function(n) util$ensemble_mcmc_quantile_est(samples1[[n]], c(0.5)))
lines(x = 1:365, y = y_sym, lwd = 6, col = 'white')
lines(x = 1:365, y = y_sym, lwd = 3, col = util$c_light_teal)
legend(x =0, y = 12, col = util$c_light_teal, lwd = 3, legend ='Symmetrical Gausian', cex = 0.7, bty = 'n')


par(mfrow = c(2,1), mar = c(4,4,1,1), cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8)
util$plot_disc_pushforward_quantiles(samples1, names, ylab = 'Residuals (GPP, gC/m²/day)', xlab = 'DOY',
                                     baseline_values = gpp_data$gpp_obs[s:e], residual = TRUE, main = 'Symmetrical Gaussian')
abline(v = c(110, 200), lty = 2)
util$plot_disc_pushforward_quantiles(samples2, names, ylab = 'Residuals (GPP, gC/m²/day)', xlab = 'DOY',
                                     baseline_values = gpp_data$gpp_obs[s:e], residual = TRUE, main = 'Asymmetrical Gaussian')
abline(v = c(110, 200), lty = 2)


## 2.4. Posterior inference
par(mfrow = c(4,1), mar = c(4,4,1,1), cex.lab = 0.8, cex.axis = 0.8)
for(s in 1:gpp_data$N_stands){
  idxs <- gpp_data$stand_start_idxs[s]:gpp_data$stand_end_idxs[s]
  names <- paste0('GSL[', idxs, ']')
  y_sym <- sapply(names, function(n) util$ensemble_mcmc_quantile_est(samples1[[n]], c(0.5)))
  util$plot_disc_pushforward_quantiles(samples2, names, ylab = 'GSL (days)', xlab = 'Year',
                                       display_ylim = c(165, 245), baseline_values = y_sym)
  text(x = 0, y = 240, labels = paste0('Stand ', s), adj = 0, cex = 0.8)
}

par(mfrow = c(4,1), mar = c(4,4,1,1), cex.lab = 0.8, cex.axis = 0.8)
for(s in 1:gpp_data$N_stands){
  idxs <- gpp_data$stand_start_idxs[s]:gpp_data$stand_end_idxs[s]
  names <- paste0('GSL[', idxs, ']')
  y_sym <- sapply(names, function(n) util$ensemble_mcmc_quantile_est(samples1[[n]], c(0.5)))
  util$plot_disc_pushforward_quantiles(samples2, names, ylab = 'Asymm model - sym. model\n(GSL, days)', xlab = 'Year',
                                       display_ylim = c(-7, 7), baseline_values = y_sym, residual = TRUE)
  text(x = 0, y = 240, labels = paste0('Stand ', s), adj = 0, cex = 0.8)
}


par(mfcol = c(4,2), mar = c(4,4,1,1), cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.85)
for(s in 1:gpp_data$N_stands){
  idxs <- gpp_data$stand_start_idxs[s]:gpp_data$stand_end_idxs[s]
  names <- paste0('SOS[', idxs, ']')
  y_sym <- sapply(names, function(n) util$ensemble_mcmc_quantile_est(samples1[[n]], c(0.5)))
  util$plot_disc_pushforward_quantiles(samples2, names, ylab = 'Asymm model - sym. model\n(SOS, doy)', xlab = 'Year',
                                       display_ylim = c(3, 20), baseline_values = y_sym, residual = TRUE,
                                       main = ifelse(s==1, 'Start of season', ''))
  text(x = 0, y = 19.5, labels = paste0('Stand ', s), adj = 0, cex = 0.8)
}

for(s in 1:gpp_data$N_stands){
  idxs <- gpp_data$stand_start_idxs[s]:gpp_data$stand_end_idxs[s]
  names <- paste0('EOS[', idxs, ']')
  y_sym <- sapply(names, function(n) util$ensemble_mcmc_quantile_est(samples1[[n]], c(0.5)))
  util$plot_disc_pushforward_quantiles(samples2, names, ylab = 'Asymm model - sym. model\n(EOS, doy)', xlab = 'Year',
                                       display_ylim = c(3, 20), baseline_values = y_sym, residual = TRUE,
                                       main = ifelse(s==1, 'End of season', ''))
  text(x = 0, y = 19.5, labels = paste0('Stand ', s), adj = 0, cex = 0.8)
}

par(mfcol = c(1,1), mar = c(4,4,1,1), cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.85)
plot(NULL, type = "l", lwd = 2, xlim = c(1, 365), ylim = c(0, 12),
     xlab = "DOY", ylab = "GPP (gC/m²/day)",  xaxt = 'n', bty = 'n')
axis(1, at = seq(0,365,60))
sy <- 58
s <- 1 + 365 * (sy - 1);
e <- s + 364
names <- paste0('gpp_pred[', s:e, ']')
y_sym <- sapply(names, function(n) util$ensemble_mcmc_quantile_est(samples1[[n]], c(0.5)))
lines(x = 1:365, y = y_sym, lwd = 3, col = 'white')
lines(x = 1:365, y = y_sym, lwd = 2, col = util$c_light_teal)
y_asym <- sapply(names, function(n) util$ensemble_mcmc_quantile_est(samples2[[n]], c(0.5)))
lines(x = 1:365, y = y_asym, lwd = 3, col = 'white')
lines(x = 1:365, y = y_asym, lwd = 2, col = util$c_mid)
abline(h = 1, lty = 2)
