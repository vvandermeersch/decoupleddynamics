rm(list = ls())
library(terra)
library(rstan)
wd <- '/home/victor/projects/decoupleddynamics'

sites_lat <- c(44.15, 65.35)
sites_lon <- c(-71.43, 19.34)
num_sites <- length(sites_lat)


NEE <- c()
days <- c()
site_idxs <- c()
for (site in 1:1){
  target_coords <- matrix(c(sites_lon[site], sites_lat[site]), ncol=2)
  
  for (year in 2001:2001){
      
      NEE_fname <- file.path(wd, 'data', 'xbase', paste0("NEE_", year, "_025_daily.nc"))
      
      NEE_raster <- rast(NEE_fname)
      cell_NEE <- terra::extract(NEE_raster, target_coords)
      NEE_vals <- as.numeric(cell_NEE[1, ])
      NEE_vals <- NEE_vals[1:365]
      
      NEE <- c(NEE, NEE_vals)
      days <- c(days, (1:365)+(365*(year-2001)))
      site_idxs <- c(site_idxs, rep(site, 365))

  }

}

plot(NEE ~ days)
abline(v = 365)
abline(h = 0)

N <- length(NEE)

mdl.data <- list(
  N = N,
  days = days, 
  nee_obs = NEE
)


modelstan <- stan_model(file.path(wd, 'analysis/stan/gsl', 'model1.stan'))
fit <- sampling(modelstan, mdl.data, chains = 4, cores = 4, 
                seed = 12345, iter = 2000, warmup = 1000)

samples <- util$extract_expectand_vals(fit)

util$plot_conn_pushforward_quantiles(samples, paste0('nee_pred[', 1:365, ']'), 1:365)
points(x = mdl.data$days, y = mdl.data$nee_obs,
       pch = 20, col = "white", cex = 2)
points(x = mdl.data$days, y = mdl.data$nee_obs,
       pch = 20, col = "black", cex = 1)
abline(h = 0, lwd = 2, col = util $c_mid_teal, lty = 2)
abline(h = -1, lwd = 2, col = util $c_mid_teal, lty = 2)

par(mfrow = c(2,2))
util$plot_expectand_pushforward(samples[['alpha']], 20, flim = c(-3,3))
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = util$c_mid_teal)
util$plot_expectand_pushforward(samples[['rho']], 20, flim = c(0,100))
prior <- rnorm(1e6, 50, 10)
lines(density(prior), col = util$c_mid_teal)
util$plot_expectand_pushforward(samples[['gamma']], 20, flim = c(0,3))
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = util$c_mid_teal)
util$plot_expectand_pushforward(samples[['sigma']], 20, flim = c(0,3))
prior <- rnorm(1e6, 0, 1)
lines(density(prior), col = util$c_mid_teal)


week_idxs <- ceiling(mdl.data$days / 7)
week_idxs <- pmin(week_idxs, 52)
N_weeks <- 52
weeks <- 1:52

mdl.data <- list(
  N = N,
  N_weeks = N_weeks,
  days = days, 
  weeks = weeks,
  week_idxs = week_idxs,
  nee_obs = NEE
)

modelstan <- stan_model(file.path(wd, 'analysis/stan/gsl', 'model2.stan'))
fit <- sampling(modelstan, mdl.data, chains = 4, cores = 4, 
                seed = 12345, iter = 2000, warmup = 1000)



samples <- util$extract_expectand_vals(fit)

util$plot_conn_pushforward_quantiles(samples, paste0('nee_pred[', 1:365, ']'), 1:365)
points(x = mdl.data$days, y = mdl.data$nee_obs,
       pch = 20, col = "white", cex = 2)
points(x = mdl.data$days, y = mdl.data$nee_obs,
       pch = 20, col = "black", cex = 1)
abline(h = 0, lwd = 2, col = util $c_mid_teal, lty = 2)
abline(h = -1, lwd = 2, col = util $c_mid_teal, lty = 2)




