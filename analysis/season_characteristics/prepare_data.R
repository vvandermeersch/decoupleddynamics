rm(list = ls());gc()
wd <- "/home/victor/projects/decoupleddynamics/analysis"
setwd(file.path(wd, 'stan'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local = util)
source('mcmc_visualization_tools.R', local = util)
setwd(wd)

library(rstan)
library(terra)
# Load treering data
datasets <- readRDS(file.path(wd, 'input', 'itrdb', 'datasets_summary_usonly.rds'))

# Same species
datasets[datasets$species_code == 'ABBI', c('species_name', 'species_code')] <- 
  unique(datasets[datasets$species_code == 'ABLA', c('species_name', 'species_code')])

# Create a grouped_stand key to gather stands if they have the same latitude and longitude (rounded to 1e-2 degree, ie ~1.1km)
group_keys <- interaction(
  round(datasets$north_lat,2),
  round(datasets$south_lat,2),
  round(datasets$east_lon,2),
  round(datasets$west_lon,2),
  drop = TRUE,
  sep = '/'
)
datasets$grouped_stand <- paste0("S", as.integer(factor(group_keys)))
group_coordinates <- data.frame(group_keys, grouped_stand = paste0("S", as.integer(factor(group_keys))))

# Temporary
datasets <- datasets[datasets$last_year >= 2016 & datasets$species_code == 'CHTH',]
uniq_stand_ids <- unique(datasets$grouped_stand)
all_years <-  2001:2021

N_years_perstand <- c()

nee_obs <- c()
gpp_obs <- c()
stand_start_idxs <- c()
stand_end_idxs <- c()
stand_years <- c()

stand_idx <- 1

for(s in uniq_stand_ids){
  
  stand_lat <- datasets[datasets$grouped_stand == s, 'north_lat']
  stand_lon <- datasets[datasets$grouped_stand == s, 'east_lon']
  target_coords <- vect(data.frame(x = stand_lon, y = stand_lat))
  
  min_year <- 1
  max_year <- 21
  
  N_years_perstand <- c(N_years_perstand, length(min_year:max_year))
  
  stand_start_idxs <- c(stand_start_idxs, stand_idx)
  stand_idx <- stand_idx + length(min_year:max_year) - 1
  stand_end_idxs <-  c(stand_end_idxs, stand_idx)
  stand_idx <- stand_idx + 1
  
  for(y in min_year:max_year){
    
    stand_years <- c(stand_years, y)
    
    ndays <- 365 # we ignore leap years
    NEE_fname <- file.path(wd, 'data', 'xbase', paste0("NEE_", all_years[y], "_025_daily.nc"))
    NEE_raster <- rast(NEE_fname, subds = "NEE")
    cell_NEE <- terra::extract(NEE_raster, target_coords, ID = FALSE)
    NEE_vals <- as.numeric(cell_NEE[1, 1:ndays])
    nee_obs <- c(nee_obs, NEE_vals)
    
    
    GPP_fname <- file.path(wd, 'data', 'xbase', paste0("GPP_", all_years[y], "_025_daily.nc"))
    GPP_raster <- rast(GPP_fname, subds = "GPP")
    cell_GPP <- terra::extract(GPP_raster, target_coords, ID = FALSE)
    GPP_vals <- as.numeric(cell_GPP[1, 1:ndays])
    gpp_obs <- c(gpp_obs, GPP_vals)
      
  }
  
}

N_stands <- length(uniq_stand_ids)
N <- length(gpp_obs)
N/365/N_stands
N_days <- 365

gpp_data <-
  list(
    N = N,
    N_stands = N_stands,
    N_years_perstand = N_years_perstand,
    N_years = sum(N_years_perstand),
    N_days = N_days,
    stand_years = stand_years,
    stand_start_idxs = stand_start_idxs, 
    stand_end_idxs = stand_end_idxs,
    gpp_obs = gpp_obs
  )

saveRDS(gpp_data, file.path('gpp_data.rds'))
