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
# datasets <- datasets[datasets$dataset %in% c("ny052"),]
datasets <- datasets[datasets$species_code == 'ACSH' & datasets$last_year >= 2015,]


# Prepare tree ring data!
ringwidth_series <- readRDS(file.path(wd, 'input', 'itrdb', 'ringwidth_series_usonly_from1896.rds'))

all_years <-  2001:2021
# all_years <- min(ringwidth_series$year): max(ringwidth_series$year)

raw_data <- data.frame()
for(d in 1:nrow(datasets)){
  
  raw_data_d <- ringwidth_series[ringwidth_series$dataset %in% datasets[d, 'dataset'], ]
  
  raw_data_d <- raw_data_d[raw_data_d$year %in% all_years, ]
  
  raw_data_d$species_code <- datasets[d, 'species_code']
  
  # temporary
  # raw_data_d <- raw_data_d[raw_data_d$year >= 1950 & raw_data_d$year <= 2025,]
  
  # create a unique tree id (across all datasets)
  raw_data_d$original_tree_id <- raw_data_d$tree_id
  raw_data_d$tree_id_uniq <- paste0(raw_data_d$dataset, "_", raw_data_d$tree_id)
  
  raw_data_d <- aggregate(
    rw_mm ~ dataset + species_code + original_tree_id + tree_id_uniq + year,
    data = raw_data_d,
    FUN = function(x) mean(x, na.rm = TRUE),
    na.action = na.pass
  )
  names(raw_data_d)[names(raw_data_d) == "rw_mm"] <- "rw_avg_mm"
  
  raw_data_d <- na.omit(raw_data_d)
  raw_data_d$year <- as.numeric(raw_data_d$year)
  
  # Here, we check that there is no missing year in the individual time series
  count_bytrees <- aggregate(
    rw_avg_mm ~ tree_id_uniq,
    data = raw_data_d,
    FUN = function(x) length(x)
  )
  years_bytrees <- aggregate(
    year ~ tree_id_uniq,
    data = raw_data_d,
    FUN = function(x) max(x)-min(x)+1
  )
  check_length <- merge(count_bytrees, years_bytrees)
  trees_to_remove <- check_length[check_length$rw_avg_mm != check_length$year, 'tree_id_uniq']
  raw_data_d <- raw_data_d[!(raw_data_d$tree_id_uniq %in% trees_to_remove),]
  
  # remove less than 10 years observed
  trees_to_remove <- check_length[check_length$rw_avg_mm < 10, 'tree_id_uniq']
  raw_data_d <- raw_data_d[!(raw_data_d$tree_id_uniq %in% trees_to_remove),]
  
  # remove trees with 0mm observations - NOT ANYMORE
  # trees_to_remove <- unique(raw_data_d[raw_data_d$rw_avg_mm == 0 , 'tree_id_uniq'])
  # raw_data_d <- raw_data_d[!(raw_data_d$tree_id_uniq %in% trees_to_remove),]
  
  raw_data <- rbind(raw_data, raw_data_d)
}
raw_data <- merge(raw_data,  datasets[, c("dataset", "grouped_stand")], by.x = 'dataset', by.y = 'dataset')
length( unique(datasets$grouped_stand))

all_years <- min(raw_data$year): max(raw_data$year)

# Deal with potential duplicates
# If several trees from different ITRDB datasets are on the same site and have the same ID, the same number of years
# and the same mean ringwidth, we keep only one of them
potential_duplicates <- aggregate(dataset ~ original_tree_id + species_code, data = raw_data, FUN = function(x) length(unique(x)))
potential_duplicates_groupedstand <- aggregate(grouped_stand ~ original_tree_id + species_code, data = raw_data, FUN = function(x) length(unique(x)))
potential_duplicates <- merge(potential_duplicates, potential_duplicates_groupedstand)
potential_duplicates <- potential_duplicates[potential_duplicates$dataset > 1 & potential_duplicates$grouped_stand == 1,]
cat(paste0(nrow(potential_duplicates), ' trees are potential duplicates!\n'))
for(t in potential_duplicates$original_tree_id){
  sp_t <- potential_duplicates[potential_duplicates$original_tree_id == t, 'species_code']
  raw_data_t <- raw_data[raw_data$original_tree_id == t & raw_data$species_code == sp_t, ] 
  count_years <- aggregate(year ~ tree_id_uniq, data = raw_data_t,  FUN = length)
  if(length(unique(count_years$year)) == 1){
    mean_rw <- aggregate(rw_avg_mm ~ tree_id_uniq, data = raw_data_t,  FUN = mean)
    if(length(unique(mean_rw$rw_avg_mm)) == 1){
      tree_to_keep <- sample(unique(raw_data_t$tree_id_uniq), 1)
      raw_data <- raw_data[-which(raw_data$original_tree_id == t & raw_data$tree_id_uniq != tree_to_keep),]
    }else{
      if(t == 'pic382'){
        # keep both trees because they are different
      }else{
        print(t)
        stop()
      }
    }
  }else if(length(unique(count_years$year)) > 1){
    if(t == 'at2022'){
      tree_to_keep <- 'co589_at2022' # longer time series
      raw_data <- raw_data[-which(raw_data$original_tree_id == t & raw_data$tree_id_uniq != tree_to_keep),]
    }else{
      print(t)
      stop()
    }
    
  }
}
datasets <- datasets[datasets$dataset %in% unique(raw_data$dataset),]
length( unique(datasets$grouped_stand))


# Sizes
uniq_tree_ids <- unique(raw_data$tree_id_uniq)
N_trees <- length(uniq_tree_ids)

uniq_stand_ids <- unique(raw_data$grouped_stand)
N_stands <- length(uniq_stand_ids)

uniq_stand_species_ids <- unique(paste(raw_data$grouped_stand, raw_data$species_code, sep = '_'))
N_stand_species <- length(uniq_stand_species_ids)

uniq_species_ids <- unique(raw_data$species_code)
N_species <- length(uniq_species_ids)
clade_idxs <- 1

# all_years <-  min(raw_data$year):max(raw_data$year)
N_all_years <- length(all_years)

# Format data into ragged arrays
rw_obs <- c()
years <- c()
all_years_idxs <- c()
N_years <- c()
stand_idxs <- c()
stand_species_idxs <- c()
species_idxs <- c()

idx <- 1
tree_start_idxs <- c()
tree_end_idxs <- c()

for(tid in uniq_tree_ids) {
  print(tid)
  
  raw_data_tree <- raw_data[raw_data$tree_id_uniq == tid & raw_data$year %in% all_years,]
  
  years_tree <- raw_data_tree$year
  all_years_idxs_tree <- sapply(years_tree, function(y) which(all_years == y))
  N_years_tree <- length(years_tree)
  # if(N_years_tree > 45 | N_years_tree < 20){stop()}
  
  rw_obs_tree <- sapply(years_tree, 
                        function(y) 
                          raw_data_tree$rw_avg_mm[raw_data_tree$year == y][1])
  rw_obs <- c(rw_obs, rw_obs_tree)
  
  years <- c(years, years_tree)
  all_years_idxs <- c(all_years_idxs, all_years_idxs_tree)
  N_years <- c(N_years, N_years_tree)
  
  stand_tree <- which(uniq_stand_ids == raw_data_tree$grouped_stand[1])
  stand_idxs <- c(stand_idxs, stand_tree)
  
  stand_species_tree <- which(uniq_stand_species_ids == paste0(raw_data_tree$grouped_stand[1], '_', raw_data_tree$species_code[1]))
  stand_species_idxs <- c(stand_species_idxs, stand_species_tree)
  
  species_tree <- which(uniq_species_ids == raw_data_tree$species_code[1])
  species_idxs <- c(species_idxs, species_tree)
  
  tree_start_idxs <- c(tree_start_idxs, idx)
  idx <- idx + N_years_tree
  tree_end_idxs <- c(tree_end_idxs, idx - 1)
}



N_stand_trees <- c()
N_stand_years <- c()
N_years_perstand <- c()

stand_tree_idxs <- c()

nee_obs <- c()
stand_start_idxs <- c()
stand_end_idxs <- c()

stand_years <- c()

year_padded_start_idxs <- c()
year_padded_end_idxs <- c()

year_true_start_idxs <- c()
year_true_end_idxs <- c()

stand_idx <- 1
padded_year_idx <- 1
true_year_idx <- 0
padding <- 15

for(s in uniq_stand_ids){
  
  # if(s == 'S116'){stop()}
  
  padding_counter <- 0
  
  stand_lat <- datasets[datasets$grouped_stand == s, 'north_lat']
  stand_lon <- datasets[datasets$grouped_stand == s, 'east_lon']
  target_coords <- vect(data.frame(x = stand_lon, y = stand_lat))

  N_stand_trees_s <- sum(stand_idxs == which(uniq_stand_ids == s))
  N_stand_trees <- c(N_stand_trees, N_stand_trees_s)
  
  stand_tree_idxs_here <- which(stand_idxs == which(uniq_stand_ids == s))
  stand_tree_idxs <- c(stand_tree_idxs, stand_tree_idxs_here)
  min_year <- min(all_years_idxs[tree_start_idxs[stand_tree_idxs_here]])
  max_year <- max(all_years_idxs[tree_end_idxs[stand_tree_idxs_here]])
  
  N_years_perstand <- c(N_years_perstand, length(min_year:max_year))

  stand_start_idxs <- c(stand_start_idxs, stand_idx)
  stand_idx <- stand_idx + length(min_year:max_year) - 1
  stand_end_idxs <-  c(stand_end_idxs, stand_idx)
  stand_idx <- stand_idx + 1
  
  for(y in min_year:max_year){
    
    stand_years <- c(stand_years, y)
    
    if(y == 1){
      true_year_idx <- true_year_idx + 1
      
      ndays <- 365 # we ignore leap years
      NEE_fname <- file.path(wd, 'data', 'xbase', paste0("NEE_", all_years[y], "_025_daily.nc"))
      NEE_raster <- rast(NEE_fname, subds = "NEE")
      cell_NEE <- terra::extract(NEE_raster, target_coords, ID = FALSE)
      NEE_vals <- as.numeric(cell_NEE[1, 1:ndays])
      nee_obs <- c(nee_obs, NEE_vals)
  
      NEE_fname <- file.path(wd, 'data', 'xbase', paste0("NEE_", all_years[y+1], "_025_daily.nc"))
      NEE_raster <- rast(NEE_fname, subds = "NEE")
      cell_NEE <- terra::extract(NEE_raster, target_coords, ID = FALSE)
      nee_obs_next <- as.numeric(cell_NEE[1, 1:padding])
      nee_obs <- c(nee_obs, nee_obs_next)
      
      ndays <- 365 + padding
    }else if(y == length(all_years)){
      true_year_idx <- true_year_idx + padding_counter + 1
      
      NEE_fname <- file.path(wd, 'data', 'xbase', paste0("NEE_", all_years[y-1], "_025_daily.nc"))
      NEE_raster <- rast(NEE_fname, subds = "NEE")
      cell_NEE <- terra::extract(NEE_raster, target_coords, ID = FALSE)
      nee_obs_prev <- as.numeric(cell_NEE[1,  seq(365-padding+1,365,1)])
      nee_obs <- c(nee_obs, nee_obs_prev)

      ndays <- 365 # we ignore leap years
      NEE_fname <- file.path(wd, 'data', 'xbase', paste0("NEE_", all_years[y], "_025_daily.nc"))
      NEE_raster <- rast(NEE_fname, subds = "NEE")
      cell_NEE <- terra::extract(NEE_raster, target_coords, ID = FALSE)
      NEE_vals <- as.numeric(cell_NEE[1, 1:ndays])
      nee_obs <- c(nee_obs, NEE_vals)
      
      ndays <- 365 + padding
    }else{
      true_year_idx <- true_year_idx + padding_counter + 1
      
      NEE_fname <- file.path(wd, 'data', 'xbase', paste0("NEE_", all_years[y-1], "_025_daily.nc"))
      NEE_raster <- rast(NEE_fname, subds = "NEE")
      cell_NEE <- terra::extract(NEE_raster, target_coords, ID = FALSE)
      nee_obs_prev <- as.numeric(cell_NEE[1,  seq(365-padding+1,365,1)])
      nee_obs <- c(nee_obs, nee_obs_prev)
      
      ndays <- 365 # we ignore leap years
      NEE_fname <- file.path(wd, 'data', 'xbase', paste0("NEE_", all_years[y], "_025_daily.nc"))
      NEE_raster <- rast(NEE_fname, subds = "NEE")
      cell_NEE <- terra::extract(NEE_raster, target_coords, ID = FALSE)
      NEE_vals <- as.numeric(cell_NEE[1, 1:ndays])
      nee_obs <- c(nee_obs, NEE_vals)
      
      NEE_fname <- file.path(wd, 'data', 'xbase', paste0("NEE_", all_years[y+1], "_025_daily.nc"))
      NEE_raster <- rast(NEE_fname, subds = "NEE")
      cell_NEE <- terra::extract(NEE_raster, target_coords, ID = FALSE)
      nee_obs_next <- as.numeric(cell_NEE[1, 1:padding])
      nee_obs <- c(nee_obs, nee_obs_next)
      
      ndays <- 365 + 2*padding
    }
    
    year_padded_start_idxs <- c(year_padded_start_idxs, padded_year_idx)
    padded_year_idx <- padded_year_idx + ndays - 1
    year_padded_end_idxs <- c(year_padded_end_idxs, padded_year_idx)
    padded_year_idx <- padded_year_idx + 1 
    
    
    year_true_start_idxs <- c(year_true_start_idxs, true_year_idx)
    true_year_idx <- true_year_idx + 364 
    year_true_end_idxs <- c(year_true_end_idxs, true_year_idx)
    
    padding_counter <- 2 * padding

  }
  
  true_year_idx <- padded_year_idx-1
  
}



N_stand_years <- sum(N_years_perstand)
N_days_peryear <- 365


# Cross check sizes
N_trees
length(rw_obs)
length(years)
length(all_years_idxs)
length(N_years)
length(tree_start_idxs)
length(tree_end_idxs)


length(N_stands)
length(N_stand_years)
length(nee_obs)
max(year_padded_end_idxs)
max(year_true_end_idxs)
length(nee_obs)/ndays





