# Generates GPP or NEE Flux Plots (Raw data + 10-day moving average) for multiple
# sites over the course of several years

library(terra)
library(ggplot2)
library(zoo)
library(patchwork)

# Storage of files in computer: Replace with path 
# to folders where .nc files are stored (GPP and NEE in separate folders)
GPP_path <- "C:/Temporal Ecology Lab/GPP/"
NEE_path <- "C:/Temporal Ecology Lab/NEE/"

# *** Plot Configuration Parameters (MODIFY) ***

# Set to TRUE to plot specific flux, FALSE otherwise
plot_GPP <- TRUE
plot_NEE <- FALSE

# Set latitude and longitude for specific sites (maps to closest cell)
site_lat <- c(44.15, 42.5, 49.28)
site_lon <- c(-71.43, 72.2, 123.12)

# Set range of years (between 2001 and 2021) to plot over
year_start <- 2001
year_end <- 2009
# *** End of Plot Configuration Parameters ***

# *** Program Variables (DO NOT MODIFY) ***

# Vectors for time frame and location
data_years <- year_start:year_end
num_sites <- length(site_lat)

# Data frame for flux values
flux_dataset <- data.frame(
  Lat = numeric(),
  Lon = numeric(),
  Date = as.Date(character(0)),
  flux_val = numeric(),
  stringsAsFactors = FALSE
)

list_subsets <- list()
plots_list <- list()

# *** MAIN PROGRAM: Extracts flux values, stores in data frame, then plots ***
# Flux values are cleared from data frame for each site iteration
if (plot_GPP == TRUE) {
  for (n in 1:num_sites){
    for (year in year_start:year_end){
      
      target_coords <- matrix(c(site_lon[n], site_lat[n]), ncol=2)
      
      GPP_fname <- paste(GPP_path, "GPP_", year, "_025_daily.nc", sep="")
      GPP_raster <- rast(GPP_fname)
      cell_GPP <- terra::extract(GPP_raster, target_coords)
      GPP_vals <- as.numeric(cell_GPP[1, ])
      GPP_vals <- GPP_vals[1:length(GPP_vals)-1]
      
      calendar_days <- as.Date(1:length(GPP_vals)-1, origin = paste(year, "-01-01", sep=""))
      
      temp_flux <- data.frame(
        Lat = rep(site_lat[n], length(GPP_vals)),
        Lon = rep(site_lon[n], length(GPP_vals)),
        Date = calendar_days,
        flux_val = GPP_vals,
        stringsAsFactors = FALSE 
      )
      
      list_subsets[[length(list_subsets) + 1]] <- temp_flux
    }
    flux_dataset <- data.table::rbindlist(list_subsets)
    
    for (k in 0:(length(data_years)-1)%/%3){
      annual_set <- subset(flux_dataset, Date >= paste(year_start+3*k, "-01-01", sep="") & Date < paste(min(year_end+1,year_start+3*k+3), "-01-01", sep=""))
      annual_set$mvavg <- rollmean(annual_set$flux_val, k = 10, fill = NA, align = "center")
      p <- ggplot()+
        geom_line(data = annual_set, mapping = aes(x = Date, y = flux_val), color = "orange", alpha = 0.5)+
        geom_line(data = annual_set, mapping = aes(x = Date, y = mvavg), color = "brown", linewidth = 1.125)+
        labs(x = "Date", y = "(gC m-2 d-1)", title = paste("GPP Flux: Site (", site_lat[n], ",", site_lon[n], ")",sep=""))
      plots_list[[k+1]] <- p
    }
    
    print(wrap_plots(plots_list, ncol = 1))
    
    flux_dataset <- flux_dataset[0, ]
    list_subsets <- list() 
    plots_list <- list() 
  }
}

if (plot_NEE == TRUE) {
  for (n in 1:num_sites){
    for (year in year_start:year_end){
      
      target_coords <- matrix(c(site_lon[n], site_lat[n]), ncol=2)
      
      NEE_fname <- paste(NEE_path, "NEE_", year, "_025_daily.nc", sep="")
      NEE_raster <- rast(NEE_fname)
      cell_NEE <- terra::extract(NEE_raster, target_coords)
      NEE_vals <- as.numeric(cell_NEE[1, ])
      NEE_vals <- NEE_vals[1:length(NEE_vals)-1]
      
      calendar_days <- as.Date(1:length(NEE_vals)-1, origin = paste(year, "-01-01", sep=""))
      
      temp_flux <- data.frame(
        Lat = rep(site_lat[n], length(NEE_vals)),
        Lon = rep(site_lon[n], length(NEE_vals)),
        Date = calendar_days,
        flux_val = NEE_vals,
        stringsAsFactors = FALSE 
      )
      
      list_subsets[[length(list_subsets) + 1]] <- temp_flux
    }
    flux_dataset <- data.table::rbindlist(list_subsets)
    
    for (k in 0:(length(data_years)-1)%/%3){
      annual_set <- subset(flux_dataset, Date >= paste(year_start+3*k, "-01-01", sep="") & Date < paste(min(year_end+1,year_start+3*k+3), "-01-01", sep=""))
      annual_set$mvavg <- rollmean(annual_set$flux_val, k = 10, fill = NA, align = "center")
      p <- ggplot()+
        geom_line(data = annual_set, mapping = aes(x = Date, y = flux_val), color = "orange", alpha = 0.5)+
        geom_line(data = annual_set, mapping = aes(x = Date, y = mvavg), color = "purple", linewidth = 1.125)+
        labs(x = "Date", y = "(gC m-2 d-1)", title = paste("NEE Flux: Site (", site_lat[n], ",", site_lon[n], ")",sep=""))
      plots_list[[k+1]] <- p
    }
    
    print(wrap_plots(plots_list, ncol = 1))
    
    flux_dataset <- flux_dataset[0, ]
    list_subsets <- list() 
    plots_list <- list() 
  }
}


