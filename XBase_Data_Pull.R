library(terra)
library(ncdf4)
library(ggplot2)

# Storage of files in computer: Replace with path 
# to folders where .nc files are stored (GPP and NEE in separate folders)
GPP_path <- "C:/Temporal Ecology Lab/GPP/"
NEE_path <- "C:/Temporal Ecology Lab/NEE/"

# Plot Configuration Parameters (Modify these)

# Set to TRUE to plot specific flux, FALSE otherwise
plot_GPP <- FALSE
plot_NEE <- TRUE

# Set latitude and longitude for specific sites (maps to closest cell)
sites_lat <- c(44.15, 65.35)
sites_lon <- c(-71.43, 19.34)
num_sites <- length(sites_lat)

# Set range of days (since Jan 1st of each year) to plot over during the year
# Min: 0    Max: 364
day_start <- 0
day_end <- 364

# Set range of years (between 2001 and 2021) to plot over
year_start <- 2001
year_end <- 2002
# End of Plot Configuration Parameters


# Vectors for time frame
data_years <- year_start:year_end
time_frame <- day_start:day_end

# Function for generating simple plot of fluxes
plot_data <- function(days_vector, flux_vector, year_num, flux_name, site_loc) {
  plot(days_vector, flux_vector,
       xlab = paste("Time (days since ",year_num, "-01-01)", sep=""),
       ylab = flux_name,
       type = 'p',
       col = 'coral1')
  
  mtext(site_loc, side = 3, line = 2, adj = 0.5)
  
  grid(col = 'lightgray',
       lty = 'solid',
       lwd = 1
  )
  return()
}

# Main Program: Loops through locations, year range specified in the parameters
# and generates plots for GPP and NEE fluxes
for (site in 1:num_sites){
  target_coords <- matrix(c(sites_lon[site], sites_lat[site]), ncol=2)

  if (plot_GPP == TRUE) {
    for (year in year_start:year_end){
    
        GPP_fname <- paste(GPP_path, "GPP_", year, "_025_daily.nc", sep="")
        
        GPP_raster <- rast(GPP_fname)
        cell_GPP <- terra::extract(GPP_raster, target_coords)
        GPP_vals <- as.numeric(cell_GPP[1, ])
        GPP_vals <- GPP_vals[1:365]
        GPP_inRange <- GPP_vals[day_start:day_end+1]
        
        site_loc <- paste("Site: (", sites_lat[site], ",", sites_lon[site], ")", sep="")
        
        plot_data(time_frame, GPP_inRange, year, "GPP (gC m-2 d-1)", site_loc)
      }
  } 
}
    
for (site in 1:num_sites){
  target_coords <- matrix(c(sites_lon[site], sites_lat[site]), ncol=2)
  
  if (plot_NEE == TRUE){
    for (year in year_start:year_end){
    
        NEE_fname <- paste(NEE_path, "NEE_", year, "_025_daily.nc", sep="")
        
        NEE_raster <- rast(NEE_fname)
        cell_NEE <- terra::extract(NEE_raster, target_coords)
        NEE_vals <- as.numeric(cell_NEE[1, ])
        NEE_vals <- NEE_vals[1:365]
        NEE_inRange <- NEE_vals[day_start:day_end+1]
        
        site_loc <- paste("Site: (", sites_lat[site], ",", sites_lon[site], ")", sep="")
        
        plot_data(time_frame, NEE_inRange, year, "NEE (gC m-2 d-1)", site_loc)
      }
  } 
}





