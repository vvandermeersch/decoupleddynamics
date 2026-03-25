# Generates interactive GPP or NEE Flux Plot for ONE specific site over
# the course of several years in the 'Viewer' tab
# - Raw data points - Zoom in to observe specific values for each day

library(terra)
library(ggplot2)
library(plotly)

# Storage of files in computer: Replace with path 
# to folders where .nc files are stored (GPP and NEE in separate folders)
GPP_path <- "C:/Temporal Ecology Lab/GPP/"
NEE_path <- "C:/Temporal Ecology Lab/NEE/"

# *** Plot Configuration Parameters (MODIFY these) ***
# Set to TRUE to plot specific flux, FALSE otherwise 
# (Only set one to TRUE at a time)
plot_GPP <- TRUE
plot_NEE <- FALSE

# Set latitude and longitude for a specific site (maps to closest cell)
site_lat <- 49.28
site_lon <- 123.12

# Set range of years (between 2001 and 2021) to plot over
year_start <- 2003
year_end <- 2005

# *** End of Plot Configuration Parameters ***

# *** Program Variables (DO NOT MODIFY) ***
# Vectors for time frame and location
data_years <- year_start:year_end

# Data frame to hold NEE and GPP raw values
flux_dataset <- data.frame(
  Lat = numeric(),
  Lon = numeric(),
  Date = as.Date(character(0)),
  GPP_val = numeric(),
  NEE_val = numeric(),
  stringsAsFactors = FALSE
)

# List for yearly subsets of data, later binded in data frame
list_subsets <- list()
# *** End of Program Variables ***

# *** Main Program Start ***
# PUSH ALL DATA TO DATAFRAME: Columns are latitude, longitude, date, GPP, NEE)
target_coords <- matrix(c(site_lon, site_lat), ncol=2)

for (year in year_start:year_end){
  GPP_fname <- paste(GPP_path, "GPP_", year, "_025_daily.nc", sep="")
  GPP_raster <- rast(GPP_fname)
  cell_GPP <- terra::extract(GPP_raster, target_coords)
  GPP_vals <- as.numeric(cell_GPP[1, ])
  GPP_vals <- GPP_vals[1:length(GPP_vals)-1]
  
  NEE_fname <- paste(NEE_path, "NEE_", year, "_025_daily.nc", sep="")
  NEE_raster <- rast(NEE_fname)
  cell_NEE <- terra::extract(NEE_raster, target_coords)
  NEE_vals <- as.numeric(cell_NEE[1, ])
  NEE_vals <- NEE_vals[1:length(NEE_vals)-1]
  
  calendar_days <- as.Date(1:length(GPP_vals)-1, origin = paste(year, "-01-01", sep=""))
  
  temp_flux <- data.frame(
    Lat = rep(site_lat, length(GPP_vals)),
    Lon = rep(site_lon, length(GPP_vals)),
    Date = calendar_days,
    GPP_val = GPP_vals,
    NEE_val = NEE_vals,
    stringsAsFactors = FALSE
  )
  
  list_subsets[[length(list_subsets) + 1]] <- temp_flux
  
}

flux_dataset <- data.table::rbindlist(list_subsets)

# PLOT GPP OR NEE FLUX VALUES
if (plot_GPP == TRUE){
  annual_set <- subset(flux_dataset, Date >= paste(year_start, "-01-01", sep="") & Date < paste(year_end+1, "-01-01", sep=""))
  p <- ggplot()+
          geom_line(data = annual_set, mapping = aes(x = Date, y = GPP_val), color = "blue")+
          labs(x = "Date", y = "GPP (gC m-2 d-1)", title = paste("GPP Flux: Site (", site_lat, ",", site_lon, ")",sep=""))
  ggplotly(p)
}

if (plot_NEE == TRUE){
    annual_set <- subset(flux_dataset, Date >= paste(year_start, "-01-01", sep="") & Date < paste(year_end+1, "-01-01", sep=""))
    p <- ggplot()+
            geom_line(data = annual_set, mapping = aes(x = Date, y = NEE_val), color = "orange")+
            labs(x = "Date", y = "NEE (gC m-2 d-1)", title = paste("NEE Flux: Site (", site_lat, ",", site_lon, ")",sep=""))
    ggplotly(p)
}
