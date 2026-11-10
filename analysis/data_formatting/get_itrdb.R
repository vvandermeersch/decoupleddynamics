## Started 09 November 2026
## By Ken

## Pulling ITRDB data

rm(list = ls())
options(stringsAsFactors = FALSE)

#setwd
if(length(grep("victor", getwd())) > 0){
  setwd("/home/victor/projects/decoupleddynamics/analysis")
} else if(length(grep("Ken", getwd()) > 0)){
  setwd("/Users/Ken Michiko Samson/Documents/Temporal Ecology Lab/decoupleddynamics/analysis")
}

library(xml2)

min_long <- -90
max_long <- -45
min_lat <- 20
max_lat <- 60

itrdb_path <- paste0(dirname(getwd()), "/data/ITRDB")
itrdb_url <- paste0("https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica/", c("canada", "usa"))

download_until_success <- function(url, destfile, ..., maxcount = 5) {
  count <- 0
  repeat{
    Sys.sleep(0.5)
    try(download.file(url, destfile, ...))
    count <- count + 1
    if (file.exists(destfile) || count >= maxcount){break}
  }
}

dref <- list()
counter <- 0

for(i in 1:length(itrdb_url)){
  filenames <- unlist(xml2::as_list(httr::content(httr::GET(itrdb_url[i]))))
  filenames <- filenames[grepl(filenames, pattern = '*.rwl$') & !grepl(filenames, pattern = '*-noaa.rwl$')]
  
  # filenames <- filenames[which(!(filenames %in% c("cana.rwl", "ar.rwl", "co.rwl", 
  #                                                 "mo.rwl", "ny.rwl", "ok.rwl",
  #                                                 "pa.rwl", "tn.rwl")))]
  
  for(j in 1201:length(filenames)){
    template_file <- gsub("\\.rwl", "-rwl-noaa.txt", filenames[j])
    
    id <- gsub("\\.rwl", "", filenames[j])
    type <- regmatches(id, regexec("[0-9]+(.*)", id))[[1]][2]
    
    if(type %in% c("", "rw", "w")){
      download_until_success(file.path(itrdb_url[i], template_file), file.path(itrdb_path, template_file))
      
      if(file.exists(file.path(itrdb_path, template_file))){
        file_lines <- readLines(file.path(itrdb_path, template_file))
      } else {
        next
      }
      
      for(k in 1:length(file_lines)){
        if(length(strsplit(file_lines[k], "Northernmost_Latitude: ", fixed = TRUE)[[1]]) == 2){
          north <- as.numeric(strsplit(file_lines[k], "Northernmost_Latitude: ", fixed = TRUE)[[1]][2])
          south <- as.numeric(strsplit(file_lines[k + 1], "Southernmost_Latitude: ", fixed = TRUE)[[1]][2])
          east <- as.numeric(strsplit(file_lines[k + 2], "Easternmost_Longitude: ", fixed = TRUE)[[1]][2])
          west <- as.numeric(strsplit(file_lines[k + 3], "Westernmost_Longitude: ", fixed = TRUE)[[1]][2])
          
          if(min_long < west & max_long > east & min_lat < south & max_lat > north){
            site <- strsplit(file_lines[k - 2], "Site_Name: ", fixed = TRUE)[[1]][2]
            location <- strsplit(file_lines[k - 1], "Location: ", fixed = TRUE)[[1]][2]
            alt <- strsplit(file_lines[k + 4], "Elevation_m: ", fixed = TRUE)[[1]][2]
            collection <- strsplit(file_lines[k + 7], "Collection_Name: ", fixed = TRUE)[[1]][2]
            first_year <- strsplit(file_lines[k + 8], "First_Year: ", fixed = TRUE)[[1]][2]
            last_year <- strsplit(file_lines[k + 9], "Last_Year: ", fixed = TRUE)[[1]][2]
            time_unit <- strsplit(file_lines[k + 10], "Time_Unit: ", fixed = TRUE)[[1]][2]
            species_name <- strsplit(file_lines[k + 16], "Species_Name: ", fixed = TRUE)[[1]][2]
            common_name <- strsplit(file_lines[k + 17], "Common_Name: ", fixed = TRUE)[[1]][2]
            species_code <- strsplit(file_lines[k + 18], "Tree_Species_Code: ", fixed = TRUE)[[1]][2]
            
            ref <- data.frame(
              "collection_name" = collection,
              "file_name" = id,
              "species_name" = species_name,
              "tree_species_code" = species_code,
              "common_name" = common_name,
              "northernmost_latitude" = north,
              "southermost_latitude" = south,
              "westernmost_longitude" = west,
              "easternmost_longitude" = east,
              "elevation (m)" = alt,
              "first_year" = first_year,
              "last_year" = last_year,
              "time_unit" = time_unit
            )
            
            download_until_success(file.path(itrdb_url[i], filenames[j]), file.path(itrdb_path, filenames[j]))
            
            counter <- counter + 1
            dref[[counter]] <- ref
            
            break
          } else {
            file.remove(file.path(itrdb_path, template_file))
            break
          }
        }
      }
    }
  }
}

dref <- do.call(rbind, dref)
write.csv(dref, file.path(itrdb_path, "_summary.csv"), row.names = FALSE)
