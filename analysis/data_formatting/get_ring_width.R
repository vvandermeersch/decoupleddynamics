## Started 10 November 2026
## By Ken

## Formatting ring width data

rm(list = ls())
options(stringsAsFactors = FALSE)

#setwd
if(length(grep("victor", getwd())) > 0){
  setwd("/home/victor/projects/decoupleddynamics/analysis")
} else if(length(grep("Ken", getwd()) > 0)){
  setwd("/Users/Ken Michiko Samson/Documents/Temporal Ecology Lab/decoupleddynamics/analysis")
}

dref <- read.csv("../data/ITRDB/_summary.csv")
rw <- read.csv("../data/ITRDB/_rw_summary.csv")

rw <- cbind(rw, data.frame(matrix(NA, nrow = nrow(rw), ncol = 3)))
colnames(rw)[5:7] <- c('first_year', 'last_year', 'raw_units')

yrs <- min(dref[["first_year"]]):max(dref[["last_year"]])
yr_data <- data.frame(matrix(NA, nrow = nrow(rw), ncol = length(yrs)))
rw <- cbind(rw, data.frame(matrix(NA, nrow = nrow(rw), ncol = length(yrs))))
colnames(rw)[(1:length(yrs))+7] <- yrs

files <- paste0("../data/ITRDB/", rw[["file_name"]], ".rwl")

for(i in 1:nrow(rw)){
  file_lines <- readLines(files[i])
  
  ids <- sapply(1:length(file_lines), function(x){return(substr(file_lines[x], 1, 8))})
  ids_nospace <- sapply(ids, function(x){paste0(strsplit(x, " ")[[1]], collapse = "")})
  
  core_line_idxs <- which(ids_nospace == rw[["tree_core_id"]][i])
  rw_vect <- paste(substr(file_lines[core_line_idxs], 13, nchar(file_lines[core_line_idxs])), collapse = " ")
  
  rw_data <- strsplit(rw_vect, " ")[[1]]
  rw_data <- as.numeric(rw_data[which(nchar(rw_data) != 0)])
  
  first_year <- as.numeric(substr(file_lines[core_line_idxs[1]], 9, 12))
  last_year <- length(rw_data) + first_year - 2
  raw_units <- rw_data[length(rw_data)]
  
  if(raw_units == 999){
    raw_units <- "0.01 mm"
    rw_data <- rw_data[1:(length(rw_data) - 1)] * 0.01
  } else if(raw_units == -9999){
    raw_units <- "0.001 mm"
    rw_data <- rw_data[1:(length(rw_data) - 1)] * 0.001
  } else {
    raw_units <- paste0(raw_units, " (unknown marker)")
    rw_data <- paste0(rw_data[1:(length(rw_data) - 1)], "?")
  }
  
  rw[["first_year"]][i] <- first_year
  rw[["last_year"]][i] <- last_year
  rw[["raw_units"]][i] <- raw_units
  which(colnames(rw) %in% first_year:last_year)
  rw[i, which(colnames(rw) %in% first_year:last_year)] <- rw_data
}

write.csv(rw, "input/rw_itrdb.csv", row.names = FALSE)
