## Started 09 November 2026
## By Ken

## Checking for multiple cores

rm(list = ls())
options(stringsAsFactors = FALSE)

#setwd
if(length(grep("victor", getwd())) > 0){
  setwd("/home/victor/projects/decoupleddynamics/analysis")
} else if(length(grep("Ken", getwd()) > 0)){
  setwd("/Users/Ken Michiko Samson/Documents/Temporal Ecology Lab/decoupleddynamics/analysis")
}

files <- list.files("../data/ITRDB", full.names = TRUE)
files <- files[grepl(files, pattern = '*.rwl$')]
file_name <- unname(sapply(files, function(x){return(regmatches(x, regexec("/([^/\\.]*)\\.rwl", x))[[1]][2])}))

rwref <- list()
counter <- 0

for(i in 1:length(files)){
  file_lines <- readLines(files[i])
  
  ids <- sapply(4:length(file_lines), function(x){return(substr(file_lines[x], 1, 8))})
  unique_ids <- unique(ids)
  unique_ids <- unique_ids[which(unique_ids != "        ")]
  unique_ids_nospace <- unique(sapply(ids, function(x){paste0(strsplit(x, " ")[[1]], collapse = "")}))
  unique_ids_nospace <- unique_ids_nospace[which(unique_ids != "        ")]
  unique_char <- sapply(unique_ids_nospace, nchar)
  
  id_sizes <- sort(unique(unique_char), decreasing = TRUE)
  
  for(j in 1:length(id_sizes)){
    same_size_idxs <- which(unique_char == id_sizes[j])
    trunc_ids <- sapply(unique_ids_nospace[same_size_idxs], function(x){return(substr(x, 1, id_sizes[j] - 1))})
    maybe_core_ids <- sapply(unique_ids_nospace[same_size_idxs], function(x){return(substr(x, id_sizes[j], id_sizes[j]))})
    
    maybe_tree_ids <- unique(trunc_ids)
    for(k in 1:length(maybe_tree_ids)){
      cores_same_tree <- which(trunc_ids == maybe_tree_ids[k])
      tree_core_id <- unique_ids_nospace[same_size_idxs[cores_same_tree]]
      if(length(cores_same_tree) == 1){
        tree_id <- tree_core_id
        core_id <- NA
      } else {
        tree_id <- maybe_tree_ids[k]
        core_id <- maybe_core_ids[cores_same_tree]
      }
      
      ref <- data.frame(
        "file_name" = file_name[i],
        "tree_core_id" = tree_core_id,
        "tree_id" = tree_id,
        "core_id" = unname(core_id)
      )
      
      counter <- counter + 1
      rwref[[counter]] <- ref
    }
  }
}

rwref <- do.call(rbind, rwref)
write.csv(rwref, "../data/ITRDB/_rw_summary.csv", row.names = FALSE)
