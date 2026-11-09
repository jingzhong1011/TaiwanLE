#################
### e0 subset ###
#################
# Keep e0 (agp == 0) from the per-draw life tables of 03_le_draws_server.R
# and combine them into one file per sex.
# Run from the repository root.

###################
### Environment ###
###################
library(data.table)

#################
### e0 subset ###
#################
e0subset <- function(sex_eng){
  file_pattern <- file.path("output/posterior/LEdraws", paste0("LEsmoothed_", sex_eng, "_chain*_draw*.csv"))
  file_list <- Sys.glob(file_pattern)

  filtered_data_list <- list()
  for(i in seq_along(file_list)){
    file_path <- file_list[i]
    cat("Data processing: (", i, "/", length(file_list), "): ", basename(file_path), "\n")

    tryCatch({
      dt <- fread(file_path)
      if("agp" %in% names(dt)){
        filtered_data_list[[length(filtered_data_list) + 1]] <- dt[agp == 0]
      }else{
        cat("Warning: File '", basename(file_path), "' did not find agp, skip this file.\n")
      }
    }, error = function(e){
      cat("Error on file '", basename(file_path), e$message, "\n")
    })
  }

  if(length(filtered_data_list) > 0){
    combined_dt <- rbindlist(filtered_data_list, fill = T)
    cat("dim of data (row, column):", dim(combined_dt), "\n")
  }else{
    combined_dt <- data.table()
    cat("\nNo data with agp == 0 was found.\n")
  }

  write.csv(combined_dt, gzfile(paste0("output/posterior/LEsmoothed_", sex_eng, ".csv.gz")), row.names = F)
}

e0subset("male")
e0subset("female")
