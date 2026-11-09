###################################
### LE for every posterior draw ###
###################################
# Life tables for 4 chains x 1000 draws x 358 townships x 22 years, by sex.
# Run on a Linux server (plan(multicore)), from the repository root.
# Input: PPd{male,female}_{chain}.csv.gz from 01_icar_poisson.ipynb
# Output: output/posterior/LEdraws/LEsmoothed_{sex}_chain{c}_draw{d}.csv

###################
### Environment ###
###################
library(data.table)
library(dplyr)
library(future.apply)
library(parallel)

source("R/lifeTable.R")

######################
### LE calculation ###
######################
# Keep this function at top level: if defined inside doLECalc, future would
# also export the enclosing frame (the PPd draws, several GB) to every worker.
compute_for_row <- function(idx, combo_dt, data_to_use, chain_id, draw_id){
  row <- combo_dt[idx, ]
  result <- leCompute_dt(data_to_use, row$YEAR, row$SEX, row$MOI) %>%
    mutate(CHAIN = chain_id, DRAW = draw_id)
  return(result)
}

doLECalc <- function(no.chain, sex){
  sex_eng = case_when(sex == 1 ~ "male", sex == 2 ~ "female")
  no_cores <- detectCores()
  if(no_cores < 1) no_cores <- 1

  plan(multicore, workers = min(no_cores, 120))
  cat("Using", no_cores, "cores for parallel processing.\n")
  on.exit(plan(sequential), add = T)

  DeathCount <- fread("data/restricted/DeathCount_v5.csv") %>%
    filter(! Urbanicity_AST2023 == 9) %>%
    rename(COUNT_RAW = COUNT) %>% filter(SEX == sex) %>% as.data.table()

  year <- as.numeric(levels(as.factor(DeathCount$YEAR)))
  sex <- as.numeric(levels(as.factor(DeathCount$SEX)))
  moi <- levels(as.factor(DeathCount$MOI))
  Combinations <- data.table(expand.grid(YEAR = year, SEX = sex, MOI = moi))

  n_rows <- dim(DeathCount)[1]
  n_draws <- 1000

  dir.create("output/posterior/LEdraws", showWarnings = F, recursive = T)

  PPd <- fread(paste0("output/posterior/PPd", sex_eng, "_", no.chain, ".csv.gz"))
  stopifnot(all(PPd$chain == no.chain))
  setindex(PPd, draw)

  for(draw in (1: n_draws) - 1){
    cat("SEX = ", sex_eng, "CHAIN =", no.chain, "DRAW =", draw, "\n")
    startTime <- Sys.time()
    dr = draw

    # Replace observed counts with the counts of this draw
    current_counts <- PPd[draw == dr, d]
    stopifnot(length(current_counts) == n_rows)

    test_dt <- copy(DeathCount)
    test_dt[, COUNT := current_counts]

    # Compute every YEAR x MOI life table in parallel
    results_list <- future_lapply(1: nrow(Combinations),
                                  compute_for_row,
                                  combo_dt = Combinations,
                                  data_to_use = test_dt,
                                  chain_id = no.chain,
                                  draw_id = draw)

    valid_results <- Filter(Negate(is.null), results_list)
    if(length(valid_results) > 0){
      LEsmoothed <- rbindlist(valid_results)
      LEsmoothed[, CHAIN := no.chain]
      LEsmoothed[, DRAW := draw]

      write.csv(LEsmoothed, paste0("output/posterior/LEdraws/LEsmoothed_", sex_eng, "_chain", no.chain, "_draw", draw, ".csv"))
    }else{
      cat("  WARN: No valid results produced for chain", no.chain, "draw", draw, "\n")
    }
    LEsmoothed <- NULL
    results_list <- NULL
    valid_results <- NULL
    costTime <- Sys.time() - startTime
    cat("  Time taken:", costTime, "\n")
    gc()
  }
}

###########
### Run ###
###########
doLECalc(0, 1)
doLECalc(1, 1)
doLECalc(2, 1)
doLECalc(3, 1)

doLECalc(0, 2)
doLECalc(1, 2)
doLECalc(2, 2)
doLECalc(3, 2)
