####################################
### Decomposition for every draw ###
####################################
# Arriaga decomposition (age and cause) of the top - bottom gap repeated for
# every posterior draw to get 95% CrI.
# Run on a Linux server (plan(multicore)), from the repository root.

###################
### Environment ###
###################
library(data.table)
library(dplyr)
library(future.apply)
library(parallel)

source("R/lifeTable.R")
source("R/arriaga.R")

YEARS   <- c(2000, 2021)
N_DRAWS <- 1000
CHAINS  <- 0: 3

plan(multicore, workers = min(detectCores(), 120))

############
### Data ###
############
DeathCount <- fread("data/restricted/DeathCount_v5.csv") %>%
  filter(! Urbanicity_AST2023 == 9) %>%
  rename(COUNT_RAW = COUNT) %>% as.data.table()

MOIsubgroup <- fread("output/MOIsubgroup.csv")

# Cause fraction of the smoothed mortality in each age group
DCc <- fread("output/DCcauseSmoothed.csv") %>% rename(MOI = subgroup)
DCc <- DCc %>%
  left_join(DCc %>% aggregate(PROP_SMOOTHED ~ YEAR + SEX + MOI + AGP, sum) %>%
              rename(TOTAL_PROP = PROP_SMOOTHED),
            by = c("YEAR", "SEX", "MOI", "AGP")) %>%
  mutate(PROP_STANDARDIZED = PROP_SMOOTHED / TOTAL_PROP) %>% as.data.table()

#####################
### Decomposition ###
#####################
runChainSex <- function(no.chain, sex){
  sex_eng <- if(sex == 1) "male" else "female"
  DC_sex  <- DeathCount[SEX == sex]
  n_rows  <- nrow(DC_sex)
  PPd <- fread(paste0("output/posterior/PPd", sex_eng, "_", no.chain, ".csv.gz"))
  stopifnot(all(PPd$chain == no.chain))
  PPd[, row_id := seq_len(.N), by = draw]
  stopifnot(PPd[, .N, by = draw][, all(N == n_rows)])

  # Sum the draws of the counts within the top / bottom group
  key <- DC_sex[, .(row_id = .I, YEAR, AGP, MOI, POP)]
  key <- merge(key, MOIsubgroup[SEX == sex, .(MOI, subgroup)], by = "MOI")
  key <- key[YEAR %in% YEARS]
  sub <- merge(PPd[row_id %in% key$row_id, .(draw, row_id, d)], key, by = "row_id")
  rm(PPd); gc()

  agg <- sub[, .(POP = sum(POP), COUNT = sum(d)),
             keyby = .(draw, YEAR, AGP, MOI = subgroup)][, SEX := sex][]
  rm(sub); gc()

  draws <- sort(unique(agg$draw))
  cb    <- as.data.table(expand.grid(YEAR = YEARS, MOI = c("top", "bottom"),
                                     stringsAsFactors = F))

  out <- future_lapply(draws, function(dr){
    a  <- agg[draw == dr]
    lt <- rbindlist(lapply(seq_len(nrow(cb)), function(i)
      leCompute_dt(a, cb$YEAR[i], sex, cb$MOI[i])))
    if(nrow(lt) == 0) return(NULL)

    age <- rbindlist(lapply(YEARS, function(y)
      as.data.table(totalContribution(lt, y, sex))))
    cau <- rbindlist(lapply(YEARS, function(y)
      as.data.table(causeContribution(lt, DCc, y, sex))))

    e0 <- as.data.table(lt)[agp == 0, .(YEAR, MOI, ex)]
    list(age = age[, `:=`(CHAIN = no.chain, DRAW = dr)],
         cau = cau[, `:=`(CHAIN = no.chain, DRAW = dr)],
         e0  = e0[,  `:=`(CHAIN = no.chain, DRAW = dr, SEX = sex)])
  }, future.seed = T)

  out <- Filter(Negate(is.null), out)
  cat(sprintf("SEX %s CHAIN %d: %d / %d draws\n", sex_eng, no.chain, length(out), length(draws)))
  list(age = rbindlist(lapply(out, `[[`, "age")),
       cau = rbindlist(lapply(out, `[[`, "cau")),
       e0  = rbindlist(lapply(out, `[[`, "e0")))
}

res <- list()
for(s in 1: 2) for(ch in CHAINS){
  st <- Sys.time()
  res[[paste(s, ch)]] <- runChainSex(ch, s)
  cat("  elapsed:", round(difftime(Sys.time(), st, units = "mins"), 2), "min\n")
}

AGE <- rbindlist(lapply(res, `[[`, "age"))
CAU <- rbindlist(lapply(res, `[[`, "cau"))
E0  <- rbindlist(lapply(res, `[[`, "e0"))
plan(sequential)

dir.create("output/posterior", showWarnings = F, recursive = T)
fwrite(AGE, "output/posterior/decom_age_draws.csv.gz")
fwrite(CAU, "output/posterior/decom_cause_draws.csv.gz")
fwrite(E0,  "output/posterior/decom_e0_draws.csv.gz")

###############
### Summary ###
###############
qs <- function(x) as.list(setNames(
  c(mean(x, na.rm = T), median(x, na.rm = T),
    quantile(x, c(0.025, 0.975), na.rm = T, names = F), sum(!is.na(x))),
  c("mean", "median", "lower", "upper", "n_draw")))

# Top - bottom gap
gap <- dcast(E0, CHAIN + DRAW + SEX + YEAR ~ MOI, value.var = "ex")[, gap := top - bottom]
fwrite(gap[, qs(gap), keyby = .(SEX, YEAR)], "output/decom_gap_summary.csv")

# Age-specific contribution
fwrite(AGE[, qs(contribution), keyby = .(SEX, YEAR, agp)], "output/decom_age_summary.csv")

# Cause-specific contribution (Figure 3), summed over age within each draw
cau_tot <- CAU[, .(contribution = sum(contribution, na.rm = T)),
               keyby = .(CHAIN, DRAW, SEX, YEAR, CAUSE)]
fwrite(cau_tot[, qs(contribution), keyby = .(SEX, YEAR, CAUSE)], "output/decom_cause_summary.csv")

# Change 2000 -> 2021 (Figure 4), difference taken within each draw
chg <- dcast(cau_tot, CHAIN + DRAW + SEX + CAUSE ~ YEAR, value.var = "contribution")
setnames(chg, c("2000", "2021"), c("y2000", "y2021"))
chg[, change := y2021 - y2000]
fwrite(chg[, qs(change), keyby = .(SEX, CAUSE)], "output/decom_cause_change_summary.csv")
