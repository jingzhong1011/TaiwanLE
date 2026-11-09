#########################
### Posterior summary ###
#########################
# Posterior summaries of township e0 and posterior contrasts. Every contrast
# is computed within each draw first and then summarised across the draws.
# Input: LEsmoothed_{male,female}.csv.gz from 04_e0_subset.R
# Run from the repository root.

###################
### Environment ###
###################
library(data.table)
library(dplyr)

#################
### Functions ###
#################
# Posterior mean, median, SD and 95% CrI of e0
calcPosteriorLEmean <- function(LEsmoothed, MOI = NULL, probs = c(0.025, 0.975)){
  dt <- as.data.table(LEsmoothed)

  out <- dt[, {
    q <- quantile(ex, probs = probs, na.rm = T, names = F)
    .(e0_mean   = mean(ex, na.rm = T),
      e0_median = as.numeric(stats::median(ex, na.rm = T)),
      e0_sd     = stats::sd(ex, na.rm = T),
      e0_lower  = q[1],
      e0_upper  = q[2],
      e0_min    = min(ex, na.rm = T),
      e0_max    = max(ex, na.rm = T),
      n_draw    = .N)
  }, keyby = .(YEAR, SEX, MOI, agp)]

  out[, max.diff := e0_max - e0_min]

  if(!is.null(MOI)) out <- merge(out, as.data.table(MOI), by = "MOI", all.x = T)

  return(out[])
}

# Difference of mean e0 between two sets of (MOI, YEAR), computed within draws
posteriorDiff <- function(LEsmoothed,
                          specA, specB,
                          agg    = c("mean", "none"),
                          label  = "contrast",
                          probs  = c(0.025, 0.975)){
  agg <- match.arg(agg)
  dt  <- as.data.table(LEsmoothed)

  specA <- as.data.table(specA)[, .(MOI, YEAR)]
  specB <- as.data.table(specB)[, .(MOI, YEAR)]

  pull <- function(spec){
    spec   <- unique(spec)
    n_unit <- nrow(spec)
    sub    <- merge(dt, spec, by = c("MOI", "YEAR"))
    agg    <- sub[, .(ex = mean(ex, na.rm = T), n_obs = .N), keyby = .(CHAIN, DRAW, SEX)]
    n_drop <- agg[n_obs < n_unit, .N]
    if(n_drop > 0L)
      message(sprintf("  [%s] %d / %d draws exclude (%.2f%%) from draws",
                      label, n_drop, nrow(agg), 100 * n_drop / nrow(agg)))
    agg[n_obs == n_unit][, n_obs := NULL][]
  }

  A <- pull(specA); setnames(A, "ex", "exA")
  B <- pull(specB); setnames(B, "ex", "exB")

  D <- merge(A, B, by = c("CHAIN", "DRAW", "SEX"))
  D[, diff := exA - exB]

  out <- D[, {
    q <- quantile(diff, probs = probs, na.rm = T, names = F)
    .(label     = label,
      diff_mean = mean(diff, na.rm = T),
      diff_med  = as.numeric(stats::median(diff, na.rm = T)),
      diff_sd   = stats::sd(diff, na.rm = T),
      diff_lower = q[1],
      diff_upper = q[2],
      p_gt0     = mean(diff > 0, na.rm = T),
      n_draw    = .N)
  }, keyby = .(SEX)]

  return(out[])
}

# Highest vs lowest township (by posterior mean) in one year
gapExtremeTownships <- function(LEsmoothed, LEmean, year, sex, ...){
  m <- as.data.table(LEmean)[YEAR == year & SEX == sex]
  hi <- m[which.max(e0_mean), .(MOI, YEAR)]
  lo <- m[which.min(e0_mean), .(MOI, YEAR)]
  message(sprintf("SEX %s, %s: highest MOI = %s (%.1f), lowest MOI = %s (%.1f)",
                  sex, year, hi$MOI, m[which.max(e0_mean), e0_mean],
                  lo$MOI, m[which.min(e0_mean), e0_mean]))
  posteriorDiff(LEsmoothed[SEX == sex], hi, lo,
                agg = "none", label = paste0("extreme_gap_", year), ...)
}

# Urban (cluster 1-3) vs rural (cluster 6-8) townships, unweighted mean
gapUrbanRural <- function(LEsmoothed, Urbanicity, year, sex, ...){
  u <- as.data.table(Urbanicity)
  urban <- u[Urbanicity_AST2023 %in% 1: 3, .(MOI)][, YEAR := year]
  rural <- u[Urbanicity_AST2023 %in% 6: 8, .(MOI)][, YEAR := year]
  posteriorDiff(LEsmoothed[SEX == sex], urban, rural,
                agg = "mean", label = paste0("urban_rural_", year), ...)
}

# Change of the township-average e0 between two years
changeOverTime <- function(LEsmoothed, moi_set, sex, year1 = 2000, year2 = 2021, ...){
  moi_set <- as.data.table(moi_set)[, .(MOI)]
  A <- copy(moi_set)[, YEAR := year2]
  B <- copy(moi_set)[, YEAR := year1]
  posteriorDiff(LEsmoothed[SEX == sex], A, B,
                agg = "mean", label = paste0("change_", year1, "_", year2), ...)
}

# Change of e0 between two years for every township
townshipChange <- function(LEsmoothed, sex, year1 = 2000, year2 = 2021, probs = c(0.025, 0.975)){
  dt <- as.data.table(LEsmoothed)[SEX == sex & YEAR %in% c(year1, year2)]

  w <- dcast(dt, CHAIN + DRAW + MOI ~ YEAR, value.var = "ex")
  setnames(w, c(as.character(year1), as.character(year2)), c("y1", "y2"))
  w <- w[!is.na(y1) & !is.na(y2)]
  w[, diff := y2 - y1]

  out <- w[, {
    q <- quantile(diff, probs = probs, na.rm = T, names = F)
    .(SEX        = sex,
      diff_mean  = mean(diff, na.rm = T),
      diff_sd    = stats::sd(diff, na.rm = T),
      diff_lower = q[1],
      diff_upper = q[2],
      p_gt0      = mean(diff > 0, na.rm = T),
      n_draw     = .N)
  }, keyby = .(MOI)]

  return(out[])
}

############
### Data ###
############
cols <- c("YEAR", "SEX", "MOI", "agp", "ex", "CHAIN", "DRAW")

LE <- rbindlist(list(
  fread("output/posterior/LEsmoothed_male.csv.gz",   select = cols),
  fread("output/posterior/LEsmoothed_female.csv.gz", select = cols)
))
gc()

LE[, MOI := as.numeric(MOI)]

# Check: 358 townships x 22 years in every draw, 4000 draws per sex
LE[, .N, keyby = .(SEX, CHAIN, DRAW)][N != 7876]
LE[, .(n_draw = uniqueN(paste(CHAIN, DRAW))), keyby = SEX]

DeathCount <- fread("data/restricted/DeathCount_v5.csv", header = T) %>%
  filter(! Urbanicity_AST2023 == 9)
MOI <- read.csv("data/MOI_ZH.csv", fileEncoding = "BIG-5")
Urbanicity <- unique(as.data.table(DeathCount)[, .(MOI = as.numeric(MOI), Urbanicity_AST2023)])
MOI_key    <- copy(as.data.table(MOI))[, MOI := as.numeric(MOI)][]

##########################
### Township LE (mean) ###
##########################
LEsmoothed_mean <- calcPosteriorLEmean(LE, MOI = MOI_key)
fwrite(LEsmoothed_mean, "output/LEsmoothed_mean.csv")

#################
### Contrasts ###
#################
gap_extreme <- rbindlist(list(
  gapExtremeTownships(LE, LEsmoothed_mean, 2000, 1),
  gapExtremeTownships(LE, LEsmoothed_mean, 2000, 2),
  gapExtremeTownships(LE, LEsmoothed_mean, 2021, 1),
  gapExtremeTownships(LE, LEsmoothed_mean, 2021, 2)
))
fwrite(gap_extreme, "output/gap_extreme_townships.csv")

gap_urbanrural <- rbindlist(list(
  gapUrbanRural(LE, Urbanicity, 2000, 1),
  gapUrbanRural(LE, Urbanicity, 2000, 2),
  gapUrbanRural(LE, Urbanicity, 2021, 1),
  gapUrbanRural(LE, Urbanicity, 2021, 2)
))
fwrite(gap_urbanrural, "output/gap_urban_rural.csv")

all_moi <- unique(LE[, .(MOI)])
change_overall <- rbindlist(list(
  changeOverTime(LE, all_moi, sex = 1),
  changeOverTime(LE, all_moi, sex = 2)
))
fwrite(change_overall, "output/change_2000_2021_overall.csv")

LEdiff <- rbindlist(list(
  townshipChange(LE, sex = 1),
  townshipChange(LE, sex = 2)
))
LEdiff <- merge(LEdiff, MOI_key, by = "MOI", all.x = T)
fwrite(LEdiff, "output/township_change_2000_2021.csv")

##########################
### Manuscript numbers ###
##########################
fmt <- function(m, l, u) sprintf("%.1f (95%% CrI, %.1f-%.1f)", m, l, u)

manuscript_numbers <- rbindlist(list(
  gap_extreme[,    .(quantity = label, SEX, value = fmt(diff_mean, diff_lower, diff_upper))],
  gap_urbanrural[, .(quantity = label, SEX, value = fmt(diff_mean, diff_lower, diff_upper))],
  change_overall[, .(quantity = label, SEX, value = fmt(diff_mean, diff_lower, diff_upper))]
))
fwrite(manuscript_numbers, "output/manuscript_numbers.csv")
print(manuscript_numbers)
