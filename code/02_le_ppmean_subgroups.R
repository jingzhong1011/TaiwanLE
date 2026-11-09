#######################################
### LE (PPmean) and township groups ###
#######################################
# Life tables from the posterior mean death counts (PPmean), the top / bottom
# 10% township subgroups, subgroup life tables and the age decomposition of
# the top - bottom gap.
# Input: PPmean_male.csv, PPmean_female.csv from 01_icar_poisson.ipynb
# Run from the repository root.

###################
### Environment ###
###################
library(data.table)
library(dplyr)
library(future.apply)
library(parallel)

source("R/lifeTable.R")
source("R/arriaga.R")

############
### Data ###
############
DeathCount <- fread("data/restricted/DeathCount_v5.csv") %>%
  filter(! Urbanicity_AST2023 == 9) %>%
  rename(COUNT_RAW = COUNT)

PPmean_male <- fread("output/posterior/PPmean_male.csv")
PPmean_female <- fread("output/posterior/PPmean_female.csv")

MOI <- read.csv("data/MOI_ZH.csv", fileEncoding = "BIG5") %>% select(MOI, County_Township_ZH)

# Replace observed counts with posterior mean counts
DeathCount_input <- bind_rows(
  DeathCount %>%
    filter(SEX == 1) %>%
    mutate(COUNT = PPmean_male$d),
  DeathCount %>%
    filter(SEX == 2) %>%
    mutate(COUNT = PPmean_female$d)
)

########################
### LE (PPmean) calc ###
########################
# Set parallel environment
no_cores <- detectCores() - 1
if(no_cores < 1) no_cores <- 1
plan(multisession, workers = no_cores) # use plan(multicore) on Linux server

DeathCount_input_dt <- as.data.table(DeathCount_input)

year <- 2000: 2021
sex <- 1: 2
moi <- levels(as.factor(DeathCount_input_dt$MOI))
Combinations <- data.table(expand.grid(YEAR = year, SEX = sex, MOI = moi))

results_list <- future_lapply(1: nrow(Combinations), function(i){
  row <- Combinations[i, ]
  leCompute_dt(DeathCount_input_dt, row$YEAR, row$SEX, row$MOI)
}, future.seed = T)

LEsmoothed_PPmean <- rbindlist(results_list)
plan(sequential)
rm(results_list)

write.csv(LEsmoothed_PPmean, "output/LEsmoothed_PPmean.csv", row.names = F)

LEsmoothed_PPmean <- fread("output/LEsmoothed_PPmean.csv") %>% left_join(MOI, by = "MOI")

########################
### Top / bottom 10% ###
########################
# Rank townships by e0 in each year; keep townships that stay in the same
# group in 2000 and 2021
calcTopBottom <- function(data, year, sex, percent){
  data1 <- data %>% filter(SEX == sex & agp == 0) %>%
    filter(YEAR %in% year) %>%
    group_by(YEAR) %>%
    mutate(rank = rank(ex, ties.method = "average")) %>%
    ungroup() %>%
    mutate(subgroup = case_when(rank > length(unique(data$MOI)) - ceiling(length(unique(data$MOI)) * percent) ~ "top",
                                rank <= ceiling(358 * percent) ~ "bottom",
                                T ~ "normal")) %>% select(MOI, subgroup) %>%
    filter(! subgroup == "normal") %>%
    unique() %>%
    mutate(SEX = sex)
  return(data1)
}

MOIsubgroup <- bind_rows(
  calcTopBottom(LEsmoothed_PPmean, c(2000, 2021), 1, 0.1) %>%
    add_count(MOI) %>%
    filter(n == 1) %>%
    select(-n),
  calcTopBottom(LEsmoothed_PPmean, c(2000, 2021), 2, 0.1) %>%
    add_count(MOI) %>%
    filter(n == 1) %>%
    select(-n)
)
write.csv(MOIsubgroup, "output/MOIsubgroup.csv", row.names = F)

#####################
### LE (subgroup) ###
#####################
DeathCount_topbottom <- DeathCount_input %>%
  filter(YEAR %in% c(2000, 2021)) %>%
  left_join(MOIsubgroup, by = c("MOI", "SEX")) %>%
  filter(! c(subgroup == "normal")) %>%
  group_by(YEAR, SEX, subgroup, AGP) %>%
  summarise(POP = sum(POP),
            COUNT = sum(COUNT)) %>%
  rename(MOI = subgroup)

Combinations <- data.table(expand.grid(YEAR = 2000: 2021, SEX = 1: 2,
                                       MOI = levels(as.factor(DeathCount_topbottom$MOI))))

plan(multisession, workers = no_cores)
results_list <- future_lapply(1: nrow(Combinations), function(i){
  row <- Combinations[i, ]
  leCompute_dt(DeathCount_topbottom, row$YEAR, row$SEX, row$MOI)
}, future.seed = T)

LEsmoothed_topbottom <- rbindlist(results_list)
plan(sequential)
rm(results_list)

write.csv(LEsmoothed_topbottom, "output/LEsmoothed_topbottom.csv", row.names = F)

LEsmoothed_topbottom <- fread("output/LEsmoothed_topbottom.csv")

#########################
### Age decomposition ###
#########################
LEsmoothed_topbottom_ageDecom <- bind_rows(
  totalContribution(LEsmoothed_topbottom, 2000, 1),
  totalContribution(LEsmoothed_topbottom, 2000, 2),
  totalContribution(LEsmoothed_topbottom, 2021, 1),
  totalContribution(LEsmoothed_topbottom, 2021, 2)
)

write.csv(LEsmoothed_topbottom_ageDecom, "output/LEsmoothed_topbottom_ageDecom.csv", row.names = F)
