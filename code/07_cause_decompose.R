###########################
### Cause decomposition ###
###########################
# Arriaga decomposition by age and cause (point estimates from PPmean):
#   1. top - bottom gap in each year
#   2. 2000 -> 2021 change within the top and the bottom group
# Input: LEsmoothed_topbottom.csv (02_le_ppmean_subgroups.R),
#        DCcauseSmoothed.csv (06_cause_model.R)
# Run from the repository root.

###################
### Environment ###
###################
library(data.table)
library(dplyr)

source("R/arriaga.R")

############
### Data ###
############
LEsmoothed_topbottom <- fread("output/LEsmoothed_topbottom.csv")

# Cause fraction of the smoothed mortality in each age group
DeathCount_cause <- fread("output/DCcauseSmoothed.csv") %>% rename(MOI = subgroup)
DeathCount_cause <- DeathCount_cause %>%
  left_join(DeathCount_cause %>% aggregate(PROP_SMOOTHED ~ YEAR + SEX + MOI + AGP, sum) %>% rename(TOTAL_PROP = PROP_SMOOTHED),
            by = c("YEAR", "SEX", "MOI", "AGP")) %>%
  mutate(PROP_STANDARDIZED = PROP_SMOOTHED / TOTAL_PROP)

########################
### Top - bottom gap ###
########################
LEsmoothed_topbottom_causeDecom <- do.call(rbind, apply(expand.grid(YEAR = 2000: 2021, SEX = c(1: 2)), 1,
                                             function(row) causeContribution(LEsmoothed_topbottom, DeathCount_cause, row["YEAR"], row["SEX"])))

write.csv(LEsmoothed_topbottom_causeDecom, "output/LEsmoothed_topbottom_causeDecom.csv")

##########################
### Change 2000 - 2021 ###
##########################
totalContribution_byYear <- function(data, moi, sex){
  data1 <- data %>% filter(MOI == moi & YEAR == 2000 & SEX == sex)
  data2 <- data %>% filter(MOI == moi & YEAR == 2021 & SEX == sex)
  record <- c()
  agp = c(0, 1, seq(5, 85, by = 5))
  for(i in 1: 19){
    # Direct and indirect (interaction) effect
    DE = data1$lx[i] / 100000 * (data2$Lx[i] / data2$lx[i] - data1$Lx[i] / data1$lx[i])
    IE = case_when(i == 19 ~ 0,
                   T ~ ((data1$lx[i] * (data2$lx[i+1] / data2$lx[i])) - data1$lx[i+1]) / 100000 * data2$ex[i+1])
    record[i] = DE + IE
  }
  contribution = data.frame(MOI = moi, SEX = sex, agp = agp, contribution = record)
  return(contribution)
}

causeContribution_byYear <- function(data, data_cause, moi, sex){
  C = totalContribution_byYear(data, moi, sex)
  agp = c(0, 1, seq(5, 85, by = 5))
  cause_factor <- levels(as.factor(data_cause$CAUSE))

  result <- matrix(nrow = length(cause_factor), ncol = length(agp))
  for(i in 1: length(agp)){
    dataA = data %>% filter(YEAR == 2000 & SEX == sex & MOI == moi)
    dataB = data %>% filter(YEAR == 2021 & SEX == sex & MOI == moi)

    for(j in 1: length(cause_factor)){
      data_causeA <- data_cause %>%
        filter(YEAR == 2000 & SEX == sex & MOI == moi & AGP == agp[i] & CAUSE == cause_factor[j])
      data_causeB <- data_cause %>%
        filter(YEAR == 2021 & SEX == sex & MOI == moi & AGP == agp[i] & CAUSE == cause_factor[j])

      mA = dataA$mx[i]
      mB = dataB$mx[i]

      RA = data_causeA$PROP_STANDARDIZED
      RB = data_causeB$PROP_STANDARDIZED

      # Share of the age-specific contribution attributable to this cause
      Prop = (RB * mB - RA * mA) / (mB - mA)

      C_i = C$contribution[i] * Prop
      result[j, i] = C_i
    }
  }
  output <- data.frame(result) %>% mutate(CAUSE = cause_factor)
  colnames(output) <-  c(0, 1, seq(5, 85, by = 5), "CAUSE")
  output <- melt(setDT(output), id.vars = "CAUSE", variable.name = "agp", value.name = "contribution") %>%
    mutate(MOI = moi, SEX = sex) %>%
    select(MOI, SEX, agp, CAUSE, contribution)
  return(output)
}

LEsmoothed_topbottom_byYear_ageDecom <- bind_rows(
  totalContribution_byYear(LEsmoothed_topbottom, "top", 1),
  totalContribution_byYear(LEsmoothed_topbottom, "top", 2),
  totalContribution_byYear(LEsmoothed_topbottom, "bottom", 1),
  totalContribution_byYear(LEsmoothed_topbottom, "bottom", 2)
)

write.csv(LEsmoothed_topbottom_byYear_ageDecom, "output/LEsmoothed_topbottom_byYear_ageDecom.csv", row.names = F)

LEsmoothed_topbottom_byYear_causeDecom <- bind_rows(
  causeContribution_byYear(LEsmoothed_topbottom, DeathCount_cause, "top", 1),
  causeContribution_byYear(LEsmoothed_topbottom, DeathCount_cause, "top", 2),
  causeContribution_byYear(LEsmoothed_topbottom, DeathCount_cause, "bottom", 1),
  causeContribution_byYear(LEsmoothed_topbottom, DeathCount_cause, "bottom", 2)
)

write.csv(LEsmoothed_topbottom_byYear_causeDecom, "output/LEsmoothed_topbottom_byYear_causeDecom.csv", row.names = F)

# Total contribution by cause
LEsmoothed_topbottom_byYear_causeDecom %>%
  aggregate(contribution ~ MOI + SEX + CAUSE, sum) %>%
  mutate(contribution = round(contribution, 3)) %>%
  arrange(SEX, MOI, desc(contribution))
