########################
### Spatial analysis ###
########################
# Age-standardized mortality (WHO standard population) by township and
# Moran's I of its spatial autocorrelation.
# Run from the repository root.

###################
### Environment ###
###################
library(data.table)
library(dplyr)
library(sf)
library(spdep)
library(stringr)

############
### Data ###
############
DeathCount <- fread("data/restricted/DeathCount_v5.csv", header = T) %>%
  filter(! Urbanicity_AST2023 == 9) %>%
  mutate(MORT = COUNT / POP)

TWbord <- TWmap::bord_town %>%
  select(town, geometry) %>%
  filter(! str_sub(town, 1, 2) == "09") %>%
  rename(MOI = town) %>%
  mutate(MOI = as.numeric(MOI))

# WHO standard population
StandPop <- data.frame(AGP = c(0, 1, seq(5, 85, by = 5)),
                       Percent = c(17917, 70652, 86870, 85970, 84670,
                                   82171, 79272, 76073, 71475, 65877,
                                   60379, 53681, 45484, 37187, 29590,
                                   22092, 15195, 9097, 6348)) %>% mutate(Percent = Percent / 1e06)

###########
### MSR ###
###########
# Age-standardized mortality (MOI)
MortStand <- DeathCount %>% left_join(StandPop, by = "AGP") %>%
  mutate(MP = MORT * Percent * 100000) %>%
  aggregate(MP ~ MOI + YEAR + SEX, FUN = sum) %>%
  rename(MORT_STAND = MP)

# Age-standardized mortality (National)
MortStand_nat <- DeathCount %>%
  group_by(YEAR, SEX, AGP) %>%
  summarise(COUNT = sum(COUNT),
            POP = sum(POP)) %>%
  mutate(MORT = COUNT / POP) %>%
  left_join(StandPop, by = "AGP") %>%
  mutate(MP = MORT * Percent * 100000) %>%
  aggregate(MP ~ YEAR + SEX, FUN = sum) %>%
  rename(MORT_STAND = MP)

# Age-standardized mortality (Both sex, coded as SEX = 3)
MortStand_bothsex <- DeathCount %>%
  aggregate(COUNT ~ YEAR + MOI + AGP, sum) %>%
  left_join(DeathCount %>% aggregate(POP ~ YEAR + MOI + AGP, sum), by = c("YEAR", "MOI", "AGP")) %>%
  left_join(StandPop, by = "AGP") %>%
  mutate(MORT = COUNT / POP,
         MP = MORT * Percent * 100000) %>%
  aggregate(MP ~ MOI + YEAR, FUN = sum) %>%
  rename(MORT_STAND = MP) %>%
  mutate(SEX = 3)

MortStand_bySex <- bind_rows(MortStand, MortStand_bothsex)

write.csv(MortStand_bySex, "output/MortStand_bySex.csv", row.names = F)

#################
### Moran's I ###
#################
moran.matrix <- matrix(NA, nrow = 22, ncol = 4)
colnames(moran.matrix) = c("YEAR", "Male", "Female", "Both")
for(i in 1: 22){
  moran.matrix[i, 1] = i + 1999
  for(j in 1: 3){
    data <- MortStand_bySex %>%
      left_join(TWbord, by = "MOI") %>%
      filter(YEAR == i + 1999 & SEX == j) %>% st_as_sf()
    NB <- poly2nb(data, queen = T)
    NBlist <- nb2listw(NB, zero.policy = T)
    moran <- moran.test(data$MORT_STAND, listw = NBlist, zero.policy = T)
    moran.matrix[i, j + 1] = moran$estimate[1]
  }
}

Moran <- data.frame(moran.matrix)
Moran
