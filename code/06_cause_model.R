###################
### Cause model ###
###################
# Smooth cause-specific mortality of the top / bottom / other township groups
# (from 02_le_ppmean_subgroups.R) with a Poisson model in INLA.
# Run from the repository root.

###################
### Environment ###
###################
library(INLA)
library(dplyr)
library(data.table)
library(ggplot2)

############
### Data ###
############
DeathCount_cause <- fread("data/restricted/DeathCount_cause.csv")
MOIsubgroup <- fread("output/MOIsubgroup.csv")

DC_cause <- DeathCount_cause %>% left_join(MOIsubgroup, by = c("MOI", "SEX")) %>%
  mutate(subgroup = case_when(is.na(subgroup) == T ~ "normal", T ~ subgroup)) %>%
  group_by(YEAR, SEX, AGP, CAUSE, subgroup) %>%
  summarise(COUNT = sum(COUNT),
            POP = sum(POP)) %>%
  ungroup() %>%
  mutate(PROP = COUNT / POP)

#############
### Model ###
#############
causeModeling <- function(data){
  results <- list()
  sex_idx = c(1: 2)
  cause_idx <- unique(data$CAUSE)

  for(i in sex_idx){
    for(j in cause_idx){
      cat("Processing: SEX =", i, "CAUSE =", j, "\n")
      dt <- data %>% filter(SEX == i & CAUSE == j)
      dt2 <- dt %>%
        mutate(AGP2 = AGP,
               AGP3 = AGP,
               YEAR2 = YEAR,
               YEAR3 = YEAR,
               subgroup2 = subgroup,
               subgroup3 = subgroup)

      # No death of this cause
      if(sum(dt2$COUNT) == 0){
        dt_fitted <- dt %>% mutate(PROP_SMOOTHED = 0)
        results[[paste(i, j, sep = "_")]] <- dt_fitted
        next
      }

      fit <- inla(COUNT ~
                    f(AGP, model = "rw1", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))) +
                    f(YEAR, model = "rw1", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))) +
                    f(subgroup, model = "iid", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))) +
                    f(YEAR2, AGP2, model = "rw1", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))) +
                    f(subgroup2, AGP3, model = "iid", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))) +
                    f(subgroup3, YEAR3, model = "rw1", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001)))),
                  family = "poisson", data = dt2, E = POP)

      dt_fitted <- dt %>% mutate(PROP_SMOOTHED = fit$summary.fitted.values$mean)
      results[[paste(i, j, sep = "_")]] <- dt_fitted
    }
  }
  return(results)
}

output <- causeModeling(DC_cause)
output2 <- bind_rows(output)

write.csv(output2, "output/DCcauseSmoothed.csv", row.names = F)

#########################
### Check (ASMR plot) ###
#########################
output2 <- fread("output/DCcauseSmoothed.csv")

# WHO standard population
stand_pop <- data.frame(AGP = c(0, 1, seq(5, 85, by = 5)),
                        Percent = c(17917, 70652, 86870, 85970, 84670,
                                    82171, 79272, 76073, 71475, 65877,
                                    60379, 53681, 45484, 37187, 29590,
                                    22092, 15195, 9097, 6348)) %>% mutate(Percent = Percent / 1e06)

output2_stand <- output2 %>% left_join(stand_pop, by = "AGP") %>%
  mutate(PROP_SMOOTHED_STAND = PROP_SMOOTHED * Percent * 100000,
         PROP_STAND = PROP * Percent * 100000) %>%
  group_by(YEAR, SEX, subgroup, CAUSE) %>%
  summarise(PROP_SMOOTHED_STAND = sum(PROP_SMOOTHED_STAND),
            PROP_STAND = sum(PROP_STAND)) %>% ungroup()

# Line: smoothed, point: observed
ggplot(data = output2_stand %>% filter(SEX == 1)) +
  geom_line(aes(x = YEAR, y = PROP_SMOOTHED_STAND, group = subgroup, colour = subgroup)) +
  geom_point(aes(x = YEAR, y = PROP_STAND, group = subgroup, colour = subgroup), pch = 1) +
  theme_bw() +
  facet_wrap(CAUSE ~ ., scales = "free") +
  xlab("Year") +
  ylab("Age-standardized Mortality Rate (per 100,000 population)") +
  ggtitle("Male") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10, family = "Times New Roman"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.spacing.y = unit(0, "cm"))

ggplot(data = output2_stand %>% filter(SEX == 2)) +
  geom_line(aes(x = YEAR, y = PROP_SMOOTHED_STAND, group = subgroup, colour = subgroup)) +
  geom_point(aes(x = YEAR, y = PROP_STAND, group = subgroup, colour = subgroup), pch = 1) +
  theme_bw() +
  facet_wrap(CAUSE ~ ., scales = "free") +
  xlab("Year") +
  ylab("Age-standardized Mortality Rate (per 100,000 population)") +
  ggtitle("Female") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10, family = "Times New Roman"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.spacing.y = unit(0, "cm"))
