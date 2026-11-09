## Arriaga decomposition of the life-expectancy gap between the "top" and
## "bottom" township subgroups, by age (totalContribution) and by cause
## (causeContribution). Shared by 02, 07 and 08.
library(data.table); library(dplyr)   ## melt resolves to data.table::melt

totalContribution <- function(data, year, sex){
  data1 <- data %>% subset(MOI == "bottom" & YEAR == year & SEX == sex)
  data2 <- data %>% subset(MOI == "top" & YEAR == year & SEX == sex)
  record <- c()
  agp = c(0, 1, seq(5, 85, by = 5))
  for(i in 1: 19){
    DE = data1$lx[i] / 100000 * (data2$Lx[i] / data2$lx[i] - data1$Lx[i] / data1$lx[i])
    IE = case_when(i == 19 ~ 0,
                   T ~ ((data1$lx[i] * (data2$lx[i+1] / data2$lx[i])) - data1$lx[i+1]) / 100000 * data2$ex[i+1])
    
    record[i] = DE + IE
    
  }
  contribution = data.frame(YEAR = year, SEX = sex, agp = agp, contribution = record)
  return(contribution)
}

causeContribution <- function(data, data_cause, year, sex){
  C = totalContribution(data, year, sex)
  agp = c(0, 1, seq(5, 85, by = 5))
  cause_factor <- levels(as.factor(data_cause$CAUSE))
  
  result <- matrix(nrow = length(cause_factor), ncol = length(agp))
  for(i in 1: length(agp)){

    dataA = data %>% filter(YEAR == year & SEX == sex & MOI == "top")
    dataB = data %>% filter(YEAR == year & SEX == sex & MOI == "bottom")
    
     for(j in 1: length(cause_factor)){

      data_causeA <- data_cause %>% 
        filter(YEAR == year & SEX == sex & MOI == "top" & AGP == agp[i] & CAUSE == cause_factor[j])
      data_causeB <- data_cause %>% 
        filter(YEAR == year & SEX == sex & MOI == "bottom" & AGP == agp[i] & CAUSE == cause_factor[j])
      
      mA = dataA$mx[i]
      mB = dataB$mx[i]
      
      RA = data_causeA$PROP_STANDARDIZED
      RB = data_causeB$PROP_STANDARDIZED
      
      Prop = (RB * mB - RA * mA) / (mB - mA)
      
      
      C_i = C$contribution[i] * Prop
      result[j, i] = C_i
    }
  }
  output <- data.frame(result) %>% mutate(CAUSE = cause_factor)
  colnames(output) <-  c(0, 1, seq(5, 85, by = 5), "CAUSE")
  output <- melt(setDT(output), id.vars = "CAUSE", variable.name = "agp", value.name = "contribution") %>% 
    mutate(YEAR = year, SEX = sex) %>% 
    select(YEAR, SEX, agp, CAUSE, contribution)
  return(output)
}
