##############
### LE raw ###
##############
# Life tables from the observed (unsmoothed) death counts, for the
# comparison of modelled and observed LE.
# Run from the repository root.

###################
### Environment ###
###################
library(data.table)
library(dplyr)
library(parallel)

source("R/lifeTable.R")

############
### Data ###
############
DeathCount <- fread("data/restricted/DeathCount_v5.csv", header = T) %>%
  filter(! Urbanicity_AST2023 == 9) %>%
  as.data.table()

#####################
### LE (raw) calc ###
#####################
combos <- as.data.table(expand.grid(YEAR = sort(unique(DeathCount$YEAR)),
                                    SEX  = sort(unique(DeathCount$SEX)),
                                    MOI  = sort(unique(DeathCount$MOI))))

no_cores <- max(1, detectCores() - 1)
res <- mclapply(seq_len(nrow(combos)), function(i)
  leCompute_dt(DeathCount, combos$YEAR[i], combos$SEX[i], combos$MOI[i]),
  mc.cores = no_cores)

n_null <- sum(!sapply(res, is.data.table))
if(n_null > 0) warning(sprintf("%d / %d life tables returned NULL", n_null, nrow(combos)))

LEraw <- rbindlist(Filter(is.data.table, res))

# Check
e0 <- LEraw[agp == 0]
cat("Rows:", nrow(LEraw), "(expected", nrow(combos) * 19, ")\n")
cat("e0 range:", paste(round(range(e0$ex, na.rm = T), 1), collapse = " - "), "\n")
cat("e0 > 100:", e0[ex > 100, .N], "  e0 NA:", e0[is.na(ex), .N], "\n")

fwrite(LEraw, "output/LEraw.csv")

###############################
### Modelled vs observed LE ###
###############################
# Pearson r in the Discussion (Supplementary Figure S7)
if(file.exists("output/LEsmoothed_mean.csv")){
  cmp <- merge(fread("output/LEsmoothed_mean.csv")[, .(YEAR, SEX, MOI, model = e0_mean)],
               e0[, .(YEAR, SEX, MOI, observed = ex)],
               by = c("YEAR", "SEX", "MOI"))
  for(s in 1: 2)
    cat(sprintf("Pearson r (%s): %.3f, n = %d\n",
                c("male", "female")[s],
                cmp[SEX == s & is.finite(observed), cor(model, observed)],
                cmp[SEX == s & is.finite(observed), .N]))
}
