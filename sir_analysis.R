# Libraries --------------------------
library("netmeta")
library("sir")
library("dplyr")


# helper functions for database analysis -----------------------------------------
get_ranks <- function(obj, sel) {
  nma <- obj$netobj
  recid <- obj$recid
  #
  ok <- netrank(nma, method = "P-score",
                small.values =
                  sel$small.values[sel$recid == recid])
  #
  ok$Pscore.random
}


sir_func <- function(x)
  sir(x, trts = seq_len(length(x)))$sir

ns_func <- function(obj) {
  nma <- obj$netob
  ns <- obj$netobj$k
  #
  ns
}

prop_func <- function(obj, alpha = 0.05) {

  mean(obj$netobj$pval.random < alpha, na.rm = T)

}

range_func <- function(pscores) {

  max(pscores) - min(pscores)

}


# Calculate SIR from the analysed datasets -------------------------------------

# First, use anlayse.R to analyse all datasets yourself using nmadb
# Or, use download the .rds object provided

# Read in all nmas analysis object and study information
allNMAs <- readRDS("allNMAs.rds") # first, set working directory to location where you saved the object
study_info <- readRDS("studyList.rds") %>% select(-contains("..choice"))

# find only the ones without errors
ran <- unlist(lapply(allNMAs, function(x) length(x) > 1))
ran_NMAs <- allNMAs[which(ran)] # 268 ran
error_NMAs <- allNMAs[-c(which(ran))] # 28 of them had errors

# Get the record IDs of the NMAs which ran
ran_NMAs_record <- unlist(lapply(ran_NMAs, function(x) x$recid))

# figure out if effect is beneficial or harmful
studies_ran <- subset(study_info,
                      Record.ID %in% ran_NMAs_record) %>%
  select(Record.ID, Harmful.Beneficial.Outcome) %>%
  rename(recid = Record.ID) %>%
  mutate(small.values =
           if_else(Harmful.Beneficial.Outcome == "Beneficial",
                   "undesirable", "desirable"))

# Calculate SIRs for everything
pscores <- lapply(ran_NMAs, get_ranks, sel = studies_ran)
sirs <- sapply(pscores, function(x) sir(x)$sir)
summary(sirs)
IQR(sirs)

hist(sirs, main = "",
     xlab= "SIR", col = "grey95", breaks = 20, xlim = c(0,1))
abline(v = quantile(sirs, c(0.5)), col = "navy", lwd = 2)
legend("topleft", text.col = c("navy", "black", "black"),
       legend = c("Median = 0.67706", "Min = 0.096408", "Max = 0.999762"), bty = "n")

# look at relationship with number of treatments
ntrts <- sapply(pscores, length)
cor(sirs, ntrts)
plot(ntrts, sirs, ylim = c(0,1), xlim = c(4, 46),
     xlab = "Number of Treatments", ylab = "SIR")

# Relationship with proportion of significantly different treatment comparisons
alphas <- rev(c(0.1, 0.05, 0.01, 0.005))
proportions <- matrix(nrow = length(alphas), ncol = length(sirs))

for (i in 1:length(alphas))
  proportions[i, ] <- sapply(ran_NMAs, prop_func, alpha = alphas[i])

plot(proportions[3,], sirs,
     xlab = "Proportion of treatment comparisons that are significant at 5% level",
     ylab = 'SIR', ylim = c(0,1))

# Relationship with range of p-scores
prange <- sapply(pscores, range_func)
plot(prange, sirs, xlim = c(0,1), ylim = c(0,1),
     xlab = "max(P-scores) - min(p-scores)", ylab = "SIR")
