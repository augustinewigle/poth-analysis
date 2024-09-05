# Libraries --------------------------
library("netmeta")
library("poth")
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


poth_func <- function(x)
  poth(x, trts = seq_len(length(x)))$poth

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

q_func <- function(obj) {

  nma <- obj$netob

  q_overall <- obj$netobj$pval.Q
  q_het <- obj$netobj$pval.Q.heterogeneity
  q_inc <- obj$netobj$pval.Q.inconsistency

  return(c(q_overall, q_het, q_inc))

}

get_sm <- function(obj) {

  obj$netobj$sm

}

tau_func <- function(obj) {

  obj$netobj$tau

}

# Calculate POTH from the analysed datasets -------------------------------------

# First, use analyse.R to analyse all datasets yourself using nmadb
# Or, use download the .rds object provided

# Read in all NMAs analysis object and study information
allNMAs <- readRDS("allNMAs.rds") # first, set working directory to location where you saved the object
study_info <- readRDS("studyList.rds") %>% select(-contains("..choice"))

# find only the ones without errors
ran <- unlist(lapply(allNMAs, function(x) length(x) > 1))
ran_NMAs <- allNMAs[which(ran)] # 266 ran

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

# Calculate poths for everything
pscores <- lapply(ran_NMAs, get_ranks, sel = studies_ran)
poths <- sapply(pscores, function(x) poth(x)$poth)
summary(poths)


# pdf("~/PhD/Freiburg/Project/SIR/results/hist_sir.pdf", width = 6, height = 4.5)
hist(poths, main = "",
     xlab= "POTH", col = "grey95", breaks = 32, xlim = c(0,1))
abline(v = quantile(poths, c(0.5)), col = "black", lwd = 2, lty = 2)
legend("topleft", text.col = c("black", "black", "black"), lty = c(2, NA, NA), lwd =c(2, NA, NA),
       legend = c("Median = 0.672", "Min = 0.096", "Max > 0.999"), bty = "n")
# dev.off()

# number of studoes
nstudies <- sapply(ran_NMAs,ns_func)


# look at relationship with number of treatments
ntrts <- sapply(pscores, length)
cor(poths, ntrts)
plot(ntrts, poths, ylim = c(0,1), xlim = c(4, 46),
     xlab = "Number of Treatments", ylab = "POTH")

# Relationship with proportion of significantly different treatment comparisons
alphas <- rev(c(0.1, 0.05, 0.01, 0.005))
proportions <- matrix(nrow = length(alphas), ncol = length(poths))

for (i in 1:length(alphas))
  proportions[i, ] <- sapply(ran_NMAs, prop_func, alpha = alphas[i])

range(poths[which(proportions[3,] == 0)])

plot(proportions[3,], poths,
     xlab = "Proportion of treatment comparisons that are stat. signif. at 5% level",
     ylab = 'POTH', ylim = c(0,1))

# plots
# pdf("~/PhD/Freiburg/Project/SIR/results/ntrt_alpha_sir.pdf", width = 8*1.5, height = 4*1.5)
par(mfrow = c(1,2))
plot(ntrts, poths, ylim = c(0,1), xlim = c(4, 46),
     xlab = "Number of Treatments", ylab = "POTH")
plot(proportions[3,], poths,
     xlab = "Proportion of treatment comparisons that are stat. signif. at 5% level",
     ylab = 'POTH', ylim = c(0,1))
# dev.off()

# Relationship with range of p-scores
prange <- sapply(pscores, range_func)
plot(prange, poths, xlim = c(0,1), ylim = c(0,1),
     xlab = "max(P-scores) - min(p-scores)", ylab = "POTH")

# Relationship with heterogeneity.inconsistency

q_pvals <- sapply(ran_NMAs, q_func)
row.names(q_pvals) <- c("overall", "heterogeneity", "inconsistency")


par(mfrow = c(1,3))

for(i in 1:3) {

  plot(q_pvals[i,], poths, ylim = c(0,1),
       xlab = paste0("P-value for Q (", row.names(q_pvals)[i], ")"), ylab = "POTH")
  abline(v = 0.05, lty = 3, col = "grey")
  print(cor(poths, q_pvals[i,]))

}

par(mfrow = c(1,1))

# Look at tau for different effect measures

taus <- sapply(ran_NMAs, tau_func)

effs <- sapply(ran_NMAs, get_sm)
table(effs)

thresh <- 0 # can use to remove numerically zero taus if desired

ors <- which(effs == "OR" & taus > thresh)
mds <- which(effs == "MD" & taus > thresh)
rrs <- which(effs == "RR" & taus > thresh)

mdsno <- which(effs == "MD" & taus > thresh & taus < 30)

summary(lm(poths[rrs] ~ taus[rrs]))
summary(lm(poths[ors] ~ taus[ors]))
summary(lm(poths[mds] ~ taus[mds]))

cor(taus[ors], poths[ors])
cor(taus[rrs], poths[rrs])
cor(taus[mds], poths[mds])


# pdf("~/PhD/Freiburg/Project/SIR/results/tau_sir.pdf", width = 8*1.2, height = 2.5*1.2)
par(mfrow = c(1,3))
plot(taus[ors], poths[ors], main = "Effect Measure: Log Odds Ratio\n119 NMAs",
     xlab = expression(tau), ylab = "POTH",
     ylim = c(0,1))
plot(taus[rrs], poths[rrs], main = "Effect Measure: Log Relative Risk\n69 NMAs",
     xlab = expression(tau), ylab = "POTH",
     ylim = c(0,1))
plot(taus[mdsno], poths[mdsno], main = "Effect Measure: Mean Difference\n44 NMAs",
     xlab = expression(tau), ylab = "POTH",
     ylim = c(0,1))

# dev.off()

# Comparisons to examples --------------------------------------------
mean(poths>0.326) # antifungals for transplants example

mean(poths > 0.838) # ICIs for cancer

mean(poths > 0.559) # depression

# max and min poth analyses ------------------------------------------
maxid <- which(poths > 0.999)
minid <- which(poths < 0.1)


max_nma <- ran_NMAs[[maxid]]$netobj
min_nma <- ran_NMAs[[minid]]$netobj
# png("largest_poth.png", width = 6, height = 2.5, units = "in", res = 320)
forest(max_nma)
# dev.off()
# png("smallest_poth.png", width = 6, height = 2.5, units = "in", res = 320)
forest(min_nma)
# dev.off()

minid2 <- ran_NMAs[[minid]]$recid
subset(study_info, Record.ID == minid2)

maxid2 <- ran_NMAs[[maxid]]$recid
subset(study_info, Record.ID == maxid2)
