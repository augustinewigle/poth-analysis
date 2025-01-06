#
# Libraries
#

library("tidyr")
library("dplyr")
library("netmeta")
library("poth")
library("ggplot2")
library("ggpubr")

#
# Load nmadb NMAs
#

load("nmadb_NMAs.rda")

study_info <- studyList %>%
  rename(recordid = Record.ID, author = First.Author,
         title = Title, outcome = Description.of.the.outcome.s.,
         small.values = Harmful.Beneficial.Outcome,
         type = Type.of.Outcome., sm = Effect.Measure,
         verified = Verified, error = Error) %>%
  mutate(small.values =
           recode(small.values,
                  "Beneficial" = "undesirable",
                  "Harmful" = "desirable",
                  "Unclear" = "unclear"),
         sm =
           recode(sm,
                  "hazard ratio" = "HR",
                  "mean difference" = "MD",
                  "odds ratio" = "OR",
                  "other" = "other",
                  "rate ratio" = "IRR",
                  "risk difference" = "RD",
                  "risk ratio" = "RR",
                  "standardized mean difference" = "SMD")) %>%
  select(recordid, author, title, outcome,
         small.values, type, sm, verified, error)


ran <- sapply(allNMAs, function(x) length(x) > 1)
ranNMAs <- allNMAs[ran]
recIDs <- sapply(ranNMAs, function(x) x$recid)


#
#
# Example 1:
# Antifungal Treatments to Prevent Mortality for Solid Organ Transplant
# Recipients
#
#

recid1 <- 501370
#
trts1 <-
  c("Control", "Fluconazole", "Itraconazole", "Liposomal amphotericin B",
    "Ketoconazole")

nma1 <- ranNMAs[[which(recIDs == recid1)]]$netobj
#
info1 <- study_info %>% filter(recordid == recid1)
info1
#
dat1 <- nma1$data %>%
  mutate(trt1 = factor(treat1, levels = seq_along(trts1), labels = trts1),
         trt2 = factor(treat2, levels = seq_along(trts1), labels = trts1))
#
dat1 %>% select(studlab, trt1, trt2, TE, event1, n1, event2, n2) %>% head()

# It is important to define argument 'small.values' which is used to calculate
# P-scores and POTHs
net1 <- netmeta(TE, seTE, trt1, trt2, studlab, data = dat1,
  sm = info1$sm, small.values = info1$small.values,
  common = FALSE, method.tau = "REML")
#
poth1 <- poth(net1)
loo1 <- loo(poth1)
bestk1 <- bestk(poth1)

# Forest plot
#
add1 <- as.data.frame(loo1) %>% arrange(trt)
all(net1$trts == add1$trt)
#
png("forest-transp.png", width = 8, height =2.5, res = 320, units = "in")
forest(net1, reference.group = "Control",
       sortvar = "-Pscore",
       add.data = add1,
       leftcols = c("studlab", "Pscore", "resid"),
       leftlabs = c(NA, NA, "rPOTH"),
       #rightcols = c("Pscore", "resid", "effect", "ci"),
       #rightlabs = c(NA, "rPOTH", NA, NA),
       digits = 2, digits.prop = 3, digits.addcols = 3,
       just.addcols = "right",
       addrows = 1,
       text.addline1 = paste("POTH =", round(poth1$poth, 3)),
       ff.addline = "bold")
#       smlab = paste("Comparison: other vs control\nPOTH =",
#                     round(poth1$poth, 3)))
dev.off()

# Bar plots
#
plt1.1 <- plot(loo1)
plt1.2 <- plot(bestk1)
#plt1.2 <- plot(bestk1, labels = TRUE)
#
ggarrange(plt1.2, plt1.1, widths = c(1.2, 2), nrow = 1)
# ggsave("diag_transp.png", width = 9, height = 3, units = "in")

# Subset POTH
#
sub1 <- subset(poth1,
  c("Fluconazole", "Itraconazole",
    "Liposomal amphotericin B", "Ketoconazole"))
sub1 # poth = 0.354


#
#
# Example 2:
# Efficacy of Pharmacological Treatments for Persistent Depressive Disorder
#
#

recid2 <- 480666
#
trts2 <-
  c("Placebo", "Fluoxetine", "Sertraline", "Paroxetine", "Reboxetine",
    "Viloxazine", "Moclobemide", "Phenelzine", "Clomipramine",
    "Imipramine", "Amineptine", "Amitriptyline", "Desipramine",
    "Trazodone", "Minaprine", "St. John's wort", "Ritanserin",
    "Amisulpride", "Flupenthixol", "Lorazepam", "Diazepam",
    "Acetyl-l-carnitine", "Fluoxetine + Bentazepam",
    "Paroxetine + Amisulpride", "Escitalopram", "Dothiepine",
    "Duloxetine", "Buprion + Escitalopram", "Venlafaxine + Mirtazapine")

nma2 <- ranNMAs[[which(recIDs == recid2)]]$netobj
#
info2 <- study_info %>% filter(recordid == recid2)
info2
#
dat2 <- nma2$data %>%
  mutate(trt1 = factor(treat1, levels = seq_along(trts2), labels = trts2),
         trt2 = factor(treat2, levels = seq_along(trts2), labels = trts2))
#
dat2 %>% select(studlab, trt1, trt2, TE, n1, n2) %>% head()

# It is important to define argument 'small.values' which is used to calculate
# P-scores and POTHs
net2 <- netmeta(TE, seTE, trt1, trt2, studlab, data = dat2,
  sm = info2$sm, small.values = info2$small.values,
  common = FALSE, method.tau = "REML")
#
poth2 <- poth(net2)
loo2 <- loo(poth2)
bestk2 <- bestk(poth2)

# Forest plot
#
add2 <- as.data.frame(loo2) %>% arrange(trt)
all(net2$trts == add2$trt)
#
png("forest-depr.png", width = 8, height =7, res = 320, units = "in")
forest(net2, reference.group = "placebo",
       sortvar = "-Pscore",
       add.data = add2,
       leftcols = c("studlab", "Pscore", "resid"),
       leftlabs = c(NA, NA, "rPOTH"),
       #rightcols = c("Pscore", "resid", "effect", "ci"),
       #rightlabs = c(NA, "rPOTH", NA, NA),
       digits = 2, digits.prop = 3, digits.addcols = 3,
       just.addcols = "right",
       addrows = 1,
       text.addline1 = paste("POTH =", round(poth2$poth, 3)),
       ff.addline = "bold")
#       smlab = paste("Comparison: other vs placebo\nPOTH =",
#                    round(poth2$poth, 3)))
dev.off()

# Bar plots
#
plt2.1 <- plot(loo2)
plt2.2 <- plot(bestk2)
ggarrange(plt2.1, plt2.2, heights = c(2, 1.1), nrow = 2)
ggsave("diag-depr.png", width = 9, height = 5, units = "in")


#
#
# Example 3:
# Safety of Immune checkpoint inhibitors in cancer
#
#

data(Xu2018, package = "poth")
#
pw <- pairwise(treat = treatment, event = adverse, n = n,
  studlab = studyID, data = Xu2018, sm = "OR")
#
net3 <- netmeta(pw, small.values = "desirable", method.tau = "REML",
  common = FALSE)
#
poth3 <- poth(net3)
loo3 <- loo(poth3)
bestk3 <- bestk(poth3)

# Forest plot
#
add3 <- as.data.frame(loo3) %>% arrange(trt)
all(sort(net3$trts) == sort(add3$trt))
add3 <- add3[net3$trts, ]
#
png("forest-ici.png",width = 8, height = 2.8, units = "in", res = 320)
forest(net3,
       reference.group = "Conventional therapy",
       sortvar = "-Pscore",
       add.data = add3,
       leftcols = c("studlab", "Pscore", "resid"),
       leftlabs = c(NA, NA, "rPOTH"),
       #rightcols = c("Pscore", "resid", "effect", "ci"),
       #rightlabs = c(NA, "rPOTH", NA, NA),
       digits = 2, digits.prop = 3, digits.addcols = 3,
       just.addcols = "right",
       addrows = 1,
       text.addline1 = paste("POTH =", round(poth3$poth, 3)),
       ff.addline = "bold")
# png("ici_example.png", width = 8, height = 2.8, units = "in", res = 320)
# forest(net3,
#        reference.group = "Conventional therapy",
#        sortvar = "-Pscore",
#        # add.data = add3,
#        leftcols = c("studlab", "Pscore"),
#        # leftlabs = c(NA, NA, "rPOTH"),
#        #rightcols = c("Pscore", "resid", "effect", "ci"),
#        #rightlabs = c(NA, "rPOTH", NA, NA),
#        digits = 2, digits.prop = 3, digits.addcols = 3,
#        # just.addcols = "right",
#        # addrows = 1,
#        # text.addline1 = paste("POTH =", round(poth3$poth, 3)),
#        ff.addline = "bold")
# #       smlab =
#        paste("Comparison: other vs Conventional therapy\nPOTH = ",
#             round(poth3$poth, 3)))
dev.off()

# Bar plots
#
plt3.1 <- plot(loo3)
plt3.2 <- plot(bestk3)
ggarrange(plt3.2, plt3.1, nrow = 1, widths = c(1.2, 2))
# ggsave("diag-ici.png", width = 9, height = 3, units = "in")

# Subset POTH
#
sub3 <- subset(poth3,
  subset = c("Nivolumab", "Pembrolizumab", "Atezolizumab", "Ipilimumab"))
#
sub3
