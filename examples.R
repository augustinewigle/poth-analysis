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
study_info <- studyList %>% select(-contains("..choice"))

ran <- sapply(allNMAs, function(x) length(x) > 1)
ranNMAs <- allNMAs[ran]
recIDs <- sapply(ranNMAs, function(x) x$recid)


#
#
# Safety of Immune checkpoint inhibitors in cancer
#
#

xu <- read.csv("page_3_table.csv")
trtnames <- c("conventional therapy",
              "nivolumab",
              "pembrolizumab",
              "two ICIs",
              "ICI and conventional",
              "atezolizumab",
              "ipilimumab")

xu <- mutate(xu, TreatmentID = trtnames[xu$TreatmentID])
pw <- pairwise(treat = TreatmentID, event = Eventcount, n = Samplesize,
  studlab = StudyID, data = xu, sm = "OR")
net.ici <- netmeta(pw, small.values = "desirable", method.tau = "REML",
  common = FALSE)
poth.ici <- poth(net.ici)
poth.ici
loo.ici <- loo(poth.ici)
loo.ici
cum.ici <- cumul(poth.ici)
cum.ici

l <- plot(loo.ici)
# l
# ggsave("loo-ici.png", width = 4, height = 4, units = "in")

subset.ici <- subset(poth.ici,
                     subset = c("nivolumab", "pembrolizumab",
                                "atezolizumab", "ipilimumab"))

subset.ici

p <- plot(cum.ici, labels = FALSE)

ggarrange(p, l, nrow = 1, widths = c(1.2, 2))
# ggsave("diag-ici.png", width = 9, height = 3, units = "in")

ici.add <- as.data.frame(loo.ici) %>% arrange(trt)
all(sort(net.ici$trts) == sort(ici.add$trt))
ici.add <- ici.add[net.ici$trts, ]

# png("forest-ici.png",width = 8, height = 2.8, units = "in", res = 320)
forest(net.ici, small.values = "desirable",
       reference.group = "conventional therapy",
       sortvar = "-Pscore",
       add.data = ici.add,
       leftcols = c("studlab", "Pscore", "resid"),
       leftlabs = c(NA, NA, "rPOTH"),
       #rightcols = c("Pscore", "resid", "effect", "ci"),
       #rightlabs = c(NA, "rPOTH", NA, NA),
       digits = 2, digits.prop = 3, digits.addcols = 3,
       just.addcols = "right",
       smlab =
         paste("Comparison: other vs conventional therapy\nPOTH = ",
               round(poth.ici$poth, 3)))

# dev.off()


#
#
# Depression example - many treatments
#
#

recid_depression <- 480655
ex_depression <- ranNMAs[[which(recIDs == recid_depression)]]$netobj

# Find out info about study, like the study name - this is how treatment names are identified
depression_info <- study_info %>% filter(Record.ID == recid_depression)

trt_names <-
  c("placebo", "fluoxetine", "sertraline", "paroxetine", "reboxetine",
    "viloxazine", "moclobemide", "phenelzine", "clomipramine",
    "imipramine", "amineptine", "amitriptyline", "desipramine",
    "trazodone", "minaprine", "st. john's wort", "ritanserin",
    "amisulpride", "flupenthixol", "lorazepam", "diazepam",
    "acetyl-l-carnitine", "fluoxetine + bentazepam",
    "paroxetine + amisulpride", "escitalopram", "dothiepine",
    "duloxetine", "buprion + escitalopram", "venlafaxine + mirtazapine")

newdata_depression <-
  ex_depression$data %>%
  mutate(trt1_new = trt_names[as.numeric(treat1)],
         trt2_new = trt_names[as.numeric(treat2)])

net.depr <-
  netmeta(TE, seTE, trt1_new, trt2_new, studlab, data = newdata_depression,
          sm = ex_depression$sm, common = FALSE, method.tau = "REML")

# Calculate POTH and do LOO POTH

poth.depr <- poth(net.depr, small.values = "undesirable")
poth.depr

loo.depr <- loo(poth.depr)
cum.depr <- cumul(poth.depr)


p1 <- plot(loo.depr)
p2 <- plot(cum.depr)
ggarrange(p1, p2, heights = c(2, 1.1), nrow = 2)
# ggsave("diag-depr.png", width = 9, height = 5, units = "in")

# plot(cumul(poth.depr), labels = FALSE)
# ggsave("cum-depr.png",width = 8, height = 5, units = "in")
# forest plot

loo.add <- as.data.frame(loo.depr) %>% arrange(trt)
all(net.depr$trts == loo.add$trt)

# png("forest-depr.png", width = 8, height =7, res = 320, units = "in")
forest(net.depr, small.values = "undesirable", reference.group = "placebo",
       sortvar = "-Pscore",
       add.data = loo.add,
       leftcols = c("studlab", "Pscore", "resid"),
       leftlabs = c(NA, NA, "rPOTH"),
       #rightcols = c("Pscore", "resid", "effect", "ci"),
       #rightlabs = c(NA, "rPOTH", NA, NA),
       digits = 2, digits.prop = 3, digits.addcols = 3,
       just.addcols = "right",
       smlab = paste("Comparison: other vs placebo\nPOTH =",
                    round(poth.depr$poth, 3)))

# dev.off()


#
#
# Transplants - low POTH
#
#

recid.transp <- 501370
nma_transp <- ranNMAs[[which(recIDs == recid.transp)]]$netobj

info_transp <- subset(study_info, Record.ID == recid.transp) %>%
  select(Record.ID, First.Author, Title, Description.of.the.outcome.s.,
         Harmful.Beneficial.Outcome,
         Type.of.Outcome., Effect.Measure, Verified, Error)
info_transp
nma_transp$data %>% select(study, treat1, treat2, TE, event1, n1, event2, n2)

trt_names <-
  c("Control", "Fluconazole", "Itraconazole", "Liposomal amphotericin B",
    "Ketoconazole")

data_transp <- nma_transp$data %>%
  mutate(trt1_new = trt_names[as.numeric(treat1)],
         trt2_new = trt_names[as.numeric(treat2)])

net.transp <-
  netmeta(TE, seTE, trt1_new, trt2_new, studlab, data = data_transp,
          sm = nma_transp$sm, common = FALSE, method.tau = "REML")

poth.transp <- poth(net.transp)
poth.transp

loo.transp <- loo(poth.transp)

cum.transp <- cumul(poth.transp)

p1 <- plot(loo.transp)
p2 <- plot(cum.transp)

ggarrange(p2, p1, widths = c(1.2, 2), nrow = 1)
# ggsave("diag_transp.png", width = 9, height = 3, units = "in")
# forest plot

transp.add <- as.data.frame(loo.transp) %>% arrange(trt)
all(net.transp$trts == transp.add$trt)

# png("forest-transp.png", width = 8, height =2.5, res = 320, units = "in")
forest(net.transp, small.values = "desirable", reference.group = "Control",
       sortvar = "-Pscore",
       add.data = transp.add,
       leftcols = c("studlab", "Pscore", "resid"),
       leftlabs = c(NA, NA, "rPOTH"),
       #rightcols = c("Pscore", "resid", "effect", "ci"),
       #rightlabs = c(NA, "rPOTH", NA, NA),
       digits = 2, digits.prop = 3, digits.addcols = 3,
       just.addcols = "right",
       smlab = paste("Comparison: other vs control\nPOTH =",
                     round(poth.transp$poth, 3)))

# dev.off()

# cumulative POTH
cum.transp <- cumul(poth.transp)
plot(cum.transp)
# ggsave("cum_transp.png", width = 8, height = 4, units = "in")

#subset POTH
sub.transp <- subset(poth.transp,
                     c("Fluconazole", "Itraconazole",
                       "Liposomal amphotericin B", "Ketoconazole"))

sub.transp # poth = 0.354
