# Libraries
library("tidyr")
library("dplyr")
library("netmeta")
library("poth")
library("ggplot2")
library("ggpubr")
# Safety of Immune checkpoint inhibitors in cancer -------------------------------

xu <- read.csv("page_3_table.csv")
trtnames <- c("conventional therapy",
              "nivolumab",
              "pembrolizumab",
              "two ICIs",
              "ICI and conventional",
              "atezolizumab",
              "ipilimumab")

xu <- mutate(xu,
             TreatmentID = trtnames[xu$TreatmentID])
pw <- pairwise(treat = TreatmentID, event = Eventcount, n = Samplesize, studlab = StudyID, data = xu, sm = "OR")
netw <- netmeta(pw, small.values = "desirable", method.tau = "REML", common = F)
poth.ici <- poth(netw)
poth.ici
indiv.ici <- loo(poth.ici)
indiv.ici
cum.ici <- cumul(poth.ici)
cum.ici

l <-plot(indiv.ici)
# ggsave("loo-ici.png", width = 4, height = 4, units = "in")

subset.ici <- subset(poth.ici, subset = c("nivolumab",
                                          "pembrolizumab",
                                          "atezolizumab",
                                          "ipilimumab"))

subset.ici

p <- plot(cum.ici, labels = F)

ggarrange(p,l, nrow = 1, widths = c(1.2,2))
# ggsave("diag-ici.png", width = 9, height = 3, units = "in")

adddata <- indiv.ici$resid
names(adddata) <- (indiv.ici$trt)
ad <- data.frame(residuals = adddata[sort(trtnames)])
# png("forest-ici.png",width = 8, height = 2.8, units = "in", res = 320)
forest(netw, small.values = "desirable", reference.group = "conventional therapy",
       sortvar = "-Pscore", rightcols = c("Pscore","residuals", "ci"), add.data = ad,
       rightlabs = c(NA, "rPOTH", NA),
       digits = 3,
       smlab= "Comparison: other vs conventional therapy\nPOTH = 0.838")
# dev.off()


# Depression example - many treatments ---------------------------------------
# Load the object of all NMA analyses
allNMAs <- readRDS("allNMAs.rds")
study_info <- readRDS("studyList.rds") %>%
  select(-contains("..choice"))

ran <- unlist(lapply(allNMAs, function(x) length(x) > 1))
ranNMAs <- allNMAs[which(ran)]

# Depression example
id_depression <- 48

recid_depression <- ranNMAs[[id_depression]]$recid
ex_depression <- ranNMAs[[id_depression]]$netobj

# Find out info about study, like the study name - this is how treatment names are identified
depression_info <- subset(study_info, Record.ID == recid_depression) %>%
  select(-contains("..choice"))

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

net.depr <- netmeta(TE, seTE, trt1_new, trt2_new, studlab,
                    method.tau = "REML",
                    data = newdata_depression,
                    sm = ex_depression$sm,
                    details.chkmultiarm = TRUE,
                    common = FALSE)

# Calculate poth and do LOO poth
poth.depr <- poth(net.depr, small.values = "undesirable")
poth.depr


indiv.depr <- loo(poth.depr)
cum.depr <- cumul(poth.depr)


p1 <- plot(indiv.depr)
p2 <- plot(cum.depr)
ggarrange(p1,p2, heights = c(2, 1.1), nrow = 2)
# ggsave("diag-depr.png", width = 9, height = 5, units = "in")

# plot(cumul(poth.depr), labels = F)
# ggsave("cum-depr.png",width = 8, height = 5, units = "in")
# forest plot
res <- indiv.depr$resid
names(res) <- (indiv.depr$trt)
dat <- data.frame(residuals = res[sort(trt_names)])

# png("forest-depr.png", width = 8, height =7, res = 320, units = "in")
forest(net.depr, small.values = "undesirable", reference.group = "placebo",
       sortvar = "-Pscore", rightcols = c("Pscore","residuals", "ci"),
       rightlabs = c(NA,"rPOTH",NA),
       add.data = dat,
       digits = 3,
       smlab= "Comparison: other vs placebo\nPOTH = 0.559")

# dev.off()

# Transplants - low poth ---------------------------------------------------
num <- 227

low_recid <- ranNMAs[[num]]$recid
low_nma <- ranNMAs[[num]]$netobj

low_info <- subset(study_info, Record.ID == low_recid) %>%
  select(Record.ID, First.Author, Title, Description.of.the.outcome.s.,
         Harmful.Beneficial.Outcome,
         Type.of.Outcome., Effect.Measure, Verified, Error)
low_nma$data %>% select(study, treat1, treat2, TE, event1, n1, event2, n2)

trt_names <-
  c("Control", "Fluconazole", "Itraconazole", "Liposomal amphotericin B",
    "Ketoconazole")

lowdata2 <- low_nma$data %>%
  mutate(trt1_new = trt_names[as.numeric(treat1)],
         trt2_new = trt_names[as.numeric(treat2)])

net.transp <- netmeta(TE, seTE, trt1_new, trt2_new, studlab,
                      data = lowdata2,
                      sm = low_nma$sm,
                      method.tau = "REML",
                      common = FALSE)

poth.transp <- poth(net.transp)
poth.transp

indiv.transp <- loo(poth.transp)

cum.transp <- cumul(poth.transp)

a <- plot(indiv.transp)

b <- plot(cum.transp)

ggarrange(b, a, widths = c(1.2, 2), nrow = 1)
# ggsave("diag_transp.png", width = 9, height = 3, units = "in")
# forest plot
res <- indiv.transp$resid
names(res) <- (indiv.transp$trt)
dat <- data.frame(residuals = res[sort(trt_names)])

# png("forest-transp.png", width = 8, height =2.5, res = 320, units = "in")
forest(net.transp, small.values = "desirable", reference.group = "Control",
       sortvar = "-Pscore", rightcols = c("Pscore","residuals", "ci"),
       rightlabs = c(NA,"rPOTH",NA),
       add.data = dat,
       digits = 3,
       smlab= "Comparison: other vs control\nPOTH = 0.326")

# dev.off()

# cumulative POTH
cum.transp <- cumul(poth.transp)
plot(cum.transp)
# ggsave("cum_transp.png", width = 8, height = 4, units = "in")

#subset POTH
sub.transp <- subset(poth.transp, c("Fluconazole",
                                   "Itraconazole",
                                   "Liposomal amphotericin B",
                                   "Ketoconazole"))

sub.transp # poth = 0.354
