# Libraries
library("tidyverse")
library("netmeta")
#remotes::install_github("audrey-b/BUGSnet")
library("BUGSnet")
library("sir")


# Diabetes prevention example ------------------------------------------------
data("diabetes")

dat <- data.prep(arm.data = diabetes,
                 varname.t = "Treatment",
                 varname.s = "Study")

random_effects_model <- nma.model(data = dat,
                                  outcome = "diabetes",
                                  N = "n",
                                  reference = "Diuretic",
                                  family = "binomial",
                                  link = "cloglog",
                                  time = "followup",
                                  effects = "random")

set.seed(2024)
random_effects_results <- nma.run(random_effects_model,
                                  n.adapt = 1000,
                                  n.burnin = 5000,
                                  n.iter = 10000)

cat(random_effects_results$model$bugs)

# Matrix of relative effect samples
all_samps <- as.matrix(random_effects_results$samples)
effect_samps <-
  all_samps[, startsWith(colnames(all_samps),prefix = "d[")]

# SIR and LOO SIR

sir1 <- sir(effect_samps,
            trts = c("Diuretic",
                     "ACE inhibitor",
                     "ARB",
                     "blocker",
                     "CCB",
                     "Placebo"),
            small.values = "desirable")

loo1 <- loo(sir1)

plot(loo1)

loo1

nma.forest(random_effects_results,
           comparator = "Placebo",
           order = rev(loo1$info$trt_name))
nma.league(random_effects_results)

# local SIR - the top 3 treatments
loc1 <- loo(subset(sir1, top = 3))

plot(loc1)





# Load the object of all NMA analyses
allNMAs <- readRDS("allNMAs.rds")
study_info <- readRDS("studyList.rds") %>%
  select(-contains("..choice"))
#
ran <- unlist(lapply(allNMAs, function(x) length(x) > 1))
ranNMAs <- allNMAs[which(ran)]


# Depression example - many treatments ---------------------------------------
id_depression <- 49
recid_depression <- ranNMAs
ex_depression <- ranNMAs[[id_depression]]$netobj

# Find out info about study, like the study name - this is how treatment names are identified
depression_info <- subset(study_info, Record.ID == recid_depression) %>%
  select(-contains("..choice"))

trt_names <- c("placebo",
               "fluoxetine",
               "sertraline",
               "paroxetine",
               "reboxetine",
               "viloxazine",
               "moclobemide",
               "phenelzine",
               "clomipramine",
               "imipramine",
               "amineptine",
               "amitriptyline",
               "desipramine",
               "trazodone",
               "minaprine",
               "st. john's wort",
               "ritanserin",
               "amisulpride",
               "flupenthixol",
               "lorazepam",
               "diazepam",
               "acetyl-l-carnitine",
               "fluoxetine + bentazepam",
               "paroxetine + amisulpride",
               "escitalopram",
               "dothiepine",
               "duloxetine",
               "buprion + escitalopram",
               "venlafaxine + mirtazapine")

newdata_depression <-
  ex_depression$data %>%
  mutate(trt1_new = trt_names[as.numeric(treat1)],
         trt2_new = trt_names[as.numeric(treat2)])

net.depr <- netmeta(TE, seTE, trt1_new, trt2_new, studlab,
                    data = newdata_depression,
                    sm = ex_depression$sm,
                    details.chkmultiarm = TRUE,
                    common = FALSE)

# Calculate SIR and do LOO SIR

sir.depr <- sir(net.depr, small.values = "undesirable")
sir.depr

loo.depr <- loo(sir.depr)

plot(loo.depr)

forest(net.depr, small.values = "undesirable", reference.group = "placebo",
       sortvar = "-Pscore", rightcols = c("Pscore","effect", "ci"),
       smlab= "Comparison: other vs placebo\nSIR = 0.538")


# Local SIR - Look at top 5 treatments
top5.depr <- subset(sir.depr, top = 5)
loo.top5.depr <- loo(top5.depr)

top5
plot(loo.top5.depr)

# Local Sir - look at top and bottom 5

top5.bottom5.depr <- subset(sir.depr, top = 5, bottom = 5)
loo.top5.bottom5.depr <- loo(top5.bottom5.depr)

top5.bottom5.
plot(loo.top5.bottom5.)


# Transplants - low SIR ---------------------------------------------------
num <- 229
low_recid <- ranNMAs[[num]]$recid
low_nma <- ranNMAs[[num]]$netobj

low_info <- subset(study_info, Record.ID == low_recid) %>%
  select(Record.ID, First.Author, Title, Description.of.the.outcome.s.,
         Harmful.Beneficial.Outcome,
         Type.of.Outcome., Effect.Measure, Verified, Error)
low_nma$data %>% select(study, treat1, treat2, TE, event1, n1, event2, n2)

trt_names <- c("Control",
               "Fluconazole",
               "Itraconazole",
               "Liposomal amphotericin B",
               "Ketoconazole")

lowdata2 <- low_nma$data %>%
  mutate(trt1_new = trt_names[as.numeric(treat1)],
         trt2_new = trt_names[as.numeric(treat2)])

net.transp <- netmeta(TE, seTE, trt1_new, trt2_new, studlab,
                      data = lowdata2,
                      sm = low_nma$sm,
                      common = FALSE)

forest(net.transp, small.values = "good", reference.group = "Control",
       sortvar = "-Pscore", rightcols = c("Pscore","effect", "ci"),
       smlab= "Comparison: other vs control\nSIR = 0.241")

sir.transp <- sir(net.transp)
sir.transp

loo3 <- loo(sir.transp)
plot(loo3)

# local SIR - compare all non-control treatments
loc3 <- subset(sir.transp,
               subset = c("Fluconazole",
                          "Itraconazole",
                          "Liposomal amphotericin B",
                          "Ketoconazole"))
loc3

plot(loc3)
