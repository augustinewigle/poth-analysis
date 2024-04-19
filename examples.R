# Libraries
library(netmeta)
library(BUGSnet)
library(sir)

# Load the object of all NMA analyses
allNMAs <- readRDS("allNMAs.rds") # first, set working directory to location where you saved the object
study_info <- readRDS("studyList.rds") %>% select(-contains("..choice"))

# Diabetes prevention example ------------------------------------------------
data("diabetes")

dat <- data.prep(arm.data = diabetes,
                 varname.t = "Treatment",
                 varname.s = "Study")

random_effects_model <- nma.model(data=dat,
                                  outcome="diabetes",
                                  N="n",
                                  reference="Diuretic",
                                  family="binomial",
                                  link="cloglog",
                                  time = "followup",
                                  effects= "random")

set.seed(2024)
random_effects_results <- nma.run(random_effects_model,
                                  n.adapt=1000,
                                  n.burnin=5000,
                                  n.iter=10000)

cat(random_effects_results$model$bugs)

# Matrix of relative effect samples
all_samps <- as.matrix(random_effects_results$samples)
effect_samps <- all_samps[,startsWith(colnames(all_samps), prefix = "d[")]

# SIR and LOO SIR
loo1 <- loo_sir(samples= effect_samps,
                trt_names = c("Diuretic",
                              "ACE inhibitor",
                              "ARB",
                              "blocker",
                              "CCB",
                              "Placebo"),
                largerbetter = F)

loo1$plot

loo1$info

nma.forest(random_effects_results,
           comparator = "Placebo",
           order = rev(loo1$info$trt_name))
nma.league(random_effects_results)

# local SIR - the top 3 treatments
loc1 <- local_sir(samples= effect_samps,
                  trt_names = c("Diuretic",
                                "ACE inhibitor",
                                "ARB",
                                "blocker",
                                "CCB",
                                "Placebo"),
                  largerbetter = F,
                  topx = 3,
                  do_loo = T)

loc1$loo$plot

# Depression example - many treatments ---------------------------------------
id_depression <- 49
recid_depression <- ran_NMAs[[id_depression]]$recid
ex_depression <- ran_NMAs[[id_depression]]$netobj

# Find out info about study, like the study name - this is how treatment names are identified
depression_info <- subset(study_info, Record.ID == recid_depression) %>% select(-contains("..choice"))

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

newdata_depression <- ex_depression$data %>% mutate(trt1_new = trt_names[as.numeric(treat1)],
                                trt2_new = trt_names[as.numeric(treat2)])

depression <- netmeta(TE = TE,
                        seTE = seTE,
                        treat1 = trt1_new,
                        treat2 = trt2_new,
                        studlab = studlab,
                        sm = ex_depression$sm,
                        details.chkmultiarm = T,
                        fixed = F,
                        random = T,
                        data = newdata_depression)

# Calculate SIR and do LOO SIR
loo2 <- loo_sir(diffs = depression$TE.random,
               ses = depression$seTE.random,
               largerbetter = T)
loo2$plot

forest(depression, small.values = "bad", reference.group = "placebo",
       sortvar = "-Pscore", rightcols = c("Pscore","effect", "ci"),
       smlab= "Comparison: other vs placebo\nSIR = 0.538")


# Local SIR - Look at top 5 treatments
top5 <- local_sir(diffs = depression$TE.random,
                  ses = depression$seTE.random,
                  largerbetter = T,
                  topx = 5,
                  do_loo = T)



top5$sir
top5$loo$plot

# Local Sir - look at top and bottom 5
subs <- local_sir(diffs = depression$TE.random,
                  ses = depression$seTE.random,
                  largerbetter = T,
                  subset = loo2$info$trt_name[c(1:5, 25:29)],
                  do_loo = T)

subs$sir

# Transplants - low SIR ---------------------------------------------------
num <- 229
low_recid <- ran_NMAs[[num]]$recid
low_nma <- ran_NMAs[[num]]$netobj

low_info <- subset(study_info, Record.ID == low_recid) %>%
  select(Record.ID, First.Author, Title, Description.of.the.outcome.s., Harmful.Beneficial.Outcome,
         Type.of.Outcome., Effect.Measure, Verified, Error)
low_nma$data %>% select(study, treat1, treat2, TE, event1, n1, event2, n2)

trt_names <- c("Control",
               "Fluconazole",
               "Itraconazole",
               "Liposomal amphotericin B",
               "Ketoconazole")

lowdata2 <- low_nma$data %>% mutate(trt1_new = trt_names[as.numeric(treat1)],
                                    trt2_new = trt_names[as.numeric(treat2)])

transplant <- netmeta(TE = TE,
                       seTE = seTE,
                       treat1 = trt1_new,
                       treat2 = trt2_new,
                       studlab = studlab,
                       sm = low_nma$sm,
                       details.chkmultiarm = T,
                       fixed = F,
                       random = T,
                       data = lowdata2)

forest(transplant, small.values = "good", reference.group = "Control",
       sortvar = "-Pscore", rightcols = c("Pscore","effect", "ci"),
       smlab= "Comparison: other vs control\nSIR = 0.241")

loo3 <- loo_sir(diffs = transplant$TE.random,
                se = transplant$seTE.random,
                largerbetter = F)
loo3$plot

# local SIR - compare all non-control treatments
loc3 <- local_sir(diffs = transplant$TE.random,
                 se = transplant$seTE.random,
                 largerbetter = F,
                 subset = c("Fluconazole",
                            "Itraconazole",
                            "Liposomal amphotericin B",
                            "Ketoconazole"),
                 do_loo = T)
loc3$loo$plot
