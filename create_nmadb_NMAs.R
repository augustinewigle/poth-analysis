#
# Libraries
#

remotes::install_github("guido-s/meta")
remotes::install_github("guido-s/netmeta")
remotes::install_github("guido-s/nmadb")
#
library("nmadb")
library("netmeta")

#
# Analyse the database
#

studyList <- getNMADB()
allStudyIDs <- subset(studyList, dataset != "")$Record.ID
#
allNMAs <-
  lapply(allStudyIDs,
         function(id) {
           tryCatch({
             net <- runnetmeta(recid = id, method.tau = "REML")
             return(list(recid = id, netobj = net))
           },
           error = function(cond) {
             res <- paste(id, conditionMessage(cond))
             res
           })
         })
#
warnings()

#
# Save the objects in rda file
#

save(allNMAs, studyList, file = "nmadb_NMAs.rda")
