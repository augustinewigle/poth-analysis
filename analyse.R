# Libraries ---------------------------------
library("nmadb")
library("netmeta")

# Analyse the database
studyList <- getNMADB()
binaryStudiesIDs <-
  subset(studyList, Effect.Measure == "odds ratio" & dataset != "")$Record.ID
allStudyIDs <- subset(studyList, dataset != "")$Record.ID

allNMAs <-
  lapply(allStudyIDs,
         #head(allStudyIDs, 3),
         function(id) {
           tryCatch({
             net <- runnetmeta(recid=id)
             return(list(recid = id, netobj = net))
           },
           error = function(cond) {
             res <- paste(id,conditionMessage(cond))
             res
           }
           )
         }
  )

# Save the objects as RDS
saveRDS(allNMAs, file = "allNMAs.rds")
saveRDS(studyList, file = "studyList.rds")
