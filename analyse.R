# Libraries ---------------------------------
library(nmadb)
library(netmeta)

# Analyse the database
studyList <- getNMADB()
binaryStudiesIDs <- subset(studyList, Effect.Measure=="odds ratio" & dataset != "")$Record.ID
allStudyIDs <- subset(studyList, dataset != "")$Record.ID

#first5Nmas <-
  #lapply( head(binaryStudiesIDs,5)
        #, function(id){
           #net <- runnetmeta(recid=id)
           #netgraph(net)
           #return(net)
  #})
allNMAs <-
  lapply( allStudyIDs
        , function(id){
            tryCatch({
              net <- runnetmeta(recid=id)
              return(list(recid = id
                         ,netobj = net))
              }, error = function(cond) {
                res <- paste(id,conditionMessage(cond))
                res
              }
            )
        })

# Save the objects as RDS
saveRDS(allNMAs, allNMAs.rds)
saveRDS(studyList, studyList.rds)
