# Libraries ---------------------------------
remotes::install_github("augustinewigle/nmadb")
library(nmadb)
library(netmeta)

# Analyse the database
studyList <- getNMADB()

allStudyIDs <- subset(studyList, dataset != "")$Record.ID


allNMAs <-
  lapply( allStudyIDs
        , function(id){
            tryCatch({
              net <- runnetmeta(recid=id, method.tau = "REML")
              return(list(recid = id
                         ,netobj = net))
              }, error = function(cond) {
                res <- paste(id,conditionMessage(cond))
                res
              }
            )
        })

# Remove unnecessary returned information due to method.tau = REML
ix <- c(156, 157, 158, 159, 160)

listfun <- function(x) {

  if(length(x) == 1) {

    return(x)

  }
  recid <- x$recid
  netobj <- x$netobj[-c(ix)]

  class(netobj) <- "netmeta"
  return(list(recid = recid, netobj = netobj))
}


allNMAsnew <- lapply(allNMAs, listfun)


saveRDS(allNMAsnew, file = "allNMAs.rds")
saveRDS(studyList, file = "studyList.rds")
