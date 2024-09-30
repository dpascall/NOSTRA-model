rm(list = ls())
source("~/Documents/NOSTRA-model/NOSTRA.R")
load("~/Documents/NOSTRA-model/finaldataformodel.Rdata")

##function generate time difference matrix from sample collection dates
timedifmat <- function(x) {
  dateformat <- as.Date(x, format = "%Y-%m-%d")
  return(do.call("cbind", lapply(seq_along(dateformat), function(x) {return(abs(as.numeric(difftime(dateformat[x],dateformat, units = "days"))))})))
}

results <- vector("list", length = length(finaldataformodel))
for (i in seq_len(length(finaldataformodel))) {
  print(paste0("Dataset ", i, " of ", length(finaldataformodel)))
  ##made data format appropriate
  finaldataformodel[[i]][[1]]$detection <- as.numeric(difftime(as.Date(finaldataformodel[[i]][[1]]$detectiondate, format = "%Y-%m-%d"),
                                           as.Date(finaldataformodel[[i]][[5]], format = "%Y-%m-%d")), units = "days") + 1
  finaldataformodel[[i]][[1]]$admission <- as.numeric(difftime(as.Date(finaldataformodel[[i]][[1]]$admission, format = "%Y-%m-%d"),
                                           as.Date(finaldataformodel[[i]][[5]], format = "%Y-%m-%d")), units = "days") + 1
  finaldataformodel[[i]][[1]]$sequencedate <- finaldataformodel[[i]][[1]]$detection
  finaldataformodel[[i]][[1]]$patient_study_id <- finaldataformodel[[i]][[1]]$sequence_id <- as.character(finaldataformodel[[i]][[1]]$ID)
  finaldataformodel[[i]][[1]]$onset_date <- format(as.Date(finaldataformodel[[i]][[1]]$detectiondate, format = "%Y-%m-%d"), "%d/%m/%Y")
  ##calculate proportion of nosocomial infections post 50 days
  prevalence <- sum(finaldataformodel[[i]][[1]]$infectionsource == "C" & finaldataformodel[[i]][[1]]$admission > 50)/
                length(finaldataformodel[[i]][[1]]$infectionsource[finaldataformodel[[i]][[1]]$admission > 50])
  ##iterate over wards
  wardresults <- vector("list", length = length(unique(finaldataformodel[[i]][[1]]$ward)))
  for (j in seq_len(length(unique(finaldataformodel[[i]][[1]]$ward)))) {
    print(paste0("Ward ", j, " of ", length(unique(finaldataformodel[[i]][[1]]$ward))))
    epidata <- finaldataformodel[[i]][[1]][finaldataformodel[[i]][[1]]$ward %in% unique(finaldataformodel[[i]][[1]]$ward)[j],]
    SNPs <- finaldataformodel[[i]][[2]][rownames(finaldataformodel[[i]][[2]]) %in% epidata$ID, colnames(finaldataformodel[[i]][[2]]) %in% epidata$ID]
    alignmentlength <- SNPs
    alignmentlength[,] <- 29811
    timedifferences <- timedifmat(epidata$detectiondate)
    rownames(timedifferences) <- colnames(timedifferences) <- as.character(epidata$ID)
    locations <- cbind(Date = finaldataformodel[[i]][[4]]$Date, finaldataformodel[[i]][[4]][,colnames(finaldataformodel[[i]][[4]]) %in% epidata$ID])
    for (l in seq_len(ncol(locations))) {
      locations[,l] <- as.character(locations[,l])
    }
    IDsgreaterthat50 <- epidata$ID[epidata$admission > 50]
    patientresults <- vector("list", length = length(IDsgreaterthat50))
    for (k in seq_len(length(IDsgreaterthat50))) {
      print(paste0("Individual ", k, " of ", length(IDsgreaterthat50)))
      modifiedepidata <- rbind(epidata[as.character(epidata$ID) %in% IDsgreaterthat50[k],], epidata[-which(epidata$ID == IDsgreaterthat50[k]),])
      patientresults[[k]] <- list(run_model(epidata = modifiedepidata, locationdata = locations, snpdistmat = SNPs,
                             timedistmat = timedifferences, alignmentlengthmat = alignmentlength,
                             priors = c(rep(0.5*(1/nrow(epidata)), nrow(epidata)), 0.5), 
                             all_comparisons = F, startdate = finaldataformodel[[i]][[5]]), modifiedepidata$infectionsource[1])
      names(patientresults[[k]][[1]]) <- c(as.character(modifiedepidata$ID[2:length(modifiedepidata$ID)]), "H", "C")
    }
    wardresults[[j]] <- list(patientresults, c(rep(0.5*(1/nrow(epidata)), nrow(epidata)), 0.5),
                             c(rep((1-prevalence)*(1/nrow(epidata)), nrow(epidata)), prevalence), epidata)
  }
  results[[i]] <- wardresults
}

##filter out admissions before 51st day
save(results, file = "~/Documents/NOSTRA-model/results.Rdata")