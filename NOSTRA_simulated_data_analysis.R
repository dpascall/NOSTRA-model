rm(list = ls())
source("~/Documents/NOSTRA-model/NOSTRA.R")
load("~/Documents/NOSTRA-model/finaldataformodel.Rdata")
requireNamespace("mlr3measures")

##function generate time difference matrix from sample collection dates
timedifmat <- function(x) {
  dateformat <- as.Date(x, format = "%Y-%m-%d")
  return(do.call("cbind", lapply(seq_along(dateformat), function(x) {return(abs(as.numeric(difftime(dateformat[x],dateformat, units = "days"))))})))
}

##run analyses
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

#save(results, file = "~/Documents/NOSTRA-model/results.Rdata")

##summarise simulation results - total
load("~/Documents/NOSTRA-model/results.Rdata")

simulationresults <- as.data.frame(matrix(NA, nrow = length(finaldataformodel), ncol = 13))
colnames(simulationresults) <- c("Prevalence", "ParamSet", "Repeat", "NosocomialityScore",
                                 "GeneralScore", "PrevPriorNosocomialityScore", "PrevPriorGeneralScore",
                                 "DefaultPriorNosocomialityScore", "DefaultPriorGeneralScore", "RuleofThumbScore",
                                 "NOSTRATreeScore", "PrevTreeScore", "DefaultTreeScore")

##iterate over datasets
for (i in seq_len(length(results))) {
  simulationresults$Prevalence[i] <- finaldataformodel[[i]][[3]][1]
  simulationresults$ParamSet[i] <- finaldataformodel[[i]][[3]][2]
  simulationresults$Repeat[i] <- finaldataformodel[[i]][[3]][3]
  ##iterate over wards
  
  briers <- matrix(NA, nrow = length(results[[i]]), ncol = 10)
  for (j in seq_len(length(results[[i]]))) {
    ##iterate over patients in wards
    patientbriers <- matrix(NA, nrow = length(results[[i]][[j]][[1]]), ncol = 10)
    for (l in seq_len(length(results[[i]][[j]][[1]]))) {
      ##if infected by someone outside their ward set source to H
      ##as they would not be in set of identified patients
      if (!results[[i]][[j]][[1]][[l]][[2]][1] %in% results[[i]][[j]][[4]]$ID & results[[i]][[j]][[1]][[l]][[2]][1] != "C") {
        results[[i]][[j]][[1]][[l]][[2]][1] <- "H"
      }
      ##get transmission trees
      results[[i]][[j]][[1]][[l]][[3]] <- results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]]
      if (sum(results[[i]][[j]][[4]]$treeID %in% results[[i]][[j]][[1]][[l]][[3]]) == 1) {
        results[[i]][[j]][[1]][[l]][[3]] <- "H"
      }
      ##adjust priors for transmission trees and correct posterior
      if (results[[i]][[j]][[1]][[l]][[3]] == "H") {
        targettrees <- unique(results[[i]][[j]][[4]]$treeID[!results[[i]][[j]][[4]]$treeID %in% results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]]])
        treeposterior <- treeprevprior <- treedefaultprior <- rep(0, length(targettrees) + 2)
        names(treeposterior) <- names(treeprevprior) <- names(treedefaultprior) <- c(as.character(targettrees), "H", "C")
        for (k in seq_len(length(targettrees))) {
          treeprevprior[k] <- sum(results[[i]][[j]][[4]]$treeID %in% targettrees[k]) * results[[i]][[j]][[3]][1]
          treedefaultprior[k] <- sum(results[[i]][[j]][[4]]$treeID %in% targettrees[k]) * results[[i]][[j]][[2]][1]
        }
        treeprevprior[length(treeprevprior) - 1] <- results[[i]][[j]][[3]][1]
        treedefaultprior[length(treedefaultprior) - 1] <- results[[i]][[j]][[2]][1]
        treeprevprior[length(treeprevprior)] <- results[[i]][[j]][[3]][length(results[[i]][[j]][[3]])]
        treedefaultprior[length(treedefaultprior)] <- 0.5
        for(k in seq_len(length(results[[i]][[j]][[1]][[l]][[1]]) - 2)) {
          treeposterior[names(treeposterior) %in% results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])[k]]] <-
            treeposterior[names(treeposterior) %in% results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])[k]]] +
            results[[i]][[j]][[1]][[l]][[1]][k]
        }
        treeposterior[length(treeposterior) - 1] <- results[[i]][[j]][[1]][[l]][[1]][length(results[[i]][[j]][[1]][[l]][[1]]) - 1]
        treeposterior[length(treeposterior)] <- results[[i]][[j]][[1]][[l]][[1]][length(results[[i]][[j]][[1]][[l]][[1]])]
      } else {
        targettrees <- unique(results[[i]][[j]][[4]]$treeID)
        treeposterior <- treeprevprior <- treedefaultprior <- rep(0, length(targettrees) + 2)
        names(treeposterior) <- names(treeprevprior) <- names(treedefaultprior) <- c(as.character(targettrees), "H", "C")
        for (k in seq_len(length(targettrees))) {
          treeprevprior[k] <- sum(results[[i]][[j]][[4]]$treeID %in% targettrees[k]) * results[[i]][[j]][[3]][1]
          treedefaultprior[k] <- sum(results[[i]][[j]][[4]]$treeID %in% targettrees[k]) * results[[i]][[j]][[2]][1]
          if (targettrees[k] == results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]]) {
            treeprevprior[k] <- treeprevprior[k] - results[[i]][[j]][[3]][1]
            treedefaultprior[k] <- treedefaultprior[k] - results[[i]][[j]][[2]][1]
          }
        }
        treeprevprior[length(treeprevprior) - 1] <- results[[i]][[j]][[3]][1]
        treedefaultprior[length(treedefaultprior) - 1] <- results[[i]][[j]][[2]][1]
        treeprevprior[length(treeprevprior)] <- results[[i]][[j]][[3]][length(results[[i]][[j]][[3]])]
        treedefaultprior[length(treedefaultprior)] <- 0.5
        for(k in seq_len(length(results[[i]][[j]][[1]][[l]][[1]]) - 2)) {
          treeposterior[names(treeposterior) %in% results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])[k]]] <-
            treeposterior[names(treeposterior) %in% results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])[k]]] +
            results[[i]][[j]][[1]][[l]][[1]][k]
        }
        treeposterior[length(treeposterior) - 1] <- results[[i]][[j]][[1]][[l]][[1]][length(results[[i]][[j]][[1]][[l]][[1]]) - 1]
        treeposterior[length(treeposterior)] <- results[[i]][[j]][[1]][[l]][[1]][length(results[[i]][[j]][[1]][[l]][[1]])]
      }
      ##get brier scores
      prevpriors <- t(as.matrix(results[[i]][[j]][[3]]))
      defaultpriors <- t(as.matrix(results[[i]][[j]][[2]]))
      colnames(prevpriors) <- colnames(defaultpriors) <- names(results[[i]][[j]][[1]][[l]][[1]])
      NOSTRAnosocomial <- matrix(c(1-results[[i]][[j]][[1]][[l]][[1]]["C"], results[[i]][[j]][[1]][[l]][[1]]["C"]), nrow = 1, ncol = 2)
      prevnosocomial <- matrix(c(1-results[[i]][[j]][[3]][length(results[[i]][[j]][[3]])], results[[i]][[j]][[3]][length(results[[i]][[j]][[3]])]), nrow = 1, ncol = 2)
      defaultnosocomial <- matrix(c(0.5, 0.5), nrow = 1, ncol = 2)
      if(results[[i]][[j]][[4]]$detection[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]] -
         results[[i]][[j]][[4]]$admission[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]] < 3) {
        thumbnosocomial <- matrix(c(0, 1), nrow = 1, ncol = 2)
      } else {
        thumbnosocomial <- matrix(c(1, 0), nrow = 1, ncol = 2)
      }
      colnames(NOSTRAnosocomial) <- colnames(prevnosocomial) <- colnames(defaultnosocomial) <- colnames(thumbnosocomial) <- c("H", "C")
      patientbriers[l,1] <- mlr3measures::mbrier(factor(ifelse(results[[i]][[j]][[1]][[l]][[2]] == "C", "C", "H"), levels = c("H","C")),
                                                 NOSTRAnosocomial)
      patientbriers[l,2] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[2]], levels = names(results[[i]][[j]][[1]][[l]][[1]])),
                                                 t(as.matrix(results[[i]][[j]][[1]][[l]][[1]])))
      patientbriers[l,3] <- mlr3measures::mbrier(factor(ifelse(results[[i]][[j]][[1]][[l]][[2]] == "C", "C", "H"), levels = c("H","C")),
                                                 prevnosocomial)
      patientbriers[l,4] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[2]], levels = names(results[[i]][[j]][[1]][[l]][[1]])),
                                                 prevpriors)
      patientbriers[l,5] <- mlr3measures::mbrier(factor(ifelse(results[[i]][[j]][[1]][[l]][[2]] == "C", "C", "H"), levels = c("H","C")),
                                                 defaultnosocomial)
      patientbriers[l,6] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[2]], levels = names(results[[i]][[j]][[1]][[l]][[1]])),
                                                 defaultpriors)
      patientbriers[l,7] <- mlr3measures::mbrier(factor(ifelse(results[[i]][[j]][[1]][[l]][[2]] == "C", "C", "H"), levels = c("H","C")),
                                                 thumbnosocomial)
      patientbriers[l,8] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[3]], levels = names(treeposterior)),
                                                 t(as.matrix(treeposterior)))
      patientbriers[l,9] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[3]], levels = names(treeprevprior)),
                                                 t(as.matrix(treeprevprior)))
      patientbriers[l,10] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[3]], levels = names(treedefaultprior)),
                                                  t(as.matrix(treedefaultprior)))
    }
    briers[j,] <- colMeans(patientbriers)
  }
  simulationresults[i,4:13] <- colMeans(briers)
}

##generate figures
library(ggplot2)
library(tidyverse)

simulationresults$Prevalence <- unlist(simulationresults$Prevalence)
simulationresults$ParamSet <- unlist(simulationresults$ParamSet)
simulationresults$Repeat <- unlist(simulationresults$Repeat)

simulationresults$Prevalence[simulationresults$Prevalence == 1] <- "High"
simulationresults$Prevalence[simulationresults$Prevalence == 2] <- "Low"
simulationresults$Prevalence[simulationresults$Prevalence == 3] <- "Intermediate"


simulationresults <- simulationresults[complete.cases(simulationresults),]

figuredata <- pivot_longer(simulationresults, cols = c("NosocomialityScore",
                                                       "GeneralScore", "PrevPriorNosocomialityScore", "PrevPriorGeneralScore",
                                                       "DefaultPriorNosocomialityScore", "DefaultPriorGeneralScore", "RuleofThumbScore",
                                                       "NOSTRATreeScore", "PrevTreeScore", "DefaultTreeScore"))

figuredata$Type <- rep(NA, nrow(figuredata))
figuredata$Type[figuredata$name %in% c("NosocomialityScore", "PrevPriorNosocomialityScore", "DefaultPriorNosocomialityScore", "RuleofThumbScore")] <- "Nosocomiality Detection"
figuredata$Type[figuredata$name %in% c("GeneralScore", "PrevPriorGeneralScore", "DefaultPriorGeneralScore")] <- "Source Detection"
figuredata$Type[figuredata$name %in% c("NOSTRATreeScore", "PrevTreeScore", "DefaultTreeScore")] <- "Transmission Tree Detection"

##tidy up names for figure
figuredata$name[figuredata$name %in% c("NOSTRATreeScore")] <- "NOSTRA"
figuredata$name[figuredata$name %in% c("NosocomialityScore")] <- "NOSTRA"
figuredata$name[figuredata$name %in% c("GeneralScore")] <- "NOSTRA"
figuredata$name[figuredata$name %in% c("PrevPriorNosocomialityScore")] <- "Prevalence Prior"
figuredata$name[figuredata$name %in% c("PrevTreeScore")] <- "Prevalence Prior"
figuredata$name[figuredata$name %in% c("PrevPriorGeneralScore")] <- "Prevalence Prior"
figuredata$name[figuredata$name %in% c("DefaultPriorNosocomialityScore")] <- "Naïve Prior"
figuredata$name[figuredata$name %in% c("DefaultPriorGeneralScore")] <- "Naïve Prior"
figuredata$name[figuredata$name %in% c("DefaultTreeScore")] <- "Naïve Prior"
figuredata$name[figuredata$name %in% c("RuleofThumbScore")] <- "96hr Categorisation"

figuredata$name <- factor(figuredata$name, levels = c("NOSTRA", "Prevalence Prior", "Naïve Prior", "96hr Categorisation"))

figuredata_means <- figuredata %>% group_by(name, Type) %>% summarise(value=mean(value))

plot <- ggplot(figuredata, aes(name, value)) + 
  facet_grid(. ~ Type, scales = "free_x") +
  ylab("Brier Score") +
  xlab("") +
  theme_bw() +
  geom_jitter(aes(colour = Prevalence), height = 0, width = 0.1) +
  geom_point(data = figuredata_means, size = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#save(results, file = "~/Documents/NOSTRA-model/results.Rdata")

##summarise simulation results - per day
#load("~/Documents/NOSTRA-model/results.Rdata")

simulationresults_perday_list <- vector("list", length = 11)
##iterate over days
for (q in 1:11) {
  simulationresults_perday <- matrix(NA, nrow = length(finaldataformodel), ncol = 13)
  colnames(simulationresults_perday) <- c("Prevalence", "ParamSet", "Repeat", "NosocomialityScore",
                                          "GeneralScore", "PrevPriorNosocomialityScore", "PrevPriorGeneralScore",
                                          "DefaultPriorNosocomialityScore", "DefaultPriorGeneralScore", "RuleofThumbScore",
                                          "NOSTRATreeScore", "PrevTreeScore", "DefaultTreeScore")
  
  ##iterate over datasets
  for (i in 1:31) {#seq_len(length(results))) {
    simulationresults_perday[i,1] <- as.numeric(finaldataformodel[[i]][[3]][1])
    simulationresults_perday[i,2] <- as.numeric(finaldataformodel[[i]][[3]][2])
    simulationresults_perday[i,3] <- as.numeric(finaldataformodel[[i]][[3]][3])
    
    ##iterate over wards
    briers <- matrix(NA, nrow = length(results[[i]]), ncol = 10)
    for (j in seq_len(length(results[[i]]))) {
      ##iterate over patients in wards
      patientbriers <- matrix(NA, nrow = length(results[[i]][[j]][[1]]), ncol = 10)
      for (l in seq_len(length(results[[i]][[j]][[1]]))) {
        ##ensure only the correct days are taken
        if (q != 11) {
          if ((results[[i]][[j]][[4]]$detection[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]] -
               results[[i]][[j]][[4]]$admission[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]]) != (q-1)) {
            next()
          }
        } else {
          if ((results[[i]][[j]][[4]]$detection[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]] -
               results[[i]][[j]][[4]]$admission[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]]) < 10) {
            next()
          }
        }
        ##adjust priors for transmission trees and correct posterior
        if (results[[i]][[j]][[1]][[l]][[3]] == "H") {
          targettrees <- unique(results[[i]][[j]][[4]]$treeID[!results[[i]][[j]][[4]]$treeID %in% results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]]])
          treeposterior <- treeprevprior <- treedefaultprior <- rep(0, length(targettrees) + 2)
          names(treeposterior) <- names(treeprevprior) <- names(treedefaultprior) <- c(as.character(targettrees), "H", "C")
          for (k in seq_len(length(targettrees))) {
            treeprevprior[k] <- sum(results[[i]][[j]][[4]]$treeID %in% targettrees[k]) * results[[i]][[j]][[3]][1]
            treedefaultprior[k] <- sum(results[[i]][[j]][[4]]$treeID %in% targettrees[k]) * results[[i]][[j]][[2]][1]
          }
          treeprevprior[length(treeprevprior) - 1] <- results[[i]][[j]][[3]][1]
          treedefaultprior[length(treedefaultprior) - 1] <- results[[i]][[j]][[2]][1]
          treeprevprior[length(treeprevprior)] <- results[[i]][[j]][[3]][length(results[[i]][[j]][[3]])]
          treedefaultprior[length(treedefaultprior)] <- 0.5
          for(k in seq_len(length(results[[i]][[j]][[1]][[l]][[1]]) - 2)) {
            treeposterior[names(treeposterior) %in% results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])[k]]] <-
              treeposterior[names(treeposterior) %in% results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])[k]]] +
              results[[i]][[j]][[1]][[l]][[1]][k]
          }
          treeposterior[length(treeposterior) - 1] <- results[[i]][[j]][[1]][[l]][[1]][length(results[[i]][[j]][[1]][[l]][[1]]) - 1]
          treeposterior[length(treeposterior)] <- results[[i]][[j]][[1]][[l]][[1]][length(results[[i]][[j]][[1]][[l]][[1]])]
        } else {
          targettrees <- unique(results[[i]][[j]][[4]]$treeID)
          treeposterior <- treeprevprior <- treedefaultprior <- rep(0, length(targettrees) + 2)
          names(treeposterior) <- names(treeprevprior) <- names(treedefaultprior) <- c(as.character(targettrees), "H", "C")
          for (k in seq_len(length(targettrees))) {
            treeprevprior[k] <- sum(results[[i]][[j]][[4]]$treeID %in% targettrees[k]) * results[[i]][[j]][[3]][1]
            treedefaultprior[k] <- sum(results[[i]][[j]][[4]]$treeID %in% targettrees[k]) * results[[i]][[j]][[2]][1]
            if (targettrees[k] == results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]]) {
              treeprevprior[k] <- treeprevprior[k] - results[[i]][[j]][[3]][1]
              treedefaultprior[k] <- treedefaultprior[k] - results[[i]][[j]][[2]][1]
            }
          }
          treeprevprior[length(treeprevprior) - 1] <- results[[i]][[j]][[3]][1]
          treedefaultprior[length(treedefaultprior) - 1] <- results[[i]][[j]][[2]][1]
          treeprevprior[length(treeprevprior)] <- results[[i]][[j]][[3]][length(results[[i]][[j]][[3]])]
          treedefaultprior[length(treedefaultprior)] <- 0.5
          for(k in seq_len(length(results[[i]][[j]][[1]][[l]][[1]]) - 2)) {
            treeposterior[names(treeposterior) %in% results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])[k]]] <-
              treeposterior[names(treeposterior) %in% results[[i]][[j]][[4]]$treeID[results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])[k]]] +
              results[[i]][[j]][[1]][[l]][[1]][k]
          }
          treeposterior[length(treeposterior) - 1] <- results[[i]][[j]][[1]][[l]][[1]][length(results[[i]][[j]][[1]][[l]][[1]]) - 1]
          treeposterior[length(treeposterior)] <- results[[i]][[j]][[1]][[l]][[1]][length(results[[i]][[j]][[1]][[l]][[1]])]
        }
        ##get brier scores
        prevpriors <- t(as.matrix(results[[i]][[j]][[3]]))
        defaultpriors <- t(as.matrix(results[[i]][[j]][[2]]))
        colnames(prevpriors) <- colnames(defaultpriors) <- names(results[[i]][[j]][[1]][[l]][[1]])
        NOSTRAnosocomial <- matrix(c(1-results[[i]][[j]][[1]][[l]][[1]]["C"], results[[i]][[j]][[1]][[l]][[1]]["C"]), nrow = 1, ncol = 2)
        prevnosocomial <- matrix(c(1-results[[i]][[j]][[3]][length(results[[i]][[j]][[3]])], results[[i]][[j]][[3]][length(results[[i]][[j]][[3]])]), nrow = 1, ncol = 2)
        defaultnosocomial <- matrix(c(0.5, 0.5), nrow = 1, ncol = 2)
        if(results[[i]][[j]][[4]]$detection[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]] -
           results[[i]][[j]][[4]]$admission[results[[i]][[j]][[4]]$ID %in% results[[i]][[j]][[4]]$ID[!results[[i]][[j]][[4]]$ID %in% names(results[[i]][[j]][[1]][[l]][[1]])]] < 3) {
          thumbnosocomial <- matrix(c(0, 1), nrow = 1, ncol = 2)
        } else {
          thumbnosocomial <- matrix(c(1, 0), nrow = 1, ncol = 2)
        }
        colnames(NOSTRAnosocomial) <- colnames(prevnosocomial) <- colnames(defaultnosocomial) <- colnames(thumbnosocomial) <- c("H", "C")
        patientbriers[l,1] <- mlr3measures::mbrier(factor(ifelse(results[[i]][[j]][[1]][[l]][[2]] == "C", "C", "H"), levels = c("H", "C")),
                                                   NOSTRAnosocomial)
        patientbriers[l,2] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[2]], levels = names(results[[i]][[j]][[1]][[l]][[1]])),
                                                   t(as.matrix(results[[i]][[j]][[1]][[l]][[1]])))
        patientbriers[l,3] <- mlr3measures::mbrier(factor(ifelse(results[[i]][[j]][[1]][[l]][[2]] == "C", "C", "H"), levels = c("H","C")),
                                                   prevnosocomial)
        patientbriers[l,4] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[2]], levels = names(results[[i]][[j]][[1]][[l]][[1]])),
                                                   prevpriors)
        patientbriers[l,5] <- mlr3measures::mbrier(factor(ifelse(results[[i]][[j]][[1]][[l]][[2]] == "C", "C", "H"), levels = c("H","C")),
                                                   defaultnosocomial)
        patientbriers[l,6] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[2]], levels = names(results[[i]][[j]][[1]][[l]][[1]])),
                                                   defaultpriors)
        patientbriers[l,7] <- mlr3measures::mbrier(factor(ifelse(results[[i]][[j]][[1]][[l]][[2]] == "C", "C", "H"), levels = c("H","C")),
                                                   thumbnosocomial)
        patientbriers[l,8] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[3]], levels = names(treeposterior)),
                                                   t(as.matrix(treeposterior)))
        patientbriers[l,9] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[3]], levels = names(treeprevprior)),
                                                   t(as.matrix(treeprevprior)))
        patientbriers[l,10] <- mlr3measures::mbrier(factor(results[[i]][[j]][[1]][[l]][[3]], levels = names(treedefaultprior)),
                                                    t(as.matrix(treedefaultprior)))
      }
      if (all(is.na(patientbriers[,1]))) {
        briers[j,] <- rep(NA, length(briers[j,]))
      } else {
        briers[j,] <- colMeans(patientbriers, na.rm = T)
      }
    }
    if (all(is.na(briers[,1]))) {
      simulationresults_perday[i,4:13] <- rep(NA, length(briers[j,]))
    } else {
      simulationresults_perday[i,4:13] <- colMeans(briers, na.rm = T)
    }
  }
  simulationresults_perday_list[[q]] <- simulationresults_perday
}

final_simulationresults_perday <- as.data.frame(do.call("rbind", simulationresults_perday_list))
colnames(final_simulationresults_perday) <- c("Prevalence", "ParamSet", "Repeat", "NosocomialityScore",
                                              "GeneralScore", "PrevPriorNosocomialityScore", "PrevPriorGeneralScore",
                                              "DefaultPriorNosocomialityScore", "DefaultPriorGeneralScore", "RuleofThumbScore",
                                              "NOSTRATreeScore", "PrevTreeScore", "DefaultTreeScore")

final_simulationresults_perday$Day <- rep(0:10, each = 180)

##generate figures
library(ggplot2)
library(tidyverse)

final_simulationresults_perday$Prevalence[final_simulationresults_perday$Prevalence == 1] <- "High"
final_simulationresults_perday$Prevalence[final_simulationresults_perday$Prevalence == 2] <- "Low"
final_simulationresults_perday$Prevalence[final_simulationresults_perday$Prevalence == 3] <- "Intermediate"

final_simulationresults_perday <- final_simulationresults_perday[complete.cases(final_simulationresults_perday),]

figuredata <- pivot_longer(final_simulationresults_perday, cols = c("NosocomialityScore",
                                                                    "GeneralScore", "PrevPriorNosocomialityScore", "PrevPriorGeneralScore",
                                                                    "DefaultPriorNosocomialityScore", "DefaultPriorGeneralScore", "RuleofThumbScore",
                                                                    "NOSTRATreeScore", "PrevTreeScore", "DefaultTreeScore"))

figuredata$Type <- rep(NA, nrow(figuredata))
figuredata$Type[figuredata$name %in% c("NosocomialityScore", "PrevPriorNosocomialityScore", "DefaultPriorNosocomialityScore", "RuleofThumbScore")] <- "Nosocomiality Detection"
figuredata$Type[figuredata$name %in% c("GeneralScore", "PrevPriorGeneralScore", "DefaultPriorGeneralScore")] <- "Source Detection"
figuredata$Type[figuredata$name %in% c("NOSTRATreeScore", "PrevTreeScore", "DefaultTreeScore")] <- "Transmission Tree Detection"

##tidy up names for figure
figuredata$name[figuredata$name %in% c("NOSTRATreeScore")] <- "NOSTRA"
figuredata$name[figuredata$name %in% c("NosocomialityScore")] <- "NOSTRA"
figuredata$name[figuredata$name %in% c("GeneralScore")] <- "NOSTRA"
figuredata$name[figuredata$name %in% c("PrevPriorNosocomialityScore")] <- "Prevalence Prior"
figuredata$name[figuredata$name %in% c("PrevTreeScore")] <- "Prevalence Prior"
figuredata$name[figuredata$name %in% c("PrevPriorGeneralScore")] <- "Prevalence Prior"
figuredata$name[figuredata$name %in% c("DefaultPriorNosocomialityScore")] <- "Naïve Prior"
figuredata$name[figuredata$name %in% c("DefaultPriorGeneralScore")] <- "Naïve Prior"
figuredata$name[figuredata$name %in% c("DefaultTreeScore")] <- "Naïve Prior"
figuredata$name[figuredata$name %in% c("RuleofThumbScore")] <- "96hr Categorisation"

figuredata$name <- factor(figuredata$name, levels = c("NOSTRA", "Prevalence Prior", "Naïve Prior", "96hr Categorisation"))

figuredata_means <- figuredata %>% group_by(name, Type, Day) %>% summarise(value=mean(value))

plot2 <- ggplot(figuredata, aes(name, value)) + 
  facet_grid(Day ~ Type, scales = "free_x") +
  ylab("Brier Score") +
  xlab("") +
  theme_bw() +
  geom_jitter(aes(colour = Prevalence), height = 0, width = 0.1) +
  geom_point(data = figuredata_means, size = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
