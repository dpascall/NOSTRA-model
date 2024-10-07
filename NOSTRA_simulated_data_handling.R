rm(list = ls())
library(lubridate)
library(plyr)
library(igraph)
library(combinat)
set.seed(789)

infectedby <- function (ident, df) {
  infections <- unique(df$ID[df$infectionsource %in% ident])
  if (length(infections) == 0) {
    return(infections)
  } else { 
    infections <- c(infections, do.call("c", lapply(infections, function(x) {return(infectedby(x, df))})))
    return(infections)
  }
}

##set parameters for simulations
mu <- 6.677*10^(-4)/365 * 29811
Neff <- floor(0.5*(exp(2.5)+exp(4.5)))
generationtime <- 5.5

runs <- c(1:1)
prev <- c("High", "Low", "Med")
rep <- c(1:20)

##get all conditions and iterate over them
conditions <- expand.grid(prev, rep, runs)
colnames(conditions) <- c("prev", "rep", "runs")

finaldataformodel <- vector("list", length = nrow(conditions))

for (x in seq_len(nrow(conditions))) {
  print(paste0("Dataset ", x, " of ", length(finaldataformodel)))

  print("Loading data")
  patientdata <- read.csv(paste("simulated_data/run", conditions$runs[x], "/pat_data_full_stay_", conditions$prev[x], conditions$rep[x], ".csv", sep = ""))
  hcwdata <- read.csv(paste("simulated_data/run", conditions$runs[x], "/hcw_data_full_timecourse_", conditions$prev[x], conditions$rep[x], ".csv", sep = ""))
  
  patientdata$Date <- as.Date(patientdata$Date)
  hcwdata$Date <- as.Date(hcwdata$Date)
  
  dataset <- data.frame(cbind(ID = unique(patientdata$ID), 
                              admission = rep(NA, length(unique(patientdata$ID))),
                              infectiondate = rep(NA, length(unique(patientdata$ID))),
                              ward = rep(NA, length(unique(patientdata$ID))),
                              infected = rep(NA, length(unique(patientdata$ID))),
                              detected = rep(NA, length(unique(patientdata$ID))),
                              discharge = rep(NA, length(unique(patientdata$ID))),
                              infectionsource = rep(NA, length(unique(patientdata$ID))),
                              detectiondate = rep(NA, length(unique(patientdata$ID)))))
  dataset$admission <- as.Date(dataset$admission)
  dataset$discharge <- as.Date(dataset$discharge)
  dataset$infectiondate <- as.Date(dataset$infectiondate)
  dataset$detectiondate <- as.Date(dataset$detectiondate)
  
  print("Converting to patient level data")
  ##extract to model usable format
  for (i in seq_along(unique(patientdata$ID))) {
    dataset$admission[i] <- patientdata$Date[patientdata$ID == unique(patientdata$ID)[i]][1]
    dataset$infected[i] <- ifelse(sum(patientdata$Infection.Status.CT[patientdata$ID == unique(patientdata$ID)[i]] != "NONE") != 0, 1, 0)
    dataset$detected[i] <- ifelse(sum(patientdata[patientdata$ID == unique(patientdata$ID)[i],"Detected"]) == 0, 0, 1)
    dataset$discharge[i] <- patientdata$Date[patientdata$ID == unique(patientdata$ID)[i]][length(patientdata$Date[patientdata$ID == unique(patientdata$ID)[i]])] + 1
    dataset$ward[i] <- patientdata$WardNumber[patientdata$ID == unique(patientdata$ID)[i]][1]
    if (dataset$infected[i] == 1) {
      if (sum(patientdata$ID == unique(patientdata$ID)[i] & patientdata$Detected == TRUE & patientdata$Infection.Status.CT != "NONE") != 0) {
        dataset$detectiondate[i] <- patientdata$Date[patientdata$ID == unique(patientdata$ID)[i] & patientdata$Detected == TRUE & patientdata$Infection.Status.CT != "NONE"][1]
        if (patientdata$Community.Infection[patientdata$ID == unique(patientdata$ID)[i] & patientdata$Detected == TRUE & patientdata$Infection.Status.CT != "NONE"][1] == TRUE) {
          
        } else {
          dataset$infectiondate[i] <- patientdata$Date[patientdata$ID == unique(patientdata$ID)[i] & patientdata$Infection.Status.CT != "NONE"][1]
        }
      } else {
        dataset$detected[i] <- 0
        if (patientdata$Community.Infection[patientdata$ID == unique(patientdata$ID)[i] & patientdata$Infection.Status.CT != "NONE"][1] == TRUE) {
          
        } else {
          dataset$infectiondate[i] <- patientdata$Date[patientdata$ID == unique(patientdata$ID)[i] & patientdata$Infection.Status.CT != "NONE"][1]
        }
      }
      if (sum(patientdata[patientdata$ID == unique(patientdata$ID)[i],"Community.Infection"]) != 0) {
        dataset$infectionsource[i] <- "C"
      } else if ("HCW" %in% patientdata[patientdata$ID == unique(patientdata$ID)[i],"Person"]) {
        dataset$infectionsource[i] <- "H"
      } else if ("ED" %in% patientdata[patientdata$ID == unique(patientdata$ID)[i],"Person"]) {
        dataset$infectionsource[i] <- "H"
      } else {
        dataset$infectionsource[i] <- unique(patientdata[patientdata$ID == unique(patientdata$ID)[i] & !is.na(patientdata$PersonID),"PersonID"])[1]
      }
    }
  }
  
  dataset$nosocomial <- ifelse(dataset$infectionsource == "C", 0, 1)
  
  dataset <- dataset[dataset$infected == 1,]
  
  reducedpatientdata <- patientdata[patientdata$ID %in% unique(dataset$ID) & patientdata$Infection.Status.CT != "NONE", c("Date", "ID", "Infection.Status.CT", "Community.Infection", "PersonID", "Detected")]
  reducedhcwdata <- hcwdata[hcwdata$Infection.Status.CT != "NONE", c("Date", "ID", "Infection.Status.CT", "Community.Infection", "PersonID",  "Detected")]
  
  reducedpatientdata$hosttype <- "Patient"
  reducedhcwdata$hosttype <- "HCW"
  
  print("Building transmission trees")
  ##extract data required to build transmission tree
  relevanttimes <- rbind(reducedpatientdata, reducedhcwdata)
  relevanttimes_start <- ddply(relevanttimes, "ID", function(z) head(z,1))
  relevanttimes_detection <- patientdata[patientdata$ID %in% unique(dataset$ID) & patientdata$Infection.Status.CT != "NONE" & patientdata$Detected == "TRUE", c("Date", "ID", "Infection.Status.CT", "Community.Infection", "PersonID",  "Detected")]
  relevanttimes_detection$hosttype <- "Patient"
  relevanttimes_detection <- ddply(relevanttimes_detection, "ID", function(z) head(z,1))
  relevanttimes <- rbind(relevanttimes_start, relevanttimes_detection)
  relevanttimes$Date <- as.Date(relevanttimes$Date)
  relevanttimes <- relevanttimes[!duplicated(relevanttimes),]
  relevanttimes$HospitalInfection <- ifelse(relevanttimes$Community.Infection, FALSE, TRUE)
  relevanttimes <- relevanttimes[order(relevanttimes$Date, relevanttimes$HospitalInfection),]
  relevanttimes$treeID <- rep(NA, nrow(relevanttimes))
  relevanttimes$eventtype <- rep(NA, nrow(relevanttimes))
  
  ##identify linked infections
  ##take first ID
  ##mark all with that ID in same tree
  ##look for infections from first ID
  ##mark all infections from that ID in same tree
  ##recurse
  ##iterate k
  
  for (i in seq_along(unique(relevanttimes$ID))) {
    focal <- relevanttimes[relevanttimes$ID == unique(relevanttimes$ID)[i],]
    if (nrow(focal) == 2) {
      relevanttimes$eventtype[relevanttimes$ID == unique(relevanttimes$ID)[i]] <- c(1,2)
    } else if (relevanttimes$Detected[relevanttimes$ID == unique(relevanttimes$ID)[i]] %in% TRUE) {
      relevanttimes$eventtype[relevanttimes$ID == unique(relevanttimes$ID)[i]] <- 3
    } else {
      relevanttimes$eventtype[relevanttimes$ID == unique(relevanttimes$ID)[i]] <- 1
    }
  }
    
  patientrelevanttimes <- relevanttimes[relevanttimes$hosttype == "Patient",]
  hcwrelevanttimes <- relevanttimes[relevanttimes$hosttype == "HCW",]
  
  patientrelevanttimes$infectionsource <- NA
  for (i in seq_along(unique(patientrelevanttimes$ID))) {
    patientrelevanttimes$infectionsource[patientrelevanttimes$ID == unique(patientrelevanttimes$ID)[i]] <- dataset$infectionsource[dataset$ID == unique(patientrelevanttimes$ID)[i]]
  }
  hcwrelevanttimes$infectionsource <- NA
  
  relevanttimes <- rbind(patientrelevanttimes, hcwrelevanttimes)
  relevanttimes <- relevanttimes[order(relevanttimes$Date, relevanttimes$HospitalInfection),]
  
  k <- 1
  for (i in unique(relevanttimes$ID[relevanttimes$infectionsource %in% c("C", "H")])) {
    if (is.na(relevanttimes$treeID[relevanttimes$ID %in% i][1])) {
      targets <- unique(c(i, infectedby(i, dataset)))
      relevanttimes$treeID[relevanttimes$ID %in% targets] <- k
      k <- k + 1
    }
  }
  
  ##assume any patients not in a transmission tree were infected by someone 
  ##who only existed in hospital between daily reporting periods
  if (length(unique(relevanttimes$ID[is.na(relevanttimes$treeID) & relevanttimes$hosttype == "Patient"])) != 0) {
    stop("Missing infections")
    for (i in unique(relevanttimes$ID[is.na(relevanttimes$treeID) & relevanttimes$hosttype == "Patient"])) {
      if (is.na(relevanttimes$treeID[relevanttimes$ID %in% i][1])) {
        dataset$infectionsource[dataset$ID %in% i] <- "H" ##correct data to note that true infection source not in dataset
        targets <- unique(c(i, infectedby(i, dataset)))
        relevanttimes$treeID[relevanttimes$ID %in% targets] <- k
        k <- k + 1
      }
    }
  }
  
  dataset$treeID <- rep(NA, nrow(dataset))
  for (i in seq_len(nrow(dataset))) {
    dataset$treeID[i] <- relevanttimes$treeID[relevanttimes$ID %in% dataset$ID[i]][1]
  }
  
  ##event type reference
  eventtypes <- data.frame(eventID = 1:3, type = c("I", "D", "ID"))
  
  ##construct transmission trees
  transtrees <- vector("list", length = k - 1)
  
  for (i in seq_along(transtrees)) {
    treeelements <- relevanttimes[relevanttimes$treeID %in% i,]
    if (sum(treeelements$eventtype %in% c(2,3)) <= 1) {
      next()
    }
    
    ##reorder to avoid reporting time issues
    ##move community infection first
    ##find all infections directly linked to community infection
    ##sort by time
    ##repeat for new infections
    if (sum(treeelements$Community.Infection) >= 1) {
      buildingtreeelements <- treeelements[treeelements$Community.Infection == TRUE,]
    } else {
      buildingtreeelements <- treeelements[treeelements$infectionsource == "H",]
    }
    targets <- infectedby(buildingtreeelements$ID[1], treeelements)
    
    for(k in seq_along(targets)) {
      working <- treeelements[treeelements$ID %in% targets[k],]
      buildingtreeelements <- rbind(buildingtreeelements, working[order(working$Date),])
    }
    treeelements <- buildingtreeelements
    
    transtrees[[i]] <- make_empty_graph()
    for (j in seq_len(nrow(treeelements))) {
      if (j == 1) {
        if (treeelements$eventtype[j] == 1) {
          transtrees[[i]] <- add_vertices(transtrees[[i]], 1)
          V(transtrees[[i]])$Event <- "Infection"
          V(transtrees[[i]])$From <- NA
          V(transtrees[[i]])$To <- treeelements$ID[j]
          V(transtrees[[i]])$Date <- treeelements$Date[j]
        } else {
          transtrees[[i]] <- add_vertices(transtrees[[i]], 2)
          transtrees[[i]] <- add_edges(transtrees[[i]], c(1, 2))
          E(transtrees[[i]])$ID <- treeelements$ID[j]
          E(transtrees[[i]])$weight <- 0
          V(transtrees[[i]])$Event <- c("Infection", "Detection")
          V(transtrees[[i]])$From <- c(NA, treeelements$ID[j])
          V(transtrees[[i]])$To <- rep(treeelements$ID[j], 2)
          V(transtrees[[i]])$Date <- rep(treeelements$Date[j], 2)
        }
      } else {
        if (treeelements$eventtype[j] == 1) {
          transtrees[[i]] <- add_vertices(transtrees[[i]], 1)
          treedataframe <- as_data_frame(transtrees[[i]], "vertices")
          target <- max(max(which(treedataframe$From == treeelements$PersonID[j])), max(which(treedataframe$To == treeelements$PersonID[j])))
          transtrees[[i]] <- add_edges(transtrees[[i]], c(target, vcount(transtrees[[i]])))
          E(transtrees[[i]])$ID[ecount(transtrees[[i]])] <- treeelements$PersonID[j]
          V(transtrees[[i]])$Event[vcount(transtrees[[i]])] <- c("Infection")
          V(transtrees[[i]])$From[vcount(transtrees[[i]])] <- treeelements$PersonID[j]
          V(transtrees[[i]])$To[vcount(transtrees[[i]])] <- treeelements$ID[j]
          V(transtrees[[i]])$Date[vcount(transtrees[[i]])] <- treeelements$Date[j]
          E(transtrees[[i]])$weight[ecount(transtrees[[i]])] <- abs(as.numeric(difftime(as.Date(V(transtrees[[i]])$Date[target]), treeelements$Date[j], units="days")))
        } else if (treeelements$eventtype[j] == 2) {
          transtrees[[i]] <- add_vertices(transtrees[[i]], 1)
          treedataframe <- as_data_frame(transtrees[[i]], "vertices")
          target <- max(max(which(treedataframe$From == treeelements$ID[j])), max(which(treedataframe$To == treeelements$ID[j])))
          transtrees[[i]] <- add_edges(transtrees[[i]], c(target, vcount(transtrees[[i]])))
          E(transtrees[[i]])$ID[ecount(transtrees[[i]])] <- treeelements$ID[j]
          V(transtrees[[i]])$Event[vcount(transtrees[[i]])] <- c("Detection")
          V(transtrees[[i]])$From[vcount(transtrees[[i]])] <- treeelements$ID[j]
          V(transtrees[[i]])$To[vcount(transtrees[[i]])] <- treeelements$ID[j]
          V(transtrees[[i]])$Date[vcount(transtrees[[i]])] <- treeelements$Date[j]
          E(transtrees[[i]])$weight[ecount(transtrees[[i]])] <- abs(as.numeric(difftime(as.Date(V(transtrees[[i]])$Date[target]), treeelements$Date[j], units="days")))
        } else {
          transtrees[[i]] <- add_vertices(transtrees[[i]], 2)
          treedataframe <- as_data_frame(transtrees[[i]], "vertices")
          target <- max(max(which(treedataframe$From == treeelements$PersonID[j])), max(which(treedataframe$To == treeelements$PersonID[j])))
          transtrees[[i]] <- add_edges(transtrees[[i]], c(target, (vcount(transtrees[[i]])-1), (vcount(transtrees[[i]])-1), vcount(transtrees[[i]])))
          E(transtrees[[i]])$ID[c((ecount(transtrees[[i]])-1), ecount(transtrees[[i]]))] <- rep(treeelements$ID[j], 2)
          V(transtrees[[i]])$Event[c((vcount(transtrees[[i]])-1), vcount(transtrees[[i]]))] <- c("Infection", "Detection")
          V(transtrees[[i]])$From[c((vcount(transtrees[[i]])-1), vcount(transtrees[[i]]))] <- c(treeelements$PersonID[j], treeelements$ID[j])
          V(transtrees[[i]])$To[c((vcount(transtrees[[i]])-1), vcount(transtrees[[i]]))] <- rep(treeelements$ID[j], 2)
          V(transtrees[[i]])$Date[c((vcount(transtrees[[i]])-1), vcount(transtrees[[i]]))] <- rep(treeelements$Date[j], 2)
          E(transtrees[[i]])$weight[c((ecount(transtrees[[i]])-1), ecount(transtrees[[i]]))] <- c( abs(as.numeric(difftime(as.Date(V(transtrees[[i]])$Date[target]), treeelements$Date[j], units="days"))), 0)
        }
      }
    }
    ##reduce to minimal graph - subgraph induced by preserving all links between detections
    detectionverts <- as.numeric(V(transtrees[[i]])[V(transtrees[[i]])$Event %in% "Detection"])
    target <- detectionverts
    for (l in 1:(length(detectionverts)-1)) {
      for (p in (l+1):length(detectionverts)) {
        target <- c(target, as.numeric(shortest_paths(transtrees[[i]], detectionverts[l], detectionverts[p], mode = "all")$vpath[[1]]))
      }
    }
    target <- unique(target)
    transtrees[[i]] <- induced_subgraph(transtrees[[i]], target)
    ##convert time to SNPs
    E(transtrees[[i]])$weight <- rpois(1, mu * E(transtrees[[i]])$weight)
  }
  ##warnings about -Inf working as designed
  
  print("Beginning SNP calculation")
  ##work out appropriate time differences from start time for SNP calculation - if singleton time from start to detection
  ##if linked earliest relevant event in tree
  dataset$timediff <- as.numeric(abs(difftime(as.Date(patientdata$Date[1]), as.Date(dataset$detectiondate), units = "days")))
  dataset$correctedtimediff <- dataset$timediff
  
  for (i in seq_len(nrow(dataset))) {
    if (is.na(dataset$timediff[i])) {
      next()
    }
    if (!is.null(transtrees[[dataset$treeID[i]]])) {
      dataset$correctedtimediff[i] <- as.numeric(abs(difftime(as.Date(patientdata$Date[1]),min(as.Date(V(transtrees[[dataset$treeID[i]]])$Date)))))
    }
  }
  
  detectedcases <- dataset[dataset$detected == 1,]

  ##fix for when infection source is not detected and thus should be H
  detectedcases$infectionsource[!detectedcases$infectionsource %in% detectedcases$ID & !detectedcases$infectionsource %in% c("C", "H")] <- "H"
  
  SNPs <- matrix(rep(0, length(unique(detectedcases$ID))^2), nrow = length(unique(unique(detectedcases$ID))), ncol = length(unique(unique(detectedcases$ID))))
  colnames(SNPs) <- rownames(SNPs) <- unique(detectedcases$ID)
  
  treeSNPs <- matrix(rep(NA, length(unique(detectedcases$treeID))^2), nrow = length(unique(detectedcases$treeID)), ncol = length(unique(detectedcases$treeID)))
  colnames(treeSNPs) <- rownames(treeSNPs) <- sort(unique(detectedcases$treeID))
  
  print("Generating location data")
  ##generate location data
  start <- min(c(patientdata$Date, hcwdata$Date))
  end <- max(c(detectedcases$admission, detectedcases$discharge, detectedcases$detectiondate))
  locations <- matrix(NA, nrow = length(start:end), ncol = nrow(detectedcases))
  colnames(locations) <- detectedcases$ID
  rownames(locations) <- as.character(as.Date(start:end))
  
  for (g in seq_len(nrow(detectedcases))) {
    locations[which(rownames(locations) == as.character(as.Date(detectedcases$admission[g]))):(which(rownames(locations) == as.character(as.Date(detectedcases$admission[g])))+length(detectedcases$admission[g]:detectedcases$discharge[g])-1),g] <- rep(detectedcases$ward[g], length(detectedcases$admission[g]:detectedcases$discharge[g]))
  }
  
  locations <- as.data.frame(locations)
  locations <- cbind(Date = rownames(locations), locations)
  rownames(locations) <- seq_len(nrow(locations))
  
  ##remove to reduce memory usage
  rm(hcwdata, reducedhcwdata)

  print("Calculating numbers of SNPs between transmission trees")
  ##tree specific SNP calculation - to minimise compute, calculate differences only between trees on the same ward
  for (q in seq_along(unique(detectedcases$ward))) {
    combinations <- combn(sort(unique(detectedcases[detectedcases$ward %in% unique(detectedcases$ward)[q],]$treeID)), 2)
    if (is.null(ncol(combinations))) {
      combinations <- as.matrix(combinations, nrow = 2, ncol = 1)
    }
    for (i in seq_len(ncol(combinations))) {
      if (is.na(treeSNPs[rownames(treeSNPs) %in% unique(as.character(detectedcases$treeID[as.character(detectedcases$treeID) %in% as.character(combinations[1,i])])), 
                         colnames(treeSNPs) %in% unique(as.character(detectedcases$treeID[as.character(detectedcases$treeID) %in% as.character(combinations[2,i])]))])) {
        treeSNPs[rownames(treeSNPs) %in% unique(as.character(detectedcases$treeID[as.character(detectedcases$treeID) %in% as.character(combinations[1,i])])), 
                         colnames(treeSNPs) %in% unique(as.character(detectedcases$treeID[as.character(detectedcases$treeID) %in% as.character(combinations[2,i])]))] <-
          rnbinom(1, 1, 1/(1 + (2 * mu * Neff * generationtime))) +
          rpois(1, mu * detectedcases$correctedtimediff[as.character(detectedcases$treeID) %in% detectedcases$treeID[as.character(detectedcases$treeID) %in% as.character(combinations[1,i])]][1]) +
          rpois(1, mu * detectedcases$correctedtimediff[as.character(detectedcases$treeID) %in% detectedcases$treeID[as.character(detectedcases$treeID) %in% as.character(combinations[2,i])]][1])
      }
    }
  }
  
  diag(treeSNPs) <- 0
  treeSNPs[lower.tri(treeSNPs)] <- t(treeSNPs)[lower.tri(treeSNPs)]
  
  print("Calculating numbers of SNPs within transmission trees")
  ##get snps from on transmission trees - calculate SNPs between isolates within trees, and SNPs to root of tree for all isolates in different trees
  for (i in seq_along(transtrees)) {
    if (!is.null(transtrees[[i]])) {
      ##break up large transmission trees by ward to minimise compute and memory usage
      ##possible as only cases on the same ward are compared genetically
      ##does not work for any future analyses where this is not the case!!
      IDsintree <- unique(c(V(transtrees[[i]])$From, V(transtrees[[i]])$To))
      wardsintree <- unique(dataset$ward[dataset$ID %in% IDsintree])
      wardtrees <- vector("list", length = length(wardsintree))
      for (j in seq_len(length(wardsintree))) {
        if (sum(IDsintree %in% unique(dataset$ID[dataset$ward %in% wardsintree[j]])) < 2) {
          next()
        }
        ##subset tree to tree within ward
        wardverts <- as.numeric(V(transtrees[[i]])[V(transtrees[[i]])$To %in% unique(dataset$ID[dataset$ward %in% wardsintree[j]])|V(transtrees[[i]])$To %in% unique(dataset$ID[dataset$ward %in% wardsintree[j]])])
        target <- wardverts
        for (l in 1:(length(wardverts)-1)) {
          for (p in (l+1):length(wardverts)) {
            target <- c(target, as.numeric(shortest_paths(transtrees[[i]], wardverts[l], wardverts[p], mode = "all")$vpath[[1]]))
          }
        }
        target <- unique(target)
        wardtrees[[j]] <- induced_subgraph(transtrees[[i]], target)
        ##calculate difference between ward subtrees
        detectionverts <- as.numeric(V(wardtrees[[j]])[V(wardtrees[[j]])$Event %in% "Detection" & (V(wardtrees[[j]])$To %in% unique(dataset$ID[dataset$ward %in% wardsintree[j]])|V(wardtrees[[j]])$To %in% unique(dataset$ID[dataset$ward %in% wardsintree[j]]))])
        if (length(detectionverts) < 2) {
          next()
        }
        detectioncombinations <- as.matrix(combn(unique(detectionverts), 2))
        if (is.null(ncol(detectioncombinations))) {
          detectioncombinations <- as.matrix(detectioncombinations, nrow = 2, ncol = 1)
        }
        for (l in seq_len(ncol(detectioncombinations))) {
          SNPs[rownames(SNPs) %in% as.character(V(wardtrees[[j]])$To[detectioncombinations[1,l]]), colnames(SNPs) %in% as.character(V(wardtrees[[j]])$To[detectioncombinations[2,l]])] <-
            sum(E(wardtrees[[j]])$weight[shortest_paths(wardtrees[[j]], detectioncombinations[1,l], detectioncombinations[2,], mode = "all", output = "both")$epath[[1]]])
        }
        for (l in seq_along(detectionverts)) {
          SNPs[rownames(SNPs) %in% as.character(detectionverts[l]), !colnames(SNPs) %in% as.character(unique(V(wardtrees[[j]])$To[V(wardtrees[[j]])$Event %in% "Detection"]))] <- sum(E(wardtrees[[j]])$weight[shortest_paths(wardtrees[[j]], detectionverts[l], which.min(as.Date(V(wardtrees[[j]])$Date)), mode = "all", output = "both")$epath[[1]]])
        }
      }
    }
  }
  
  diag(SNPs) <- 0
  
  print("Generating final SNP matrix")
  ##to minimise compute, sum between transmission tree and within transmission tree SNPs per ward 
  for (q in seq_along(unique(detectedcases$ward))) {
    combinations <- combn(unique(detectedcases[detectedcases$ward %in% unique(detectedcases$ward)[q],]$ID), 2)
    if (is.null(ncol(combinations))) {
      combinations <- as.matrix(combinations, nrow = 2, ncol = 1)
    }
    for (i in seq_len(ncol(combinations))) {
        SNPs[rownames(SNPs) %in% as.character(combinations[1,i]), colnames(SNPs) %in% as.character(combinations[2,i])] <-
        SNPs[rownames(SNPs) %in% as.character(combinations[1,i]), colnames(SNPs) %in% as.character(combinations[2,i])] +
        treeSNPs[rownames(treeSNPs) %in% as.character(detectedcases$treeID[as.character(detectedcases$ID) %in% as.character(combinations[1,i])]), 
                 colnames(treeSNPs) %in% as.character(detectedcases$treeID[as.character(detectedcases$ID) %in% as.character(combinations[2,i])])]
    }
  }

  SNPs[lower.tri(SNPs)] <- t(SNPs)[lower.tri(SNPs)]

  print("Saving data for dataset")
  finaldataformodel[[x]] <- list(detectedcases, SNPs, conditions[x,], locations, start)
}
save(finaldataformodel, file = "finaldataformodel.Rdata")
