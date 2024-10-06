rm(list=ls())
source("NOSTRA.R")
library(lubridate)
library(seqinr)
library(pheatmap)

##NOSTRA illustrative analysis

##function generate time difference matrix from sample collection dates
timedifmat <- function(x) {
  dateformat <- as.Date(x, format = "%d/%m/%Y")
  return(do.call("cbind", lapply(1:length(dateformat), function(x) {return(abs(as.numeric(difftime(dateformat[x],dateformat, units = "days"))))})))
}

##the data is not sharable, so this analysis is not replicable 
##the code is included for completeness, but has no identifying data

##read in wards
#A <- read.csv("")

##read in sequences
#seqs <- read.fasta(file = "")

##read in location information
#location <- read.csv("")

##reduce to required sequences
neededseqs <- c(A$sequence_id)
neededseqs <- unique(neededseqs)
neededseqs <- neededseqs[!neededseqs %in% "sequence_missing"]
seqs <- seqs[names(seqs) %in% neededseqs]

##find patient data without corresponding sequence in the sequence data
missingseqs <- neededseqs[!neededseqs %in% names(seqs)]

##TEMPORARY - set them to missing
A$sequence_id[A$sequence_id %in% missingseqs] <- "sequence_missing"

##generate time differences
Atimedifference <-timedifmat(A$sample_collection_date)

colnames(Atimedifference) <- row.names(Atimedifference) <- A$patient_study_id

##generate snp differences - note that due the the symmetry between pairwise consensuses, this is also each sequences' hamming distance to its mutual consensus
SNPs <- matrix(0, nrow = length(seqs), ncol = length(seqs))
alignmentlength <- matrix(0, nrow = length(seqs), ncol = length(seqs))

for (i in 1:length(seqs)) {
  for (j in i:length(seqs)) {
    for (k in 1:length(seqs[[1]])) { 
      ##assumes sequences have been aligned and thus are of the same length
      if (seqs[[i]][k] %in% c("c","g","a","t") & seqs[[j]][k] %in% c("c","g","a","t")) {
        alignmentlength[i,j] <- alignmentlength[i,j] + 1
        if (seqs[[i]][k] != seqs[[j]][k]) {
          SNPs[i,j] <- SNPs[i,j] + 1
        }
      }
    }
    print(paste(i, j, sep = " "))
  }
}

SNPs[lower.tri(SNPs)] <- t(SNPs)[lower.tri(SNPs)]
colnames(SNPs) <- row.names(SNPs) <- names(seqs)

alignmentlength[lower.tri(alignmentlength)] <- t(alignmentlength)[lower.tri(alignmentlength)]
colnames(alignmentlength) <- row.names(alignmentlength) <- names(seqs)


##read in patient locations
locationslong <- read.csv("~/Documents/locations.csv")

positions <- matrix(NA, nrow = abs(as.numeric(difftime(min(as.Date(locationslong$StartDate, format = "%d/%m/%Y")), 
                                                       max(as.Date(locationslong$EndDate, format = "%d/%m/%Y")),units = "days"))) + 1, 
                    ncol = length(unique(locationslong$patient_study_id)) + 1)
positions <- as.data.frame(positions)
colnames(positions) <- c("Date", unique(locationslong$patient_study_id))

##convert to format appropriate for NOSTRA implementation
positions[1,1] <- as.character(min(as.Date(locationslong$StartDate, format = "%d/%m/%Y")))
k <- 2
while (k <= nrow(positions)) {
  positions[k,1] <- as.character(min(as.Date(locationslong$StartDate, format = "%d/%m/%Y")) + (k - 1))
  k <- k + 1
}

for(i in 1:length(unique(locationslong$patient_study_id))) {
  working <- locationslong[locationslong$patient_study_id %in% unique(locationslong$patient_study_id)[i],]
  k <- 1
  for (j in 1:nrow(positions)) {
    if (is.na(as.Date(working$EndDate[k], format = "%d/%m/%Y"))) {
      break()
    } else if (as.Date(positions[j,1], format = "%Y-%m-%d") < as.Date(working$EndDate[k], format = "%d/%m/%Y") & as.Date(positions[j,1], format = "%Y-%m-%d") >= as.Date(working$StartDate[1], format = "%d/%m/%Y")) {
      if (is.na(positions[j,unique(locationslong$patient_study_id)[i]])) {
        positions[j,unique(locationslong$patient_study_id)[i]] <- working$LocationName[k]
      } else {
        positions[j,unique(locationslong$patient_study_id)[i]] <- paste(positions[j,unique(locationslong$patient_study_id)[i]], 
                                                                        working$LocationName[k], sep = "/")
      }
    } else if (as.Date(positions[j,1], format = "%Y-%m-%d") == as.Date(working$EndDate[k], format = "%d/%m/%Y") & as.Date(positions[j,1], format = "%Y-%m-%d") >= as.Date(working$StartDate[1], format = "%d/%m/%Y")) {
      targets <- as.Date(working$StartDate, format = "%d/%m/%Y") %in% as.Date(working$EndDate[k], format = "%d/%m/%Y")
      print(sum(targets))
      if (is.na(positions[j,unique(locationslong$patient_study_id)[i]]) & sum(targets) == 0) {
        positions[j,unique(locationslong$patient_study_id)[i]] <- working$LocationName[k]
      } else if (is.na(positions[j,unique(locationslong$patient_study_id)[i]]) & sum(targets) != 0) {
        positions[j,unique(locationslong$patient_study_id)[i]] <- paste(unique(c(working$LocationName[k], 
                                                                        working$LocationName[targets])), collapse = "/")
      }
    } else if (as.Date(positions[j,1], format = "%Y-%m-%d") > as.Date(working$EndDate[k], format = "%d/%m/%Y") & as.Date(positions[j,1], format = "%Y-%m-%d") >= as.Date(working$StartDate[1], format = "%d/%m/%Y")) {
      k <- k + 1
      if (is.na(positions[j,unique(locationslong$patient_study_id)[i]])) {
        positions[j,unique(locationslong$patient_study_id)[i]] <- working$LocationName[k]
      } else {
        positions[j,unique(locationslong$patient_study_id)[i]] <- paste(positions[j,unique(locationslong$patient_study_id)[i]], 
                                                                        working$LocationName[k], sep = "/")
      }
    }
  }
}

##run without X_H as no admission dates
##convert epi data for correct format for NOSTRA implementation
A$detection <- as.numeric(difftime(as.Date(A$onset_date, format = "%d/%m/%Y"),
                                   as.Date("30/12/2019", format = "%d/%m/%Y")), units = "days") + 1
A$admission <- rep(NA, nrow(A))
A$sequencedate <- as.numeric(difftime(as.Date(A$sample_collection_date, format = "%d/%m/%Y"),
                                      as.Date("30/12/2019", format = "%d/%m/%Y")), units = "days") + 1
A$sequence_id[A$sequence_id == "sequence_missing"] <- NA
A$sequencedate[is.na(A$sequence_id)] <- NA

##fix incorrect date
A$onset_date[8] <- "09/04/2020"

for(i in 1:length(unique(A$patient_study_id))) {
  print(i)
  print(location$StartDate_0[location$patient_study_id %in% unique(A$patient_study_id)[i]][1])
  A$admission[A$patient_study_id == unique(A$patient_study_id)[i]] <- 
    as.numeric(difftime(as.Date(location$StartDate_0[location$patient_study_id %in% unique(A$patient_study_id)[i]][1], format = "%d/%m/%Y"),
                        as.Date("30/12/2019", format = "%d/%m/%Y")), units = "days") + 1
}

##reduce to minimal data - one sequence per person per bout, earliest taken
indexes <- c()
for (i in 1:length(unique(A$patient_study_id))) {
  temp <- A[A$patient_study_id %in% unique(A$patient_study_id)[i],]
  if (nrow(temp) == 1) {
    indexes <- c(indexes, which(A$patient_study_id == unique(A$patient_study_id)[i]))
  } else (
    indexes <- c(indexes, which(A$patient_study_id == unique(A$patient_study_id)[i] & 
                                  as.Date(A$sample_collection_date, format = "%d/%m/%Y") == 
                                  min(as.Date(temp$sample_collection_date, format = "%d/%m/%Y"))))
  )
}

##fix likely sequence error
indexes[indexes == 20] <- 21

##finalise data format
Areduced <- A[indexes,]
Areducedtimedifference <-timedifmat(Areduced$sample_collection_date)
colnames(Areducedtimedifference) <- row.names(Areducedtimedifference) <- Areduced$patient_study_id

##generate NA SNPS for genetic data impact testing
SNPsNA <- SNPs
SNPsNA[SNPsNA >= 0] <- NA

##run complete model
Aposteriorall <- run_model(epidata = Areduced, locationdata = positions, snpdistmat = SNPs, timedistmat = Areducedtimedifference, alignmentlengthmat = alignmentlength, priors = c(rep(1/(2*nrow(Areduced)), nrow(Areduced)), 0.5), all_comparisons = T, startdate = "2019-12-30")
pheatmap(Aposteriorall, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

##run models missing different data components
Aposteriornodates <- run_model(epidata = Areduced, locationdata = NA, snpdistmat = SNPs, timedistmat = Areducedtimedifference, alignmentlengthmat = alignmentlength, priors = c(rep(1/(2*nrow(Areduced)), nrow(Areduced)), 0.5), all_comparisons = T, startdate = "2019-12-30")
Aposteriornodatesnogenetics <- run_model(epidata = Areduced, locationdata = NA, snpdistmat = SNPsNA, alignmentlengthmat = alignmentlength, timedistmat = Areducedtimedifference, priors = c(rep(1/(2*nrow(Areduced)), nrow(Areduced)), 0.5), all_comparisons = T, startdate = "2019-12-30")
Aposteriornogenetics <- run_model(epidata = Areduced, locationdata = positions, snpdistmat = SNPsNA, alignmentlengthmat = alignmentlength, timedistmat = Areducedtimedifference, priors = c(rep(1/(2*nrow(Areduced)), nrow(Areduced)), 0.5), all_comparisons = T, startdate = "2019-12-30")

##calculate and add nosocomial component (1-C)
fixedcolnames <- c(colnames(Aposteriorall), "Nosocomial")
Aposterioralln <- cbind(Aposteriorall, rowSums(Aposteriorall[,1:(ncol(Aposteriorall)-1)], na.rm = T))
Aposteriornodatesn <- cbind(Aposteriornodates, rowSums(Aposteriornodates[,1:(ncol(Aposteriornodates)-1)], na.rm = T))
Aposteriornodatesnogeneticsn <- cbind(Aposteriornodatesnogenetics, rowSums(Aposteriornodatesnogenetics[,1:(ncol(Aposteriornodatesnogenetics)-1)], na.rm = T))
Aposteriornogeneticsn <- cbind(Aposteriornogenetics, rowSums(Aposteriornogenetics[,1:(ncol(Aposteriornogenetics)-1)], na.rm = T))

##set correct colnames
colnames(Aposterioralln) <- colnames(Aposteriornodatesn) <- colnames(Aposteriornodatesnogeneticsn) <- colnames(Aposteriornogeneticsn) <- fixedcolnames

##remove patients with missing location data
Aposteriorallnreduced <- Aposterioralln[!rownames(Aposterioralln) %in% c("CAMP004621", "CAMP004628"),]
Aposteriornodatesn <- Aposteriornodatesn[!rownames(Aposteriornodatesn) %in% c("CAMP004621", "CAMP004628"),]
Aposteriornodatesnogeneticsn <- Aposteriornodatesnogeneticsn[!rownames(Aposteriornodatesnogeneticsn) %in% c("CAMP004621", "CAMP004628"),]
Aposteriornogeneticsn <- Aposteriornogeneticsn[!rownames(Aposteriornogeneticsn) %in% c("CAMP004621", "CAMP004628"),]

##save plots
png("Aall.png", width = 2000, height = round(2000*419/481), res = 350)
pheatmap(Aposterioralln, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
dev.off()

png("Aallfromgenetics.png", width = 2000, height = round(2000*419/481), res = 350)
pheatmap(Aposteriorallnreduced-Aposteriornodatesn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(-0.5, 0.5, by = 0.01)))
dev.off()

png("Aallfromdates.png", width = 2000, height = round(2000*419/481), res = 350)
pheatmap(Aposteriorallnreduced-Aposteriornogeneticsn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(-0.5, 0.5, by = 0.01)))
dev.off()

png("Anodates.png", width = 2000, height = round(2000*419/481), res = 350)
pheatmap(Aposteriornodatesn-Aposteriornodatesnogeneticsn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(-0.5, 0.5, by = 0.01)))
dev.off()

png("Anogenetics.png", width = 2000, height = round(2000*419/481), res = 350)
pheatmap(Aposteriornogeneticsn-Aposteriornodatesnogeneticsn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(-0.5, 0.5, by = 0.01)))
dev.off()

png("Anodatesnogenetics.png", width = 2000, height = round(2000*419/481), res = 350)
pheatmap(Aposteriornodatesnogeneticsn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
dev.off()