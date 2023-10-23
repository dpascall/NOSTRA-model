source("HOCI_prototype.R")
library(lubridate)
library(seqinr)
library(pheatmap)

##A2B analysis

#function generate time difference matrix from sample collection dates
timedifmat <- function(x) {
  dateformat <- as.Date(x, format = "%d/%m/%Y")
  return(do.call("cbind", lapply(1:length(dateformat), function(x) {return(abs(as.numeric(difftime(dateformat[x],dateformat, units = "days"))))})))
}

##read in clusters
A <- read.csv("~/Downloads/Ward_data/cluster_A_genetic_temporal_data_NoPII_20200807.csv")
B <- read.csv("~/Downloads/Ward_data/cluster_B_genetic_temporal_data_NoPII_20200807.csv")
C <- read.csv("~/Downloads/Ward_data/cluster_C_genetic_temporal_data_NoPII_20200807.csv")
D <- read.csv("~/Downloads/Ward_data/cluster_D_genetic_temporal_data_NoPII_20200807.csv")
E <- read.csv("~/Downloads/Ward_data/cluster_E_extended_genetic_temporal_data_NoPII_20200807.csv")

##read in sequences
#seqs <- read.fasta(file = "~/Downloads/Ward_data/Seqs_editN20_manali.fa")
seqs <- read.fasta(file = "~/Downloads/Ward_data/AlignedWithConsensus.fasta")


##read in location information
location <- read.csv("~/Downloads/Ward_data/ward_movement_network_edit_anonymised_20200811_NoPII.csv")

##reduce to required sequences
neededseqs <- c(A$sequence_id, B$sequence_id, C$sequence_id, D$sequence_id, E$sequence_id)
neededseqs <- unique(neededseqs)
neededseqs <- neededseqs[!neededseqs %in% "sequence_missing"]
seqs <- seqs[names(seqs) %in% neededseqs]

##find patient data without corresponding sequence in the sequence data
missingseqs <- neededseqs[!neededseqs %in% names(seqs)]

##TEMPORARY - set them to missing
A$sequence_id[A$sequence_id %in% missingseqs] <- "sequence_missing"
B$sequence_id[B$sequence_id %in% missingseqs] <- "sequence_missing"
C$sequence_id[C$sequence_id %in% missingseqs] <- "sequence_missing"
D$sequence_id[D$sequence_id %in% missingseqs] <- "sequence_missing"
E$sequence_id[E$sequence_id %in% missingseqs] <- "sequence_missing"

##if sequence but no sample collection date, take received date
C$sample_collection_date[38] <- C$sample_received_date[38]
E$sample_collection_date[38] <- E$sample_received_date[38]

##generate time differences
Atimedifference <-timedifmat(A$sample_collection_date)
Btimedifference <-timedifmat(B$sample_collection_date)
Ctimedifference <-timedifmat(C$sample_collection_date)
Dtimedifference <-timedifmat(D$sample_collection_date)
Etimedifference <-timedifmat(E$sample_collection_date)

colnames(Atimedifference) <- row.names(Atimedifference) <- A$patient_study_id
colnames(Btimedifference) <- row.names(Btimedifference) <- B$patient_study_id
colnames(Ctimedifference) <- row.names(Ctimedifference) <- C$patient_study_id
colnames(Dtimedifference) <- row.names(Dtimedifference) <- D$patient_study_id
colnames(Etimedifference) <- row.names(Etimedifference) <- E$patient_study_id

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


##read in fixed file
locationslong <- read.csv("~/Documents/locations.csv")

positions <- matrix(NA, nrow = abs(as.numeric(difftime(min(as.Date(locationslong$StartDate, format = "%d/%m/%Y")), max(as.Date(locationslong$EndDate, format = "%d/%m/%Y")),units = "days"))) + 1, ncol = length(unique(locationslong$patient_study_id)) + 1)
positions <- as.data.frame(positions)
colnames(positions) <- c("Date", unique(locationslong$patient_study_id))

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

##convert to correct format
##run without X_H as no admission dates
##assume constant contact for initial run
A$detection <- as.numeric(difftime(as.Date(A$onset_date, format = "%d/%m/%Y"),as.Date("01/01/2020", format = "%d/%m/%Y")), units = "days") + 1
A$admission <- rep(NA, nrow(A))
A$sequencedate <- as.numeric(difftime(as.Date(A$sample_collection_date, format = "%d/%m/%Y"),as.Date("01/01/2020", format = "%d/%m/%Y")), units = "days") + 1
A$sequence_id[A$sequence_id == "sequence_missing"] <- NA
A$sequencedate[is.na(A$sequence_id)] <- NA

##fix incorrect date
A$onset_date[8] <- "09/04/2020"

##reduce to minimal data - one sequence per person per bout, earliest taken
indexes <- c()
for (i in 1:length(unique(A$patient_study_id))) {
  temp <- A[A$patient_study_id %in% unique(A$patient_study_id)[i],]
  if (nrow(temp) == 1) {
    indexes <- c(indexes, which(A$patient_study_id == unique(A$patient_study_id)[i]))
  } else (
    indexes <- c(indexes, which(A$patient_study_id == unique(A$patient_study_id)[i] & 
                                  as.Date(A$sample_collection_date, format = "%d/%m/%Y") == min(as.Date(temp$sample_collection_date, format = "%d/%m/%Y"))))
  )
}
##fix likely error
indexes[indexes == 20] <- 21
Areduced <- A[indexes,]
Areducedtimedifference <-timedifmat(Areduced$sample_collection_date)
colnames(Areducedtimedifference) <- row.names(Areducedtimedifference) <- Areduced$patient_study_id

SNPsNA <- SNPs
SNPsNA[SNPsNA >= 0] <- NA 

Aposteriorall <- run_model(epidata = Areduced, locationdata = positions, snpdistmat = SNPs, timedistmat = Areducedtimedifference, alignmentlengthmat = alignmentlength, priors = rep(1/(nrow(Areduced)+1), nrow(Areduced)+1), all_comparisons = T, startdate = "2020-01-01")
pheatmap(Aposteriorall, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

Aposteriornodates <- run_model(epidata = Areduced, locationdata = NA, snpdistmat = SNPs, timedistmat = Areducedtimedifference, alignmentlengthmat = alignmentlength, priors = rep(1/(nrow(Areduced)+1), nrow(Areduced)+1), all_comparisons = T, startdate = "2020-01-01")
Aposteriornodatesnogenetics <- run_model(epidata = Areduced, locationdata = NA, snpdistmat = SNPsNA, alignmentlengthmat = alignmentlength, timedistmat = Areducedtimedifference, priors = rep(1/(nrow(Areduced)+1), nrow(Areduced)+1), all_comparisons = T, startdate = "2020-01-01")
Aposteriornogenetics <- run_model(epidata = Areduced, locationdata = positions, snpdistmat = SNPsNA, alignmentlengthmat = alignmentlength, timedistmat = Areducedtimedifference, priors = rep(1/(nrow(Areduced)+1), nrow(Areduced)+1), all_comparisons = T, startdate = "2020-01-01")

fixedcolnames <- c(colnames(Aposteriorall), "Nosocomial")
# 
# ##set diagonals to 0 as you cannot infect yourself (aka 0 prior probability)
# #Aposterior[is.na(Aposterior)] <- 0
# #Aposterior[Aposterior == 0] <- NA
# 
Aposterioralln <- cbind(Aposteriorall, rowSums(Aposteriorall[,1:(ncol(Aposteriorall)-1)], na.rm = T))
Aposteriornodatesn <- cbind(Aposteriornodates, rowSums(Aposteriornodates[,1:(ncol(Aposteriornodates)-1)], na.rm = T))
Aposteriornodatesnogeneticsn <- cbind(Aposteriornodatesnogenetics, rowSums(Aposteriornodatesnogenetics[,1:(ncol(Aposteriornodatesnogenetics)-1)], na.rm = T))
Aposteriornogeneticsn <- cbind(Aposteriornogenetics, rowSums(Aposteriornogenetics[,1:(ncol(Aposteriornogenetics)-1)], na.rm = T))

colnames(Aposterioralln) <- colnames(Aposteriornodatesn) <- colnames(Aposteriornodatesnogeneticsn) <- colnames(Aposteriornogeneticsn) <- fixedcolnames
# 
# pheatmap(Aposteriorall, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
# pheatmap(Aposteriornodates, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
# pheatmap(Aposteriornodatesnogenetics, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
# pheatmap(Aposteriornogenetics, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
# 
png("~/Documents/HOCI_prototype/Aall.png", width = 2000, height = round(2000*419/481), res = 300)
pheatmap(Aposteriornodatesnogeneticsn-Aposterioralln, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(-0.2, 0.2, by = 0.01)))
dev.off()

png("~/Documents/HOCI_prototype/Aallnochange.png", width = 2000, height = round(2000*419/481), res = 300)
pheatmap(Aposterioralln, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
dev.off()

png("~/Documents/HOCI_prototype/Anodates.png", width = 2000, height = round(2000*419/481), res = 300)
pheatmap(Aposteriornodatesnogeneticsn-Aposteriornodatesn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(-0.2, 0.2, by = 0.01)))
dev.off()

png("~/Documents/HOCI_prototype/Anogenetics.png", width = 2000, height = round(2000*419/481), res = 300)
pheatmap(Aposteriornodatesnogeneticsn-Aposteriornogeneticsn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(-0.2, 0.2, by = 0.01)))
dev.off()

png("~/Documents/HOCI_prototype/Anodatesnogenetics.png", width = 2000, height = round(2000*419/481), res = 300)
pheatmap(Aposteriornodatesnogeneticsn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
dev.off()

##run without X_H as no admission dates
##assume constant contact for initial run
B$detection <- as.numeric(difftime(as.Date(B$onset_date, format = "%d/%m/%Y"),as.Date("01/01/2020", format = "%d/%m/%Y")), units = "days") + 1
B$admission <- rep(NA, nrow(B))
B$sequencedate <- as.numeric(difftime(as.Date(B$sample_collection_date, format = "%d/%m/%Y"),as.Date("01/01/2020", format = "%d/%m/%Y")), units = "days") + 1
B$sequence_id[B$sequence_id == "sequence_missing"] <- NA
B$sequencedate[is.na(B$sequence_id)] <- NA

##reduce B
Breduced <- B[-5,]
Breducedtimedifference <-timedifmat(Breduced$sample_collection_date)
colnames(Breducedtimedifference) <- row.names(Breducedtimedifference) <- Breduced$patient_study_id

Bposteriorall <- run_model(epidata = Breduced, locationdata = positions, snpdistmat = SNPs, timedistmat = Breducedtimedifference, priors = rep(1/(nrow(Breduced)+1), nrow(Breduced)+1), all_comparisons = T, startdate = "2020-01-01")
Bposteriornodates <- run_model(epidata = Breduced, locationdata = NA, snpdistmat = SNPs, timedistmat = Breducedtimedifference, priors = rep(1/(nrow(Areduced)+1), nrow(Areduced)+1), all_comparisons = T, startdate = "2020-01-01")
Bposteriornodatesnogenetics <- run_model(epidata = Breduced, locationdata = NA, snpdistmat = SNPsNA, timedistmat = Breducedtimedifference, priors = rep(1/(nrow(Areduced)+1), nrow(Areduced)+1), all_comparisons = T, startdate = "2020-01-01")
Bposteriornogenetics <- run_model(epidata = Breduced, locationdata = positions, snpdistmat = SNPsNA, timedistmat = Breducedtimedifference, priors = rep(1/(nrow(Areduced)+1), nrow(Areduced)+1), all_comparisons = T, startdate = "2020-01-01")

fixedcolnames <- c(colnames(Bposteriorall), "Nosocomial")
Bposterioralln <- cbind(Bposteriorall, rowSums(Bposteriorall[,1:(ncol(Bposteriorall)-1)], na.rm = T))
Bposteriornodatesn <- cbind(Bposteriornodates, rowSums(Bposteriornodates[,1:(ncol(Bposteriornodates)-1)], na.rm = T))
Bposteriornodatesnogeneticsn <- cbind(Bposteriornodatesnogenetics, rowSums(Bposteriornodatesnogenetics[,1:(ncol(Bposteriornodatesnogenetics)-1)], na.rm = T))
Bposteriornogeneticsn <- cbind(Bposteriornogenetics, rowSums(Bposteriornogenetics[,1:(ncol(Bposteriornogenetics)-1)], na.rm = T))

colnames(Bposterioralln) <- colnames(Bposteriornodatesn) <- colnames(Bposteriornodatesnogeneticsn) <- colnames(Bposteriornogeneticsn) <- fixedcolnames


#Bposterior[is.na(Bposterior)] <- 0
#Bposterior[Bposterior == 0] <- NA

pheatmap(Bposterioralln, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
pheatmap(Bposteriornodatesn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
pheatmap(Bposteriornodatesnogeneticsn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
pheatmap(Bposteriornogeneticsn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

png("~/Documents/HOCI_prototype/Ball.png", width = 2000, height = round(2000*419/481), res = 300)
pheatmap(Bposteriornodatesnogeneticsn-Bposterioralln, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(-0.2, 0.2, by = 0.01)))
dev.off()

png("~/Documents/HOCI_prototype/Bnodates.png", width = 2000, height = round(2000*419/481), res = 300)
pheatmap(Bposteriornodatesnogeneticsn-Bposteriornodatesn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(-0.2, 0.2, by = 0.01)))
dev.off()

png("~/Documents/HOCI_prototype/Bnogenetics.png", width = 2000, height = round(2000*419/481), res = 300)
pheatmap(Bposteriornodatesnogeneticsn-Bposteriornogeneticsn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(-0.2, 0.2, by = 0.01)))
dev.off()

png("~/Documents/HOCI_prototype/Bnodatesnogenetics.png", width = 2000, height = round(2000*419/481), res = 300)
pheatmap(Bposteriornodatesnogeneticsn, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
dev.off()

##run without X_H as no admission dates
##assume constant contact for initial run
C$detection <- as.numeric(difftime(as.Date(C$onset_date, format = "%d/%m/%Y"),as.Date("01/01/2020", format = "%d/%m/%Y")), units = "days") + 1
C$admission <- rep(NA, nrow(C))
C$sequencedate <- as.numeric(difftime(as.Date(C$sample_collection_date, format = "%d/%m/%Y"),as.Date("01/01/2020", format = "%d/%m/%Y")), units = "days") + 1
C$sequence_id[C$sequence_id == "sequence_missing"] <- NA
C$sequencedate[is.na(C$sequence_id)] <- NA

##reduce to minimal data - one sequence per person per bout, earliest taken
indexes <- c()
for (i in 1:length(unique(C$patient_study_id))) {
  temp <- C[C$patient_study_id %in% unique(C$patient_study_id)[i],]
  if (nrow(temp) == 1) {
    indexes <- c(indexes, which(C$patient_study_id == unique(C$patient_study_id)[i]))
  } else (
    indexes <- c(indexes, which(C$patient_study_id == unique(C$patient_study_id)[i] &
                                  as.Date(C$sample_collection_date, format = "%d/%m/%Y") == min(as.Date(temp$sample_collection_date[!is.na(temp$sequence_id)], format = "%d/%m/%Y"))))
  )
}
Creduced <- C[indexes,]
Creducedtimedifference <-timedifmat(Creduced$sample_collection_date)
colnames(Creducedtimedifference) <- row.names(Creducedtimedifference) <- Creduced$patient_study_id

Cposterior <- run_model(epidata = Creduced, locationdata = positions, snpdistmat = SNPs, timedistmat = Creducedtimedifference, priors = rep(1/(nrow(Creduced)+1), nrow(Creduced)+1), all_comparisons = T, startdate = "2020-01-01")
#Cposterior[is.na(Cposterior)] <- 0
#Cposterior[Cposterior == 0] <- NA

pheatmap(Cposterior, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

##run without X_H as no admission dates
##assume constant contact for initial run
D$detection <- as.numeric(difftime(as.Date(D$onset_date, format = "%d/%m/%Y"),as.Date("01/01/2020", format = "%d/%m/%Y")), units = "days") + 1
D$admission <- rep(NA, nrow(D))
D$sequencedate <- as.numeric(difftime(as.Date(D$sample_collection_date, format = "%d/%m/%Y"),as.Date("01/01/2020", format = "%d/%m/%Y")), units = "days") + 1
D$sequence_id[D$sequence_id == "sequence_missing"] <- NA
D$sequencedate[is.na(D$sequence_id)] <- NA

indexes <- c()
for (i in 1:length(unique(D$patient_study_id))) {
  temp <- D[D$patient_study_id %in% unique(D$patient_study_id)[i],]
  if (nrow(temp) == 1) {
    indexes <- c(indexes, which(D$patient_study_id == unique(D$patient_study_id)[i]))
  } else (
    indexes <- c(indexes, which(D$patient_study_id == unique(D$patient_study_id)[i] &
                                  as.Date(D$sample_collection_date, format = "%d/%m/%Y") == min(as.Date(temp$sample_collection_date[!is.na(temp$sequence_id)], format = "%d/%m/%Y"))))
  )
}
Dreduced <- D[indexes,]
Dreducedtimedifference <-timedifmat(Dreduced$sample_collection_date)
colnames(Dreducedtimedifference) <- row.names(Dreducedtimedifference) <- Dreduced$patient_study_id

Dposterior <- run_model(epidata = Dreduced, locationdata = positions, snpdistmat = SNPs, timedistmat = Dreducedtimedifference, priors = rep(1/(nrow(Dreduced)+1), nrow(Dreduced)+1), all_comparisons = T, startdate = "2020-01-01" )
#Dposterior[is.na(Dposterior)] <- 0
#Dposterior[Dposterior == 0] <- NA

pheatmap(Dposterior, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

##run without X_H as no admission dates
##assume constant contact for initial run
E$detection <- as.numeric(difftime(as.Date(E$onset_date, format = "%d/%m/%Y"),as.Date("01/01/2020", format = "%d/%m/%Y")), units = "days") + 1
E$admission <- rep(NA, nrow(E))
E$sequencedate <- as.numeric(difftime(as.Date(E$sample_collection_date, format = "%d/%m/%Y"),as.Date("01/01/2020", format = "%d/%m/%Y")), units = "days") + 1
E$sequence_id[E$sequence_id == "sequence_missing"] <- NA
E$sequencedate[is.na(E$sequence_id)] <- NA

indexes <- c()
for (i in 1:length(unique(E$patient_study_id))) {
  temp <- E[E$patient_study_id %in% unique(E$patient_study_id)[i],]
  if (nrow(temp) == 1) {
    indexes <- c(indexes, which(E$patient_study_id == unique(E$patient_study_id)[i]))
  } else (
    indexes <- c(indexes, which(E$patient_study_id == unique(E$patient_study_id)[i] &
                                  as.Date(E$sample_collection_date, format = "%d/%m/%Y") == min(as.Date(temp$sample_collection_date[!is.na(temp$sequence_id)], format = "%d/%m/%Y"))))
  )
}
Ereduced <- E[indexes,]
Ereducedtimedifference <-timedifmat(Ereduced$sample_collection_date)
colnames(Ereducedtimedifference) <- row.names(Ereducedtimedifference) <- Ereduced$patient_study_id

# Ereducedtest <- Ereduced
# Ereducedtest$sequencedate <- NA
# Ereducedtest$sequence_id <- NA

Eposterior <- run_model(epidata = Ereduced, locationdata = positions, snpdistmat = SNPs, timedistmat = Ereducedtimedifference, priors = rep(1/(nrow(Ereduced)+1), nrow(Ereduced)+1), all_comparisons = T, startdate = "2020-01-01")
#Eposterior[is.na(Eposterior)] <- 0
#Eposterior[Eposterior == 0] <- NA

pheatmap(Eposterior, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

##pure A2B test

seqs <- read.fasta("~/Downloads/Ward_data/AlignedWithConsensus.fasta")
seqs <- seqs[names(seqs) %in% c("Consensus", A$sequence_id, B$sequence_id, C$sequence_id, D$sequence_id, E$sequence_id)[!is.na(c("Consensus", A$sequence_id, B$sequence_id, C$sequence_id, D$sequence_id, E$sequence_id))]]

H <- array(0 ,dim = c(length(seqs) - 1, length(seqs) - 1, 2), dimnames = list(names(seqs)[1:(length(seqs)-1)], names(seqs)[1:(length(seqs)-1)], c("A", "B")))


##breaks for indels
for (i in 1:(length(seqs)-1)) {
  for (j in 1:(length(seqs)-1)) {
    for (k in 1:length(seqs[[1]])) { 
      ##assumes sequences have been aligned and thus are of the same length
      if (seqs[[i]][k] %in% c("c","g","a","t","-") & seqs[[j]][k] %in% c("c","g","a","t","-")) {
        if (seqs[[i]][k] != seqs[[j]][k]) {
          if (seqs[[i]][k] != seqs[["Consensus"]][k]) {
            H[i,j,"A"] <- H[i,j,"A"] + 1
          }
          if (seqs[[j]][k] != seqs[["Consensus"]][k]) {
            H[i,j,"B"] <- H[i,j,"B"] + 1
          }
        }
      }
    }
    print(paste(i, j, sep = " "))
  }
}

Areduced$snpdifference_A <- rep(NA, nrow(Areduced))
Areduced$snpdifference_B <- rep(NA, nrow(Areduced))

AposteriorA2B <- run_model(epidata = Areduced, locationdata = positions, snpdistmat = SNPs, timedistmat = Areducedtimedifference, priors = rep(1/(nrow(Areduced)+1), nrow(Areduced)+1), all_comparisons = T, startdate = "2020-01-01",  H = H, pureA2B = T)

AposteriorA2B[AposteriorA2B == 0] <- NA

pheatmap(AposteriorA2B, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

Breduced$snpdifference_A <- rep(NA, nrow(Breduced))
Breduced$snpdifference_B <- rep(NA, nrow(Breduced))

BposteriorA2B <- run_model(epidata = Breduced, locationdata = positions, snpdistmat = SNPs, timedistmat = Breducedtimedifference, priors = rep(1/(nrow(Breduced)+1), nrow(Breduced)+1), all_comparisons = T, startdate = "2020-01-01",  H = H, pureA2B = T)

BposteriorA2B[BposteriorA2B == 0] <- NA

pheatmap(BposteriorA2B, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

Creduced$snpdifference_A <- rep(NA, nrow(Creduced))
Creduced$snpdifference_B <- rep(NA, nrow(Creduced))

CposteriorA2B <- run_model(epidata = Creduced, locationdata = positions, snpdistmat = SNPs, timedistmat = Creducedtimedifference, priors = rep(1/(nrow(Creduced)+1), nrow(Creduced)+1), all_comparisons = T, startdate = "2020-01-01",  H = H, pureA2B = T)

CposteriorA2B[CposteriorA2B == 0] <- NA

pheatmap(CposteriorA2B, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

Dreduced$snpdifference_A <- rep(NA, nrow(Dreduced))
Dreduced$snpdifference_B <- rep(NA, nrow(Dreduced))

DposteriorA2B <- run_model(epidata = Dreduced, locationdata = positions, snpdistmat = SNPs, timedistmat = Dreducedtimedifference, priors = rep(1/(nrow(Dreduced)+1), nrow(Dreduced)+1), all_comparisons = T, startdate = "2020-01-01",  H = H, pureA2B = T)

DposteriorA2B[DposteriorA2B == 0] <- NA

pheatmap(DposteriorA2B, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

Ereduced$snpdifference_A <- rep(NA, nrow(Ereduced))
Ereduced$snpdifference_B <- rep(NA, nrow(Ereduced))

EposteriorA2B <- run_model(epidata = Ereduced, locationdata = positions, snpdistmat = SNPs, timedistmat = Ereducedtimedifference, priors = rep(1/(nrow(Ereduced)+1), nrow(Ereduced)+1), all_comparisons = T, startdate = "2020-01-01",  H = H, pureA2B = T)

EposteriorA2B[EposteriorA2B == 0] <- NA

pheatmap(EposteriorA2B, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

### Sheffield HOCI data

#function generate time difference matrix from sample collection dates
timedifmat <- function(x) {
  dateformat <- as.Date(x, format = "%m/%d/%Y")
  return(do.call("cbind", lapply(1:length(dateformat), function(x) {return(abs(as.numeric(difftime(dateformat[x],dateformat, units = "days"))))})))
}

#read in data
shefdata <- read.csv("~/Downloads/HOCI_Sheffield_Glasgow/shef_hoci_metadata_movements_withseq.csv")
shefdatafocus <- shefdata[shefdata$Sequence_Type == "Focus",]

shefseqs <- read.fasta("~/Downloads/HOCI_Sheffield_Glasgow/combined_shef_seq_set.fasta")
shefseqs <- shefseqs[names(shefseqs) %in% shefdatafocus$sequenceID]

#shefdatafocus <- read.csv("~/Downloads/HOCI_Sheffield/SRT_Sheffield_complete_data_10Oct2020.csv")

##generate snp differences - note that due the the symmetry between pairwise consensuses, this is also each sequences' hamming distance to its mutual consensus
SNPs <- matrix(0, nrow = length(shefseqs), ncol = length(shefseqs))

for (i in 1:length(shefseqs)) {
  for (j in i:length(shefseqs)) {
    for (k in 1:length(shefseqs[[1]])) { 
      ##assumes sequences have been aligned and thus are of the same length
      if (shefseqs[[i]][k] %in% c("c","g","a","t","-") & shefseqs[[j]][k] %in% c("c","g","a","t","-")) {
        if (shefseqs[[i]][k] != shefseqs[[j]][k]) {
          SNPs[i,j] <- SNPs[i,j] + 1
        }
      }
    }
    print(paste(i, j, sep = " "))
  }
}

SNPs[lower.tri(SNPs)] <- t(SNPs)[lower.tri(SNPs)]
colnames(SNPs) <- row.names(SNPs) <- names(shefseqs)

shefdatafocustimedifference <-timedifmat(shefdatafocus$SampleDate)
colnames(shefdatafocustimedifference) <- row.names(shefdatafocustimedifference) <- shefdatafocus$sequenceID

shefdatafocus$admission <- as.numeric(difftime(as.Date(shefdatafocus$admissionDate, format = "%m/%d/%Y"),as.Date("20/11/2019", format = "%d/%m/%Y")), units = "days") + 1
shefdatafocus$detection <- as.numeric(difftime(as.Date(shefdatafocus$SampleDate, format = "%m/%d/%Y"),as.Date("20/11/2019", format = "%d/%m/%Y")), units = "days") + 1

shefdatafocus$patient_study_id <- shefdatafocus$sequenceID

shefdatafocus$sequencedate <- format(as.Date(shefdatafocus$SampleDate, format = "%m/%d/%Y"), "%d/%m/%Y")

shefdatafocus$sequencedate <- as.numeric(difftime(as.Date(shefdatafocus$sequencedate, format = "%d/%m/%Y"),as.Date("20/11/2019", format = "%d/%m/%Y")), units = "days") + 1

shefdatafocus$sequence_id <- shefdatafocus$sequenceID

# shefdatafocus$sequencedate <- NA
# 
# shefdatafocus$sequence_id <- NA

shefdatafocus$onset_date <- format(as.Date(shefdatafocus$SampleDate, format = "%m/%d/%Y"), "%d/%m/%Y")

shefnulllocations <- matrix(NA, nrow = 3, ncol = 2)

colnames(shefnulllocations) <- c("Date", "NULL")

shefposterior <- list()

for(i in 1:length(unique(shefdatafocus$unitID))) {
  if (nrow(shefdatafocus[shefdatafocus$unitID == unique(shefdatafocus$unitID)[i],]) != 1) {
    shefposterior[[i]] <- run_model(epidata = shefdatafocus[shefdatafocus$unitID == unique(shefdatafocus$unitID)[i],], locationdata = shefnulllocations, snpdistmat = SNPs, timedistmat = shefdatafocustimedifference, priors = rep(1/(nrow(shefdatafocus[shefdatafocus$unitID == unique(shefdatafocus$unitID)[i],])+1), nrow(shefdatafocus[shefdatafocus$unitID == unique(shefdatafocus$unitID)[i],])+1), all_comparisons = T, startdate = "2019-11-20")
    pheatmap(shefposterior[[i]], cluster_rows = F, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
  }
}

#shefposterior[shefposterior == 0] <- NA

pheatmap(shefposterior[[1]], cluster_rows = F, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))


##all data
shefseqs <- read.fasta("~/Downloads/HOCI_Sheffield_Glasgow/combined_shef_seq_set.fasta")
shefseqs <- shefseqs[names(shefseqs) %in% shefdata$sequenceID]

##generate snp differences - note that due the the symmetry between pairwise consensuses, this is also each sequences' hamming distance to its mutual consensus
SNPs <- list()

for (q in 1:length(unique(shefdata$unitID))) {
  shefseqs_working <- shefseqs[names(shefseqs) %in% shefdata$sequenceID[shefdata$unitID == unique(shefdata$unitID)[q]]]
  SNPs[[q]] <- matrix(0, nrow = length(shefseqs_working), ncol = length(shefseqs_working))
  colnames(SNPs[[q]]) <- row.names(SNPs[[q]]) <- names(shefseqs_working)
  for (i in 1:length(shefseqs_working)) {
    for (j in i:length(shefseqs_working)) {
      for (k in 1:length(shefseqs_working[[1]])) { 
        ##assumes sequences have been aligned and thus are of the same length
        if (shefseqs_working[[i]][k] %in% c("c","g","a","t","-") & shefseqs_working[[j]][k] %in% c("c","g","a","t","-")) {
          if (shefseqs_working[[i]][k] != shefseqs_working[[j]][k]) {
            SNPs[[q]][i,j] <- SNPs[[q]][i,j] + 1
          }
        }
      }
      print(paste(i, j, sep = " "))
    }
  }
  SNPs[[q]][lower.tri(SNPs[[q]])] <- t(SNPs[[q]])[lower.tri(SNPs[[q]])]
  print(q)
}

shefdatatimedifference <-timedifmat(shefdata$SampleDate)
colnames(shefdatatimedifference) <- row.names(shefdatatimedifference) <- shefdata$sequenceID

shefdata$admission <- as.numeric(difftime(as.Date(shefdata$admissionDate, format = "%m/%d/%Y"),as.Date("20/11/2019", format = "%d/%m/%Y")), units = "days") + 1
shefdata$detection <- as.numeric(difftime(as.Date(shefdata$SampleDate, format = "%m/%d/%Y"),as.Date("20/11/2019", format = "%d/%m/%Y")), units = "days") + 1

shefdata$patient_study_id <- shefdata$sequenceID

shefdata$sequencedate <- format(as.Date(shefdata$SampleDate, format = "%m/%d/%Y"), "%d/%m/%Y")

shefdata$sequencedate <- as.numeric(difftime(as.Date(shefdata$sequencedate, format = "%d/%m/%Y"),as.Date("20/11/2019", format = "%d/%m/%Y")), units = "days") + 1

shefdata$sequence_id <- shefdata$sequenceID

# shefdata$sequencedate <- NA
# 
# shefdata$sequence_id <- NA

shefdata$onset_date <- format(as.Date(shefdata$SampleDate, format = "%m/%d/%Y"), "%d/%m/%Y")

shefnulllocations <- matrix(NA, nrow = 3, ncol = 2)

colnames(shefnulllocations) <- c("Date", "NULL")

shefposterior <- list()

for(i in 1:length(unique(shefdata$unitID))) {
  if (nrow(shefdata[shefdata$unitID == unique(shefdata$unitID)[i],]) != 1) {
    shefposterior[[i]] <- run_model(epidata = shefdata[shefdata$unitID == unique(shefdata$unitID)[i],], locationdata = shefnulllocations, snpdistmat = SNPs[[i]], timedistmat = shefdatatimedifference, priors = rep(1/(nrow(shefdata[shefdata$unitID == unique(shefdata$unitID)[i],])+1), nrow(shefdata[shefdata$unitID == unique(shefdata$unitID)[i],])+1), all_comparisons = T, startdate = "2019-11-20")
    pheatmap(shefposterior[[i]], cluster_rows = F, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))
  }
  print(i)
}

#shefposterior[shefposterior == 0] <- NA

pheatmap(shefposterior[[1]], cluster_rows = F, cluster_cols = FALSE, breaks = c(seq(0, 1, by = 0.01)))

##runtime tests

runtimes <- data.frame(cbind(rep(c("SONIIC", "A2B", "A2B-Network", "HOCI", "Outbreaker2", "TransPhylo", "TiTUS", "SCOTTI"), times = 3), rep(c("A2B", "Sheffield", "Glasgow"), each = 8), rep(NA, 24), rep(NA, 24)))
colnames(runtimes) <- c("Model", "Dataset", "Start", "End")
