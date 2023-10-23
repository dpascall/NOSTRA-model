library(Delaporte)
library(matrixStats)
library(extraDistr)

#under uniform infection time probability
x_h_likelihood <- function(hypothesis, admission, detection, MU = 1.434, SIGMA = 0.6612) {
  if (hypothesis == "C") {
    return(
      log((1/detection)*
        (plnorm(detection, meanlog = MU, sdlog = SIGMA)-plnorm((detection-admission), meanlog = MU, sdlog = SIGMA)))
    )
  } else {
    return(
      log((1/detection)*
        (plnorm(detection-admission, meanlog = MU, sdlog = SIGMA)))
    )
  }
}

#under uniform infection time probability
x_ai_likelihood_not_target <- function(detection, C, MU = 1.434, SIGMA = 0.6612) {
    return(
      log((1/detection)*(plnorm(detection, meanlog = MU, sdlog = SIGMA))) + log(0.5^C)
    )
}

#under uniform infection time probability
x_ai_likelihood_onset <- function(detection, MU = 1.434, SIGMA = 0.6612) {
  return(
    log((1/detection)*(plnorm(detection, meanlog = MU, sdlog = SIGMA)))
  )
}

#genetics under non-infection - uses empirical A2B error
x_ai_genetic_likelihood <- function(Neff = floor(0.5*(exp(2.5)+exp(4.5))), mutationrate = 6.677*10^(-4)/365, errorrate = 0.404, snpdifference, timedifference, alignmentlength, generationtime = 5.5, genomesize = 29811) {
  return(ddelap(snpdifference, alpha = 2*mutationrate*generationtime*Neff*alignmentlength, beta = 1, lambda = errorrate + mutationrate*timedifference*alignmentlength/genomesize, log = T))
}

#modified a2b likelihood - continuous and only certain locations - TOST distribution changed to scaled t per Ferretti et al
a2b_likelihood <- function(sharedlocation, symptomdate_B, symptomdate_A, sequencedate_B = NA, sequencedate_A = NA, snpdifference = NA, genomesize = 29811, mutationrate = 6.677*10^(-4)/365, E = 0.404, C, alignmentlength = NA) {
  accumulate <- rep(-Inf, C)
  if (length(sharedlocation) != 0) {
    for (i in 1:length(sharedlocation)) {
      if (!is.na(snpdifference) & !(is.na(symptomdate_A) | is.na(symptomdate_B))) {
        accumulate[sharedlocation[i]] <- log(plst((sharedlocation[i] - symptomdate_A + 1), df = 3.3454, mu = -0.0747, sigma = 1.8567) - 
                                               plst((sharedlocation[i] - symptomdate_A), df = 3.3454, mu = -0.0747, sigma = 1.8567)) + 
          log(plnorm((symptomdate_B - sharedlocation[i] + 1), meanlog = 1.434, sdlog = 0.6612) - 
                plnorm((symptomdate_B - sharedlocation[i]), meanlog = 1.434, sdlog = 0.6612)) + 
          dpois(snpdifference, lambda = E + alignmentlength*mutationrate*(abs(sequencedate_B - sharedlocation[i]) + 
                                                              abs(sequencedate_A - sharedlocation[i])), log = T) +
          log(0.5^(C-1)) +
          log(1)
      } else if (is.na(snpdifference) & !(is.na(symptomdate_A) | is.na(symptomdate_B))) {
        accumulate[sharedlocation[i]] <- log(plst((sharedlocation[i] - symptomdate_A + 1), df = 3.3454, mu = -0.0747, sigma = 1.8567) - 
                                               plst((sharedlocation[i] - symptomdate_A), df = 3.3454, mu = -0.0747, sigma = 1.8567)) + 
          log(plnorm((symptomdate_B - sharedlocation[i] + 1), meanlog = 1.434, sdlog = 0.6612) - 
                plnorm((symptomdate_B - sharedlocation[i]), meanlog = 1.434, sdlog = 0.6612)) +
          log(0.5^(C-1)) +
          log(1)
      } else {
        accumulate[sharedlocation[i]] <- log(0.5^(C-1)) +
          log(1)
      }
    }
  }
  return(logSumExp(accumulate))
}

#pure a2b likelihood
pure_a2b_likelihood <- function(sharedlocation, symptomdate_B, symptomdate_A, sequencedate_B = NA, sequencedate_A = NA, snpdifference_A = NA, snpdifference_B = NA, genomesize = 29811, mutationrate = 6.677*10^(-4)/365, E = 0.404, C) {
  accumulate <- rep(-Inf, C)
  if (length(sharedlocation) != 0) {
    for (i in 1:length(sharedlocation)) {
      if ((!is.na(snpdifference_A) | !is.na(snpdifference_B)) & (!(is.na(symptomdate_A) | is.na(symptomdate_B))) & (!(is.na(sequencedate_A) | is.na(sequencedate_B)))) {
        accumulate[sharedlocation[i]] <- log(pgamma((sharedlocation[i] - symptomdate_A + 1 + 25.625), shape = 97.1875, scale = 0.2689) - 
                                               pgamma((sharedlocation[i] - symptomdate_A + 25.625), shape = 97.1875, scale = 0.2689)) + 
          log(plnorm((symptomdate_B - sharedlocation[i] + 1), meanlog = 1.434, sdlog = 0.6612) - 
                plnorm((symptomdate_B - sharedlocation[i]), meanlog = 1.434, sdlog = 0.6612)) + 
          dpois(snpdifference_A, lambda = (E/2 + mutationrate*genomesize*max(0, sequencedate_A - sharedlocation[i])), log = T) + 
          dpois(snpdifference_B, lambda = (E/2 + mutationrate*genomesize*(sequencedate_B - min(sequencedate_A, sharedlocation[i]))), log = T) + 
          log(0.5^(C-1)) +
          log(1)
      } else if (((is.na(snpdifference_A) | is.na(snpdifference_B)) | (is.na(sequencedate_A) | is.na(sequencedate_B))) & !(is.na(symptomdate_A) | is.na(symptomdate_B))) {
        accumulate[sharedlocation[i]] <- log(pgamma((sharedlocation[i] - symptomdate_A + 1 + 25.625), shape = 97.1875, scale = 0.2689) - 
                                               pgamma((sharedlocation[i] - symptomdate_A + 25.625), shape = 97.1875, scale = 0.2689)) + 
          log(plnorm((symptomdate_B - sharedlocation[i] + 1), meanlog = 1.434, sdlog = 0.6612) - 
                plnorm((symptomdate_B - sharedlocation[i]), meanlog = 1.434, sdlog = 0.6612)) + 
          log(0.5^(C-1)) +
          log(1)
      } else {
        accumulate[sharedlocation[i]] <- log(0.5^(C-1)) +
          log(1)
      }
    }
  }
  return(logSumExp(accumulate))
}

likelihood_model <- function(epidata, locationdata, priors, number_of_hypotheses, hypothesis, startdate, pureA2B) {
  function_vector <- rep(1, number_of_hypotheses - 2)
  if (is.numeric(hypothesis)) {
    function_vector[hypothesis] <- 2
  } else {
    function_vector <- c(function_vector, 3)
  }
  l_prime <- rep(-Inf, length(function_vector))
  for (q in 1:length(function_vector)) {
    if (function_vector[q] == 1) {
      if (!is.na(epidata$snpdifference[q + 1]) & !is.na(epidata$detection[q + 1])) {
        l_prime[q] <- x_ai_likelihood_not_target(detection = epidata$detection[q + 1], 
                                                 C = as.numeric(difftime(as.Date(epidata$onset_date[1], format = "%d/%m/%Y"),
                                                                         as.Date(startdate, format = "%Y-%m-%d"),
                                                                         units = "days")) + 1) +
          x_ai_genetic_likelihood(snpdifference = epidata$snpdifference[q + 1], timedifference = epidata$timediff[q + 1], alignmentlength = epidata$alignmentlength[q + 1])
      } else if (is.na(epidata$snpdifference[q + 1]) & !is.na(epidata$detection[q + 1])) {
        l_prime[q] <- x_ai_likelihood_not_target(detection = epidata$detection[q + 1], 
                                                 C = as.numeric(difftime(as.Date(epidata$onset_date[1], format = "%d/%m/%Y"),
                                                                         as.Date(startdate, format = "%Y-%m-%d"),
                                                                         units = "days")) + 1) + log(1)
      } else if (!is.na(epidata$snpdifference[q + 1]) & is.na(epidata$detection[q + 1])) {
        l_prime[q] <- x_ai_genetic_likelihood(snpdifference = epidata$snpdifference[q + 1], timedifference = epidata$timediff[q + 1], alignmentlength = epidata$alignmentlength[q + 1]) + 
          log(0.5^(as.numeric(difftime(as.Date(startdate, format = "%Y-%m-%d"), 
                                       as.Date(epidata$onset_date[1], format = "%d/%m/%Y"), 
                                       units = "days")) + 1)) +
          log(1)
      } else {
        l_prime[q] <- log(1)
      }
    } else if (function_vector[q] == 2) {
      if (pureA2B == F) {
        if ((epidata$patient_study_id[q + 1] %in% colnames(locationdata)) & (epidata$patient_study_id[1] %in% colnames(locationdata))) {
          l_prime[q] <- a2b_likelihood(symptomdate_B = epidata$detection[1], symptomdate_A = epidata$detection[q + 1], 
                                       sequencedate_B = epidata$sequencedate[1], sequencedate_A = epidata$sequencedate[q + 1],
                                       snpdifference = epidata$snpdifference[q + 1],
                                       alignmentlength = epidata$alignmentlength[q + 1],
                                       sharedlocation = generate_location(targetA = epidata$patient_study_id[q + 1], targetB = epidata$patient_study_id[1],
                                                                          startdate = startdate, enddate = epidata$onset_date[1], locationmatrix = locationdata), 
                                       C = as.numeric(difftime(as.Date(epidata$onset_date[1], format = "%d/%m/%Y"),
                                                               as.Date(startdate, format = "%Y-%m-%d"), units = "days")) + 1) +
            x_ai_likelihood_onset(detection = epidata$detection[q + 1])
        } else if (!is.na(epidata$admission[1]) & !is.na(epidata$admission[q + 1])) {
          if (epidata$admission[q + 1] <= epidata$detection[1]) {
            l_prime[q] <- a2b_likelihood(symptomdate_B = epidata$detection[1], symptomdate_A = epidata$detection[q + 1], 
                                         sequencedate_B = epidata$sequencedate[1], sequencedate_A = epidata$sequencedate[q + 1],
                                         snpdifference = epidata$snpdifference[q + 1],
                                         alignmentlength = epidata$alignmentlength[q + 1],
                                         sharedlocation = max(epidata$admission[1], epidata$admission[q + 1]):epidata$detection[1], 
                                         C = as.numeric(difftime(as.Date(epidata$onset_date[1], format = "%d/%m/%Y"),
                                                                 as.Date(startdate, format = "%Y-%m-%d"), units = "days")) + 1) +
              x_ai_likelihood_onset(detection = epidata$detection[q + 1])
          }
        } else {
          l_prime[q] <- a2b_likelihood(symptomdate_B = epidata$detection[1], symptomdate_A = epidata$detection[q + 1], 
                                       sequencedate_B = epidata$sequencedate[1], sequencedate_A = epidata$sequencedate[q + 1],
                                       snpdifference = epidata$snpdifference[q + 1],
                                       alignmentlength = epidata$alignmentlength[q + 1],
                                       sharedlocation = 1:epidata$detection[1], 
                                       C = as.numeric(difftime(as.Date(epidata$onset_date[1], format = "%d/%m/%Y"),
                                                               as.Date(startdate, format = "%Y-%m-%d"), units = "days")) + 1)  +
            x_ai_likelihood_onset(detection = epidata$detection[q + 1])
        }
      } else {
        if ((epidata$patient_study_id[q + 1] %in% colnames(locationdata)) & (epidata$patient_study_id[1] %in% colnames(locationdata))) {
          l_prime[q] <- pure_a2b_likelihood(symptomdate_B = epidata$detection[1], symptomdate_A = epidata$detection[q + 1], 
                                            sequencedate_B = epidata$sequencedate[1], sequencedate_A = epidata$sequencedate[q + 1],
                                            snpdifference_A = epidata$snpdifference_A[q + 1], snpdifference_B = epidata$snpdifference_B[q + 1],
                                            alignmentlength = epidata$alignmentlength[q + 1],
                                            sharedlocation = generate_location(targetA = epidata$patient_study_id[q + 1], targetB = epidata$patient_study_id[1],
                                                                               startdate = startdate, enddate = epidata$onset_date[1], locationmatrix = locationdata), 
                                            C = as.numeric(difftime(as.Date(epidata$onset_date[1], format = "%d/%m/%Y"),
                                                                    as.Date(startdate, format = "%Y-%m-%d"), units = "days")) + 1) +
            x_ai_likelihood_onset(detection = epidata$detection[q + 1])
        } else if (!is.na(epidata$admission[1]) & !is.na(epidata$admission[q + 1])) {
          if (epidata$admission[q + 1] <= epidata$detection[1]) {
            l_prime[q] <- pure_a2b_likelihood(symptomdate_B = epidata$detection[1], symptomdate_A = epidata$detection[q + 1], 
                                              sequencedate_B = epidata$sequencedate[1], sequencedate_A = epidata$sequencedate[q + 1],
                                              snpdifference_A = epidata$snpdifference_A[q + 1], snpdifference_B = epidata$snpdifference_B[q + 1],
                                              alignmentlength = epidata$alignmentlength[q + 1],
                                              sharedlocation = max(epidata$admission[1], epidata$admission[q + 1]):epidata$detection[1], 
                                              C = as.numeric(difftime(as.Date(epidata$onset_date[1], format = "%d/%m/%Y"),
                                                                      as.Date(startdate, format = "%Y-%m-%d"), units = "days")) + 1) +
              x_ai_likelihood_onset(detection = epidata$detection[q + 1])
          }
        } else {
          l_prime[q] <- pure_a2b_likelihood(symptomdate_B = epidata$detection[1], symptomdate_A = epidata$detection[q + 1], 
                                            sequencedate_B = epidata$sequencedate[1], sequencedate_A = epidata$sequencedate[q + 1],
                                            snpdifference_A = epidata$snpdifference_A[q + 1], snpdifference_B = epidata$snpdifference_B[q + 1],
                                            alignmentlength = epidata$alignmentlength[q + 1],
                                            sharedlocation = 1:epidata$detection[1], 
                                            C = as.numeric(difftime(as.Date(epidata$onset_date[1], format = "%d/%m/%Y"),
                                                                    as.Date(startdate, format = "%Y-%m-%d"), units = "days")) + 1) +
            x_ai_likelihood_onset(detection = epidata$detection[q + 1])
        }
      }
    } else if (function_vector[q] == 3) {
      if (!is.na(epidata$admission[1])) {
        l_prime[q] <- x_h_likelihood(hypothesis = hypothesis, detection = epidata$detection[1], admission = epidata$admission[1])
      } else {
        l_prime[q] <- log(1)
      }
    } else {
      
    }
  }
  return(sum(l_prime))
}

generate_location <- function(targetA, targetB, locationmatrix, startdate, enddate) {
  locationmatrixreduced <- locationmatrix[as.Date(locationmatrix$Date, format = "%Y-%m-%d") >= as.Date(startdate, format = "%Y-%m-%d"),]
  locationmatrixreduced <- locationmatrixreduced[as.Date(locationmatrixreduced$Date, format = "%Y-%m-%d") <= as.Date(enddate, format = "%d/%m/%Y"),]
  locationsA <- lapply(locationmatrixreduced[,targetA], function (x) return(strsplit(x, "/")))
  locationsB <- lapply(locationmatrixreduced[,targetB], function (x) return(strsplit(x, "/")))
  coincident <- rep(FALSE, nrow(locationmatrixreduced))
  for (o in 1:nrow(locationmatrixreduced)) {
    if ((sum(locationsA[[o]][[1]] %in% locationsB[[o]][[1]]) >= 1) & sum(locationsA[[o]][[1]] %in% NA) == 0 & sum(locationsB[[o]][[1]] %in% NA) == 0) {
      coincident[o] <- TRUE
    }
  }
  return(as.numeric(difftime(as.Date(locationmatrixreduced$Date, format = "%Y-%m-%d"), as.Date(startdate, format = "%Y-%m-%d"), units = "days")[coincident] + 1))
}

normalise <- function(x) {
  return(exp(x-logSumExp(x)))
}

run_model <- function(epidata, locationdata, snpdistmat, timedistmat, priors, startdate, alignmentlengthmat, all_comparisons = F, pureA2B = F, H = NULL) {
  number_of_hypotheses <- nrow(epidata) + 1
  hypotheses <- as.list(1:number_of_hypotheses)
  hypotheses[[number_of_hypotheses]] <- "C"
  hypotheses[[number_of_hypotheses-1]] <- "H"
  if (all_comparisons == F) {
    hypliks <- rep(NA, number_of_hypotheses)
    epidata$snpdifference <- NA
    epidata$timedifference <- NA
    for (j in 1:nrow(epidata)) {
      epidata$timedifference[j] <- timedistmat[epidata$patient_study_id[1], epidata$patient_study_id[j]]
      if(pureA2B == F) {
        if (!is.na(epidata$sequence_id[1]) & !is.na(epidata$sequence_id[j])) {
          epidata$snpdifference[j] <- snpdistmat[epidata$sequence_id[1], epidata$sequence_id[j]]
          epidata$alignmentlength[j] <- alignmentlengthmat[epidata$sequence_id[1], epidata$sequence_id[j]]
        }
      } else if (pureA2B == T) {
        if (!is.na(epidata$sequence_id[1]) & !is.na(epidata$sequence_id[j])) {
          epidata$snpdifference_A[j] <- H[epidata$sequence_id[1], epidata$sequence_id[j],"A"]
          epidata$snpdifference_B[j] <- H[epidata$sequence_id[1], epidata$sequence_id[j],"B"]
        }
      }
    }
    for (h in 1:number_of_hypotheses) {
      hypliks[h] <- log(priors[h]) + likelihood_model(epidata = epidata, locationdata = locationdata, number_of_hypotheses = number_of_hypotheses, hypothesis = hypotheses[[h]], startdate = startdate, pureA2B = pureA2B)
    }
    posterior <- normalise(hypliks)
    return(posterior)
  } else {
    posterior <- matrix(NA, nrow = nrow(epidata), ncol = number_of_hypotheses + 1)
    row.names(posterior) <- epidata$patient_study_id
    colnames(posterior) <- c(epidata$patient_study_id, "Hospital", "Community")
    for (u in 1:nrow(epidata)) {
      hypliks <- rep(NA, number_of_hypotheses)
      workingepidata <- rbind(epidata[u,], epidata[-u,])
      names(hypliks) <- c(workingepidata$patient_study_id[2:length(workingepidata$patient_study_id)], "Hospital", "Community")
      workingepidata$snpdifference <- NA
      workingepidata$timedifference <- NA
      for (j in 1:nrow(workingepidata)) {
        workingepidata$timedifference[j] <- timedistmat[workingepidata$patient_study_id[1], workingepidata$patient_study_id[j]]
        if(pureA2B == F) {
          if (!is.na(workingepidata$sequence_id[1]) & !is.na(workingepidata$sequence_id[j])) {
            workingepidata$snpdifference[j] <- snpdistmat[workingepidata$sequence_id[1], workingepidata$sequence_id[j]]
            workingepidata$alignmentlength[j] <- alignmentlengthmat[workingepidata$sequence_id[1], workingepidata$sequence_id[j]]
          }
        } else if (pureA2B == T) {
          if (!is.na(workingepidata$sequence_id[1]) & !is.na(workingepidata$sequence_id[j])) {
            workingepidata$snpdifference_A[j] <- H[workingepidata$sequence_id[1], workingepidata$sequence_id[j],"A"]
            workingepidata$snpdifference_B[j] <- H[workingepidata$sequence_id[1], workingepidata$sequence_id[j],"B"]
          }
        }
      }
      for (h in 1:number_of_hypotheses) {
        hypliks[h] <- log(priors[h]) + likelihood_model(epidata = workingepidata, locationdata = locationdata, number_of_hypotheses = number_of_hypotheses, hypothesis = hypotheses[[h]], startdate = startdate, pureA2B = pureA2B)
      }
      hypliks <- normalise(hypliks)
      for (s in 1:length(hypliks)) {
        posterior[u, names(hypliks)[s]] <- hypliks[s]
      }
    }
    return(posterior)
  }
}
