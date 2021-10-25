# datasets are created using the 'Create_Dataset.R' script

############# Note: we do not know how freely we can share the DDSI function, created by Lars Kulik and Roger Mundry (Kulik 2015). Therefore, here we only provide a script that calculates the other two indices. If you are interested in the DDSI, please contact Lars Kulik

# load functions for indices
source("Scripts/index_calculations.R")
library(betareg)
library(parallel)
library(doParallel)
library(tidyverse)

iterations = 100 # how many random focal assignments per dataset
cores = 4 # how many cores for parallelisation

################# Now, take each of these data sets, create the different relationship indices for different amount of data, and correlated them with each other, their full values, and the expected probability
sim.results.q1 <- lapply(1:length(sim.data.sets), function(z) {
  ##### Create all the measures for the full data set ####
  
  # assigned 'true' values: underlying probabilities, full rates, full Sociality Indices
  data.set <- sim.data.sets[[z]]
  data.set$dyad.frame$dyad <- apply(cbind(data.set$dyad.frame$focal, data.set$dyad.frame$partner), 1, function(x) paste(sort(x), collapse = "_"))
  data.set$dyad.frame <- data.set$dyad.frame[order(data.set$dyad.frame$focal, data.set$dyad.frame$partner), ]
  expected.probs <- data.set$dyad.frame$expected.probs
  full.rates.beh1 <- data.set$dyad.frame$behaviour.1.ph
  full.rates.beh2 <- data.set$dyad.frame$behaviour.2.ph
  full.rates.beh3 <- data.set$dyad.frame$behaviour.3.ph
  full.rates.beh4 <- data.set$dyad.frame$behaviour.4.ph
  
  ##### Assign 'uninformative' behaviour from other dataset ####
  set.ids <- sapply(sim.data.sets, function(x) {
    length(unique(x[[1]]$focal))
  })
  set.ids <- which(set.ids == length(unique(data.set[[1]]$focal)))
  set.ids <- set.ids[set.ids != z]
  set.ids <- set.ids[sample(x = 1:length(set.ids), 1)]
  data.set.2 <- sim.data.sets[[set.ids]]
  
  data.set.2$dyad.frame$dyad <- apply(cbind(data.set.2$dyad.frame$focal, data.set.2$dyad.frame$partner), 1, function(x) paste(sort(x), collapse = "_"))
  data.set.2$dyad.frame <- data.set.2$dyad.frame[order(data.set.2$dyad.frame$focal, data.set.2$dyad.frame$partner), ]
  
  # integrate uninformative in data set
  data.set$dyad.frame$behaviour.5.ph <- data.set.2$dyad.frame$behaviour.1.ph
  data.set$dyad.frame$expected.probs.2 <- data.set.2$dyad.frame$expected.probs
  expected.probs.2 <- data.set$dyad.frame$expected.probs.2
  full.rates.beh5 <- data.set$dyad.frame$behaviour.5.ph
  dyads <- data.set$dyad.frame$dyad
  
  # also integrate interactions in list
  interactions.2 <- data.set.2$interaction.list[data.set.2$interaction.list$behaviour == "behaviour.1", ]
  interactions.2$behaviour <- "behaviour.5"
  
  data.set$interaction.list <- rbind(data.set$interaction.list, interactions.2)
  data.set$interaction.list <- data.set$interaction.list[order(data.set$interaction.list$rdate), ]
  
  ##### Create Indices #####
  # Dyadic Composite Sociality Index: Positive Behaviours, each divided by its mean and combined
  full.DSI <- calc_index(dyads = dyads,
                         behaviours = cbind(
                           full.rates.beh1,
                           full.rates.beh2,
                           full.rates.beh3
                         ),
                         type = "DSI"
  )
  
  # 'bad' indices throughout are uninformative - they have the same properties like a real behaviour but are not correlated with underlying relationship
  full.DSI.bad <- calc_index(dyads = dyads,
                             behaviours = cbind(
                               full.rates.beh2,
                               full.rates.beh3,
                               full.rates.beh5
                             ),
                             type = "DSI"
  )
  
  # Dyadic Sociality Index - Common: Positive Behaviours without food sharing, each divided by its mean and combined
  full.DSI.common <- calc_index(dyads = dyads,
                                behaviours = cbind(
                                  full.rates.beh1,
                                  full.rates.beh3
                                ),
                                type = "DSI"
  )
  
  
  # Dyadic Sociality Index - Interaction: Positive Behaviours without proximity, each divided by its mean and combined
  full.DSI.interaction <- calc_index(dyads = dyads,
                                     behaviours = cbind(
                                       full.rates.beh1,
                                       full.rates.beh2
                                     ),
                                     type = "DSI"
  )
  
  # Composite Relationship Index: Positive and Negative Behaviours, common positive combined
  full.CRI <-
    calc_index(dyads = dyads,
               behaviours = cbind(
                 full.rates.beh1,
                 full.rates.beh2,
                 full.rates.beh3,
                 full.rates.beh4
               ),
               valence = c("pos.c", "pos.r", "pos.c", "neg"),
               type = "CRI"
    )
  
  full.CRI.bad <-
    calc_index(dyads = dyads,
               behaviours = cbind(
                 full.rates.beh5,
                 full.rates.beh2,
                 full.rates.beh3,
                 full.rates.beh4
               ),
               valence = c("pos.c", "pos.r", "pos.c", "neg"),
               type = "CRI"
    )
  
  
  ##### Start parallelisation for randomized iterations #####
  
  ############ next, take all those things and select only one observer per day, different amounts of days of the year
  mycluster <- makeCluster(cores, type = "PSOCK")
  # export the relevant information to each core
  clusterExport(
    cl = mycluster,
    c(
      "data.set",
      "expected.probs",
      "full.CRI",
      "full.DSI",
      "full.DSI.common",
      "full.CRI.bad",
      "full.DSI.bad",
      "full.rates.beh1",
      "full.rates.beh2",
      "full.rates.beh3",
      "full.rates.beh4",
      "full.rates.beh5",
      "sim.pos.neg",
      'calc_index',
      'calc_DSI',
      'calc_CRI'
    ),
    envir = environment()
  )
  registerDoParallel(mycluster)
  # run parallel loop
  
  ##### randomly select focal days, create indices, and check their abilities
  ran.measures <- parLapply(cl = mycluster, X = 1:iterations, function(y) {
    # random dataset
    ran.data.set <- data.set
    ran.number.days <- sample(30:length(unique(data.set$interaction.list$rdate)), 1)
    
    # set likelihood of individuals for unbalanced datasets: each individual has either a 0.8, 0.5, or 0.2 chance of being focaled
    ran.llh <- rep(c(0.8, 0.5, 0.2), 1000)
    ran.llh <- sample(ran.llh, size = length(unique(data.set$interaction.list$focal)))
    
    # select focal per day of the year
    ran.focal.day <- data.frame(rdate = unique(as.Date(data.set$interaction.list$rdate)))
    ran.focal.day$focal <- sapply(ran.focal.day$rdate, function(x) sample(unique(data.set[[1]]$focal), 1))
    ran.focal.day$focal.select <- sapply(ran.focal.day$rdate, function(x) sample(unique(data.set[[1]]$focal), 1, prob = ran.llh))    
    ran.focal.day <- ran.focal.day[ran.focal.day$rdate %in% sample(ran.focal.day$rdate, ran.number.days), ]
    
    ## only select those days where the focal or partner are those focals
    ran.list.data <- ran.data.set$interaction.list
    
    # remove days that are not in
    ran.list.data <- ran.list.data[ran.list.data$rdate %in% ran.focal.day$rdate, ]
    
    # order by date
    ran.list.data <- ran.list.data[order(ran.list.data$rdate), ]
    
    # create a second dataset that will serve as the dataset with unbalanced data
    ran.list.data.select <- ran.list.data
    ran.list.data.select <- ran.list.data.select[!(ran.list.data.select$behaviour %in% c("behaviour.5")), ] # remove behaviour 5 from unbalanced dataset
    
    # remove non-focals
    # selecting focals by "focal" column - balanced sampling
    keep.list <- unlist(lapply(1:nrow(ran.focal.day), function(x) {
      xx.set <- ran.list.data[ran.list.data$rdate == ran.focal.day$rdate[x], ]
      return(as.numeric(xx.set$focal == ran.focal.day$focal[x] | xx.set$partner == ran.focal.day$focal[x]))
    }))
    ran.list.data <- ran.list.data[keep.list == 1, ]
    
    # selecting focals by "focal.select" column - unbalanced sampling
    keep.list.select <- unlist(lapply(1:nrow(ran.focal.day), function(x) {
      xx.set <- ran.list.data.select[ran.list.data.select$rdate == ran.focal.day$rdate[x], ]
      return(as.numeric(xx.set$focal == ran.focal.day$focal.select[x] | xx.set$partner == ran.focal.day$focal.select[x]))
    }))
    ran.list.data.select <- ran.list.data.select[keep.list.select == 1, ]
    
    ran.list.data$dyad <- apply(cbind(ran.list.data$focal, ran.list.data$partner), 1, function(x) paste(sort(x), collapse = "_"))
    ran.list.data.select$dyad <- apply(cbind(ran.list.data.select$focal, ran.list.data.select$partner), 1, function(x) paste(sort(x), collapse = "_"))
    
    
    ran.focal.day$observation.time <- 12
    
    ran.beh.rates <- data.frame(expand.grid(unique(c(data.set$dyad.frame$focal, data.set$dyad.frame$partner)), unique(c(data.set$dyad.frame$focal, data.set$dyad.frame$partner))))
    colnames(ran.beh.rates) <- c("focal", "partner")
    ran.beh.rates <- ran.beh.rates[ran.beh.rates$focal != ran.beh.rates$partner, ]
    ran.beh.rates$behaviour.1 <- unlist(lapply(1:nrow(ran.beh.rates), function(x) {
      return(sum(ran.list.data$duration[ran.list.data$behaviour == "behaviour.1" &
                                                      ((ran.list.data$focal == ran.beh.rates$focal[x] &
                                                          ran.list.data$partner == ran.beh.rates$partner[x]) |
                                                         (ran.list.data$partner == ran.beh.rates$focal[x] &
                                                            ran.list.data$focal == ran.beh.rates$partner[x])
                                                      )]))
    }))
    ran.beh.rates$behaviour.1[is.na(ran.beh.rates$behaviour.1)] <- 0
    
    ran.beh.rates$behaviour.2 <- unlist(lapply(1:nrow(ran.beh.rates), function(x) {
      return(sum(ran.list.data$duration[ran.list.data$behaviour == "behaviour.2" &
                                                      ((ran.list.data$focal == ran.beh.rates$focal[x] &
                                                          ran.list.data$partner == ran.beh.rates$partner[x]) |
                                                         (ran.list.data$partner == ran.beh.rates$focal[x] &
                                                            ran.list.data$focal == ran.beh.rates$partner[x])
                                                      )]))
    }))
    ran.beh.rates$behaviour.2[is.na(ran.beh.rates$behaviour.2)] <- 0
    
    ran.beh.rates$behaviour.3 <- unlist(lapply(1:nrow(ran.beh.rates), function(x) {
      return(sum(ran.list.data$duration[ran.list.data$behaviour == "behaviour.3" &
                                                      ((ran.list.data$focal == ran.beh.rates$focal[x] &
                                                          ran.list.data$partner == ran.beh.rates$partner[x]) |
                                                         (ran.list.data$partner == ran.beh.rates$focal[x] &
                                                            ran.list.data$focal == ran.beh.rates$partner[x])
                                                      )]))
    }))
    ran.beh.rates$behaviour.3[is.na(ran.beh.rates$behaviour.3)] <- 0
    
    ran.beh.rates$behaviour.4 <- unlist(lapply(1:nrow(ran.beh.rates), function(x) {
      return(sum(ran.list.data$duration[ran.list.data$behaviour == "behaviour.4" &
                                                      ((ran.list.data$focal == ran.beh.rates$focal[x] &
                                                          ran.list.data$partner == ran.beh.rates$partner[x]) |
                                                         (ran.list.data$partner == ran.beh.rates$focal[x] &
                                                            ran.list.data$focal == ran.beh.rates$partner[x])
                                                      )]))
    }))
    ran.beh.rates$behaviour.4[is.na(ran.beh.rates$behaviour.4)] <- 0
    
    ran.beh.rates$behaviour.5 <- unlist(lapply(1:nrow(ran.beh.rates), function(x) {
      return(sum(ran.list.data$duration[ran.list.data$behaviour == "behaviour.5" &
                                                 ((ran.list.data$focal == ran.beh.rates$focal[x] &
                                                     ran.list.data$partner == ran.beh.rates$partner[x]) |
                                                    (ran.list.data$partner == ran.beh.rates$focal[x] &
                                                       ran.list.data$focal == ran.beh.rates$partner[x])
                                                 )]))
    }))
    ran.beh.rates$behaviour.5[is.na(ran.beh.rates$behaviour.5)] <- 0
    
    
    ran.beh.rates$observation.time <- unlist(lapply(1:nrow(ran.beh.rates), function(x) {
      return(sum(ran.focal.day$observation.time[ran.focal.day$focal %in% c(ran.beh.rates$focal[x], ran.beh.rates$partner[x])]))
    }))
    ran.beh.rates$dyad <- apply(cbind(ran.beh.rates$focal, ran.beh.rates$partner), 1, function(x) paste(sort(x), collapse = "_"))
    ran.beh.rates <- ran.beh.rates[order(ran.beh.rates$focal, ran.beh.rates$partner), ]
    
    
    # makes rates for each of the behaviours
    ran.rates.beh1 <- ran.beh.rates$behaviour.1 / ran.beh.rates$observation.time
    ran.rates.beh2 <- ran.beh.rates$behaviour.2 / ran.beh.rates$observation.time
    ran.rates.beh3 <- ran.beh.rates$behaviour.3 / ran.beh.rates$observation.time
    ran.rates.beh4 <- ran.beh.rates$behaviour.4 / ran.beh.rates$observation.time
    ran.rates.beh5 <- ran.beh.rates$behaviour.5 / ran.beh.rates$observation.time
    
    # they can be infinite or NA, turn those to 0
    ran.rates.beh1[is.infinite(ran.rates.beh1)] <- 0
    ran.rates.beh2[is.infinite(ran.rates.beh2)] <- 0
    ran.rates.beh3[is.infinite(ran.rates.beh3)] <- 0
    ran.rates.beh4[is.infinite(ran.rates.beh4)] <- 0
    ran.rates.beh5[is.infinite(ran.rates.beh5)] <- 0
    ran.rates.beh1[is.na(ran.rates.beh1)] <- 0
    ran.rates.beh2[is.na(ran.rates.beh2)] <- 0
    ran.rates.beh3[is.na(ran.rates.beh3)] <- 0
    ran.rates.beh4[is.na(ran.rates.beh4)] <- 0
    ran.rates.beh5[is.na(ran.rates.beh5)] <- 0
    
    dyads <- data.set$dyad.frame$dyad
    
    ## creating behaviour rates - unbalanced dataset
    ran.beh.rates.unbalanced <- data.frame(expand.grid(unique(c(data.set$dyad.frame$focal, data.set$dyad.frame$partner)), unique(c(data.set$dyad.frame$focal, data.set$dyad.frame$partner))))
    colnames(ran.beh.rates.unbalanced) <- c("focal", "partner")
    ran.beh.rates.unbalanced <- ran.beh.rates.unbalanced[ran.beh.rates.unbalanced$focal != ran.beh.rates.unbalanced$partner, ]
    ran.beh.rates.unbalanced$behaviour.1 <- unlist(lapply(1:nrow(ran.beh.rates.unbalanced), function(x) {
      return(sum(ran.list.data.select$duration[ran.list.data.select$behaviour == "behaviour.1" &
                                                 ((ran.list.data.select$focal == ran.beh.rates.unbalanced$focal[x] &
                                                     ran.list.data.select$partner == ran.beh.rates.unbalanced$partner[x]) |
                                                    (ran.list.data.select$partner == ran.beh.rates.unbalanced$focal[x] &
                                                       ran.list.data.select$focal == ran.beh.rates.unbalanced$partner[x])
                                                 )]))
    }))
    ran.beh.rates.unbalanced$behaviour.1[is.na(ran.beh.rates.unbalanced$behaviour.1)] <- 0
    
    
    ran.beh.rates.unbalanced$behaviour.2 <- unlist(lapply(1:nrow(ran.beh.rates.unbalanced), function(x) {
      return(sum(ran.list.data.select$duration[ran.list.data.select$behaviour == "behaviour.2" &
                                                 ((ran.list.data.select$focal == ran.beh.rates.unbalanced$focal[x] &
                                                     ran.list.data.select$partner == ran.beh.rates.unbalanced$partner[x]) |
                                                    (ran.list.data.select$partner == ran.beh.rates.unbalanced$focal[x] &
                                                       ran.list.data.select$focal == ran.beh.rates.unbalanced$partner[x])
                                                 )]))
    }))
    ran.beh.rates.unbalanced$behaviour.2[is.na(ran.beh.rates.unbalanced$behaviour.2)] <- 0
    
    ran.beh.rates.unbalanced$behaviour.3 <- unlist(lapply(1:nrow(ran.beh.rates.unbalanced), function(x) {
      return(sum(ran.list.data.select$duration[ran.list.data.select$behaviour == "behaviour.3" &
                                                 ((ran.list.data.select$focal == ran.beh.rates.unbalanced$focal[x] &
                                                     ran.list.data.select$partner == ran.beh.rates.unbalanced$partner[x]) |
                                                    (ran.list.data.select$partner == ran.beh.rates.unbalanced$focal[x] &
                                                       ran.list.data.select$focal == ran.beh.rates.unbalanced$partner[x])
                                                 )]))
    }))
    ran.beh.rates.unbalanced$behaviour.3[is.na(ran.beh.rates.unbalanced$behaviour.3)] <- 0
    
    ran.beh.rates.unbalanced$behaviour.4 <- unlist(lapply(1:nrow(ran.beh.rates.unbalanced), function(x) {
      return(sum(ran.list.data.select$duration[ran.list.data.select$behaviour == "behaviour.4" &
                                                 ((ran.list.data.select$focal == ran.beh.rates.unbalanced$focal[x] &
                                                     ran.list.data.select$partner == ran.beh.rates.unbalanced$partner[x]) |
                                                    (ran.list.data.select$partner == ran.beh.rates.unbalanced$focal[x] &
                                                       ran.list.data.select$focal == ran.beh.rates.unbalanced$partner[x])
                                                 )]))
    }))
    ran.beh.rates.unbalanced$behaviour.4[is.na(ran.beh.rates.unbalanced$behaviour.4)] <- 0
    
    ### repeating the same with the unbalanced data
    ran.beh.rates.unbalanced$observation.time <- unlist(lapply(1:nrow(ran.beh.rates.unbalanced), function(x) {
      return(sum(ran.focal.day$observation.time[ran.focal.day$focal.select %in% c(ran.beh.rates.unbalanced$focal[x], ran.beh.rates.unbalanced$partner[x])]))
    }))
    ran.beh.rates.unbalanced$dyad <- apply(cbind(ran.beh.rates.unbalanced$focal, ran.beh.rates.unbalanced$partner), 1, function(x) paste(sort(x), collapse = "_"))
    ran.beh.rates.unbalanced <- ran.beh.rates.unbalanced[order(ran.beh.rates.unbalanced$focal, ran.beh.rates.unbalanced$partner), ]
    
    # makes rates for each of the behaviours
    ran.rates.beh1.unbalanced <- ran.beh.rates.unbalanced$behaviour.1 / ran.beh.rates.unbalanced$observation.time
    ran.rates.beh2.unbalanced <- ran.beh.rates.unbalanced$behaviour.2 / ran.beh.rates.unbalanced$observation.time
    ran.rates.beh3.unbalanced <- ran.beh.rates.unbalanced$behaviour.3 / ran.beh.rates.unbalanced$observation.time
    ran.rates.beh4.unbalanced <- ran.beh.rates.unbalanced$behaviour.4 / ran.beh.rates.unbalanced$observation.time
    
    # they can be infinite or NA, turn those to 0
    ran.rates.beh1.unbalanced[is.infinite(ran.rates.beh1.unbalanced)] <- 0
    ran.rates.beh2.unbalanced[is.infinite(ran.rates.beh2.unbalanced)] <- 0
    ran.rates.beh3.unbalanced[is.infinite(ran.rates.beh3.unbalanced)] <- 0
    ran.rates.beh4.unbalanced[is.infinite(ran.rates.beh4.unbalanced)] <- 0
    ran.rates.beh1.unbalanced[is.na(ran.rates.beh1.unbalanced)] <- 0
    ran.rates.beh2.unbalanced[is.na(ran.rates.beh2.unbalanced)] <- 0
    ran.rates.beh3.unbalanced[is.na(ran.rates.beh3.unbalanced)] <- 0
    ran.rates.beh4.unbalanced[is.na(ran.rates.beh4.unbalanced)] <- 0
    
    
    
    
    
    # Dyadic Sociality Index: Positive Behaviours, each divided by its mean and combined
    ran.DSI <- calc_index(dyads = dyads,
                          behaviours = cbind(
                            ran.rates.beh1,
                            ran.rates.beh2,
                            ran.rates.beh3
                          ),
                          type = "DSI"
    )
    ran.DSI.bad <- calc_index(dyads = dyads,
                              behaviours = cbind(
                                ran.rates.beh2,
                                ran.rates.beh3,
                                ran.rates.beh5
                              ),
                              type = "DSI"
    )
    
    # Dyadic Sociality Index - Common: Positive Behaviours without food sharing, each divided by its mean and combined
    ran.DSI.common <- calc_index(dyads = dyads,
                                 behaviours = cbind(
                                   ran.rates.beh1,
                                   ran.rates.beh3
                                 ),
                                 type = "DSI"
    )
    ran.DSI.common.bad <- calc_index(dyads = dyads,
                                     behaviours = cbind(
                                       ran.rates.beh3,
                                       ran.rates.beh5
                                     ),
                                     type = "DSI"
    )
    
    ran.DSI.unbalanced <- calc_index(dyads = dyads,
                                     behaviours = cbind(
                                       ran.rates.beh1.unbalanced,
                                       ran.rates.beh2.unbalanced,
                                       ran.rates.beh3.unbalanced
                                     ),
                                     type = "DSI"
    )
    ran.DSI.common.unbalanced <- calc_index(dyads = dyads,
                                            behaviours = cbind(
                                              ran.rates.beh1.unbalanced,
                                              ran.rates.beh3.unbalanced
                                            ),
                                            type = "DSI"
    )
    
    
    # Dyadic Sociality Index - Interaction: Positive Behaviours without proximity, each divided by its mean and combined
    ran.DSI.interaction <- calc_index(dyads = dyads,
                                      behaviours = cbind(
                                        ran.rates.beh1,
                                        ran.rates.beh2
                                      ),
                                      type = "DSI"
    )
    ran.DSI.interaction.bad <- calc_index(dyads = dyads,
                                          behaviours = cbind(
                                            ran.rates.beh5,
                                            ran.rates.beh2
                                          ),
                                          type = "DSI"
    )
    
    # Composite Relationship Index: Positive and Negative Behaviours, common positive combined
    ran.CRI <-
      calc_index(dyads = dyads,
                 behaviours = cbind(
                   ran.rates.beh1,
                   ran.rates.beh2,
                   ran.rates.beh3,
                   ran.rates.beh4
                 ),
                 valence = c("pos.c", "pos.r", "pos.c", "neg"),
                 type = "CRI"
      )
    
    ran.CRI.bad <-
      calc_index(dyads = dyads,
                 behaviours = cbind(
                   ran.rates.beh2,
                   ran.rates.beh3,
                   ran.rates.beh4,
                   ran.rates.beh5
                 ),
                 valence = c("pos.r", "pos.c", "neg", 'pos.c'),
                 type = "CRI"
      )
    
    ran.CRI.unbalanced <-
      calc_index(dyads = dyads,
                 behaviours = cbind(
                   ran.rates.beh1.unbalanced,
                   ran.rates.beh2.unbalanced,
                   ran.rates.beh3.unbalanced,
                   ran.rates.beh4.unbalanced
                 ),
                 valence = c("pos.c", "pos.r", "pos.c", "neg"),
                 type = "CRI"
      )
    
    
      
    ran.DSI <- ifelse(is.na(ran.DSI) | is.infinite(ran.DSI), 0.5, ran.DSI)
    ran.DSI.unbalanced <- ifelse(is.na(ran.DSI.unbalanced) | is.infinite(ran.DSI.unbalanced), 0.5, ran.DSI.unbalanced)
    ran.DSI.common <- ifelse(is.na(ran.DSI.common) | is.infinite(ran.DSI.common), 0.5, ran.DSI.common)
    ran.DSI.common.unbalanced <- ifelse(is.na(ran.DSI.common.unbalanced) | is.infinite(ran.DSI.common.unbalanced), 0.5, ran.DSI.common.unbalanced)
    ran.DSI.common.bad <- ifelse(is.na(ran.DSI.common.bad) | is.infinite(ran.DSI.common.bad), 0.5, ran.DSI.common.bad)
    ran.DSI.interaction <- ifelse(is.na(ran.DSI.interaction) | is.infinite(ran.DSI.interaction), 0.5, ran.DSI.interaction)
    ran.DSI.interaction.bad <- ifelse(is.na(ran.DSI.interaction.bad) | is.infinite(ran.DSI.interaction.bad), 0.5, ran.DSI.interaction.bad)
    ran.CRI <- ifelse(is.na(ran.CRI) | is.infinite(ran.CRI), 0.5, ran.CRI)
    ran.CRI.unbalanced <- ifelse(is.na(ran.CRI.unbalanced) | is.infinite(ran.CRI.unbalanced), 0.5, ran.CRI.unbalanced)
    ran.CRI.bad <- ifelse(is.na(ran.CRI.bad) | is.infinite(ran.CRI.bad), 0.5, ran.CRI.bad)
    ran.rates.beh1 <- ifelse(is.na(ran.rates.beh1) | is.infinite(ran.rates.beh1), 0.05, ran.rates.beh1)
    ran.rates.beh2 <- ifelse(is.na(ran.rates.beh2) | is.infinite(ran.rates.beh2), 0.05, ran.rates.beh2)
    ran.rates.beh3 <- ifelse(is.na(ran.rates.beh3) | is.infinite(ran.rates.beh3), 0.05, ran.rates.beh3)
    ran.rates.beh4 <- ifelse(is.na(ran.rates.beh4) | is.infinite(ran.rates.beh4), 0.05, ran.rates.beh4)
    ran.rates.beh5 <- ifelse(is.na(ran.rates.beh5) | is.infinite(ran.rates.beh5), 0.05, ran.rates.beh5)
    
    # Do a simple lm with each predictor, test which one predict the distribution best
    expected.probs1 <- expected.probs + 0.0001
    lm.DSI <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.DSI)
    lm.DSI.unbalanced <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.DSI.unbalanced)
    lm.DSI.common <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.DSI.common)
    lm.DSI.interaction <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.DSI.interaction)
    lm.DSI.bad <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.DSI.bad)
    lm.DSI.common.bad <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.DSI.common.bad)
    lm.CRI <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.CRI)
    lm.CRI.unbalanced <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.CRI.unbalanced)
    lm.CRI.bad <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.CRI.bad)
    lm.beh1 <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.rates.beh1)
    lm.beh2 <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.rates.beh2)
    lm.beh3 <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.rates.beh3)
    lm.beh4 <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.rates.beh4)
    lm.beh5 <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.rates.beh5)
    lm.combination <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.rates.beh1 + ran.rates.beh2 + ran.rates.beh3 + ran.rates.beh4)
    lm.combination.bad <- betareg::betareg(na.action=na.omit, expected.probs1 ~ ran.rates.beh1 + ran.rates.beh2 + ran.rates.beh3 + ran.rates.beh4 + ran.rates.beh5)
    
    AIClms <- AIC(
      lm.DSI, lm.DSI.common, lm.CRI, 
      lm.DSI.bad, lm.DSI.common.bad, lm.CRI.bad,
      lm.beh1, lm.beh2, lm.beh3, lm.beh4, lm.beh5, lm.combination, lm.combination.bad, lm.DSI.interaction, lm.CRI.unbalanced, lm.DSI.unbalanced
    )
    AIClms$DeltaAIC <- AIClms$AIC - min(AIClms$AIC)
    
    
    #### compile results in a big big frame
    results <- data.frame(
      number.days = ran.number.days,
      number.individuals = length(unique(data.set$interaction.list$focal)),
      number.interactions = nrow(ran.list.data),
      number.interactions.beh1 = sum(ran.list.data$behaviour == "behaviour.1"),
      number.interactions.beh2 = sum(ran.list.data$behaviour == "behaviour.2"),
      number.interactions.beh3 = sum(ran.list.data$behaviour == "behaviour.3"),
      number.interactions.beh4 = sum(ran.list.data$behaviour == "behaviour.4"),
      interactions.per.dyad = nrow(ran.list.data) /
        (((length(unique(data.set$interaction.list$focal))^2) - length(unique(data.set$interaction.list$focal))) / 2),
      cor.exp.beh1 = cor(expected.probs, ran.rates.beh1, use = "complete.obs"),
      cor.exp.beh2 = cor(expected.probs, ran.rates.beh2, use = "complete.obs"),
      cor.exp.beh3 = cor(expected.probs, ran.rates.beh3, use = "complete.obs"),
      cor.exp.beh4 = cor(expected.probs, ran.rates.beh4, use = "complete.obs"),
      cor.exp.beh5 = cor(expected.probs, ran.rates.beh5, use = "complete.obs"),
      cor.full.beh1 = cor(full.rates.beh1, ran.rates.beh1, use = "complete.obs"),
      cor.full.beh2 = cor(full.rates.beh2, ran.rates.beh2, use = "complete.obs"),
      cor.full.beh3 = cor(full.rates.beh3, ran.rates.beh3, use = "complete.obs"),
      cor.full.beh4 = cor(full.rates.beh4, ran.rates.beh4, use = "complete.obs"),
      cor.full.beh4 = cor(full.rates.beh5, ran.rates.beh5, use = "complete.obs"),
      cor.exp.DSI = cor(expected.probs, ran.DSI, use = "complete.obs"),
      cor.exp.DSI.common = cor(expected.probs, ran.DSI.common, use = "complete.obs"),
      cor.exp.DSI.interaction = cor(expected.probs, ran.DSI.interaction, use = "complete.obs"),
      cor.exp.CRI = cor(expected.probs, ran.CRI, use = "complete.obs"),
      cor.exp.DSI.unbalanced = cor(expected.probs, ran.DSI.unbalanced, use = "complete.obs"),
      cor.exp.CRI.unbalanced = cor(expected.probs, ran.CRI.unbalanced, use = "complete.obs"),
      cor.full.DSI = cor(full.DSI, ran.DSI, use = "complete.obs"),
      cor.full.DSI.common = cor(full.DSI.common, ran.DSI.common, use = "complete.obs"),
      cor.full.DSI.interaction = cor(full.DSI.interaction, ran.DSI.interaction, use = "complete.obs"),
      cor.full.CRI = cor(full.CRI, ran.CRI, use = "complete.obs"),
      cor.full.DSI.unbalanced = cor(full.DSI, ran.DSI.unbalanced, use = "complete.obs"),
      cor.full.CRI.unbalanced = cor(full.CRI, ran.CRI.unbalanced, use = "complete.obs"),
      cor.exp.DSI.bad = cor(expected.probs, ran.DSI.bad, use = "complete.obs"),
      cor.exp.DSI.common.bad = cor(expected.probs, ran.DSI.common.bad, use = "complete.obs"),
      cor.exp.CRI.bad = cor(expected.probs, ran.CRI.bad, use = "complete.obs"),
      cor.DSI.CRI = cor(ran.DSI, ran.CRI, use = "complete.obs"),
      best.model.AIC = rownames(AIClms)[AIClms$AIC == min(AIClms$AIC)],
      AIC.DSI = AIClms$DeltaAIC[rownames(AIClms) == "lm.DSI"],
      AIC.DSI.common = AIClms$DeltaAIC[rownames(AIClms) == "lm.DSI.common"],
      AIC.DSI.interaction = AIClms$DeltaAIC[rownames(AIClms) == "lm.DSI.interaction"],
      AIC.CRI = AIClms$DeltaAIC[rownames(AIClms) == "lm.CRI"],
      AIC.DSI.unbalanced = AIClms$DeltaAIC[rownames(AIClms) == "lm.DSI.unbalanced"],
      AIC.CRI.unbalanced = AIClms$DeltaAIC[rownames(AIClms) == "lm.CRI.unbalanced"],
      AIC.DSI.bad = AIClms$DeltaAIC[rownames(AIClms) == "lm.DSI.bad"],
      AIC.DSI.common.bad = AIClms$DeltaAIC[rownames(AIClms) == "lm.DSI.common.bad"],
      AIC.CRI.bad = AIClms$DeltaAIC[rownames(AIClms) == "lm.CRI.bad"],
      AIC.beh1 = AIClms$DeltaAIC[rownames(AIClms) == "lm.beh1"],
      AIC.beh2 = AIClms$DeltaAIC[rownames(AIClms) == "lm.beh2"],
      AIC.beh3 = AIClms$DeltaAIC[rownames(AIClms) == "lm.beh3"],
      AIC.beh4 = AIClms$DeltaAIC[rownames(AIClms) == "lm.beh4"],
      AIC.beh5 = AIClms$DeltaAIC[rownames(AIClms) == "lm.beh5"],
      AIC.combination = AIClms$DeltaAIC[rownames(AIClms) == "lm.combination"],
      AIC.combination.bad = AIClms$DeltaAIC[rownames(AIClms) == "lm.combination.bad"]
    )
    return(results)
  })
  stopCluster(mycluster)
  result.frame <- do.call(rbind, ran.measures)
  return(result.frame)
})
