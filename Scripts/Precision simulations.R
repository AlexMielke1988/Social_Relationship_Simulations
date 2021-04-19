source("Scripts/index_calculations.R")
library(parallel)
library(doParallel)

iterations = 100 # how many random focal assignments per dataset
cores = 4 # how many cores for parallelisation
##### this function doesn't look at the 'full' values but at the variation in random indices
##### the number of focal days is fixed to 30, 90, 180, 270, and 360 days, rather than varying randomly

sim.variation <- do.call(rbind, lapply(1:length(sim.data.sets), function(z) {
  ######## First, create all the measures for the full data set
  
  # assigned 'true' values: underlying probabilities, full rates, full Sociality Indices
  data.set <- sim.data.sets[[z]]
  data.set$dyad.frame$dyad <- apply(cbind(data.set$dyad.frame$focal, data.set$dyad.frame$partner), 1, function(x) paste(sort(x), collapse = "_"))
  data.set$dyad.frame <- data.set$dyad.frame[order(data.set$dyad.frame$focal, data.set$dyad.frame$partner), ]
  expected.probs <- data.set$dyad.frame$expected.probs
  full.rates.beh1 <- data.set$dyad.frame$behaviour.1.ph
  full.rates.beh2 <- data.set$dyad.frame$behaviour.2.ph
  full.rates.beh3 <- data.set$dyad.frame$behaviour.3.ph
  full.rates.beh4 <- data.set$dyad.frame$behaviour.4.ph
  
  # For each data set, use another one with the same sample size as data.set.2
  set.ids <- sapply(sim.data.sets, function(x) {
    length(unique(x[[1]]$focal))
  })
  set.ids <- which(set.ids == length(unique(data.set[[1]]$focal)))
  set.ids <- set.ids[set.ids != z]
  set.ids <- set.ids[sample(x = 1:length(set.ids), 1)]
  data.set.2 <- sim.data.sets[[set.ids]]
  
  data.set.2$dyad.frame$dyad <- apply(cbind(data.set.2$dyad.frame$focal, data.set.2$dyad.frame$partner), 1, function(x) paste(sort(x), collapse = "_"))
  data.set.2$dyad.frame <- data.set.2$dyad.frame[order(data.set.2$dyad.frame$focal, data.set.2$dyad.frame$partner), ]
  
  # data set 2 provides the new behaviour - this should be only weakly correlated with the others
  # integrate in the normal dataset
  data.set$dyad.frame$behaviour.5.ph <- data.set.2$dyad.frame$behaviour.1.ph
  data.set$dyad.frame$expected.probs.2 <- data.set.2$dyad.frame$expected.probs
  expected.probs.2 <- data.set$dyad.frame$expected.probs.2
  full.rates.beh5 <- data.set$dyad.frame$behaviour.5.ph
  
  # also integrate interactions in list
  interactions.21 <- data.set.2$interaction.list[data.set.2$interaction.list$behaviour == "behaviour.1", ]
  interactions.21$behaviour <- "behaviour.5"
  
  data.set$interaction.list <- rbind(data.set$interaction.list, interactions.21, interactions.22)
  data.set$interaction.list <- data.set$interaction.list[order(data.set$interaction.list$rdate), ]
  
  dyads <- data.set$dyad.frame$dyad
  
  
  # Dyadic Sociality Index: Positive Behaviours, each divided by its mean and combined
  full.DSI <- calc_index(dyads = dyads,
                         behaviours = cbind(
                           full.rates.beh1,
                           full.rates.beh2,
                           full.rates.beh3
                         ),
                         type = "DSI"
  )
  
  
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
                 full.rates.beh2,
                 full.rates.beh3,
                 full.rates.beh4,
                 full.rates.beh5
               ),
               valence = c("pos.r", "pos.c", "neg", 'pos.c'),
               type = "CRI"
    )
  
  
  # set likelihood of individuals for unbalanced datasets: each individual has either a 0.8, 0.2 chance of being focaled
  ran.llh <- rep(c(0.2, 0.8), 1000)
  ran.llh <- sample(ran.llh, size = length(unique(data.set$interaction.list$focal)))
  data.set$dyad.frame$selection1 <- ran.llh[data.set$dyad.frame$focal]
  data.set$dyad.frame$selection2 <- ran.llh[data.set$dyad.frame$partner]
  selection.llh <- apply(cbind(data.set$dyad.frame$selection1, data.set$dyad.frame$selection2), 1, function(x) paste(sort(x), collapse = "_"))
  
  ##### the number of focal days is fixed to 30, 90, 180, 270, and 360 days, rather than varying randomly
  sim.list <- lapply(c(60, 180, 360), function(subs) {
    
    ############ next, take all those things and select only one observer per day, different amounts of days of the year
    mycluster <- makeCluster(cores, type = "PSOCK")
    # export the relevant information to each core
    clusterExport(
      cl = mycluster,
      c(
        "data.set",
        "expected.probs",
        'calc_index',
        'calc_DSI',
        'calc_CRI',
        "subs",
        "sim.pos.neg",
        "selection.llh"
      ),
      envir = environment()
    )
    registerDoParallel(mycluster)
    
    # run parallel loop
    ran.measures <- parLapply(cl = mycluster, X = 1:iterations, function(y) {
      # random dataset
      ran.data.set <- data.set
      ran.number.days <- subs
      
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
      #
      # ran.list.data$dyad <- apply(cbind(ran.list.data$focal, ran.list.data$partner), 1, function(x) paste(sort(x), collapse = "_"))
      # ran.list.data.select$dyad <- apply(cbind(ran.list.data.select$focal, ran.list.data.select$partner), 1, function(x) paste(sort(x), collapse = "_"))
      #
      #
      ran.focal.day$observation.time <- 12
      
      ## behaviour rates - balanced data including a bad behaviour (behaviour.5)
      ran.beh.rates <- data.frame(expand.grid(unique(c(data.set$dyad.frame$focal, data.set$dyad.frame$partner)), unique(c(data.set$dyad.frame$focal, data.set$dyad.frame$partner))))
      colnames(ran.beh.rates) <- c("focal", "partner")
      ran.beh.rates <- ran.beh.rates[ran.beh.rates$focal != ran.beh.rates$partner, ]
      ran.beh.rates$behaviour.1 <- unlist(lapply(1:nrow(ran.beh.rates), function(x) {
        return(sum(ran.list.data.select$duration[ran.list.data$behaviour == "behaviour.1" &
                                                   ((ran.list.data$focal == ran.beh.rates$focal[x] &
                                                       ran.list.data$partner == ran.beh.rates$partner[x]) |
                                                      (ran.list.data$partner == ran.beh.rates$focal[x] &
                                                         ran.list.data$focal == ran.beh.rates$partner[x])
                                                   )]))
      }))
      ran.beh.rates$behaviour.1[is.na(ran.beh.rates$behaviour.1)] <- 0
      
      ran.beh.rates$behaviour.2 <- unlist(lapply(1:nrow(ran.beh.rates), function(x) {
        return(sum(ran.list.data.select$duration[ran.list.data$behaviour == "behaviour.2" &
                                                   ((ran.list.data$focal == ran.beh.rates$focal[x] &
                                                       ran.list.data$partner == ran.beh.rates$partner[x]) |
                                                      (ran.list.data$partner == ran.beh.rates$focal[x] &
                                                         ran.list.data$focal == ran.beh.rates$partner[x])
                                                   )]))
      }))
      ran.beh.rates$behaviour.2[is.na(ran.beh.rates$behaviour.2)] <- 0
      
      ran.beh.rates$behaviour.3 <- unlist(lapply(1:nrow(ran.beh.rates), function(x) {
        return(sum(ran.list.data.select$duration[ran.list.data$behaviour == "behaviour.3" &
                                                   ((ran.list.data$focal == ran.beh.rates$focal[x] &
                                                       ran.list.data$partner == ran.beh.rates$partner[x]) |
                                                      (ran.list.data$partner == ran.beh.rates$focal[x] &
                                                         ran.list.data$focal == ran.beh.rates$partner[x])
                                                   )]))
      }))
      ran.beh.rates$behaviour.3[is.na(ran.beh.rates$behaviour.3)] <- 0
      
      ran.beh.rates$behaviour.4 <- unlist(lapply(1:nrow(ran.beh.rates), function(x) {
        return(sum(ran.list.data.select$duration[ran.list.data$behaviour == "behaviour.4" &
                                                   ((ran.list.data$focal == ran.beh.rates$focal[x] &
                                                       ran.list.data$partner == ran.beh.rates$partner[x]) |
                                                      (ran.list.data$partner == ran.beh.rates$focal[x] &
                                                         ran.list.data$focal == ran.beh.rates$partner[x])
                                                   )]))
      }))
      ran.beh.rates$behaviour.4[is.na(ran.beh.rates$behaviour.4)] <- 0
      
      ran.beh.rates$behaviour.5 <- unlist(lapply(1:nrow(ran.beh.rates), function(x) {
        return(sum(ran.list.data.select$duration[ran.list.data$behaviour == "behaviour.5" &
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
      
      
      
      ### now calculating different sociality indices that include: 1. varying number of interactions, 2. Varying number of interactions + bad behaviour, 3. Varying number of interactions + unbalanced data
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
      
      
      # Dyadic Sociality Index - Interaction: Positive Behaviours without proximity, each divided by its mean and combined
      ran.DSI.interaction <- calc_index(dyads = dyads,
                                        behaviours = cbind(
                                          ran.rates.beh1,
                                          ran.rates.beh2
                                        ),
                                        type = "DSI"
      )
      
      
      
      # Dyadic Sociality Index - Unbalanced focal data: Positive Behaviours
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
                     ran.rates.beh5,
                     ran.rates.beh2,
                     ran.rates.beh3,
                     ran.rates.beh4
                   ),
                   valence = c("pos.c", "pos.r", "pos.c", "neg"),
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
      
      
      # ## Create DDSI positive only
      
      normalise <- function(x) {
        return((x - min(x)) / (max(x) - min(x)))
      }
      
      ###### store all the different vectors
      results <- list(
        beh1 = normalise(ran.rates.beh1),
        beh2 = normalise(ran.rates.beh2),
        beh3 = normalise(ran.rates.beh3),
        DSI = normalise(ran.DSI),
        DSI.common = normalise(ran.DSI.common),
        DSI.interaction = normalise(ran.DSI.interaction),
        DSI.bad = normalise(ran.DSI.bad),
        DSI.unbalanced = normalise(ran.DSI.unbalanced),
        
        CRI = normalise(ran.CRI),
        CRI.bad = normalise(ran.CRI.bad),
        CRI.unbalanced = normalise(ran.CRI.unbalanced)
      )
      return(results)
    })
    stopCluster(mycluster)
    
    #### take all the vectors with the same index and put them together
    sim.list.results <- list(
      beh1 = do.call(cbind, lapply(ran.measures, function(x) x$beh1)),
      beh2 = do.call(cbind, lapply(ran.measures, function(x) x$beh2)),
      beh3 = do.call(cbind, lapply(ran.measures, function(x) x$beh3)),
      DSI = do.call(cbind, lapply(ran.measures, function(x) x$DSI)),
      DSI.common = do.call(cbind, lapply(ran.measures, function(x) x$DSI.common)),
      DSI.interaction = do.call(cbind, lapply(ran.measures, function(x) x$DSI.interaction)),
      DSI.bad = do.call(cbind, lapply(ran.measures, function(x) x$DSI.bad)),
      DSI.unbalanced = do.call(cbind, lapply(ran.measures, function(x) x$DSI.unbalanced)),
      
      CRI = do.call(cbind, lapply(ran.measures, function(x) x$CRI)),
      CRI.bad = do.call(cbind, lapply(ran.measures, function(x) x$CRI.bad)),
      CRI.unbalanced = do.call(cbind, lapply(ran.measures, function(x) x$CRI.unbalanced))
    )
    return(sim.list.results)
  })
  
  names(sim.list) <- c(60, 180, 360)
  
  #### create a frame that has for each dyad for each subset size and index the mean and different variation measures
  all.variation <- do.call(rbind, lapply(1:length(sim.list), function(l1) {
    x.subs <- names(sim.list)[[l1]]
    x.index <- names(sim.list[[l1]])
    x.sum <- do.call(rbind, lapply(1:length(sim.list[[l1]]), function(l2) {
      x.names <- x.index[[l2]]
      x.means <- rowMeans(sim.list[[l1]][[l2]])
      x.sd <- apply(sim.list[[l1]][[l2]], 1, sd)
      x.min <- apply(sim.list[[l1]][[l2]], 1, min)
      x.max <- apply(sim.list[[l1]][[l2]], 1, max)
      x.range <- apply(sim.list[[l1]][[l2]], 1, quantile, 0.9, na.rm = T) - 
        apply(sim.list[[l1]][[l2]], 1, quantile, 0.1, na.rm = T)
      x.cv <- x.sd / x.means
      x.se <- x.sd / sqrt(iterations)
      
      # select correct full value
      x.full <- rep(NA, times = length(x.means))
      if (x.names == "beh1") {
        x.full <- full.rates.beh1
      }
      if (x.names == "beh2") {
        x.full <- full.rates.beh2
      }
      if (x.names == "beh3") {
        x.full <- full.rates.beh3
      }
      if (x.names == "CRI") {
        x.full <- full.CRI
      }
      if (x.names == "CRI.bad") {
        x.full <- full.CRI.bad
      }
      if (x.names == "CRI.unbalanced") {
        x.full <- full.CRI
      }
      if (x.names == "DSI") {
        x.full <- full.DSI
      }
      if (x.names == "DSI.bad") {
        x.full <- full.DSI.bad
      }
      if (x.names == "DSI.interaction") {
        x.full <- full.DSI.interaction
      }
      if (x.names == "DSI.common") {
        x.full <- full.DSI.common
      }
      if (x.names == "DSI.unbalanced") {
        x.full <- full.DSI
      }
      x.full <- normalise(x.full)
      
      x.all <- data.frame(
        dataset = z,
        subset = x.subs,
        index = x.names,
        mean = x.means,
        sd = x.sd,
        cv = x.cv,
        se = x.se,
        min = x.min,
        max = x.max,
        range = x.range,
        focal.llh = selection.llh,
        expected.probs = expected.probs,
        full.value = x.full
      )
      return(x.all)
    }))
    return(x.sum)
  }))
  save(all.variation, file = paste(c("precision_", z, ".rda"), collapse = ""))
  return(1)
}))
