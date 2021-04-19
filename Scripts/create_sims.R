##### function to simulate interaction data for social indices
### the user determines the number of individuals, different interaction types, days, certainty of partner choice, interactions per day, number of observers, and number of observation days, as well as mean party size
### each dyad is assigned an expected interaction probability (underlying 'relationship') from a beta-distribution
### The 'certainties' list determines for each behaviour how strongly choice is determined
### outcomes are different ways of presenting the date: a list of interactions (for DDSI), daily interactions per focal with all possible partners (consistency measure), and a dyadic overview that also contains the expected probability for each dyad
## speed is mainly determined by the number of observers and interactions you put


create.sims <- function(n.subj = 20, # number subjects
                        n.behaviours = 3, # number of interaction types
                        n.days = 365, # number of days
                        certainties = list(3, 3, 3), # list of same length as n.behaviours that determines the factor by which the probabilities to interact are increased; above 1.25 would be 'high' and means that every individual has one friend, between 0.75 and 1.25 is medium, below 0.75 is low
                        n.int.per.day = list(5, 5, 5), # list of same length as n.behaviours; contains number of interactions each individual has on average on a day
                        pos.neg = list("pos", "pos", "neg"), # is behaviour positive or negative? 
                        state.event = list("event", "state", "event"), # list of same length as n.behaviours, determines whether interaction type is an event(duration = 1) or gets a
                        n.observers = 1, # number of individuals focaled per day
                        n.obs.days = 0.66, # proportion of days with focal observer
                        party.size = 1 # value between 0.1 and 1; if 1, all individuals are always present, if 0.5 only half of individuals are usually present as partner etc; ID of present individuals is random; selection is with replacement, so even if you put '1', not everyone will be present
) 
{
  library(tidyverse)
  # name individuals
  ids <- n.subj:1
  
  # create frame for each focal - partner dyads
  dyad.frame = data.frame(expand.grid(list(ids, ids)))
  dyad.frame = dyad.frame[order(dyad.frame$Var1, dyad.frame$Var2),]
  dyad.frame = subset(dyad.frame, Var1 != Var2)
  colnames(dyad.frame) = c('focal', 'partner')
  
  # add expected probabilities for each dyad
  expected.probs = rbeta(n = 10000, shape1 = 0.5, shape2 = 2)
  dyad.frame$expected.probs = sample(expected.probs, nrow(dyad.frame))
  
  # standardise expected probabilities to sum to 1 within individuals 
  dyad.frame$expected.probs = sapply(1:nrow(dyad.frame), function(x){
    xx.all = dyad.frame$expected.probs[dyad.frame$focal == dyad.frame$focal[x]]
    xx.all.stan = xx.all / sum(xx.all)
    return(xx.all.stan[xx.all == dyad.frame$expected.probs[x]][1])
  })
  
  
  # create dates
  rdate <- seq.Date(from = as.Date("2020-01-01"), to = as.Date("2020-01-01") + n.days, by = 1)
  # randomly select dates based on n.obs.days
  rdate <- sample(rdate, size = round(n.obs.days * n.days))
  # assign focals to each day based on n.observers
  sim.data <- as_tibble(data.frame(do.call(rbind, lapply(rdate, function(x) {
    xx.ids <- data.frame(sample(ids, n.observers), x)
    xx.ids
  }))))
  colnames(sim.data) <- c("focal", "rdate")
  sim.data$rdate <- as.Date(sim.data$rdate, origin = "1970-01-01")
  
  # create data frame with observations times per day
  obs.times <- sim.data
  obs.times$observation.hours <- 12
  
  # create data frame with each individual and each possible partner for a day and how often they interacted
  daily.ints <- as_tibble(data.frame(do.call(rbind, lapply(1:nrow(sim.data), function(x) {
    xx.ids <- data.frame(focal = sim.data$focal[x], rdate = sim.data$rdate[x], partner = ids, observation.time = 12, interactions = 0)
    xx.ids <- xx.ids[xx.ids$focal != xx.ids$partner, ]
    xx.ids
  }))))
  
  ### here need to incorporate different time periods
  
  
  # for each interaction type, select a partner based on the certainty and number of interactions
  interaction.type <- lapply(1:n.behaviours, function(z) {
    xdata <- sim.data
    # determine how many interactions each individual has per day
    xdata <- as_tibble(data.frame(do.call(rbind, lapply(1:nrow(xdata), function(x) { # assign focals to each day based on n.observers
      nr.partners <- round(rnorm(100, mean = n.int.per.day[[z]], sd = n.int.per.day[[z]] / 3))
      nr.partners <- sample(nr.partners[nr.partners > 0], 1)
      xx.ids <- t(matrix(unlist(rep(as.vector(xdata[x, ]), nr.partners)), nrow = 2)) # duplicate each focal per day by the number of interactions it should have
      xx.ids
    }))))
    colnames(xdata) <- c("focal", "rdate")
    xdata$rdate <- as.Date(xdata$rdate, origin = "1970-01-01")
    
    # partner choice is determined by dyadic expected interaction probability
    partner.choice <- lapply(xdata$focal, function(x) {
      probs <- dyad.frame$expected.probs[dyad.frame$focal == x]
      
      probs <- probs ^ (certainties[[z]]) # change distribution based on the 'certainty' parameter for the interaction type
      probs <- (probs / sum(probs))
      
      # for negative behaviours, turn probability
      if (pos.neg[[z]] == "neg") {
        probs <- dyad.frame$expected.probs[dyad.frame$focal == x]
        probs <- 1 - sqrt(sqrt(sqrt(probs)))
        
        probs <- (probs / sum(probs))
        probs <- probs ^ (2 * certainties[[z]])
        probs <- (probs / sum(probs))
      }
      
      # randomly select visible subset of individuals by randomly selecting without replacement, based on the party.size parameter
      ran.sel <- unique(sample(1:length(probs), size = round(party.size * length(probs), 0), replace = F))
      ids.sel <- dyad.frame$partner[dyad.frame$focal == x]
      probs <- probs[ran.sel]
      ids.sel <- ids.sel[ran.sel]
      
      xx <- sample(ids.sel, size = 1, prob = probs) # select partner by high certainty
      res <- c(partner = xx, probs = as.numeric(probs[ids.sel == xx]), group = paste(sort(ran.sel), collapse = "/"))
      return(res)
    })
    partner.ch <- do.call(rbind, partner.choice)
    xdata$partner <- as.numeric(partner.ch[, 1])
    xdata$prob <- as.numeric(partner.ch[, 2])
    xdata$party <- partner.ch[, 3]
    # insert duration
    xdata$duration <- 1
    if (state.event[[z]] == "state") {
      xdata$duration <- (sample(rnbinom(1000, s = 100, m = 5), replace = T, size = nrow(xdata))) + 0.5
    }
    xdata$behaviour <- paste(c("behaviour", z), collapse = ".")
    
    return(xdata)
  })
  
  interactions.list <- do.call(rbind, interaction.type)
  
  # create number of interactions per dyad per day
  ints.to.add <- lapply(1:n.behaviours, function(z) {
    xx <- interaction.type[[z]]
    xx.ints <- unlist(lapply(1:nrow(daily.ints), function(x) {
      sum(xx$duration[xx$focal == daily.ints$focal[x] &
                        xx$partner == daily.ints$partner[x] &
                        xx$rdate == daily.ints$rdate[x]])
    }))
    xx.ints
  })
  
  names(ints.to.add) <- sapply(1:n.behaviours, function(z) {
    paste(c("behaviour", z), collapse = ".")
  })
  daily.ints <- cbind(daily.ints, do.call(cbind, ints.to.add))
  daily.ints <- daily.ints[order(daily.ints$rdate), ]
  
  # create frame with all interactions per hour per dyad and 'real' value
  dyad.frame <- dyad.frame[order(dyad.frame$focal), ]
  
  # summarize for each dyad the observation time
  dyad.frame$observation.time <- unlist(lapply(1:nrow(dyad.frame), function(x) {
    sum(obs.times$observation.hours[obs.times$focal %in% c(dyad.frame$focal[x], dyad.frame$partner[x])])
  }))
  
  # create number of interactions per dyad
  ints.to.add <- lapply(1:n.behaviours, function(z) {
    xx <- interaction.type[[z]]
    xx.ints <- unlist(lapply(1:nrow(dyad.frame), function(x) {
      sum(xx$duration[xx$focal %in% c(dyad.frame$focal[x], dyad.frame$partner[x]) &
                        xx$partner %in% c(dyad.frame$focal[x], dyad.frame$partner[x])])
    }))
    xx.ints
  })
  names(ints.to.add) <- sapply(1:n.behaviours, function(z) {
    paste(c("behaviour", z), collapse = ".")
  })
  dyad.frame <- cbind(dyad.frame, do.call(cbind, ints.to.add))
  
  # create number of interactions per dyad per observation hour
  ints.to.add <- lapply(ints.to.add, function(z) {
    z / dyad.frame$observation.time
  })
  names(ints.to.add) <- sapply(1:n.behaviours, function(z) {
    paste(c("behaviour", z, "ph"), collapse = ".")
  })
  dyad.frame <- cbind(dyad.frame, do.call(cbind, ints.to.add))
  
  return(
    list(
      interaction.list = interactions.list, # list of interactions between focal and receiver by date
      daily.interactions = daily.ints, # frame with all focals on all days and all possible partners
      dyad.frame = dyad.frame
    )
  )
}
