### load packages ####
source("Scripts/create_sims.R")
library(tidyverse)

### Create datasets with different properties for different group sizes and store them somewhere ####

nr.ids <- sort(
  c(rep(c(5), 3))) # we will create 10 groups with 25 individuals

##### all frames have 4 behaviours:
# Frequent Positive Medium-Certainty (Grooming, 'behaviour.1')
# Rare Positive High-Certainty (Food Sharing, 'behaviour.2')
# Frequent Positive Low-Certainty (Contact Sitting, 'behaviour.3')
# Medium-Rare Negative Medium-Certainty (Aggression, 'behaviour.4')

sim.certainties <- list(1.2, 3, 0.5, 2.5) # 'certainty' here refers to the selectiveness of individuals when making decisions. The higher, the more selective. Check things like coefficient of variance after creating data to check how different they are
sim.int.per.day <- list(3, 1, 8, 3) # number of interactions of this type each individual has eaech day
sim.pos.neg <- list("pos", "pos", "pos", "neg") # for each interaction type, decide whether they reflect relationships positively or negatively
sim.state.event <- list("state", "event", "event", "event") # select whether interaction type is an 'event' (i.e., 1 aggression) or a state (i.e., 65sec of grooming). Durations will be randomly attached

##### all datasets have 360 days of observation,
sim.n.days <- 360

##### all datasets have the same number of observers as they have individuals. This allows us to know ALL interactions
sim.n.observers <- nr.ids

##### all interactions take place in the presence of half of the group
sim.party.size <- 0.5


# simulations for all 'groups'. This takes a really long time. Can be parallelised
sim.data.sets <- lapply(nr.ids, function(y) {
  sim.data.set <- create.sims(
    n.subj = y, # number subjects
    n.behaviours = 4, # number of behaviours
    n.days = sim.n.days, # number of days
    certainties = sim.certainties, # list of same length as n.behaviours; 'high' means individuals choose always same partner when available, 'medium' are quite choosy, 'low' are egalitarian
    n.int.per.day = sim.int.per.day, # list of same length as n.behaviours; contains number of interactions each individual has on average on a day
    pos.neg = sim.pos.neg, # is behaviour positive or negative? positive behaviours target closeby individuals and kin; negative behaviours target faraway individuals and non-kin
    state.event = sim.state.event, # decides whether each interaction type is a state ( == has duration) or an event
    n.observers = y, # number of individuals focaled per day
    n.obs.days = 1, # proportion of days with focal observer
    party.size = sim.party.size
  )
  return(sim.data.set)
})



# save image
# save(list = c("sim.data.sets", "nr.ids", "sim.n.days", "sim.certainties", "sim.state.event", "sim.int.per.day", "sim.pos.neg", "sim.n.observers", "sim.party.size"), file = "simulated_data.RData")
