
# Load Packages -----------------------------------------------------------
library(tidyverse)

# Create Individuals ------------------------------------------------------

# Set seed for reproducibility
set.seed(123)
nr_individuals <- 20
error_ranks <- 0.05
rank_subset <- 0.3
nr_males <- 8
male_male_bias <- 0.1
female_female_bias <- 4

# Generate x individuals with hierarchy rank values between 1 and x

individuals <- data.frame(
  ID = 1:nr_individuals,
  Rank = 1:nr_individuals
)

probabilities <- function(numbers, difference) {

  #Calculate the probability of 10
  probability_max <- difference

  #Calculate the probability of all other numbers
  probability_min <- 1

  #Create a vector of probabilities
  probabilities <- seq(from = probability_min, to = probability_max, length.out = length(numbers))

  #Turn into probabilities
  probabilities <- probabilities/sum(probabilities)

  #Return the vector of probabilities
  return(probabilities)
}


# Generate Rank Interactions ----------------------------------------------


# Generate 5000 interactions between pairs of individuals
interactions <- data.frame()
for (i in 1:5000) {
  # Select two individuals with similar rank values
  sender <- individuals[sample.int(nrow(individuals), 1), ]
  receiver <- individuals[sample.int(sender$ID, 1, prob = probabilities(1:sender$ID, difference = 10)), ]

  # Swap sender and receiver for 5% of the interactions
  if (runif(1) <= error_ranks) {
    temp <- sender
    sender <- receiver
    receiver <- temp
  }

  # Swap sender and receiver for dyads with 1 rank difference in 5% of cases
  if (abs(sender$Rank - receiver$Rank) == 1 && runif(1) <= error_ranks) {
    temp <- sender
    sender <- receiver
    receiver <- temp
  }

  # Add the interaction to the data frame with probability
  interactions <- rbind(interactions, data.frame(
    Sender = sender$ID,
    Receiver = receiver$ID
  ))
}

# Randomly assign a date between 01/01/2022 and 31/12/2022
start_date <- as.Date("2022-01-01")
end_date <- as.Date("2022-12-31")
interactions$Date <- sample(seq(start_date, end_date, by = "day"), nrow(interactions), replace = TRUE)

# Order the interactions data frame by date
interactions <- interactions[order(interactions$Date), ]

# Randomly Select Subset of Rank Interactions -----------------------------

subset_interactions <- interactions %>% sample_frac(rank_subset)


# Calculate Elo Ratings ---------------------------------------------------

# Calculate Elo ratings for dominance hierarchy
elo_ratings <- rep(1000, nrow(individuals))  # Initialize Elo ratings for each individual
k_factor <- 20  # K-factor determines the impact of each interaction on Elo ratings

# Function to calculate expected probability of winning for a given Elo rating difference
calculate_expected <- function(rating_diff) {
  return(1 / (1 + 10^(-rating_diff/400)))
}

# Iterate through interactions to update Elo ratings
for (i in 1:nrow(subset_interactions)) {
  sender <- subset_interactions$Sender[i]
  receiver <- subset_interactions$Receiver[i]

  sender_rank <- individuals$Rank[individuals$ID == sender]
  receiver_rank <- individuals$Rank[individuals$ID == receiver]

  sender_elo <- elo_ratings[sender]
  receiver_elo <- elo_ratings[receiver]

  elo_diff <- sender_elo - receiver_elo

  # Calculate expected probabilities
  expected_sender <- calculate_expected(elo_diff)
  expected_receiver <- 1 - expected_sender

  # Update Elo ratings
  elo_ratings[sender] <- sender_elo + k_factor * (1 - expected_sender)
  elo_ratings[receiver] <- receiver_elo + k_factor * (0 - expected_receiver)
}

# Add Elo ratings to individuals data frame
individuals$EloRating <- elo_ratings



# Add Sex ----------------------------------------------------------------

individuals$Sex <- ifelse(individuals$ID >= nr_individuals - nr_males, "m", "f")




# Similarity Based on Elo -------------------------------------------------

# Calculate Elo similarity matrix
elo_similarity <- as.matrix(dist(individuals$EloRating, method = "euclidean"))

# Convert Elo similarity matrix to probability matrix
similarity_prob <- 1 - (elo_similarity/max(elo_similarity))
diag(similarity_prob) <- 0
similarity_prob <- similarity_prob ^ 3
similarity_prob <- similarity_prob / rowSums(similarity_prob, na.rm = TRUE)



# Rank as variable --------------------------------------------------------

# Create a matrix with the same dimensions and row/column names as similarity_prob
elo_matrix <- matrix(0, nrow = nrow(similarity_prob), ncol = ncol(similarity_prob))
rownames(elo_matrix) <- rownames(similarity_prob)
colnames(elo_matrix) <- colnames(similarity_prob)

# Attach corresponding Elo rating to each row
for (i in 1:nrow(elo_matrix)) {
  elo_matrix[i, ] <- individuals %>% arrange(ID) %>% pull(EloRating)
}
elo_matrix <- elo_matrix + abs(min(elo_matrix))
diag(elo_matrix) <- 0
elo_matrix <- elo_matrix/max(elo_matrix)

# Square the values to ensure stronger impact of higher rank
elo_matrix <- exp(elo_matrix ^ 3)

# Calculate probability of being chosen as interaction partner
elo_prob <- elo_matrix / rowSums(elo_matrix, na.rm = TRUE)



# Interactions based on Sex -----------------------------------------------

# Create an adjacency matrix with the same dimensions and row/column names as similarity_prob
adjacency_matrix <- matrix(0, nrow = nrow(similarity_prob), ncol = ncol(similarity_prob))
rownames(adjacency_matrix) <- rownames(similarity_prob)
colnames(adjacency_matrix) <- colnames(similarity_prob)

# Assign values to the adjacency matrix based on the 'sex' column
for (i in 1:nrow(adjacency_matrix)) {
  for (j in 1:ncol(adjacency_matrix)) {
    if (individuals$Sex[i] == individuals$Sex[j]) {
      if (individuals$Sex[i] == "m") {
        adjacency_matrix[i, j] <- male_male_bias  # Both individuals are male
      } else {
        adjacency_matrix[i, j] <- female_female_bias  # Both individuals are female
      }
    } else {
      adjacency_matrix[i, j] <- 1  # Different values (male and female)
    }
  }
}

# Calculate probability of being chosen as interaction partner
sex_prob <- adjacency_matrix / rowSums(adjacency_matrix, na.rm = TRUE)

# Define Rules for Interactions -------------------------------------------

# Set the seed for reproducibility
set.seed(123)

# Create the 'grooming_data' data frame
grooming_data <- data.frame(Sender = character(10000), Receiver = character(10000), Date = as.Date(character(10000), format = "%Y-%m-%d"), stringsAsFactors = FALSE)
all_matrix <- sex_prob * similarity_prob * elo_prob
all_matrix <- all_matrix ^ 3
all_matrix <- all_matrix/rowSums(all_matrix)

individuals$Sociability <- ifelse(individuals$Sex == 'm', sample(rnorm(100, mean = 20, sd = 3)), sample(rnorm(100, mean = 40, sd = 5)))
individuals$Sociability <- individuals$Sociability / sum(individuals$Sociability)

# Iterate through each row of 'grooming_data' to simulate data points
for (i in 1:10000) {
  # Randomly select a Sender
  grooming_data$Sender[i] <- sample(individuals$ID, 1, prob = individuals$Sociability)

  # Randomly select the Date between 01/01/2022 and 31/12/2022
  grooming_data$Date[i] <- as.Date(sample(ISOdate(2022, 1:12, 1:31), 1))

  # Randomly select the Receiver based on the probabilities of the Sender
  sender_index <- match(grooming_data$Sender[i], individuals$ID)

  grooming_data$Receiver[i] <-
    sample(individuals$ID, 1, prob = all_matrix[sender_index, ])
}


xx = grooming_data %>%
  group_by(Sender, Receiver) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(Sender = as.character(Sender)) %>%
  left_join(individuals %>%
              mutate(ID = as.character(ID)) , by = c('Sender' = 'ID')) %>%
  rename(c(EloSender = 'EloRating', SexSender = 'Sex')) %>%
  left_join(individuals %>%
              mutate(ID = as.character(ID)) , by = c('Receiver' = 'ID')) %>%
  rename(c(EloReceiver = 'EloRating', SexReceiver = 'Sex'))

model = lm(data = xx, count ~ EloSender * EloReceiver + SexReceiver * SexSender)

library(broom)
library(ggplot2)
library(plotly)

# Get predicted values from the model
predicted_data <- augment(model)


# Create the 3D scatter plot
plot_ly(predicted_data, x = ~EloSender, y = ~EloReceiver, z = ~.fitted, color = ~.fitted,
        colors = viridisLite::viridis(256), type = "scatter3d", mode = "markers") %>%
  layout(scene = list(xaxis = list(title = "EloSender"),
                      yaxis = list(title = "EloReceiver"),
                      zaxis = list(title = "Predicted Count"),
                      aspectmode = "manual",
                      aspectratio = list(x = 1, y = 1, z = 0.8)))

plot_ly(predicted_data, x = ~SexSender, y = ~SexReceiver, z = ~.fitted, color = ~.fitted,
        colors = viridisLite::viridis(256), type = "scatter3d", mode = "markers") %>%
  layout(scene = list(xaxis = list(title = "EloSender"),
                      yaxis = list(title = "EloReceiver"),
                      zaxis = list(title = "Predicted Count"),
                      aspectmode = "manual",
                      aspectratio = list(x = 1, y = 1, z = 0.8)))

