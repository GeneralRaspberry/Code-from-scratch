# Trying to solve parallel issue for Katherine's student
# Written by Dave S (david.schoeman@gmail.com)
# 30 November 2022


# Packages ----------------------------------------------------------------

library(tidyverse)
library(furrr)


# Setting up ----------------------------------------------------------------

# Make your parameter data frame like this
d <- expand.grid(betavalues = seq(from = 50, to = 150, length.out = 3),
                 thetavalues = seq(from = 200, to = 300, length.out = 3),
                 deltavalues = seq(from = 20, to = 140, length.out = 3)) %>% 
  mutate(Counter = 1:nrow(.)) # Add a variable called counter

# Build a function to manipulate the variables you have
do_stuff <- function(Counter, betavalues, thetavalues, deltavalues,...) { # The function has arguments that match your parameters
  output <- data.frame(i = Counter, beta = betavalues, theta = thetavalues, delta = deltavalues) # Use those argumnents to do something for each line of your data frame
  return(output) # Output the result of the function
}

# The purrr version of a straight nested for loop
result <- pmap_dfr(d, do_stuff) # pmap() would return a list by rows of d, pmap_dfr() returns the corresponding data frame
result

# The furrr version of a parallelised nested for loop
cores <- availableCores()-1 # How many free cores do I have, less 1 to be safe
plan(multisession, workers = cores) # Start paralell processing
#*** NOTE that on Windoze systems, you might use multicore instead of multisession
result <- future_pmap_dfr(d, do_stuff) # future_pmap() would return a list by rows of d, future_pmap_dfr() returns the corresponsing data frame
result
plan(sequential) # Return to single-core

