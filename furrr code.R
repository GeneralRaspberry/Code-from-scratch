library(tidyverse)
library(furrr)
library(spatstat)

# Setting up ----------------------------------------------------------------

# Make your parameter data frame like this
d <- expand.grid(betavalues = seq(from = 50, to = 150, length.out = 3),
                 thetavalues = seq(from = 200, to = 300, length.out = 3),
                 deltavalues = seq(from = 20, to = 140, length.out = 3)) %>% 
  mutate(Counter = 1:nrow(.)) # Add a variable called counter

# Build a function to manipulate the variables you have
do_stuff <- function(beta, # transmission rate
                     theta, # dispersal scale
                     b=1, # kernel shape parameter, 1 for exponential
                     sigma=0, # asymptomatic period, used for outputting the time series
                     q0=0, # starting incidence if ppp is without marks
                     q.end=1, # stopping condition 1: incidence lvl
                     t.end=100, # stopping condition 2: time after first simulated time step
                     area.host=1, # surface area occupied by one host
                     delta.t=10, # time step
                     ppp, # point pattern as a ppp object, optionally with marks 1/0 for infected/healthy
                     dist.mat=NULL){ # matrix distance if its computation is to be avoided here (for e.g. repeated calls)
  ## if the point pattern has no marks, generate some randomly that fits q0
  if (is.null(marks(ppp))){
    inf.start <- max(1, round(ppp$n * q0))
    marks(ppp) <- sample(c(rep(FALSE, ppp$n-inf.start), rep(TRUE, inf.start)))
  }
  
  ## compute distance matrix if not provided
  if (is.null(dist.mat)){ 
    ## add the kernel computation that can be added directly on the dist matrix to reduce comp time
    dist.mat <- exp(-pairdist(ppp)^b / theta^b)
    diag(dist.mat) <- NA
  }
  
  ## function that compute infection event probability, based on the dispersal kernel
  k.norm <- beta * area.host * (b/(2*pi*theta^2*gamma(2/b))) # constant part of the exponential power kernel
  infection <- function(infected, dist){
    inf <-  matrix(k.norm * dist[infected,!infected],
                   nrow=sum(infected), byrow=FALSE)
    inf[is.na(inf)] <- 0
    inf
  }
  
  ## starting time
  time <- 0
  ## inititate the heavy dataframe that will aggregate all changes
  df.big <- data.frame(time=0, who=which(ppp$marks), t(ppp$marks))
  
  
  ## computation loop
  while (any(!ppp$marks) & time <= t.end & mean(ppp$marks) < q.end){
    
    ## infection event probaility
    events <- infection(ppp$marks, dist=dist.mat)
      ## random proisson draws
      new.infected <- which(!ppp$marks)[rpois(n=sum(!ppp$marks), lambda=apply(events, 2, sum) * delta.t) > 0]
      ## change marks of newly infected
      ppp$marks[new.infected] <- TRUE
      ## increment time
      time <- time + delta.t
      ## if some infection, increment the big dataframe
      if (length(new.infected) > 0){
        df.big <- rbind(df.big, data.frame(time=time, who=new.infected, t(ppp$marks)))
      }
      ## print a dot per new infection
      # cat(paste0(rep('.', length(new.infected)), collapse = '')) ## comment for quiet
      ## increment counter and print value after each simulation
    }
    ## make compact, time only, version of the big dataframe
    times.i <- unique(df.big[,1])
    times.d <- times.i + sigma
    times <- sort(unique(c(times.i, times.d)))
    infected <- sapply(times, FUN=function(t) sum(t >= df.big[,1]))
    detectable <- sapply(times, FUN=function(t) sum(t >= df.big[,1] + sigma))
    df.small <- data.frame(time=times, infected=infected, detectable=detectable)
    
    ## out put the simplified time series, and the big one
    list(df.small[df.small$time <= max(df.big$time),], df.big) 
  } 

# The furrr version of a parallelised nested for loop
cores <- availableCores()-1 # How many free cores do I have, less 1 to be safe
plan(multisession, workers = cores) # Start parallel processing
#*** NOTE that on Windoze systems, you might use multicore instead of multisession
result <- future_pmap_dfr(d, do_stuff) # future_pmap() would return a list by rows of d, future_pmap_dfr() returns the corresponsing data frame
plan(sequential) # Return to single-core