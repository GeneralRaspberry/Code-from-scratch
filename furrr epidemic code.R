#define beta and theta for later concatenation...

be<-100
th<-100


tauLeapG <- function(beta=be, # transmission rate
                     theta=th, # dispersal scale
                     inf.start=1, #number of initial infected
                     b=1, # kernel shape parameter, 1 for exponential
                     sigma=0, # asymptomatic period, used for outputting the time series
                     q0=0, # starting incidence if ppp is without marks
                     q.end=.01, # stopping condition 1: incidence lvl
                     t.end=100, # stopping condition 2: time after first simulated time step
                     area.host=1, # surface area occupied by one host
                     delta.t=100, # time step
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
  #add counter
  
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
    
    ## make compact, time only, version of the big dataframe
    times.i <- unique(df.big[,1])
    times.d <- times.i + sigma
    times <- sort(unique(c(times.i, times.d)))
    infected <- sapply(times, FUN=function(t) sum(t >= df.big[,1]))
    detectable <- sapply(times, FUN=function(t) sum(t >= df.big[,1] + sigma))
    df.small <- data.frame(time=times, infected=infected, detectable=detectable)
    
    ## out put the simplified time series, and the big one
    return(list(df.small[df.small$time <= max(df.big$time),], df.big)) 
  } 
}


#create parameter combinations
#d <- expand.grid(betavalues = seq(from = 50, to = 150, length.out = 3),
                # thetavalues = seq(from = 200, to = 300, length.out = 3))#%>% 
 # mutate(Counter = 1:nrow(.)) # Add a variable called counter
cores <- availableCores()-2
plan(multisession, workers=cores)
output <- future_map(landscapes, function(x) tauLeapG(ppp = x))
plan(sequential)
#isolate df.big
df.bigtemp <- lapply(output,function(x) x[[2]])

#isolate time stamp of infection and who
temp<-lapply(df.bigtemp,function(x) x[,1:2][order(x[,2]),])

#add co-ordinates to output

for (i in 1:length(temp)) {
  # Get the current output and corresponding landscape
  current_output<-temp[[i]]
  current_landscape<-landscapes[[i]]
  
  # Add the coordinates to the current output
  current_output$x <- current_landscape$x[current_output$who]
  current_output$y <- current_landscape$y[current_output$who]
  
  # Update the output list with the modified output
  temp[[i]] <- current_output
}


final_output <- list(c(beta = be, theta = th, randmod = rlf), temp)

