############################################################################################
# Author: Hana Fu
# Date: Jan 23, 2022
#
# This code performs the basic microsimulation model for the PHO Practicum 2021-22.
# Adapted to performed a three-stage microsimulation model.
#
# Codes are largely based on those by this tutorial paper:
# Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22.
#
# GitHub link for the original codes based upon is found here:
# https://github.com/DARTH-git/Microsimulation-tutorial/blob/master/Appendix%20B_online_supp-samplev.R
# 
############################################################################################


##################################### Functions ###########################################

# THE NEW samplev() FUNCTION
# efficient implementation of the rMultinom() function of the Hmisc package #### 

samplev <- function (probs, m) { # probs is a n x k dimenion matrix, with n = n.i (number of individuals) and k = n.s (number of health states), that is it contains the transition probabilities to various health state for each individual. m is the number of health state(s) to return per person for the next cycle, not entirely sure why this option is offered, it has a default of m = 1, maybe it gives an option to sample again from the possible m health states? 
  d <- dim(probs) # gives the dimension of the probs matrix
  n <- d[1] # return n.i
  k <- d[2] # return n.s
  lev <- dimnames(probs)[[2]] # return names of health states
  if (!length(lev)) # if the probs matrix doesn't have column names 
    lev <- 1:k # set lev as 1:k to represent names of health states 
  ran <- matrix(lev[1], ncol = m, nrow = n) # create a matrix with the name of first health state, to store health state in next cycle
  un <- matrix(NA, ncol = m, nrow = n) # create a matrix to store randomly generated uniform numbers
  U <- t(probs) # U is a transpose of probs matrix, thus it is a k x n matrix
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ] # This converts the transition probabilities to "cumulative probabilities", it applies to all individual at once. 
    # For example, if the current health state for individual j is "healthy", the transition probabilities is c(0.845, 0.15, 0.005), then U[,j] becomes c(0.845,0.995,1).
  }
  if (any((U[k, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) { # looping through 
    un[, j] <- runif(n) 
    un.k <- rep(un[, j], rep(k, n)) # returns a n x k length vector, where random uniform number for each person is repeated k times (correspond to the number of health states)
    ran[, j] <- lev[1 + colSums(un.k > U)] # colSums(un > U) obtain the sum for each person that the "un" is greater than the transition probabilities of the health states. 
    # For example, if the transition probabilies are c(0.845,0.995,1) and "un" for the person is c(0.85, 0.85, 0.85), then colSUms(un > U) equals to 1. 
    # Including "1 +" so to select the health state for the next cycle, that is lev[1+1] is "S", the person's health state in the next cycle is "S". 
  }
  return(list(ran=ran, un=un))
}
# The MicroSim function for the simple microsimulation of the 'Sick-Sicker' model keeps track of what happens to each individual during each cycle. 

MicroSim <- function(v.M_1, n.i, n.t, v.n, samplev.m, X = NULL, d.e, TR.out = TRUE, TS.out = TRUE, seed = 1) {
  # Arguments:  
  # v.M_1:   vector of initial states for individuals
  # n.i:     number of individuals
  # n.t:     total number of cycles to run the model
  # v.n:     vector of health state names
  # X:       vector or matrix of individual characteristics
  # d.e:     discount rate for health outcome (QALYs)
  # TR.out:  should the output include a Microsimulation trace? (default is TRUE)
  # TS.out:  should the output include a matrix of transitions between states? (default is TRUE)
  # seed:    starting seed number for random number generator (default is 1)
  # Makes use of:
  # Probs:   function for the estimation of transition probabilities
  # Effs:    function for the estimation of state specific health outcomes (QALYs)
  
  v.dwe <- 1 / (1 + d.e) ^ (0:n.t)   # calculate the QALY discount weight based on the discount rate d.e
  
  # Create the matrix capturing the state name/health outcomes for all individuals at each time point 
  m.M <- m.E <- matrix(nrow = n.i, ncol = n.t + 1, 
                       dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                       paste("cycle", 0:n.t, sep = " ")))  
  
  m.M[, 1] <- v.M_1                     # indicate the initial health state
  
  # Create matrix to store randomly generated uniform numbers
  m.un <- matrix(rep(NA, n.i*n.t), nrow = n.i, ncol = n.t,
                  dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                  paste("cycle", 1:n.t, sep = " "))) 
  
  # Create array to store transition probabilities
  arr.p <- array(NA, dim = c(n.i, n.s, n.t + 1),
                 dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                 paste("state", 1:n.s, sep = " "),
                                 paste("cycle", 0:n.t, sep = " "))) 
  
  arr.p[,, 1] <- rep(c(1 - p.HS - p.HD, p.HS, p.HD), each = n.i)
  
  set.seed(seed)                  # set the seed for every individual for the random number generator
  
  # create the dur variable that stores the number of consecutive cycles the individual occupies either when sick or sicker
  dur <- rep(0, n.i)                # the individual start without history        
  m.E[, 1] <- Effs (m.M[, 1], dur, X = X)  # estimate QALYs per individual for the initial health state
  
  for (t in 1:n.t) { # t <- 3
    m.p <- Probs(m.M[, t], dur)           # calculate the transition probabilities at cycle t 
    
    arr.p[,, t+1] <-  m.p # store the transition probabilities
    
    ran.un <- samplev(prob = m.p, m = samplev.m) # returns a matrix of the health states and a matrix of randomly generated uniform numbers
  
    # this is Hana's guess of the function of the m argument in samplev() function
    if (samplev.m > 1) { 
      which.m <- sample(1:samplev.m, size = 1) # pick which m index to use 
      
      m.un[, t] <- ran.un[['un']][, which.m] # store the matrix of randomly generated uniform numbers for each cycle
      m.M[, t + 1] <- ran.un[['ran']][, which.m]  # sample the next health state and store that state in matrix m.M
    } else {
      m.un[, t] <- ran.un[['un']]
      m.M[, t + 1] <- ran.un[['ran']]
    }
    
    m.E[, t + 1] <- Effs( m.M[, t + 1], dur, X = X)   # estimate QALYs per individual during cycle t + 1 conditional on treatment
    
    dur <- ifelse(m.M[, t + 1] == "S", 
                  dur + 1, 
                  0) # updates duration (add 1) if person is in "Sick" state, else assign duration as zero, can be updated to consider previous duration in "Sick" state even if the person returns back to "Healthy" states
    
    cat('\r', paste(round(t/n.t * 100), "% done", sep = " "))       # display the progress of the simulation
    
  } # close the loop for the time points 
  
  
  te <- m.E %*% v.dwe       # total (discounted) QALYs per individual 
  
  te_hat <- mean(te)        # average (discounted) QALYs
  
  if (TS.out == TRUE) {  # create a matrix of transitions across states
    TS <- paste(m.M, cbind(m.M[, -1], NA), sep = "->") # transitions from one state to the other
    TS <- matrix(TS, nrow = n.i)
    rownames(TS) <- paste("Ind",   1:n.i, sep = " ")   # name the rows 
    colnames(TS) <- paste("Cycle", 0:n.t, sep = " ")   # name the columns 
  } else {
    TS <- NULL
  }
  
  if (TR.out == TRUE) {
    TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE))))
    TR <- TR / n.i                                       # create a distribution trace
    rownames(TR) <- paste("Cycle", 0:n.t, sep = " ")     # name the rows 
    colnames(TR) <- v.n                                  # name the columns 
  } else {
    TR <- NULL
  }
  results <- list(m.p = m.p, arr.p = arr.p, m.M = m.M, m.E = m.E, te = te, te_hat = te_hat, m.un = m.un, TS = TS, TR = TR) # store the results from the simulation in a list  
  return(results)  # return the results
}  # end of the MicroSim function  


#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below.

Probs <- function(M_it, dur) { 
  # M_it:   health state occupied by individual i at cycle t (character variable)
  # dur:    the duration of being sick (sick/sicker)
  
  m.p.it <- matrix(NA, n.s, n.i)     # create vector of state transition probabilities
  rownames(m.p.it) <- v.n            # assign names to the vector
  
  # update probabilities of death after first converting them to rates and applying the rate ratio
  r.SD <-  - log(1 - p.SD)
  p.SD <- 1 - exp(- r.SD * (1 + dur[M_it == "S"] * rp.S)) # calculate p.SD conditional on duration of being sick
  
  # update the v.p with the appropriate probabilities   
  m.p.it[,M_it == "H"]  <- c(1 - p.HS - p.HD, p.HS, p.HD)                  # transition probabilities when healthy
  m.p.it[,M_it == "S"] <- rbind(p.SH, 1- p.SH - p.SD, p.SD)   # transition probabilities when sick
  m.p.it[,M_it == "D"]  <- c(0, 0, 1)                                        # transition probabilities when dead   
  ifelse(colSums(m.p.it) == 1, return(t(m.p.it)), print("Probabilities do not sum to 1")) # return the transition probabilities or produce an error
}


### Health outcome function 
# The Effs function to update the utilities at every cycle.

Effs <- function (M_it, dur, cl = 1, X = NULL) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # dur:  the duration of being sick
  # cl:   cycle length (default is 1)
  # X:    the vector or matrix of individual characteristics (optional)
  
  u.it               <- 0        # by default the utility for everyone is zero
  u.it[M_it == "H"]  <- u.H      # update the utility if healthy
  u.it[M_it == "S"] <- X[which(M_it == "S")] * ( u.S + dur[which(M_it == "S")] * ru.S ) # update the utility if sick conditional on treatment and duration of being sick
  u.it[M_it == "D"]  <- 0        # update the utility if dead
  
  QALYs <-  u.it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
}


##################################### Run the simulation ##################################

#sim.list <- list()
comp.time <- list()

data_dir <- 'C:\\Users\\Hana.Fu\\Documents'

for (i in 1:nrow(crudeTab)) {

  # Data subset
  cchs2013_14_Ontario_strata_sub <- cchs2013_14_Ontario_strata %>% filter(stratum == i)
  crudeTabSub <- crudeTab[i,variable.names]
  
  # Model input
  n.i   <- sum(round(cchs2013_14_Ontario_strata_sub$WTS_S))   # number of simulated individuals, since it is weighted sum it needs to be rounded to join back to CCHS data
  
  #round(sum(cchs2013_14_Ontario_strata_sub$WTS_S)) # previously

  print(n.i)

  n.t   <- 30                    # time horizon, 30 cycles
  v.n   <- c("H","S","D")        # the model states: Healthy (H), Sick (S), Dead (D)
  n.s   <- length(v.n)           # the number of health states
  v.M_1 <- rep("H", n.i)         # everyone begins in the healthy state 
  d.c   <- d.e <- 0.03           # equal discounting of costs and QALYs by 3%
  samplev.m <- 1                 # argument use for samplev() function 
  
  # Transition probabilities (per cycle)
  p.HD   <- transProbMat[i,'p.HD'] %>% unlist                # probability to die when healthy
  p.HS   <- transProbMat[i,'p.HS'] %>% unlist           	   # probability to become sick when healthy
  p.SH   <- transProbMat[i,'p.SH'] %>% unlist          	     # probability to become healthy when sick
  rr.S   <- transProbMat[i,'rr.S'] %>% unlist             	 # rate ratio of death when sick vs healthy
  r.HD   <- -log(1 - p.HD) 	     # rate of death when healthy 
  r.SD   <- rr.S * r.HD  	       # rate of death when sick
  p.SD   <- 1 - exp(- r.SD)      # probability to die when sick
  rp.S   <- transProbMat[i,'rp.S'] %>% unlist                # increase of the mortality rate with every additional year being sick
  
  # QALY inputs
  u.H     <- 1                   # utility when healthy 
  u.S     <- 0.75                # increase in utility of treated sick individuals with every additional year being sick
  v.x     <- runif(n.i, 0.95, 1.05) # vector capturing individuals' effect modifier at baseline
  
  # START SIMULATION
  cat('\n# Stratum:', paste(apply(crudeTabSub,2,as.character), collapse=':'),', n =',n.i, '\n')      # display the stratum
  
  p = Sys.time()
  sim <- MicroSim(v.M_1, n.i, n.t, v.n, samplev.m, X = v.x, d.e, seed = 100)  # run for treatment
  comp.time[[i]] = Sys.time() - p
  
  comp.time[[i]]
  
  # SAVE SIMULATION RESULTS TO LIST (DOES NOT WORK NOW, SINCE RUNNING OUT OF MEMORY)
  #sim.list[[i]] <- sim
  
  filefullpath <- file.path(data_dir,paste0('sim_',i,'.rds'))
  saveRDS(sim, filefullpath)
  
  rm(sim)
}



################################# Summary of microsimulation model #############################

library(tidyverse)

# Load results
i <- 1 # select stratum

filefullpath <- file.path(data_dir,paste0('sim_',i,'.rds'))
sim <- readRDS(filefullpath)

m.M <- sim$m.M
m.E <- sim$m.E
m.un <- sim$m.un
arr.p <- sim$arr.p

m.M.df <- data.frame(m.M)
m.M.df$ind <- as.numeric(gsub("ind ","",row.names(m.M.df)))
m.M.df.melt <- reshape2::melt(m.M.df, id.vars = "ind", variable.name = "cycle")
m.M.df.melt$cycle <- as.numeric(gsub("cycle.","",m.M.df.melt$cycle))
m.M.df.melt$state <- factor(m.M.df.melt$value, levels=c("H","S","D"))
m.M.df.melt.sub <- m.M.df.melt %>% filter(ind %in% 1:50)

m.un.df <- data.frame(m.un)
m.un.df$ind <- as.numeric(gsub("ind ","",row.names(m.un.df)))
m.un.df.melt <- reshape2::melt(m.un.df, id.vars = "ind", variable.name = "cycle")
m.un.df.melt$cycle <- as.numeric(gsub("cycle.","",m.un.df.melt$cycle))
m.un.df.melt.sub <- m.un.df.melt %>% filter(ind %in% 1:50)

m.E.df <- data.frame(m.E)
m.E.df$ind <- as.numeric(gsub("ind ","",row.names(m.E.df)))
m.E.df.melt <- reshape2::melt(m.E.df, id.vars = "ind", variable.name = "cycle")
m.E.df.melt$cycle <- as.numeric(gsub("cycle.","",m.E.df.melt$cycle))
m.E.df.melt.sub <- m.E.df.melt %>% filter(ind %in% 1:50)

# Plot heat maps
library(ggplot2)
ggplot(m.M.df.melt.sub, aes(cycle, ind)) + 
  geom_tile(aes(fill = state), colour = "white")  +   
  scale_fill_manual(values= c("darkgreen","purple","red"))

ggplot(m.E.df.melt.sub, aes(cycle, ind)) + 
  geom_tile(aes(fill = value)) 

ggplot(m.un.df.melt.sub, aes(cycle, ind)) + 
  geom_tile(aes(fill = value)) 

# Count by health state by cycle
m.M.df.agg <- m.M.df.melt %>%
  group_by(state, cycle) %>%
  summarise(count = n()) %>% 
  pivot_wider(names_from = state, values_from = count, values_fill = 0) 

kable(m.M.df.agg) %>%
  kable_styling("striped")


######################## Join results back to CCHS data ########################
############## Using first stratum as an example, due to data size #############

# Load results
i <- 1 # select stratum

filefullpath <- file.path(data_dir,paste0('sim_',i,'.rds'))
sim <- readRDS(filefullpath)

m.M <- sim$m.M
m.E <- sim$m.E
m.un <- sim$m.un
arr.p <- sim$arr.p

# Select data from the first stratum
cchs2013_14_Ontario_strata_sub <- cchs2013_14_Ontario_strata %>% filter(stratum == i)

# Survey weights for expanding CCHS data
cchs_wgts <- round(cchs2013_14_Ontario_strata_sub$WTS_S)
cchs_nrow <- length(cchs_wgts)
expand_wgts_list <- sapply(1:cchs_nrow, function (i) {
  rep(i, times=cchs_wgts[i])
})
expand_wgts <- unlist(expand_wgts_list)

# Expand data according to survey weights
cchs_sub_weighted <- cchs2013_14_Ontario_strata_sub[expand_wgts,
                         c('agecat','sex','educat','incpercat','inchhcat',
                           'ethncat','maritalcat','geourcat')]

# Calculate number of healthy and sick health states
num_hlth_sick_states <- apply(m.M, 1, 
      function(i) {
        table(i)[c("H","S")]
      }) 
num_hlth_sick_states <- data.frame(num_hlth_sick_states)
cchs_sub_weighted[c("num_hlth_state","num_sick_state")] <- t(num_hlth_sick_states)

# Time of death and Assign health state at the end of cycles
death_time_cycles_end_list <- apply(m.M, 1, 
                               function(i) {
                                 
                                 d_times <- grep("D",i)
                                 
                                 if (length(d_times) != 0) {
                                   
                                   min_pos <- min(grep("D",i))
                                   time_of_death <- min_pos 
                                 } else {
                                   time_of_death <- NA 
                                 }
                                 
                                 end_of_cycles <- tail(i, 1)
                                 
                                 return(
                                   list(time_of_death = time_of_death, end_of_cycles=end_of_cycles))
                               }) 
death_time_cycles_end <- do.call(rbind, death_time_cycles_end_list)
death_time_cycles_end <- data.frame(death_time_cycles_end)
cchs_sub_weighted[c('time_of_death','end_of_cycles')] <- death_time_cycles_end

# QALYs per individual
cchs_sub_weighted$qalys <- rowSums(m.E) 

# Show the results
DT::datatable(cchs_sub_weighted[1:100,], 
              colnames=c("Age group","Sex","Education","Income (personal)","Income (household)",
                         "Ethnicity","Marital Status","Urban/Rural", 
                         "Number of healthy state","Number of sick state",
                         'Time of death','End of cycles',
                         "QALYs"))

################################# Summary of QALYs #############################
# store the mean QALYs (and the MCSE) of each strategy in a new variable E (vector effects)
v.E  <-  sim.list[[1]]$te_hat
sd.E <-  sd(sim.list[[1]]$te) / sqrt(nrow(sim.list[[1]]$m.p))

