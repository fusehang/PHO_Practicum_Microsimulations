##################################### Functions ###########################################

#### THE NEW samplev() FUNCTION
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

#### The MicroSim functions 
# for the simple two-stage microsimulation model keeps track of what happens to each individual during each cycle. 

## Version 1, "select" either alcoholic & non-alcoholic transition probabilities
MicroSimV1 <- function(v.M_1, n.i, n.t, v.n, samplev.m, X = NULL, d.e, TR.out = TRUE, TS.out = TRUE, seed = 1) {
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
  
  # Create matrix to store transition type
  m.type <- matrix(rep(NA, n.i*n.t), nrow = n.i, ncol = n.t,
                   dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                   paste("cycle", 1:n.t, sep = " "))) 
  
  # Create array to store transition probabilities
  arr.p <- array(NA, dim = c(n.i, n.s, n.t),
                 dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                 paste("state", 1:n.s, sep = " "),
                                 paste("cycle", 1:n.t, sep = " "))) 
  
  set.seed(seed)                  # set the seed for every individual for the random number generator
  
  m.E[, 1] <- Effs(m.M[, 1], X = X)  # estimate QALYs per individual for the initial health state
  
  # Loop through time points with transition probabilites and state updated at each time point
  for (t in 1:n.t) { # t <- 3
    
    p.type <- sample(c('p.HD1','p.HD2'), prob = c(p.AA, 1-p.AA), size = n.i, replace = TRUE) # Sample either to use alcoholic or non-alcoholic attributable probabilities
    p.HD <- sapply(p.type, get)
    
    m.p <- Probs(m.M[, t])         # calculate the transition probabilities at cycle t 
    
    arr.p[,, t] <-  m.p            # store the transition probabilities
    
    ran.un <- samplev(prob = m.p, m = samplev.m) # returns a matrix of the health states and a matrix of randomly generated uniform numbers
    
    m.type[,t] <- p.type             # store type of HD transition probability
    
    # this is Hana's guess of the function of the m argument in samplev() function
    if (samplev.m > 1) { 
      which.m <- sample(1:samplev.m, size = 1) # pick which m index to use 
      
      m.un[, t] <- ran.un[['un']][, which.m] # store the matrix of randomly generated uniform numbers for each cycle
      m.M[, t + 1] <- ran.un[['ran']][, which.m]  # sample the next health state and store that state in matrix m.M
    } else {
      m.un[, t] <- ran.un[['un']]
      m.M[, t + 1] <- ran.un[['ran']]
    }
    
    m.E[, t + 1] <- Effs( m.M[, t + 1], X = X)   # estimate QALYs per individual during cycle t + 1
    
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
  results <- list(m.p = m.p, arr.p = arr.p, m.M = m.M, m.E = m.E, te = te, te_hat = te_hat, 
                  m.un = m.un, m.type = m.type, TS = TS, TR = TR) # store the results from the simulation in a list  
  return(results)  # return the results
}  # end of the MicroSim function  


## Version 2, "combine" alcoholic & non-alcoholic transition probabilities

MicroSimV2 <- function(v.M_1, n.i, n.t, v.n, samplev.m, X = NULL, d.e, TR.out = TRUE, TS.out = TRUE, seed = 1) {
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
  
  # Create the matrix capturing the state name/health outcomes for all individuals at each time point; m.M1 is for alcoholic attributable mortality, m.M2 for non-alcoholic attr. mort., m.M combines m.M1 and m.M2
  m.M1 <- m.M2 <- m.M <- m.E <- matrix(nrow = n.i, ncol = n.t + 1, 
                                       dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                                       paste("cycle", 0:n.t, sep = " ")))  
  
  m.M1[, 1] <- m.M2[, 1] <- m.M[, 1] <- v.M_1                     # indicate the initial health state
  
  # Create matrix to store randomly generated uniform numbers, no.1 for alcoholic-attributable mortality, no.2 for non-alcoholic-attributable mortality
  m.un2 <- m.un1 <- matrix(rep(NA, n.i*n.t), nrow = n.i, ncol = n.t,
                           dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                           paste("cycle", 1:n.t, sep = " "))) # for alcoholic-attributable mortality
  
  # Create matrix to indicate the cause of death
  m.dtype <- matrix(NA, nrow = n.i, ncol = n.t+1,
                    dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                    paste("cycle", 0:n.t, sep = " "))) 
  
  # Create array to store transition probabilities, no.1 for alcoholic-attributable mortality, no.2 for non-alcoholic-attributable mortality
  arr.p2 <- arr.p1 <- array(NA, dim = c(n.i, n.s, n.t),
                            dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                            paste("state", 1:n.s, sep = " "),
                                            paste("cycle", 1:n.t, sep = " "))) 
  
  set.seed(seed)                  # set the seed for every individual for the random number generator
  
  m.E[, 1] <- Effs(m.M[, 1], X = X)  # estimate QALYs per individual for the initial health state
  
  # Loop through time points with transition probabilites and state updated at each time point
  for (t in 1:n.t) { # t <- 3
    #for (t in 1:5) {
    m.p1 <- ProbsV1(m.M[, t])      # calculate the transition probabilities at cycle t, for alcoholic attr mort
    m.p2 <- ProbsV2(m.M[, t])      # calculate the transition probabilities at cycle t, for non-alcoholic attr mort
    
    ran.un1 <- samplev(prob = m.p1, m = samplev.m) # returns a matrix of the health states and a matrix of randomly generated uniform numbers, for alcoholic attr mort
    ran.un2 <- samplev(prob = m.p2, m = samplev.m) # returns a matrix of the health states and a matrix of randomly generated uniform numbers, for non-alcoholic attr mort
    
    arr.p1[,, t] <-  m.p1            # store the transition probabilities, for alcoholic attr mort
    arr.p2[,, t] <-  m.p2            # store the transition probabilities, for non-alcoholic attr mort
    
    # some checks
    if (FALSE) {
      m.p1[56,]; ran.un1[['un']] [56]; ran.un1[['ran']][56]
      m.p2[56,]; ran.un2[['un']] [56]; ran.un2[['ran']][56]
    }
    
    # this is Hana's guess of the function of the m argument in samplev() function
    if (samplev.m > 1) { 
      which.m <- sample(1:samplev.m, size = 1) # pick which m index to use 
      
      m.un1[, t] <- ran.un1[['un']][, which.m] # store the matrix of randomly generated uniform numbers for each cycle, for alcoholic attr mort
      m.un2[, t] <- ran.un2[['un']][, which.m] # store the matrix of randomly generated uniform numbers for each cycle, for non-alcoholic attr mort
      
      m.M1[, t + 1] <- ran.un1[['ran']][, which.m]  # sample the next health state, for alcoholic attr mort
      m.M2[, t + 1] <- ran.un2[['ran']][, which.m]  # sample the next health state, for non-alcoholic attr mort
    } else {
      m.un1[, t] <- ran.un1[['un']] # store the matrix of randomly generated uniform numbers for each cycle, for alcoholic attr mort
      m.un2[, t] <- ran.un2[['un']] # store the matrix of randomly generated uniform numbers for each cycle, for non-alcoholic attr mort
      
      m.M1[, t + 1] <- ran.un1[['ran']]  # sample the next health state, for alcoholic attr mort
      m.M2[, t + 1] <- ran.un2[['ran']]  # sample the next health state, for non-alcoholic attr mort
    }
    
    m.M[, t + 1] <- ifelse(m.M1[, t + 1] == "H" & m.M2[, t + 1] == "H", "H", "D") # if both alcoholic and non-alcoholic are healthy, then return "H", else when either one is deceased, then return "D" # combine alcoholic and non-alcoholic states and store that state in matrix m.M
    
    m.dtype[, t+1] <- ifelse(!is.na(m.dtype[, t]), m.dtype[, t],
                             ifelse(m.M1[, t + 1] == "H" & m.M2[, t + 1] == "H", NA, 
                                    ifelse(m.M1[, t + 1] == "D" & m.M2[, t + 1] == "H", "Alcoholic",
                                           ifelse(m.M1[, t + 1] == "H" & m.M2[, t + 1] == "D", "Non-alcoholic", 
                                                  ifelse(m.M1[, t + 1] == "D" & m.M2[, t + 1] == "D", "Both", NA)))))
    
    m.E[, t + 1] <- Effs( m.M[, t + 1], X = X)   # estimate QALYs per individual during cycle t + 1
    
    cat('\r', paste(round(t/n.t * 100), "% done", sep = " "))       # display the progress of the simulation
    
  } # close the loop for the time points 
  
  # check m.M, m.M1, m.M2
  if (FALSE) {
    cbind(m.M[56,],m.M1[56,],m.M2[56,])
  }
  
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
  results <- list(arr.p1 = arr.p1, arr.p2 = arr.p2, m.M1 = m.M1, m.M2 = m.M2, m.M = m.M, 
                  m.dtype = m.dtype, m.E = m.E, te = te, te_hat = te_hat, 
                  m.un1 = m.un1, m.un2 = m.un2, TS = TS, TR = TR) # store the results from the simulation in a list  
  return(results)  # return the results
}  # end of the MicroSim function  


#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below.

Probs <- function(M_it) { 
  # M_it:   health state occupied by individual i at cycle t (character variable)
  
  m.p.it <- matrix(NA, n.s, n.i)     # create vector of state transition probabilities
  rownames(m.p.it) <- v.n            # assign names to the vector
  
  # update the v.p with the appropriate probabilities   
  m.p.it[,M_it == "H"]  <- c(1 - p.HD, p.HD)                  # transition probabilities when healthy
  m.p.it[,M_it == "D"]  <- c(0, 1)                            # transition probabilities when dead   
  #ifelse(colSums(m.p.it) == 1, return(t(m.p.it)), print("Probabilities do not sum to 1")) # return the transition probabilities or produce an error # not using this as == leads to error for some runs, while all.equal() solves the problem
  ifelse(sapply(colSums(m.p.it), function(x) all.equal(x,1)), return(t(m.p.it)), print("Probabilities do not sum to 1")) # return the transition probabilities or produce an error
}

# Probs function V1 to store alcohol attributable mortality transition probabilities
ProbsV1 <- function(M_it) { 
  # M_it:   health state occupied by individual i at cycle t (character variable)
  
  m.p1.it <- matrix(NA, n.s, n.i)     # create vector of state transition probabilities
  rownames(m.p1.it) <- v.n            # assign names to the vector
  
  # update the v.p with the appropriate probabilities   
  m.p1.it[,M_it == "H"]  <- c(1 - p.HD1, p.HD1)                # transition probabilities when healthy
  m.p1.it[,M_it == "D"]  <- c(0, 1)                            # transition probabilities when dead   
  
  ifelse(sapply(colSums(m.p1.it), function(x) all.equal(x,1)), return(t(m.p1.it)), print("Probabilities do not sum to 1")) # return the transition probabilities or produce an error
}

# Probs function V2 to store non-alcohol attributable mortality transition probabilities
ProbsV2 <- function(M_it) { 
  # M_it:   health state occupied by individual i at cycle t (character variable)
  
  m.p2.it <- matrix(NA, n.s, n.i)     # create vector of state transition probabilities
  rownames(m.p2.it) <- v.n            # assign names to the vector
  
  # update the v.p with the appropriate probabilities   
  m.p2.it[,M_it == "H"]  <- c(1 - p.HD2, p.HD2)                  # transition probabilities when healthy
  m.p2.it[,M_it == "D"]  <- c(0, 1)                            # transition probabilities when dead   
  ifelse(sapply(colSums(m.p2.it), function(x) all.equal(x,1)), return(t(m.p2.it)), print("Probabilities do not sum to 1")) # return the transition probabilities or produce an error
}


### Health outcome function 
# The Effs function to update the utilities at every cycle.

Effs <- function (M_it, cl = 1, X = NULL) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # cl:   cycle length (default is 1)
  # X:    the vector or matrix of individual characteristics (optional)
  
  u.it               <- 0        # by default the utility for everyone is zero
  u.it[M_it == "H"]  <- u.H      # update the utility if healthy
  u.it[M_it == "D"]  <- 0        # update the utility if dead
  
  QALYs <-  u.it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
}

