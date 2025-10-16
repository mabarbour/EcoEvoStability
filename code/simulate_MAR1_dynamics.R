
EcoEvo_MAR1_Dynamics <- function(
    Additions, # matrix of additions to system (usually just initial conditions)
    X.vector, # vector of intrinsic growth rates for abundances and baseline selection gradients on traits
    J.matrix, # matrix of eco-evolutionary effects. Important to add identity matrix before running simulation!
    Noise.Additions = NULL, # Flexibility to add stochasticity in intrinsic growth rates or selection, as inferred from the model
    Abund.Positions, # which positions of the X.vector correspond to abundances
    Freq.Positions, # which positions of the X.vector correspond to traits (frequency)
    Duration) # duration of simulation
  {
  
  # specify positions for abundance and trait values
  Npos <- Abund.Positions
  Zpos <- Freq.Positions
  
  # generate noise matrix
  if (!is.null(Noise.Additions)) {
    Noise.Matrix <- Noise.Additions 
  } else {
    Noise.Matrix <- matrix(0, nrow = length(X.vector), ncol = Duration)
  }
  
  # create empty list simulation
  Xt <- list()
  fixed <- rep(FALSE, length(X.vector))  # Track fixation status
  
  # for-loop, which also applies constraints on abundance (if < 0, set to zero) and frequencies (if < 0, set to zero, if > 1, set to 1)
  for (i in 1:Duration) {
    if (i == 1) {
      Xt[[i]] <- Additions[, 1]
    } else {
      prev <- Xt[[i - 1]]
      
      # Apply extinction and fixation constraints to previous state
      prev[Npos][prev[Npos] <= 0] <- 0
      prev[Zpos][prev[Zpos] <= 0] <- 0
      prev[Zpos][prev[Zpos] >= 1] <- 1
      
      # Update fixation status
      fixed[Zpos] <- fixed[Zpos] | (prev[Zpos] == 0 | prev[Zpos] == 1)
      
      # Update state
      updated <- X.vector + J.matrix %*% prev + Additions[, i] + Noise.Matrix[, i]
      updated <- as.vector(updated)
      
      # Enforce extinction and fixation
      updated[Npos][prev[Npos] == 0] <- 0
      updated[Zpos][fixed[Zpos]] <- prev[Zpos][fixed[Zpos]]
      
      # Apply bounds again after update
      updated[Npos][updated[Npos] <= 0] <- 0
      updated[Zpos][updated[Zpos] <= 0] <- 0
      updated[Zpos][updated[Zpos] >= 1] <- 1
      
      Xt[[i]] <- updated
      
      # Stop simulation if any species in Npos goes extinct
      if (any(updated[Npos] == 0)) {
        break
      }
    }
  }
  
  # organize and return data after simulation
  Xt.df <- t(as.data.frame(Xt)) %>%
    as.data.frame()
  colnames(Xt.df) <- rownames(X.vector)
  Xt.df <- Xt.df %>%
    tibble::as_tibble() %>%
    mutate(Week = 0:(nrow(.) - 1))
  return(Xt.df)
}
