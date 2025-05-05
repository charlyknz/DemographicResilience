# Model simulations ------------------------------------

simulate_disturbed_pop <- function(mat, rlwr = NULL, rupr = NULL,
                                   disturb_time = 50, total_time = 200,
                                   case = "random") {
  A_original <- mat$matrix_A
  n_stages <- nrow(A_original)
  init <- rep(1, n_stages)
  
  # Set up population matrix
  N <- matrix(NA, nrow = n_stages, ncol = total_time + 1)
  N[, 1] <- init
  
  for (t in 1:total_time) {
    A <- A_original  # reset A each timestep (to undo compensation effects)
    
    # Apply disturbance only at the defined time step (both to population and matrix)
    if (t == disturb_time && case != "control") {
      # Modify population vector
      if (case == "random") {
       # drop_frac <- runif(1, 0.1, 0.9)
        drop_frac<-0.5
      } else if (case == "resilience" && !is.null(rlwr)) {
        drop_frac <- (1 - rlwr)
        drop_frac <- min(drop_frac, 1) #  keeps 10%  of individuals
      } else {
        drop_frac <- 0.5
      }
      
      # Apply population crash
      N[, t] <- N[, t] * drop_frac
      
      # Apply structural damage to the matrix
      A <- A * drop_frac  # reduce all transitions to 50%
    }
    
    # Compensation phase (short-term boost in survival/recovery)
    if (case == "resilience" && !is.null(rupr) && t %in% (disturb_time + 1):(disturb_time + 5)) {
      A <- A + diag(rep(rupr / 10, n_stages))
    }
    
    # Project next time step
    N[, t + 1] <- A %*% N[, t]
  }
  
  tibble(
    time = rep(0:total_time, each = n_stages),
    stage = rep(1:n_stages, times = total_time + 1),
    abundance = as.vector(N),
    total = rep(colSums(N), each = n_stages)
  )
}
