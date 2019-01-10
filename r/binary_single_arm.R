# Simulate Bayesian single-arm non-inferiority trial
# with early termination for inferiority/non-inferiority
# Model:
#  p         ~ beta(a, b)
#  x_k|p     ~ binomial(n_k, p) k = 1...K
#  p|x_{1:k} ~ beta(a + sum(x_{1:k}), b + n - sum(x_{1:k}))
#
# Hypothesis:
#   H_0: p > p_0 + d
#   H_1: p < p_0 + d
#
# Terimnal decision rule:
#   Pr(p < p_0 + d | x_{(1:k)}) > k_non => non-inferior
#   Pr(p < p_0 + d | x_{(1:k)}) < k_inf => inferior
#   Otherwise => inconclusive
#
# Interim decision rule:
#   Pr(p < p_0 + d | x_{(1:k)}) > k_non => non-inferior
#   Pr(p < p_0 + d | x_{(1:k)}) < k_inf => inferior
#   Otherwise => continue to k + 1
sim_trial <- function(
  sim_id = 1,
  ptru   = 0.08,
  p0     = 0.05,
  d      = 0.03,
  nint   = seq(100, 500, 50),
  a      = 1,
  b      = 1
) {
  
  n_analyses <- length(nint)

  # Generate data
  y <- cumsum(rbinom(length(nint), diff(c(0, nint)), ptru))
  
  # Simulate interim analyses
  sim_interim <- function(interim) {
    # Posterior parameters
    ycurr <- y[interim]
    ncurr <- nint[interim]
    apost <- a + ycurr
    bpost <- b + ncurr - ycurr
    # Posterior probability of H_1
    ptail <- pbeta(p0 + d, apost, bpost)
    return(list(n = ncurr, y = ycurr, a = apost, b = bpost, ptail = ptail))
  }
  # Iterate through all the interims
  intr <- lapply(1:n_analyses, sim_interim)
  # Collect the results
  do.call(rbind, intr)
  sim_trial <- data.frame(sim_id = sim_id, stage = 1:n_analyses, p0, d, ptru, do.call(rbind.data.frame, intr))
  return(sim_trial)
}


# Simulate a scenario - a collection of independent trials
# with the same set of parameters
sim_scenario <- function(sims, ...) {
  require(doParallel)
  res <- foreach(i = 1:sims, .packages = "rmutil", .export = "sim_trial") %dopar% {
    sim_trial(i, ...)
  }
  return(res)
  # sim_trials <- do.call(rbind.data.frame, res)
  # return(sim_trials, sim_results))
}

# Apply a decision rule to a simulated Bayesian single-arm trial
# Args:
#   sim_trial - the data.frame output from sim_trial()
#   k_non     - the non-inferior cut-point testing P_k > k_non
#   k_inf     - the inferior cut-point testing P_k < k_inf
trial_decision <- function(sim_trial, k_non = 0.975, k_inf = 0.025) {
  n_analyses <- nrow(sim_trial)
  
  # Perform input checks...
  if((length(k_non) > 1 && length(k_non) != n_analyses) ||
     (length(k_inf) > 1 && length(k_inf) != n_analyses))
    stop("Interim and cut-off dimensions do not match")
  
  if(length(k_non) == 1) k_non <- rep(k_non, n_analyses)
  if(length(k_inf) == 1) k_inf <- rep(k_inf, n_analyses)
  
  ptail <- sim_trial$ptail
  
  for(stage in 1:n_analyses) {
    est <- sim_trial[stage, "a"] / (sim_trial[stage, "a"] + sim_trial[stage, "b"])
    if(stage < n_analyses) {
      # Check non-inferiority
      if(ptail[stage] > k_non[stage]) {
        return(list(action = "stop", decision = "non-inferior", stage = stage, est = est))
      }
      # Check inferiority
      if(ptail[stage] < k_inf[stage]) {
        return(list(action = "stop", decision = "inferior", stage = stage, est = est))
      } 
    } else {
      # Check non-inferiority
      if(ptail[stage] > k_non[stage]) {
        return(list(action = "none", decision = "non-inferior", stage = stage, est = est))
      }
      # Check inferiority
      if(ptail[stage] <  k_inf[stage]) {
        return(list(action = "none", decision = "inferior", stage = stage, est = est))
      }
      return(list(action = "none", decision = "inconclusive", stage = stage, est = est))
    }
  }
}