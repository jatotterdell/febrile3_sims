dbetabinom <- function(x, n, a = 1, b = 1){
  num <- lgamma(a + b) + lgamma(n + 1) + lgamma(x + a) + lgamma(n - x + b)
  den <- lgamma(a) + lgamma(b) + lgamma(x + 1) + lgamma(n - x + 1) + lgamma(n + a + b)
  prob <- exp(num - den)
  prob
}

hdiBeta <- function(p = 0.95, a, b) {
  if (a == b) {
    lower <- qbeta((1-p)/2, a, b)
    upper <- qbeta((1+p)/2, a, b)
  } else if ((a <= 1) & (b > 1)) {
    lower <- 0
    upper <- qbeta(p, a, b)
  } else if ((a > 1) & (b <= 1)) {
    lower <- qbeta(1-p, a, b)
    upper <- 1
  } else {
    zerofn <- function(x){(dbeta(qbeta(p+pbeta(x,a,b),a,b),a,b)-dbeta(x,a,b))^2}
    maxl <- qbeta(1-p,a,b)
    res <- optimize(zerofn,interval=c(0,maxl))
    lower <- res$minimum; upper=qbeta(p+pbeta(lower,a,b),a,b)
  }
  c(lower,upper)
}


fixed_trial <- function(
  ptru   = 0.08,
  p0     = 0.05,
  d      = 0.03,
  nmax   = 500,
  k_non  = 0.975,
  k_inf  = 0.025,
  a      = 1,
  b      = 1
) {
  n <- 0:nmax
  i_non <- which(pbeta(p0 + d, a + n, b + nmax - n) >= k_non) - 1
  i_inf <- which(pbeta(p0 + d, a + n, b + nmax - n) <= k_inf) - 1
  p_non <- sum(dbinom(i_non, nmax, ptru))
  p_inf <- sum(dbinom(i_inf, nmax, ptru))
  return(c(p_non, p_inf))
}


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
  x <- rbinom(length(nint), diff(c(0, nint)), ptru)
  y <- cumsum(x)
  
  # Simulate interim analyses
  sim_interim <- function(interim) {
    # Posterior parameters
    ycurr <- y[interim]
    ncurr <- nint[interim]
    apost <- a + ycurr
    bpost <- b + ncurr - ycurr
    # Posterior probability of H_1
    ptail <- pbeta(p0 + d, apost, bpost)
    return(list(n = ncurr, x = x[interim], y = ycurr, a = apost, b = bpost, ptail = ptail))
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
  sim_id <- sim_trial$sim_id[1]
  
  for(stage in 1:n_analyses) {
    n <- sim_trial[stage, "n"]
    est <- sim_trial[stage, "a"] / (sim_trial[stage, "a"] + sim_trial[stage, "b"])
    if(stage < n_analyses) {
      # Check non-inferiority
      if(ptail[stage] > k_non[stage]) {
        return(list(sim_id = sim_id, action = "stop", decision = "non-inferior", stage = stage, n = n, est = est, ptail = ptail[stage]))
      }
      # Check inferiority
      if(ptail[stage] < k_inf[stage]) {
        return(list(sim_id = sim_id, action = "stop", decision = "inferior", stage = stage, n = n, est = est, ptail = ptail[stage]))
      } 
    } else {
      # Check non-inferiority
      if(ptail[stage] > k_non[stage]) {
        return(list(sim_id = sim_id, action = "none", decision = "non-inferior", stage = stage, n = n, est = est, ptail = ptail[stage]))
      }
      # Check inferiority
      if(ptail[stage] <  k_inf[stage]) {
        return(list(sim_id = sim_id, action = "none", decision = "inferior", stage = stage, n = n, est = est, ptail = ptail[stage]))
      }
      return(list(sim_id = sim_id, action = "none", decision = "inconclusive", stage = stage, n = n, est = est, ptail = ptail[stage]))
    }
  }
}

trial_decision_pp <- function(sim_trial, k_non = 0.975, k_inf = 0.025, k_fut = 0.01, k_suc = 0.99) {
  n_analyses <- nrow(sim_trial)
  sim_id <- sim_trial$sim_id[1]
  ptail <- sim_trial$ptail
  p0 <- sim_trial$p0
  d <- sim_trial$d
  a <- sim_trial$a
  b <- sim_trial$b
  n <- sim_trial$n
  n_max <- max(n)
  
  for(stage in 1:n_analyses) {
    est <- a[stage] / (a[stage] + b[stage])
    if(stage < n_analyses) {
      n_rem <- n_max - n[stage]
      y_suc <- which(pbeta(p0 + d, a[stage] + 0:n_rem, b[stage] + n_rem - 0:n_rem) > k_non) - 1
      ppos <- sum(dbetabinom(y_suc, n_rem, a[stage], b[stage]))
      if(ppos < k_fut) {
        return(list(sim_id = sim_id, action = "stop", decision = "futile", stage = stage, n = n[stage], est = est, ptail = ptail[stage], ppos = ppos))
      }
      if(ppos > k_suc) {
        return(list(sim_id = sim_id, action = "stop", decision = "predict success", stage = stage, n = n[stage], est = est, ptail = ptail[stage], ppos = ppos))
      }
    } else {
      # Check non-inferiority
      if(ptail[stage] > k_non) {
        return(list(sim_id = sim_id, action = "none", decision = "non-inferior", stage = stage, n = n[stage], est = est, ptail = ptail[stage], ppos = NA))
      }
      # Check inferiority
      if(ptail[stage] <  k_inf) {
        return(list(sim_id = sim_id, action = "none", decision = "inferior", stage = stage, n = n[stage], est = est, ptail = ptail[stage], ppos = NA))
      }
      return(list(sim_id = sim_id, action = "none", decision = "inconclusive", stage = stage, n = n[stage], est = est, ptail = ptail[stage], ppos = NA))
    }
  }
}

decision_final <- function(sim_scenario, k_non = 0.975, k_inf = 0.025) {
  df <- do.call(rbind.data.frame, lapply(sims, "[", 9, ))
  df$final_decision <- ifelse(df$ptail > k_non, "non-inferior",
                        ifelse(df$ptail < k_inf, "inferior", "inconclusive"))
  return(df[, c("sim_id", "ptail", "final_decision")])
}



### TESTING ###


# sims <- sim_scenario(1e4, ptru = 0.05)
# dec1 <- do.call(rbind.data.frame, lapply(sims, trial_decision))
# dec2 <- do.call(rbind.data.frame, lapply(sims, trial_decision_pp))
# dec3 <- decision_final(sims)
# dec <- left_join(left_join(dec1, dec2, by = "sim_id"), dec3, by = "sim_id")

# library(dplyr)
# 
# sce <- sim_scenario(1000, ptru = 0.08)
# 
# # Apply decision rules and record frequentist results
# 
# get_results <- function(case, k_non, k_inf) {
#   dec <- cbind(
#     trial_id = 1:length(sce),
#     do.call(rbind.data.frame,
#             lapply(sce, 
#                    trial_decision, 
#                    k_non = k_non,
#                    k_inf = k_inf)))
#   dec_df <- dec %>% 
#     group_by(stage) %>% 
#     summarise(y_non = sum(decision == "non-inferior"),
#               y_inf = sum(decision == "inferior")) %>% 
#     ungroup() %>% 
#     mutate(case = case,
#            y_non_c = cumsum(y_non),
#            y_inf_c = cumsum(y_inf), 
#            p_non_c = y_non_c / 1000,
#            p_inf_c = y_inf_c / 1000) %>%
#     select(case, stage, y_non_c, y_inf_c, p_non_c, p_inf_c)
#   return(bind_rows(tibble(case = case, stage = 0, y_non_c = 0, y_inf_c = 0, p_non_c = 0, p_inf_c = 0), dec_df))
# }
# 
# k_non <- seq(0.90, 0.99, 0.01)
# out <- do.call(rbind.data.frame, lapply(1:length(k_non), function(i) get_results(i, k_non[i], 0.025)))
# ggplot(out, aes(stage, p_non_c, colour = factor(case), group = case)) + 
#   geom_line()
# 
# sce <- sim_scenario(1000, ptru = 0.1)
# 
# 
# dec1 <- cbind(
#   trial_id = 1:length(sce),
#   do.call(rbind.data.frame,
#           lapply(sce, 
#                  trial_decision, 
#                  k_non = seq(0.99, 0.975, length.out = 9),
#                  k_inf = seq(0.01, 0.025, length.out = 9))))
# dec2 <- cbind(
#   trial_id = 1:length(sce),
#   do.call(rbind.data.frame,
#           lapply(sce, 
#                  trial_decision, 
#                  k_non = rep(0.975, 9),
#                  k_inf = rep(0.05, 9))))
# 
# prop.table(table(dec1$decision))
# prop.table(table(dec2$decision))
# 
# cp <- dec2 %>% 
#   group_by(stage) %>% 
#   summarise(y_non = sum(decision == "non-inferior"),
#             y_inf = sum(decision == "inferior")) %>% 
#   ungroup() %>% 
#   mutate(y_non_c = cumsum(y_non),
#          y_inf_c = cumsum(y_inf), 
#          p_non_c = y_non_c / 1000,
#          p_inf_c = y_inf_c / 1000)
# 
# 
# dec1 <- cbind(trial_id = 1:1000, do.call(rbind.data.frame, lapply(sce, trial_decision)))
# dec2 <- cbind(trial_id = 1:1000, do.call(rbind.data.frame, lapply(sce, trial_decision_pp)))
# 
# # Average sample size
# mean(do.call(c, lapply(1:1000, function(i) sce[[i]][dec1$stage[i], "n"])))
# mean(do.call(c, lapply(1:1000, function(i) sce[[i]][dec2$stage[i], "n"])))
# 
# # How many predicited success trials did not end in success?
# mean(do.call(rbind, lapply(sce[dec2[dec2$decision == "predict success", "trial_id"]], "[", 9, "ptail", drop = F)) > 0.975)
# # How many trials stopped for futility would have ended in success?
# mean(do.call(rbind, lapply(sce[dec2[dec2$decision == "futile", "trial_id"]], "[", 9, "ptail", drop = F)) > 0.975)
