# helpers_descriptives.R
# Kernel z-scoring + attrition tracking.

# ---- Kernel standardization --------------------------------------------
# Smooth kernel-weighted z-scoring. For each birth year t0, the local mean
# and SD of x are estimated with a Gaussian kernel K((t - t0) / h) where
# h = bandwidth.

kernel_z_score <- function(x, t, bandwidth = KERNEL_BW, min_eff_n = 10) {
  # x:         variable to standardize (e.g. parental_education)
  # t:         continuous index  (e.g. birth_year)
  # bandwidth: Gaussian kernel bandwidth (in years)
  # min_eff_n: minimum effective sample size for valid z-score

  cc <- !is.na(x) & !is.na(t)
  if (sum(cc) < min_eff_n) return(rep(NA_real_, length(x)))

  t_obs <- t[cc]
  x_obs <- x[cc]

  unique_t <- sort(unique(t[!is.na(t)]))

  mu_vec   <- setNames(rep(NA_real_, length(unique_t)), as.character(unique_t))
  sd_vec   <- mu_vec
  neff_vec <- mu_vec

  for (i in seq_along(unique_t)) {
    ti <- unique_t[i]
    w  <- dnorm((t_obs - ti) / bandwidth)
    w_sum <- sum(w)
    if (w_sum < 1e-10) next

    w_norm <- w / w_sum
    n_eff  <- 1 / sum(w_norm^2)           # Kish effective sample size
    if (n_eff < min_eff_n) next

    mu      <- sum(w_norm * x_obs)
    var_est <- sum(w_norm * (x_obs - mu)^2)
    sd_est  <- sqrt(var_est)

    mu_vec[i]   <- mu
    sd_vec[i]   <- sd_est
    neff_vec[i] <- n_eff
  }

  result <- rep(NA_real_, length(x))
  for (i in seq_along(unique_t)) {
    ti <- unique_t[i]
    if (is.na(sd_vec[i]) || sd_vec[i] < 1e-10) next
    idx <- which(t == ti & !is.na(x))
    if (length(idx) > 0) {
      result[idx] <- (x[idx] - mu_vec[i]) / sd_vec[i]
    }
  }
  result
}

# Diagnostic: kernel-smoothed mean, SD and n_eff at each unique t
kernel_moments <- function(x, t, bandwidth = KERNEL_BW) {
  cc <- !is.na(x) & !is.na(t)
  t_obs <- t[cc]; x_obs <- x[cc]
  unique_t <- sort(unique(t[!is.na(t)]))

  out <- data.frame(t = unique_t, mu = NA_real_, sd = NA_real_,
                    n_eff = NA_real_, n_raw = NA_integer_)
  for (i in seq_along(unique_t)) {
    ti <- unique_t[i]
    w  <- dnorm((t_obs - ti) / bandwidth)
    w_sum <- sum(w)
    if (w_sum < 1e-10) next
    w_norm <- w / w_sum
    mu <- sum(w_norm * x_obs)
    out$mu[i]    <- mu
    out$sd[i]    <- sqrt(sum(w_norm * (x_obs - mu)^2))
    out$n_eff[i] <- 1 / sum(w_norm^2)
    out$n_raw[i] <- sum(t[!is.na(x)] == ti)
  }
  out
}

compute_attrition_flow <- function(dat, steps, cohort_var = "cohort", pipeline = "Core") {
  cohorts <- sort(unique(as.character(dat[[cohort_var]])))
  cohorts <- cohorts[!is.na(cohorts)]

  bind_rows(lapply(cohorts, function(coh) {
    current <- dat[as.character(dat[[cohort_var]]) == coh, , drop = FALSE]
    bind_rows(lapply(seq_along(steps), function(i) {
      step <- steps[[i]]
      keep <- step$keep(current)
      keep[is.na(keep)] <- FALSE
      n_before <- nrow(current)
      n_after  <- sum(keep)
      out <- data.frame(
        pipeline = pipeline,
        cohort = coh,
        step_order = i,
        step = step$label,
        n_before = n_before,
        n_after = n_after,
        n_dropped = n_before - n_after,
        pct_dropped = if (n_before > 0) 100 * (n_before - n_after) / n_before else NA_real_,
        stringsAsFactors = FALSE
      )
      current <<- current[keep, , drop = FALSE]
      out
    }))
  }))
}
