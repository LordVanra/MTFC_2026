################################################################################
# ADVANCED EXTENSIONS FOR AV ADOPTION MODEL
#
# This file contains:
# 1. Parameter estimation from data
# 2. State-space matrix/Jacobian formulation
# 3. Extended Kalman Filter (EKF) implementation
# 4. Sensitivity analysis
# 5. Policy optimization
################################################################################

# Load main model first
source("av_adoption_model.R")


# ==============================================================================
# 1. PARAMETER ESTIMATION FROM DATA
# ==============================================================================

#' Estimate adoption model parameters from historical data
#' 
#' Uses logistic regression to estimate beta coefficients from
#' observed adoption rates, prices, and infrastructure
#' 
#' @param data data.frame with columns: adoption, price_adv, infrastructure
#' @return list of estimated parameters
estimate_parameters <- function(data) {
  # Fit logistic regression
  # logit(a) = beta_0 + beta_r * price_adv + beta_i * I
  
  model <- glm(
    adoption ~ price_adv + infrastructure,
    data = data,
    family = binomial(link = "logit")
  )
  
  coefs <- coef(model)
  
  list(
    beta_0 = coefs[1],
    beta_r = coefs[2],
    beta_i = coefs[3],
    model = model
  )
}

#' Generate synthetic data for testing parameter estimation
generate_synthetic_data <- function(n = 100, params) {
  set.seed(123)
  
  # Random policy values
  price_adv <- runif(n, 0, 0.3)
  infrastructure <- runif(n, 0, 1)
  
  # True adoption
  z <- params$beta_0 + params$beta_r * price_adv + params$beta_i * infrastructure
  a_true <- sigmoid(z)
  
  # Add noise
  adoption <- pmin(pmax(a_true + rnorm(n, 0, 0.05), 0), 1)
  
  data.frame(
    adoption = adoption,
    price_adv = price_adv,
    infrastructure = infrastructure
  )
}

# Example usage
cat("\n=== PARAMETER ESTIMATION EXAMPLE ===\n")
synthetic_data <- generate_synthetic_data(200, params)
est_params <- estimate_parameters(synthetic_data)

cat("\nTrue vs Estimated Parameters:\n")
cat(sprintf("beta_0: True = %.3f, Est = %.3f\n", params$beta_0, est_params$beta_0))
cat(sprintf("beta_r: True = %.3f, Est = %.3f\n", params$beta_r, est_params$beta_r))
cat(sprintf("beta_i: True = %.3f, Est = %.3f\n", params$beta_i, est_params$beta_i))


# ==============================================================================
# 2. STATE-SPACE MATRIX FORMULATION
# ==============================================================================

#' Compute Jacobian matrix of state transition function
#' 
#' Linearization around current state for EKF
#' 
#' @param state current state vector
#' @param policy current policy
#' @param params model parameters
#' @return 3x3 Jacobian matrix
compute_jacobian <- function(state, policy, params) {
  A <- state["A"]
  C <- state["C"]
  I <- state["I"]
  
  r <- policy["r"]
  s <- policy["s"]
  i <- policy["i"]
  
  # Compute intermediate values
  price_adv <- r / params$P_AV
  z <- params$beta_0 + params$beta_r * price_adv + params$beta_i * I
  sig <- sigmoid(z)
  cap <- params$cap_0 + params$gamma * s
  a <- min(sig, cap)
  
  # Jacobian elements
  J <- matrix(0, nrow = 3, ncol = 3)
  
  # ∂A_next/∂A
  J[1, 1] <- 1 - params$delta
  
  # ∂A_next/∂I (if adoption not supply-constrained)
  if (sig < cap) {
    dsig_dz <- sig * (1 - sig)
    dz_dI <- params$beta_i
    J[1, 3] <- params$S * dsig_dz * dz_dI
  }
  
  # ∂C_next/∂C
  J[2, 2] <- 1 - params$delta
  
  # ∂I_next/∂I
  J[3, 3] <- 1 - params$rho
  
  return(J)
}

#' State transition in matrix form (linear approximation)
#' 
#' x_{t+1} ≈ F_t * x_t + B_t * u_t + c_t
#' 
#' @return list with matrices F, B, and vector c
linearize_model <- function(state, policy, params) {
  # Jacobian (state transition matrix)
  F <- compute_jacobian(state, policy, params)
  
  # Control matrix (simplified - assumes local linearity)
  B <- matrix(0, nrow = 3, ncol = 3)
  B[3, 3] <- 1  # infrastructure input
  
  # Constant term
  next_state <- step_model(state, policy, params)
  c <- next_state[1:3] - F %*% state - B %*% policy
  
  list(F = F, B = B, c = as.vector(c))
}


# ==============================================================================
# 3. EXTENDED KALMAN FILTER
# ==============================================================================

#' Extended Kalman Filter for state estimation with noisy observations
#' 
#' @param observations matrix of noisy state observations (3 x T)
#' @param policy policy trajectory
#' @param params model parameters
#' @param Q process noise covariance
#' @param R observation noise covariance
#' @return filtered state estimates
extended_kalman_filter <- function(observations, policy, params, 
                                   Q = diag(c(1e8, 1e8, 0.01)),
                                   R = diag(c(5e8, 5e8, 0.05))) {
  T <- ncol(observations)
  n <- 3  # state dimension
  
  # Initialize
  x_est <- matrix(NA, nrow = n, ncol = T)
  P <- diag(c(1e10, 1e10, 0.1))  # initial covariance
  
  x_est[, 1] <- observations[, 1]
  
  # Kalman filter loop
  for (t in 2:T) {
    # Get policy for this step
    u_t <- policy[t - 1, ]
    
    # Predict step
    x_pred <- step_model(x_est[, t-1], u_t, params)[1:3]
    F_t <- compute_jacobian(x_est[, t-1], u_t, params)
    P_pred <- F_t %*% P %*% t(F_t) + Q
    
    # Update step
    H <- diag(3)  # observation matrix (directly observe states)
    y <- observations[, t] - H %*% x_pred  # innovation
    S <- H %*% P_pred %*% t(H) + R  # innovation covariance
    K <- P_pred %*% t(H) %*% solve(S)  # Kalman gain
    
    x_est[, t] <- x_pred + K %*% y
    P <- (diag(n) - K %*% H) %*% P_pred
  }
  
  return(t(x_est))
}


# ==============================================================================
# 4. SENSITIVITY ANALYSIS
# ==============================================================================

#' One-at-a-time sensitivity analysis
#' 
#' @param state_0 initial state
#' @param policy baseline policy
#' @param years simulation horizon
#' @param params baseline parameters
#' @param param_name parameter to vary
#' @param range variation range (fraction)
#' @return data.frame with sensitivity results
sensitivity_analysis <- function(state_0, policy, years, params, 
                                param_name = "beta_r", range = c(0.5, 1.5)) {
  n_points <- 20
  multipliers <- seq(range[1], range[2], length.out = n_points)
  
  results <- data.frame(
    multiplier = multipliers,
    final_share = NA_real_
  )
  
  for (i in seq_along(multipliers)) {
    params_temp <- params
    params_temp[[param_name]] <- params[[param_name]] * multipliers[i]
    
    sim <- simulate_model(state_0, policy, years, params_temp)
    results$final_share[i] <- sim$pct_AV[years + 1]
  }
  
  return(results)
}

#' Tornado diagram for multi-parameter sensitivity
plot_tornado <- function(state_0, policy, years, params, 
                         param_names = c("beta_0", "beta_r", "beta_i", "delta")) {
  
  impacts <- data.frame(
    parameter = param_names,
    low = NA_real_,
    high = NA_real_
  )
  
  # Baseline
  baseline <- simulate_model(state_0, policy, years, params)
  baseline_share <- baseline$pct_AV[years + 1]
  
  for (i in seq_along(param_names)) {
    # Low case (0.8x)
    params_low <- params
    params_low[[param_names[i]]] <- params[[param_names[i]]] * 0.8
    sim_low <- simulate_model(state_0, policy, years, params_low)
    impacts$low[i] <- sim_low$pct_AV[years + 1] - baseline_share
    
    # High case (1.2x)
    params_high <- params
    params_high[[param_names[i]]] <- params[[param_names[i]]] * 1.2
    sim_high <- simulate_model(state_0, policy, years, params_high)
    impacts$high[i] <- sim_high$pct_AV[years + 1] - baseline_share
  }
  
  # Sort by total range
  impacts$range <- abs(impacts$high - impacts$low)
  impacts <- impacts[order(impacts$range, decreasing = TRUE), ]
  
  # Plot
  par(mar = c(4, 6, 3, 2))
  n <- nrow(impacts)
  plot(NULL, xlim = c(min(impacts$low), max(impacts$high)) * 100,
       ylim = c(0.5, n + 0.5), xlab = "Change in Final AV Share (pp)",
       ylab = "", yaxt = "n", main = "Tornado Diagram: Parameter Sensitivity")
  
  for (i in 1:n) {
    segments(impacts$low[i] * 100, n - i + 1, 
             impacts$high[i] * 100, n - i + 1, 
             lwd = 6, col = "#2E86AB")
  }
  
  abline(v = 0, lty = 2, col = "red")
  axis(2, at = 1:n, labels = rev(impacts$parameter), las = 1)
  grid(col = "gray80", lty = 2)
}


# ==============================================================================
# 5. POLICY OPTIMIZATION
# ==============================================================================

#' Simple policy optimizer (grid search)
#' 
#' Finds optimal constant policy within budget constraint
#' 
#' @param state_0 initial state
#' @param years horizon
#' @param params model parameters
#' @param budget total annual budget (billion $)
#' @param objective "share" or "speed"
#' @return optimal policy
optimize_policy <- function(state_0, years, params, 
                           budget = 100e9, objective = "share") {
  
  # Grid search over policy space
  r_vals <- seq(0, 15000, by = 2500)
  s_vals <- seq(0, 8000, by = 2000)
  i_vals <- seq(0, 0.20, by = 0.05)
  
  best_obj <- -Inf
  best_policy <- NULL
  
  for (r in r_vals) {
    for (s in s_vals) {
      for (i in i_vals) {
        # Check budget constraint (simplified)
        cost <- r * params$S + s * params$S + i * 1e11
        if (cost > budget) next
        
        policy <- c(r = r, s = s, i = i)
        sim <- simulate_model(state_0, policy, years, params)
        
        if (objective == "share") {
          obj <- sim$pct_AV[years + 1]
        } else if (objective == "speed") {
          obj <- mean(diff(sim$pct_AV))
        }
        
        if (obj > best_obj) {
          best_obj <- obj
          best_policy <- policy
        }
      }
    }
  }
  
  list(policy = best_policy, objective = best_obj)
}


# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

cat("\n=== SENSITIVITY ANALYSIS ===\n")
state_0 <- c(A = 10e6, C = 250e6, I = 0.30)
policy <- c(r = 5000, s = 2000, i = 0.05)

sens_beta_r <- sensitivity_analysis(state_0, policy, 30, params, "beta_r")

plot(sens_beta_r$multiplier, sens_beta_r$final_share * 100,
     type = "l", lwd = 2, col = "#2E86AB",
     xlab = "beta_r Multiplier", ylab = "Final AV Share (%)",
     main = "Sensitivity to Rebate Effectiveness")
grid(col = "gray80", lty = 2)
abline(v = 1, lty = 2, col = "red")

cat("\n=== TORNADO DIAGRAM ===\n")
plot_tornado(state_0, policy, 30, params)

cat("\n=== POLICY OPTIMIZATION ===\n")
opt <- optimize_policy(state_0, 30, params, budget = 150e9)
cat("\nOptimal policy:\n")
cat(sprintf("  Rebate:      $%,.0f\n", opt$policy["r"]))
cat(sprintf("  Subsidy:     $%,.0f\n", opt$policy["s"]))
cat(sprintf("  Infra:       %.1f%%\n", opt$policy["i"] * 100))
cat(sprintf("  Final share: %.1f%%\n", opt$objective * 100))

cat("\n=== EXTENSIONS COMPLETE ===\n")