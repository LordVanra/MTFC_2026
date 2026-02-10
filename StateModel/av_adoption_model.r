################################################################################
# AV ADOPTION STATE-SPACE MODEL WITH POLICY LEVERS
# 
# This script implements a discrete-time state-space model for autonomous
# vehicle (AV) adoption with three policy instruments:
#   - r: consumer rebates ($)
#   - s: manufacturer subsidies ($)
#   - i: infrastructure investment (fraction)
#
# State variables:
#   - A: AV fleet size
#   - C: total vehicle fleet size
#   - I: infrastructure readiness (0-1)
#
# SETUP INSTRUCTIONS:
# 1. Save this file as av_adoption_model.R
# 2. Open R or RStudio
# 3. Run: source("av_adoption_model.R")
# 4. No packages required - uses base R only
#
# The script will:
# - Run baseline simulation (30 years)
# - Compare 5 policy scenarios
# - Show time-varying policy phase-out
# - Run Monte Carlo uncertainty analysis
# - Export results to CSV files
################################################################################

# --- 1. MODEL PARAMETERS (dummy values for demonstration) ---

params <- list(
  # Fleet dynamics
  delta = 0.07,        # annual vehicle retirement rate
  rho = 0.05,          # infrastructure decay rate
  S = 15e6,            # annual new car sales
  P_AV = 40000,        # average AV price ($)
  
  # Adoption model (logistic regression coefficients)
  beta_0 = -2.0,       # baseline (intercept)
  beta_r = 3.0,        # rebate sensitivity
  beta_i = 4.0,        # infrastructure sensitivity
  
  # Supply-side constraints
  cap_0 = 0.20,        # baseline production capacity (20% of sales)
  gamma = 0.00002      # subsidy effectiveness on capacity
)


# --- 2. HELPER FUNCTIONS ---

#' Sigmoid (logistic) function
#' @param z numeric vector
#' @return values in (0, 1)
sigmoid <- function(z) {
  1 / (1 + exp(-z))
}


# --- 3. CORE STATE-SPACE TRANSITION FUNCTION ---

#' One-step state update: x_{t+1} = f(x_t, u_t)
#' 
#' @param state named vector with A, C, I
#' @param policy named vector with r, s, i
#' @param params list of model parameters
#' @return named vector: next state + adoption fraction
step_model <- function(state, policy, params) {
  # Unpack state
  A <- state["A"]  # AV count
  C <- state["C"]  # total cars
  I <- state["I"]  # infrastructure readiness
  
  # Unpack policy
  r <- policy["r"]  # rebates
  s <- policy["s"]  # subsidies
  i <- policy["i"]  # infrastructure spending
  
  # --- Infrastructure dynamics ---
  I_next <- (1 - params$rho) * I + i
  I_next <- min(max(I_next, 0), 1)  # clamp to [0, 1]
  
  # --- Adoption fraction (demand-side) ---
  price_adv <- r / params$P_AV
  z <- params$beta_0 + params$beta_r * price_adv + params$beta_i * I
  a_uncapped <- sigmoid(z)
  
  # --- Supply constraint ---
  cap <- params$cap_0 + params$gamma * s
  a <- min(a_uncapped, cap)
  
  # --- Fleet dynamics ---
  A_next <- (1 - params$delta) * A + a * params$S
  C_next <- (1 - params$delta) * C + params$S
  
  # Return next state and diagnostics
  c(
    A = A_next,
    C = C_next,
    I = I_next,
    adoption = a,
    a_uncapped = a_uncapped,
    cap = cap
  )
}


# --- 4. SIMULATION FUNCTION ---

#' Simulate AV adoption over multiple years
#' 
#' @param state_0 initial state (named vector)
#' @param policy policy trajectory (data.frame or named vector for constant policy)
#' @param years number of years to simulate
#' @param params model parameters
#' @return data.frame with year-by-year results
simulate_model <- function(state_0, policy, years, params) {
  # Initialize results storage
  results <- data.frame(
    year = 0:years,
    A = NA_real_,
    C = NA_real_,
    I = NA_real_,
    pct_AV = NA_real_,
    adoption = NA_real_,
    a_uncapped = NA_real_,
    cap = NA_real_,
    r = NA_real_,
    s = NA_real_,
    i = NA_real_
  )
  
  # Store initial state
  results[1, c("A", "C", "I")] <- state_0[c("A", "C", "I")]
  results$pct_AV[1] <- state_0["A"] / state_0["C"]
  
  # If policy is constant, replicate it
  if (is.vector(policy)) {
    policy_df <- as.data.frame(t(replicate(years, policy)))
  } else {
    policy_df <- policy
  }
  
  # Store initial policy
  results[1, c("r", "s", "i")] <- policy_df[1, c("r", "s", "i")]
  
  # Current state
  state <- state_0
  
  # Simulation loop
  for (t in 1:years) {
    # Get policy for this period
    u_t <- as.numeric(policy_df[t, c("r", "s", "i")])
    names(u_t) <- c("r", "s", "i")
    
    # Update state
    next_state <- step_model(state, u_t, params)
    
    # Store results
    results[t + 1, c("A", "C", "I")] <- next_state[c("A", "C", "I")]
    results[t + 1, "pct_AV"] <- next_state["A"] / next_state["C"]
    results[t + 1, "adoption"] <- next_state["adoption"]
    results[t + 1, "a_uncapped"] <- next_state["a_uncapped"]
    results[t + 1, "cap"] <- next_state["cap"]
    results[t + 1, c("r", "s", "i")] <- u_t
    
    # Advance state
    state[c("A", "C", "I")] <- next_state[c("A", "C", "I")]
  }
  
  return(results)
}


# --- 5. VISUALIZATION FUNCTIONS ---

# Plot AV penetration over time
plot_av_share <- function(results, title = "AV Market Share Over Time") {
  plot(
    results$year,
    results$pct_AV * 100,
    type = "l",
    lwd = 2.5,
    col = "#2E86AB",
    xlab = "Year",
    ylab = "AV Share (%)",
    main = title,
    las = 1,
    ylim = c(0, max(results$pct_AV * 100) * 1.1)
  )
  grid(col = "gray80", lty = 2)
}

#' Plot all state variables
plot_states <- function(results) {
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  # AV penetration
  plot(results$year, results$pct_AV * 100, type = "l", lwd = 2, col = "#2E86AB",
       xlab = "Year", ylab = "AV Share (%)", main = "Market Penetration")
  grid(col = "gray80", lty = 2)
  
  # Infrastructure
  plot(results$year, results$I * 100, type = "l", lwd = 2, col = "#A23B72",
       xlab = "Year", ylab = "Infrastructure (%)", main = "Infrastructure Readiness")
  grid(col = "gray80", lty = 2)
  
  # Adoption rate
  plot(results$year, results$adoption * 100, type = "l", lwd = 2, col = "#F18F01",
       xlab = "Year", ylab = "Adoption Rate (%)", main = "Annual AV Adoption Rate")
  lines(results$year, results$cap * 100, lty = 2, col = "red")
  legend("topright", c("Actual", "Capacity"), lty = c(1, 2), col = c("#F18F01", "red"), bty = "n")
  grid(col = "gray80", lty = 2)
  
  # Fleet size
  plot(results$year, results$A / 1e6, type = "l", lwd = 2, col = "#6A994E",
       xlab = "Year", ylab = "AVs (millions)", main = "AV Fleet Size")
  grid(col = "gray80", lty = 2)
  
  par(mfrow = c(1, 1))
}

#' Compare multiple policy scenarios
plot_scenarios <- function(scenarios, labels = NULL) {
  if (is.null(labels)) {
    labels <- paste("Scenario", seq_along(scenarios))
  }
  
  colors <- c("#2E86AB", "#A23B72", "#F18F01", "#6A994E", "#C73E1D", "#7209B7")
  
  plot(
    scenarios[[1]]$year,
    scenarios[[1]]$pct_AV * 100,
    type = "l",
    lwd = 2.5,
    col = colors[1],
    xlab = "Year",
    ylab = "AV Share (%)",
    main = "Policy Scenario Comparison",
    las = 1,
    ylim = c(0, max(sapply(scenarios, function(s) max(s$pct_AV))) * 100 * 1.1)
  )
  
  for (i in 2:length(scenarios)) {
    lines(scenarios[[i]]$year, scenarios[[i]]$pct_AV * 100, 
          lwd = 2.5, col = colors[i %% length(colors) + 1])
  }
  
  legend("topleft", labels, lwd = 2.5, col = colors[1:length(scenarios)], bty = "n")
  grid(col = "gray80", lty = 2)
}


# --- 6. EXAMPLE: BASELINE SIMULATION ---

cat("=== BASELINE SCENARIO ===\n")

# Initial state
state_0 <- c(
  A = 10e6,      # 10 million AVs
  C = 250e6,     # 250 million total vehicles
  I = 0.30       # 30% infrastructure readiness
)

# Constant policy
policy_baseline <- c(
  r = 5000,      # $5,000 rebate
  s = 2000,      # $2,000 subsidy
  i = 0.05       # 5% annual infrastructure investment
)

# Run simulation
years <- 30
results_baseline <- simulate_model(state_0, policy_baseline, years, params)

# Display final state
cat(sprintf("Year %d: AV share = %.1f%%\n", 
            years, 
            results_baseline$pct_AV[years + 1] * 100))

# Visualize
plot_states(results_baseline)


# --- 7. EXAMPLE: POLICY EXPERIMENTS ---

cat("\n=== POLICY EXPERIMENTS ===\n")

# Scenario 1: No policy
policy_none <- c(r = 0, s = 0, i = 0.01)
results_none <- simulate_model(state_0, policy_none, years, params)

# Scenario 2: High rebates
policy_high_rebate <- c(r = 10000, s = 2000, i = 0.05)
results_high_rebate <- simulate_model(state_0, policy_high_rebate, years, params)

# Scenario 3: Infrastructure focus
policy_infra <- c(r = 5000, s = 2000, i = 0.15)
results_infra <- simulate_model(state_0, policy_infra, years, params)

# Scenario 4: Aggressive all-around
policy_aggressive <- c(r = 12000, s = 5000, i = 0.20)
results_aggressive <- simulate_model(state_0, policy_aggressive, years, params)

# Compare scenarios
scenarios <- list(
  results_none,
  results_baseline,
  results_high_rebate,
  results_infra,
  results_aggressive
)

labels <- c(
  "No Policy",
  "Baseline",
  "High Rebates",
  "Infra Focus",
  "Aggressive"
)

plot_scenarios(scenarios, labels)

# Print final outcomes
cat("\nFinal AV Share (Year 30):\n")
for (i in seq_along(scenarios)) {
  cat(sprintf("  %-15s: %5.1f%%\n", 
              labels[i], 
              scenarios[[i]]$pct_AV[years + 1] * 100))
}


# --- 8. EXAMPLE: TIME-VARYING POLICY (PHASE-OUT) ---

cat("\n=== TIME-VARYING POLICY: GRADUAL PHASE-OUT ===\n")

# Create declining rebate schedule
policy_phaseout <- data.frame(
  r = seq(10000, 0, length.out = years),  # linear decline
  s = rep(3000, years),                    # constant subsidy
  i = seq(0.15, 0.05, length.out = years)  # declining infra investment
)

results_phaseout <- simulate_model(state_0, policy_phaseout, years, params)

par(mfrow = c(1, 2))
plot_av_share(results_phaseout, "Phase-Out Policy")

# Show policy trajectory
plot(results_phaseout$year[-1], results_phaseout$r[-1] / 1000, 
     type = "l", lwd = 2, col = "#2E86AB",
     xlab = "Year", ylab = "Rebate ($1000s)", 
     main = "Policy Schedule")
grid(col = "gray80", lty = 2)

par(mfrow = c(1, 1))


# --- 9. MONTE CARLO UNCERTAINTY ANALYSIS ---

cat("\n=== MONTE CARLO SIMULATION ===\n")

#' Run Monte Carlo with parameter uncertainty
monte_carlo <- function(state_0, policy, years, params, n_sim = 100) {
  results_mc <- vector("list", n_sim)
  
  for (i in 1:n_sim) {
    # Perturb parameters (±20% uniform)
    params_i <- params
    params_i$beta_0 <- params$beta_0 * runif(1, 0.8, 1.2)
    params_i$beta_r <- params$beta_r * runif(1, 0.8, 1.2)
    params_i$beta_i <- params$beta_i * runif(1, 0.8, 1.2)
    params_i$delta <- params$delta * runif(1, 0.8, 1.2)
    
    results_mc[[i]] <- simulate_model(state_0, policy, years, params_i)
  }
  
  return(results_mc)
}

# Run 100 simulations
set.seed(42)
mc_sims <- monte_carlo(state_0, policy_baseline, years, params, n_sim = 100)

# Extract AV shares
av_shares <- sapply(mc_sims, function(r) r$pct_AV)

# Plot uncertainty bands
plot(
  results_baseline$year,
  results_baseline$pct_AV * 100,
  type = "n",
  xlab = "Year",
  ylab = "AV Share (%)",
  main = "Monte Carlo Uncertainty Bands",
  ylim = c(0, max(av_shares) * 100)
)

# 90% confidence interval
q05 <- apply(av_shares, 1, quantile, probs = 0.05) * 100
q95 <- apply(av_shares, 1, quantile, probs = 0.95) * 100
polygon(
  c(results_baseline$year, rev(results_baseline$year)),
  c(q05, rev(q95)),
  col = rgb(0.18, 0.52, 0.67, 0.3),
  border = NA
)

# Median
q50 <- apply(av_shares, 1, quantile, probs = 0.5) * 100
lines(results_baseline$year, q50, lwd = 2.5, col = "#2E86AB")

legend("topleft", 
       c("Median", "90% CI"), 
       lwd = c(2.5, 10), 
       col = c("#2E86AB", rgb(0.18, 0.52, 0.67, 0.3)),
       bty = "n")
grid(col = "gray80", lty = 2)


# --- 10. EXPORT RESULTS ---

cat("\n=== SAVING RESULTS ===\n")

# Save baseline results to CSV
write.csv(results_baseline, "av_baseline_results.csv", row.names = FALSE)
cat("Baseline results saved to av_baseline_results.csv\n")

# Save all scenarios
scenarios_combined <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
  df <- scenarios[[i]]
  df$scenario <- labels[i]
  df
}))

write.csv(scenarios_combined, "av_all_scenarios.csv", row.names = FALSE)
cat("All scenarios saved to av_all_scenarios.csv\n")

cat("\n=== MODEL RUN COMPLETE ===\n")
cat("To customize:\n")
cat("  1. Edit params list for your calibration\n")
cat("  2. Modify policy vectors for different interventions\n")
cat("  3. Use simulate_model() for custom runs\n")
cat("  4. Extend step_model() for additional dynamics\n")