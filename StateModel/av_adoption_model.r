# Complete implementation with:
# - Core state-space dynamics
# - Policy scenario analysis
# - Monte Carlo uncertainty
# - Parameter estimation
# - Visualization

cat("  AV ADOPTION STATE-SPACE MODEL\n")

# PART 1: MODEL PARAMETERS

params <- list(
  # Fleet dynamics
  delta = 0.07,        # annual vehicle retirement rate (7%)
  rho = 0.05,          # infrastructure decay rate (5%)
  S = 15e6,            # annual new car sales (15 million)
  P_AV = 40000,        # average AV price ($40k)
  
  # Adoption model (logistic regression coefficients)
  beta_0 = -2.0,       # baseline (intercept)
  beta_r = 3.0,        # rebate sensitivity
  beta_i = 4.0,        # infrastructure sensitivity
  
  # Supply-side constraints
  cap_0 = 0.20,        # baseline production capacity (20% of sales)
  gamma = 0.00002      # subsidy effectiveness on capacity
)

# PART 2: CORE FUNCTIONS

# Sigmoid (logistic) function
sigmoid <- function(z) {
  1 / (1 + exp(-z))
}

# One-step state transition: x_{t+1} = f(x_t, u_t)
step_model <- function(state, policy, params) {
  A <- as.numeric(state["A"])  # AV count
  C <- as.numeric(state["C"])  # total cars
  I <- as.numeric(state["I"])  # infrastructure
  
  r <- as.numeric(policy["r"])  # rebates
  s <- as.numeric(policy["s"])  # subsidies
  i <- as.numeric(policy["i"])  # infrastructure spending
  
  # Infrastructure dynamics
  I_next <- (1 - params$rho) * I + i
  I_next <- max(0, min(1, I_next))
  
  # Adoption fraction (demand-side)
  price_adv <- r / params$P_AV
  z <- params$beta_0 + params$beta_r * price_adv + params$beta_i * I
  a_uncapped <- sigmoid(z)
  
  # Supply constraint
  cap <- params$cap_0 + params$gamma * s
  a <- min(a_uncapped, cap)
  
  # Fleet dynamics
  A_next <- (1 - params$delta) * A + a * params$S
  C_next <- (1 - params$delta) * C + params$S
  
  c(A = A_next, C = C_next, I = I_next, adoption = a, cap = cap)
}

# Multi-year simulation
simulate_model <- function(state_0, policy, years, params) {
  # Initialize output
  out <- data.frame(
    year = 0:years,
    A = numeric(years + 1),
    C = numeric(years + 1),
    I = numeric(years + 1),
    pct_AV = numeric(years + 1),
    adoption = numeric(years + 1),
    r = numeric(years + 1),
    s = numeric(years + 1),
    i = numeric(years + 1)
  )
  
  # Initial state
  out$A[1] <- state_0["A"]
  out$C[1] <- state_0["C"]
  out$I[1] <- state_0["I"]
  out$pct_AV[1] <- state_0["A"] / state_0["C"]
  
  # Handle constant vs time-varying policy
  if (is.vector(policy) && !is.data.frame(policy)) {
    # Constant policy
    policy_r <- rep(policy["r"], years)
    policy_s <- rep(policy["s"], years)
    policy_i <- rep(policy["i"], years)
  } else {
    # Time-varying policy
    policy <- as.data.frame(policy)
    policy_r <- policy$r[1:years]
    policy_s <- policy$s[1:years]
    policy_i <- policy$i[1:years]
  }
  
  # Store initial policy
  out$r[1] <- policy_r[1]
  out$s[1] <- policy_s[1]
  out$i[1] <- policy_i[1]
  
  # Current state
  current <- state_0
  
  # Time loop
  for (t in 1:years) {
    # Current policy
    u_t <- c(r = policy_r[t], s = policy_s[t], i = policy_i[t])
    
    # Step forward
    next_state <- step_model(current, u_t, params)
    
    # Store
    out$A[t + 1] <- next_state["A"]
    out$C[t + 1] <- next_state["C"]
    out$I[t + 1] <- next_state["I"]
    out$pct_AV[t + 1] <- next_state["A"] / next_state["C"]
    out$adoption[t + 1] <- next_state["adoption"]
    out$r[t + 1] <- u_t["r"]
    out$s[t + 1] <- u_t["s"]
    out$i[t + 1] <- u_t["i"]
    
    # Update state
    current <- c(A = next_state["A"], C = next_state["C"], I = next_state["I"])
  }
  
  return(out)
}

# PART 3: VISUALIZATION FUNCTIONS

plot_av_share <- function(results, title = "AV Market Share Over Time") {
  plot(results$year, results$pct_AV * 100,
       type = "l", lwd = 2.5, col = "#2E86AB",
       xlab = "Year", ylab = "AV Share (%)", main = title, las = 1)
  grid(col = "gray80", lty = 2)
}

plot_all_states <- function(results) {
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  plot(results$year, results$pct_AV * 100, type = "l", lwd = 2, col = "#2E86AB",
       xlab = "Year", ylab = "AV Share (%)", main = "Market Penetration")
  grid()
  
  plot(results$year, results$I * 100, type = "l", lwd = 2, col = "#A23B72",
       xlab = "Year", ylab = "Infrastructure (%)", main = "Infrastructure Readiness")
  grid()
  
  plot(results$year, results$adoption * 100, type = "l", lwd = 2, col = "#F18F01",
       xlab = "Year", ylab = "Adoption Rate (%)", main = "Annual AV Adoption Rate")
  grid()
  
  plot(results$year, results$A / 1e6, type = "l", lwd = 2, col = "#6A994E",
       xlab = "Year", ylab = "AVs (millions)", main = "AV Fleet Size")
  grid()
  
  par(mfrow = c(1, 1))
}

plot_scenarios <- function(scenario_list, labels) {
  colors <- c("#2E86AB", "#A23B72", "#F18F01", "#6A994E", "#C73E1D", "#7209B7")
  
  max_share <- max(sapply(scenario_list, function(s) max(s$pct_AV, na.rm = TRUE)))
  
  plot(scenario_list[[1]]$year, scenario_list[[1]]$pct_AV * 100,
       type = "l", lwd = 2.5, col = colors[1],
       xlab = "Year", ylab = "AV Share (%)", main = "Policy Scenario Comparison",
       ylim = c(0, max_share * 100 * 1.1), las = 1)
  
  for (i in 2:length(scenario_list)) {
    lines(scenario_list[[i]]$year, scenario_list[[i]]$pct_AV * 100,
          lwd = 2.5, col = colors[i])
  }
  
  legend("topleft", labels, lwd = 2.5, col = colors[1:length(labels)], bty = "n")
  grid()
}

# PART 4: MONTE CARLO UNCERTAINTY

monte_carlo <- function(state_0, policy, years, params, n_sim = 100) {
  results_mc <- vector("list", n_sim)
  
  for (i in 1:n_sim) {
    params_i <- params
    params_i$beta_0 <- params$beta_0 * runif(1, 0.8, 1.2)
    params_i$beta_r <- params$beta_r * runif(1, 0.8, 1.2)
    params_i$beta_i <- params$beta_i * runif(1, 0.8, 1.2)
    params_i$delta <- params$delta * runif(1, 0.8, 1.2)
    
    results_mc[[i]] <- simulate_model(state_0, policy, years, params_i)
  }
  
  return(results_mc)
}

plot_monte_carlo <- function(mc_sims, baseline) {
  av_shares <- sapply(mc_sims, function(r) r$pct_AV)
  
  q05 <- apply(av_shares, 1, quantile, probs = 0.05, na.rm = TRUE) * 100
  q50 <- apply(av_shares, 1, quantile, probs = 0.50, na.rm = TRUE) * 100
  q95 <- apply(av_shares, 1, quantile, probs = 0.95, na.rm = TRUE) * 100
  
  plot(baseline$year, q50, type = "n",
       xlab = "Year", ylab = "AV Share (%)", main = "Monte Carlo Uncertainty (90% CI)",
       ylim = c(0, max(q95, na.rm = TRUE) * 1.1))
  
  polygon(c(baseline$year, rev(baseline$year)), c(q05, rev(q95)),
          col = rgb(0.18, 0.52, 0.67, 0.3), border = NA)
  
  lines(baseline$year, q50, lwd = 2.5, col = "#2E86AB")
  legend("topleft", c("Median", "90% CI"), 
         lwd = c(2.5, 10), col = c("#2E86AB", rgb(0.18, 0.52, 0.67, 0.3)), bty = "n")
  grid()
}

# PART 5: PARAMETER ESTIMATION

generate_synthetic_data <- function(n = 100, params) {
  set.seed(123)
  price_adv <- runif(n, 0, 0.3)
  infrastructure <- runif(n, 0, 1)
  z <- params$beta_0 + params$beta_r * price_adv + params$beta_i * infrastructure
  a_true <- sigmoid(z)
  adoption <- pmin(pmax(a_true + rnorm(n, 0, 0.05), 0), 1)
  data.frame(adoption = adoption, price_adv = price_adv, infrastructure = infrastructure)
}

estimate_parameters <- function(data) {
  model <- glm(adoption ~ price_adv + infrastructure, data = data,
               family = binomial(link = "logit"))
  coefs <- coef(model)
  list(beta_0 = coefs[1], beta_r = coefs[2], beta_i = coefs[3], model = model)
}

# PART 6: SENSITIVITY ANALYSIS

sensitivity_analysis <- function(state_0, policy, years, params, 
                                param_name = "beta_r", range = c(0.5, 1.5)) {
  multipliers <- seq(range[1], range[2], length.out = 20)
  results <- data.frame(multiplier = multipliers, final_share = NA_real_)
  
  for (i in seq_along(multipliers)) {
    params_temp <- params
    params_temp[[param_name]] <- params[[param_name]] * multipliers[i]
    sim <- simulate_model(state_0, policy, years, params_temp)
    results$final_share[i] <- sim$pct_AV[years + 1]
  }
  
  return(results)
}

plot_sensitivity <- function(sens_results, param_name) {
  plot(sens_results$multiplier, sens_results$final_share * 100,
       type = "l", lwd = 2, col = "#2E86AB",
       xlab = paste(param_name, "Multiplier"), ylab = "Final AV Share (%)",
       main = paste("Sensitivity to", param_name))
  abline(v = 1, lty = 2, col = "red")
  grid()
}

# PART 7: RUN EXAMPLES

# Initial state
state_0 <- c(A = 10e6, C = 250e6, I = 0.30)

# BASELINE SCENARIO
cat("1. BASELINE SCENARIO\n")
policy_baseline <- c(r = 5000, s = 2000, i = 0.05)
results_baseline <- simulate_model(state_0, policy_baseline, 30, params)
cat(sprintf("   Year 30: AV share = %.1f%%\n\n", results_baseline$pct_AV[31] * 100))

# POLICY SCENARIOS
cat("2. POLICY SCENARIOS\n")

policy_none <- c(r = 0, s = 0, i = 0.01)
results_none <- simulate_model(state_0, policy_none, 30, params)

policy_high_rebate <- c(r = 10000, s = 2000, i = 0.05)
results_high_rebate <- simulate_model(state_0, policy_high_rebate, 30, params)

policy_infra <- c(r = 5000, s = 2000, i = 0.15)
results_infra <- simulate_model(state_0, policy_infra, 30, params)

policy_aggressive <- c(r = 12000, s = 5000, i = 0.20)
results_aggressive <- simulate_model(state_0, policy_aggressive, 30, params)

scenarios <- list(results_none, results_baseline, results_high_rebate, 
                 results_infra, results_aggressive)
labels <- c("No Policy", "Baseline", "High Rebates", "Infra Focus", "Aggressive")

cat("   Final AV Share (Year 30):\n")
for (i in seq_along(scenarios)) {
  cat(sprintf("   %-15s: %5.1f%%\n", labels[i], scenarios[[i]]$pct_AV[31] * 100))
}
cat("\n")

# TIME-VARYING POLICY
cat("3. TIME-VARYING POLICY (Phase-Out)\n")
policy_phaseout <- data.frame(
  r = seq(10000, 0, length.out = 30),
  s = rep(3000, 30),
  i = seq(0.15, 0.05, length.out = 30)
)
results_phaseout <- simulate_model(state_0, policy_phaseout, 30, params)
cat(sprintf("   Year 30: AV share = %.1f%%\n\n", results_phaseout$pct_AV[31] * 100))

# MONTE CARLO
cat("4. MONTE CARLO UNCERTAINTY (100 simulations)\n")
set.seed(42)
mc_sims <- monte_carlo(state_0, policy_baseline, 30, params, n_sim = 100)
av_shares_mc <- sapply(mc_sims, function(r) r$pct_AV[31])
cat(sprintf("   Year 30: %.1f%% (median), 90%% CI [%.1f%%, %.1f%%]\n\n",
            median(av_shares_mc, na.rm = TRUE) * 100,
            quantile(av_shares_mc, 0.05, na.rm = TRUE) * 100,
            quantile(av_shares_mc, 0.95, na.rm = TRUE) * 100))

# SENSITIVITY ANALYSIS
cat("5. SENSITIVITY ANALYSIS\n")
sens_beta_r <- sensitivity_analysis(state_0, policy_baseline, 30, params, "beta_r")
cat(sprintf("   beta_r sensitivity: %.1f%% to %.1f%% (0.5x to 1.5x)\n\n",
            min(sens_beta_r$final_share, na.rm = TRUE) * 100,
            max(sens_beta_r$final_share, na.rm = TRUE) * 100))

# PARAMETER ESTIMATION
cat("6. PARAMETER ESTIMATION (Synthetic Data)\n")
synthetic_data <- generate_synthetic_data(200, params)
est_params <- estimate_parameters(synthetic_data)
cat(sprintf("   beta_0: True = %.3f, Est = %.3f\n", params$beta_0, est_params$beta_0))
cat(sprintf("   beta_r: True = %.3f, Est = %.3f\n", params$beta_r, est_params$beta_r))
cat(sprintf("   beta_i: True = %.3f, Est = %.3f\n\n", params$beta_i, est_params$beta_i))

# PART 8: GENERATE PLOTS

cat("GENERATING PLOTS\n")

# # Plot 1: Baseline trajectory
# cat("Plot 1: Baseline state trajectories\n")
# plot_all_states(results_baseline)

# # Plot 2: Scenario comparison
# cat("Plot 2: Policy scenario comparison\n")
# dev.new()
# plot_scenarios(scenarios, labels)

# # Plot 3: Monte Carlo uncertainty
# cat("Plot 3: Monte Carlo uncertainty bands\n")
# dev.new()
# plot_monte_carlo(mc_sims, results_baseline)

# # Plot 4: Sensitivity analysis
# cat("Plot 4: Sensitivity to rebate effectiveness\n")
# dev.new()
# plot_sensitivity(sens_beta_r, "beta_r")

# PART 9: EXPORT RESULTS

cat("EXPORTING RESULTS\n")

write.csv(results_baseline, "av_baseline_results.csv", row.names = FALSE)
cat("Baseline results saved to av_baseline_results.csv\n")

scenarios_combined <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
  df <- scenarios[[i]]
  df$scenario <- labels[i]
  df
}))
write.csv(scenarios_combined, "av_all_scenarios.csv", row.names = FALSE)
cat("All scenarios saved to av_all_scenarios.csv\n")
