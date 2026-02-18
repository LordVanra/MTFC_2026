params <- list(
  delta_1 = 0.07,        # annual vehicle retirement rate (7%)
  delta_2 = 0.05,          # infrastructure depreciation rate (5%)
  S = 15e6,            # annual new car sales (15 million)
  Price_AV = 40000,        # average AV price ($40k)
  
  beta_0 = -2.0,       # baseline (intercept)
  beta_1 = 3.0,        # rebate sensitivity
  beta_2 = 4.0,        # infrastructure sensitivity
  
  cap_0 = 0.20,        # baseline production capacity (20% of sales)
  gamma = 0.00002      # subsidy effectiveness on capacity
)

# This might need to be some form of static param because the same thing is used in 3 different places , but thats up to Tanay
sigmoid <- function(z) {
  1 / (1 + exp(-z))
}

# One-step state transition: x_{t+1} = f(x_t, u_t)
step_model <- function(state, policy, params) {
  A_i <- as.numeric(state["A"])  # AV count
  C_i <- as.numeric(state["C"])  # total cars
  I_i <- as.numeric(state["I"])  # infrastructure
  k_i <- as.numeric(state["K"])  # cap as percent
  
  r <- as.numeric(policy["r"])  # rebates
  s <- as.numeric(policy["s"])  # subsidies
  i <- as.numeric(policy["i"])  # infrastructure spending
  
  # Infrastructure dynamics
  I_iPlus1 <- (1 - params$delta_2) * I_i + i
  I_scaled <- sigmoid(I_iPlus1)
  
  # Adoption fraction (demand-side)
  price_adv <- r / params$Price_AV
  a_uncapped <- sigmoid(params$beta_0 + params$beta_1 * price_adv + params$beta_2 * I_i)
  
  # Supply constraint
  cap <- params$cap_0 + params$gamma * s + k_i
  K_iPlus1 <- max(0, min(1, (k_i + params$gamma * s + sigmoid(k_i))))
  a <- min(a_uncapped, cap)
  
  # Fleet dynamics
  A_iPlus1 <- (1 - params$delta_1) * A_i + a * params$S
  C_iPlus1 <- (1 - params$delta_1) * C_i + params$S
  
  c(A = A_iPlus1, C = C_iPlus1, I = I_iPlus1, K = K_iPlus1, adoption = a, cap = cap)
}

# Multi-year simulation
simulate_model <- function(state_0, policy, years, params) {
  # Initialize output
  out <- data.frame(
    year = 0:years,
    A = numeric(years + 1),
    C = numeric(years + 1),
    I = numeric(years + 1),
    K = numeric(years + 1),
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
  out$K[1] <- state_0["K"]
  out$pct_AV[1] <- state_0["A"] / state_0["C"]
  
  # Handle constant vs time-varying policy
  if (is.vector(policy) && !is.data.frame(policy)) {
    # Constant policy
    policy_r <- rep(as.numeric(policy["r"]), years)
    policy_s <- rep(as.numeric(policy["s"]), years)
    policy_i <- rep(as.numeric(policy["i"]), years)
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
    u_t <- c(r = as.numeric(policy_r[t]), s = as.numeric(policy_s[t]), i = as.numeric(policy_i[t]))
    
    # Step forward
    next_state <- step_model(current, u_t, params)
    
    # Store
    out$A[t + 1] <- next_state["A"]
    out$C[t + 1] <- next_state["C"]
    out$I[t + 1] <- next_state["I"]
    out$K[t + 1] <- next_state["K"]
    out$pct_AV[t + 1] <- next_state["A"] / next_state["C"]
    out$adoption[t + 1] <- next_state["adoption"]
    out$r[t + 1] <- u_t["r"]
    out$s[t + 1] <- u_t["s"]
    out$i[t + 1] <- u_t["i"]
    
    # Update state
    current <- next_state[c("A", "C", "I", "K")]
  }
  
  return(out)
}

# Visualization 
plot_all_states <- function(results, scenario_name) {
  
  # Start clean page
  plot.new()
  title(main = paste("Scenario:", scenario_name), cex.main = 1.8)
  
  # Layout grid
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  plot(results$year, results$pct_AV * 100, type = "l", lwd = 2, col = "#2E86AB",
       xlab = "Year", ylab = "AV Share (%)",
       main = "AV Market Penetration")
  grid()
  
  plot(results$year, results$I * 100, type = "l", lwd = 2, col = "#A23B72",
       xlab = "Year", ylab = "Infrastructure (%)",
       main = "Infrastructure")
  grid()
  
  plot(results$year, results$adoption * 100, type = "l", lwd = 2, col = "#F18F01",
       xlab = "Year", ylab = "Adoption Rate (%)",
       main = "Adoption Rate")
  grid()
  
  plot(results$year, results$A / 1e6, type = "l", lwd = 2, col = "#6A994E",
       xlab = "Year", ylab = "AV Fleet (millions)",
       main = "Fleet Size")
  grid()
  
  # Reset layout
  par(mfrow = c(1, 1))
}



plot_scenarios <- function(result_list, labels, var = "pct_AV", scale = 100) {
  
  colors <- c("#2E86AB", "#A23B72", "#F18F01",
              "#6A994E", "#C73E1D", "#7209B7")  

  max_val <- max(sapply(result_list,
                          function(s) max(s[[var]], na.rm = TRUE)))
    
    plot(result_list[[1]]$year, result_list[[1]][[var]] * scale, type = "l",
         lwd = 2.5, col = colors[1], xlab = "Year", ylab = "AV Share (%)",
         main = "Comparison Plot", ylim = c(0, max_val * scale * 1.1), las = 1)
    
    if (length(result_list) > 1) {
      for (i in 2:length(result_list)) {
        lines(result_list[[i]]$year,
              result_list[[i]][[var]] * scale,
              lwd = 2.5,
              col = colors[i])
      }
    }
    
    legend("topleft",
           labels,
           lwd = 2.5,
           col = colors[1:length(labels)],
           bty = "n")
    
    grid()
  
  vars  <- c("pct_AV", "I", "adoption", "A")
  ylabs <- c("AV Share (%)", "Infrastructure (%)", "Adoption Rate (%)", "Fleet (millions)")
  scales <- c(100, 100, 100, 1e-6)
  
  par(mfrow = c(2, 2), mar = c(4,4,2,1))
  
  for (v in 1:4) {
    
    max_val <- max(sapply(result_list,
                          function(s) max(s[[vars[v]]], na.rm = TRUE)))
    
    plot(result_list[[1]]$year,
        result_list[[1]][[vars[v]]] * scales[v],
        type = "l",
        lwd = 2,
        col = colors[1],
        xlab = "Year",
        ylab = ylabs[v],
        main = ylabs[v],
        ylim = c(0, max_val * scales[v] * 1.1),
        las = 1)
    
    if (length(result_list) > 1) {
      for (i in 2:length(result_list)) {
        lines(result_list[[i]]$year,
              result_list[[i]][[vars[v]]] * scales[v],
              lwd = 2,
              col = colors[i])
      }
    }
    
    grid()
  }
  
  legend("bottomleft",
        labels,
        lwd = 2,
        col = colors[1:length(labels)],
        bty = "n")
  
  par(mfrow = c(1,1))
}


# Monte Carlo for Uncertainty
# TS DEF WRONG REWRITE WITH SAME STATS PER SIM
monte_carlo <- function(state_0, policy, years, params, n_sim = 1000) {
  n_sim <- 1
  results_mc <- vector("list", n_sim)
  
  for (i in 1:n_sim) {
    params_i <- params
    params_i$beta_0 <- params$beta_0 * runif(1, 0.8, 1.2)
    params_i$beta_1 <- params$beta_1 * runif(1, 0.8, 1.2)
    params_i$beta_2 <- params$beta_2 * runif(1, 0.8, 1.2)
    params_i$delta_1 <- params$delta_1 * runif(1, 0.8, 1.2)
    params_i$delta_2 <- params$delta_2 * runif(1, 0.8, 1.2)
    
    results_mc[[i]] <- simulate_model(state_0, policy, years, params_i)
  }

  av_shares_mc <- sapply(results_mc, function(r) r$pct_AV[31])
  cat(sprintf("   Year 30: %.1f%% (median), 90%% CI [%.1f%%, %.1f%%]\n\n",
            median(av_shares_mc, na.rm = TRUE) * 100,
            quantile(av_shares_mc, 0.05, na.rm = TRUE) * 100,
            quantile(av_shares_mc, 0.95, na.rm = TRUE) * 100))
  
  return(results_mc)
}

plot_monte_carlo <- function(mc_sims, baseline) {
  av_shares <- sapply(mc_sims, function(r) r$pct_AV)
  
  q05 <- apply(av_shares, 1, quantile, probs = 0.05, na.rm = TRUE) * 100
  q50 <- apply(av_shares, 1, quantile, probs = 0.50, na.rm = TRUE) * 100
  q95 <- apply(av_shares, 1, quantile, probs = 0.95, na.rm = TRUE) * 100
  
  plot(baseline$year, q50, type = "n",
       xlab = "Year", ylab = "AV Share (%)",
       main = paste("Monte Carlo (90% CI)"),
       ylim = c(0, max(q95, na.rm = TRUE) * 1.1))
  
  polygon(c(baseline$year, rev(baseline$year)), c(q05, rev(q95)),
          col = rgb(0.18, 0.52, 0.67, 0.3), border = NA)
  
  lines(baseline$year, q50, lwd = 2.5, col = "#2E86AB")
  grid()
}


# Sensitivity Analysis
sensitivity_analysis <- function(state_0, policy, years, params, param_name, range = c(0.5, 1.5)) {
  multipliers <- seq(range[1], range[2], length.out = 20)
  results <- data.frame(multiplier = multipliers, final_share = NA_real_)
  
  for (i in seq_along(multipliers)) {
    params_temp <- params
    params_temp[[param_name]] <- params[[param_name]] * multipliers[i]
    sim <- simulate_model(state_0, policy, years, params_temp)
    results$final_share[i] <- sim$pct_AV[years + 1]
  }

  cat(sprintf("   %s sensitivity: %.1f%% to %.1f%% (0.5x to 1.5x)\n\n",
            param_name,
            min(results$final_share, na.rm = TRUE) * 100,
            max(results$final_share, na.rm = TRUE) * 100))
  
  return(results)
}

plot_sensitivity <- function(sens_results, param_name) {
  plot(sens_results$multiplier, sens_results$final_share * 100,
       type = "l", lwd = 2, col = "#2E86AB",
       xlab = paste(param_name, "Multiplier"),
       ylab = "Final AV Share (%)",
       main = paste("Sensitivity to", param_name))
  abline(v = 1, lty = 2, col = "red")
  grid()
}

plot_all_sensitivies <- function(result_list,
                                 labels,
                                 title,
                                 var = "final_share",
                                 scale = 100) {
  colors <- c("#2E86AB", "#A23B72", "#F18F01",
              "#6A994E", "#C73E1D", "#7209B7")
  
  max_val <- max(sapply(result_list, function(s) max(s[[var]], na.rm = TRUE)))
  
  # Use 'multiplier' as x-axis
  plot(result_list[[1]]$multiplier, result_list[[1]][[var]] * scale, type = "l", lwd = 2.5,
       col = colors[1], xlab = "Parameter Multiplier", ylab = "Final AV Share (%)",
       main = paste("Sensitivity Comparison for ", title), ylim = c(0, max_val * scale * 1.1), las = 1)
  
  if (length(result_list) > 1) {
    for (i in 2:length(result_list)) {
      lines(result_list[[i]]$multiplier, result_list[[i]][[var]] * scale, lwd = 2.5, col = colors[i])
    }
  }
  
  legend("topleft", labels, lwd = 2.5, col = colors[1:length(labels)], bty = "n")
  grid()
}

state_0 <- c(A = 10e6, C = 250e6, I = 0.30, K = 0.20)

# Policy scenarios
cat("POLICY SCENARIOS")

policy_none <- c(r = 0, s = 0, i = 0)
results_none <- simulate_model(state_0, policy_none, 30, params)

policy_high_rebate <- c(r = 10000, s = 500, i = 0.05)
results_high_rebate <- simulate_model(state_0, policy_high_rebate, 30, params)

policy_infra <- c(r = 1000, s = 500, i = 0.15)
results_infra <- simulate_model(state_0, policy_infra, 30, params)

policy_aggressive <- c(r = 12000, s = 5000, i = 0.20)
results_aggressive <- simulate_model(state_0, policy_aggressive, 30, params)

scenarios <- list(results_none, results_high_rebate, 
                  results_infra, results_aggressive)
labels <- c("No Policy", "High Rebates", "Infra Focus", "Aggressive")

cat("   Final AV Share (Year 30):\n")
for (i in seq_along(scenarios)) {
  cat(sprintf("   %-15s: %5.1f%%\n", labels[i], scenarios[[i]]$pct_AV[31] * 100))
}

# TS DEF UNTRUE BUT WHATEVER
cat("3. TIME-VARYING POLICY (Phase-Out)\n")
policy_phaseout <- data.frame(
  r = seq(10000, 0, length.out = 30),
  s = rep(3000, 30),
  i = seq(0.15, 0.05, length.out = 30)
)
results_phaseout <- simulate_model(state_0, policy_phaseout, 30, params)
cat(sprintf("   Year 30: AV share = %.1f%%\n\n", results_phaseout$pct_AV[31] * 100))

# TS DEF UNTRUE BUT WHATEVER
cat("4. MONTE CARLO UNCERTAINTY (100 simulations)\n")
set.seed(42)

mc_sims <- lapply(
  scenarios,
  function(s) monte_carlo(state_0, s, 30, params, n_sim = 100)
)

cat("5. SENSITIVITY ANALYSIS\n")

sens_results <- lapply(
  scenarios,
  function(s) {
    list(
      beta_0 = sensitivity_analysis(state_0, s, 30, params, "beta_0"),
      beta_1 = sensitivity_analysis(state_0, s, 30, params, "beta_1"),
      beta_2 = sensitivity_analysis(state_0, s, 30, params, "beta_2")
    )
  }
)

pdf("StateModel/av_model_plots.pdf", width = 10, height = 8)

for (i in seq_along(scenarios)) {
    plot_all_states(scenarios[[i]], labels[i])
    plot_monte_carlo(mc_sims[[i]], scenarios[[i]])
    plot_sensitivity(sens_results[[i]]$beta_0, "beta_0")
    plot_sensitivity(sens_results[[i]]$beta_1, "beta_1")
    plot_sensitivity(sens_results[[i]]$beta_2, "beta_2")
}

plot_scenarios(scenarios, labels)
plot_all_sensitivies(lapply(sens_results, function(x) x$beta_0), labels, "Beta 0", var = "final_share", scale = 100)
plot_all_sensitivies(lapply(sens_results, function(x) x$beta_1), labels, "Beta 1", var = "final_share", scale = 100)
plot_all_sensitivies(lapply(sens_results, function(x) x$beta_2), labels, "Beta 2", var = "final_share", scale = 100)

dev.off()

# Export Results

scenarios_combined <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
  df <- scenarios[[i]]
  df$scenario <- labels[i]
  df
}))
write.csv(scenarios_combined, "StateModel/av_all_scenarios.csv", row.names = FALSE)