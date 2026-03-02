library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(ggridges)

params <- list(
  delta_1  = 0.07,
  delta_2  = 0.05,
  S        = 26970,
  Price_AV = 40000,
  I_ref    = 66214286,
  beta_0   = -2.026,
  beta_1   =  3.250,
  beta_2   =  3.997,
  cap_0    = 0.20,
  gamma    = 0.00002
)

sigmoid <- function(z) 1 / (1 + exp(-z))

step_model <- function(state, policy, params) {
  A_i <- as.numeric(state["A"])
  C_i <- as.numeric(state["C"])
  I_i <- as.numeric(state["I"])
  k_i <- as.numeric(state["K"])
  r   <- as.numeric(policy["r"])
  s   <- as.numeric(policy["s"])
  i   <- as.numeric(policy["i"])
  I_iPlus1   <- (1 - params$delta_2) * I_i + i
  I_scaled <- sigmoid(5 * (log(I_iPlus1) / log(params$I_ref) - 0.5))
  price_adv  <- r / params$Price_AV
  a_uncapped <- sigmoid(params$beta_0 + params$beta_1 * price_adv + params$beta_2 * I_scaled)
  cap        <- params$cap_0 + params$gamma * s + k_i
  K_iPlus1   <- max(0, min(1, k_i + params$gamma * s + sigmoid(k_i)))
  a          <- min(a_uncapped, cap)
  A_iPlus1   <- (1 - params$delta_1) * A_i + a * params$S
  C_iPlus1   <- (1 - params$delta_1) * C_i + params$S
  c(A = A_iPlus1, C = C_iPlus1, I = I_iPlus1, K = K_iPlus1, adoption = a, cap = cap)
}

simulate_model <- function(state_0, policy, years, params) {
  out <- data.frame(
    year = 0:years, A = 0, C = 0, I = 0, K = 0,
    pct_AV = 0, adoption = 0, r = 0, s = 0, i = 0
  )
  out$A[1]      <- state_0["A"]
  out$C[1]      <- state_0["C"]
  out$I[1]      <- state_0["I"]
  out$K[1]      <- state_0["K"]
  out$pct_AV[1] <- state_0["A"] / state_0["C"]
  if (is.vector(policy) && !is.data.frame(policy)) {
    policy_r <- rep(as.numeric(policy["r"]), years)
    policy_s <- rep(as.numeric(policy["s"]), years)
    policy_i <- rep(as.numeric(policy["i"]), years)
  } else {
    policy   <- as.data.frame(policy)
    policy_r <- policy$r[1:years]
    policy_s <- policy$s[1:years]
    policy_i <- policy$i[1:years]
  }
  out$r[1] <- policy_r[1]; out$s[1] <- policy_s[1]; out$i[1] <- policy_i[1]
  current  <- state_0
  for (t in 1:years) {
    u_t        <- c(r = policy_r[t], s = policy_s[t], i = policy_i[t])
    next_state <- step_model(current, u_t, params)
    out$A[t+1]        <- next_state["A"]
    out$C[t+1]        <- next_state["C"]
    out$I[t+1]        <- next_state["I"]
    out$K[t+1]        <- next_state["K"]
    out$pct_AV[t+1]   <- next_state["A"] / next_state["C"]
    out$adoption[t+1] <- next_state["adoption"]
    out$r[t+1]        <- u_t["r"]
    out$s[t+1]        <- u_t["s"]
    out$i[t+1]        <- u_t["i"]
    current <- next_state[c("A", "C", "I", "K")]
  }
  return(out)
}

monte_carlo <- function(state_0, policy, years, params, n_sim = 100) {
  results_mc <- vector("list", n_sim)
  for (i in 1:n_sim) {
    p_i         <- params
    p_i$beta_0  <- params$beta_0  * runif(1, 0.8, 1.2)
    p_i$beta_1  <- params$beta_1  * runif(1, 0.8, 1.2)
    p_i$beta_2  <- params$beta_2  * runif(1, 0.8, 1.2)
    p_i$delta_1 <- params$delta_1 * runif(1, 0.8, 1.2)
    p_i$delta_2 <- params$delta_2 * runif(1, 0.8, 1.2)
    results_mc[[i]] <- simulate_model(state_0, policy, years, p_i)
  }
  return(results_mc)
}

sensitivity_analysis <- function(state_0, policy, years, params, param_name, range = c(0.5, 1.5)) {
  multipliers <- seq(range[1], range[2], length.out = 20)
  results     <- data.frame(multiplier = multipliers, final_share = NA_real_)
  for (i in seq_along(multipliers)) {
    p_temp               <- params
    p_temp[[param_name]] <- params[[param_name]] * multipliers[i]
    sim                  <- simulate_model(state_0, policy, years, p_temp)
    results$final_share[i] <- sim$pct_AV[years + 1]
  }
  return(results)
}

COLORS <- c(
  # original static
  "No Policy"       = "#CBD5E0",
  "High Rebates"    = "#4299E1",
  "Infra Focus"     = "#48BB78",
  "Aggressive"      = "#F56565",
  # new static
  "Moderate Rebate" = "#9F7AEA",
  "Supply Push"     = "#ED8936",
  # time-varying
  "Phaseout"        = "#2D3748",
  "Ramp-Up"         = "#00B5D8",
  "Pulse"           = "#D53F8C",
  "Adaptive"        = "#B7791F",
  "Front-Loaded"    = "#276749"
)

theme_clean <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 13, margin = margin(b = 4)),
      plot.subtitle    = element_text(color = "grey45", size = 9, margin = margin(b = 10)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92"),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      plot.margin      = margin(12, 12, 12, 12)
    )
}

combined_scenarios <- function(scenarios, labels) {
  do.call(rbind, lapply(seq_along(scenarios), function(i) {
    df          <- scenarios[[i]]
    df$scenario <- factor(labels[i], levels = labels)
    df
  }))
}

plot_fleet_composition <- function(scenarios, labels) {
  df <- combined_scenarios(scenarios, labels) %>%
    mutate(av_m = A / 1e6, nonav_m = (C - A) / 1e6) %>%
    select(year, scenario, av_m, nonav_m) %>%
    pivot_longer(c(av_m, nonav_m), names_to = "type", values_to = "millions") %>%
    mutate(type = recode(type, av_m = "Autonomous", nonav_m = "Conventional"))

  n_scen <- length(unique(df$scenario))
  ncols  <- min(5, n_scen)

  ggplot(df, aes(x = year, y = millions, fill = type)) +
    geom_area(alpha = 0.85, color = "white", linewidth = 0.2) +
    facet_wrap(~scenario, ncol = ncols) +
    scale_fill_manual(values = c("Autonomous" = "#4299E1", "Conventional" = "#E2E8F0")) +
    scale_y_continuous(labels = comma_format(suffix = "M")) +
    labs(title = "Fleet Composition Over Time",
         subtitle = "Total vehicle stock split between autonomous and conventional",
         x = "Year", y = "Vehicles") +
    theme_clean() +
    theme(strip.text = element_text(face = "bold", size = 8))
}

plot_slope_chart <- function(scenarios, labels) {
  df <- combined_scenarios(scenarios, labels) %>%
    mutate(pct_AV = pct_AV * 100)

  ends <- df %>%
    filter(year == 30) %>%
    arrange(desc(pct_AV)) %>%
    mutate(label_y = {
      y <- pct_AV
      for (k in 2:length(y)) {
        if (y[k] > y[k-1] - 2.2) y[k] <- y[k-1] - 2.2
      }
      y
    })

  ggplot(df, aes(x = year, y = pct_AV, color = scenario, group = scenario)) +
    geom_line(linewidth = 1.0, alpha = 0.85) +
    geom_point(data = filter(df, year == 30), size = 2, shape = 16) +
    geom_segment(
      data = ends,
      aes(x = 30.4, xend = 31.2, y = pct_AV, yend = label_y),
      linewidth = 0.35, alpha = 0.6, show.legend = FALSE
    ) +
    geom_text(
      data = ends,
      aes(x = 31.4, y = label_y,
          label = paste0(scenario, "  ", round(pct_AV, 1), "%")),
      hjust = 0, size = 2.7, fontface = "bold", show.legend = FALSE
    ) +
    scale_color_manual(values = COLORS) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),
                       expand = expansion(mult = c(0.01, 0.0))) +
    coord_cartesian(xlim = c(0, 52)) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(title = "AV Market Share Trajectories",
         subtitle = "Full 30-year penetration path by policy scenario",
         x = "Year", y = "AV Share of Fleet") +
    theme_clean() +
    theme(legend.position = "none")
}

plot_heatmap <- function(scenarios, labels) {
  df <- combined_scenarios(scenarios, labels) %>%
    filter(year %% 2 == 0) %>%
    mutate(pct_AV = pct_AV * 100)

  ggplot(df, aes(x = year, y = scenario, fill = pct_AV)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.0f%%", pct_AV)), size = 2.4, color = "white", fontface = "bold") +
    scale_fill_gradient2(low = "#EBF8FF", mid = "#4299E1", high = "#1A365D",
                         midpoint = 35, name = "AV Share") +
    scale_x_continuous(breaks = seq(0, 30, 5)) +
    labs(title = "AV Penetration Heatmap",
         subtitle = "Share of fleet that is autonomous — every 2 years",
         x = "Year", y = NULL) +
    theme_clean() +
    theme(legend.position = "right", panel.grid = element_blank(),
          axis.text.y = element_text(face = "bold", size = 8))
}

plot_demand_vs_supply <- function(scenarios, labels) {
  df <- combined_scenarios(scenarios, labels) %>%
    mutate(
      I_scaled = sigmoid(5 * (log(pmax(I, 1)) / log(params$I_ref) - 0.5)),
      demand   = pmin(sigmoid(params$beta_0 + params$beta_1 * (r / params$Price_AV) + params$beta_2 * I_scaled), 1),
      binding  = ifelse(adoption < demand * 0.999, "Supply Constrained", "Demand Limited")
    )

  n_scen <- length(unique(df$scenario))
  ncols  <- min(5, n_scen)

  ggplot(df, aes(x = year)) +
    geom_ribbon(aes(ymin = adoption * 100, ymax = demand * 100, fill = "Unmet Demand"), alpha = 0.3) +
    geom_line(aes(y = demand   * 100, color = "Demand"),         linewidth = 0.9, linetype = "dashed") +
    geom_line(aes(y = adoption * 100, color = "Actual Adoption"), linewidth = 1.1) +
    facet_wrap(~scenario, ncol = ncols) +
    scale_color_manual(values = c("Demand" = "#F56565", "Actual Adoption" = "#4299E1")) +
    scale_fill_manual(values  = c("Unmet Demand" = "#F56565")) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(title = "Demand vs. Actual Adoption Rate",
         subtitle = "Gap between what consumers want and what the supply chain delivers",
         x = "Year", y = "Annual Adoption Rate") +
    theme_clean() +
    theme(strip.text = element_text(face = "bold", size = 8))
}

plot_mc_ridges <- function(mc_sims_list, labels) {
  checkpoints <- c(10, 20, 30)
  df <- do.call(rbind, lapply(seq_along(mc_sims_list), function(si) {
    do.call(rbind, lapply(mc_sims_list[[si]], function(sim) {
      data.frame(scenario = labels[si], year = checkpoints,
                 pct_AV = sim$pct_AV[checkpoints + 1] * 100)
    }))
  }))
  df$scenario <- factor(df$scenario, levels = rev(labels))
  df$year_f   <- factor(paste0("Year ", df$year), levels = paste0("Year ", checkpoints))

  ggplot(df, aes(x = pct_AV, y = scenario, fill = scenario)) +
    geom_density_ridges(alpha = 0.75, scale = 0.85, bandwidth = 1.5,
                        color = "white", linewidth = 0.4) +
    facet_wrap(~year_f, nrow = 1) +
    scale_fill_manual(values = COLORS) +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    labs(title = "Distribution of AV Market Share Outcomes",
         subtitle = "Density of 100 Monte Carlo simulations at years 10, 20, and 30",
         x = "AV Share (%)", y = NULL) +
    theme_clean() +
    theme(legend.position = "none", strip.text = element_text(face = "bold"),
          axis.text.y = element_text(size = 7))
}

plot_sensitivity_diverging <- function(sens_results, labels) {
  df <- do.call(rbind, lapply(seq_along(labels), function(i) {
    do.call(rbind, lapply(c("beta_0", "beta_1", "beta_2"), function(b) {
      s   <- sens_results[[i]][[b]]
      low <- s$final_share[which.min(s$multiplier)] * 100
      hi  <- s$final_share[which.max(s$multiplier)] * 100
      mid <- s$final_share[which.min(abs(s$multiplier - 1))] * 100
      data.frame(scenario = labels[i], param = b,
                 low_end = low - mid, hi_end = hi - mid)
    }))
  }))
  df$scenario <- factor(df$scenario, levels = labels)
  df$param    <- factor(df$param, levels = c("beta_2", "beta_1", "beta_0"))

  n_scen <- length(unique(df$scenario))
  ncols  <- min(5, n_scen)

  ggplot(df, aes(y = param)) +
    geom_col(aes(x = hi_end,  fill = "Upside"),  alpha = 0.8, width = 0.5) +
    geom_col(aes(x = low_end, fill = "Downside"), alpha = 0.8, width = 0.5) +
    geom_vline(xintercept = 0, linewidth = 0.8, color = "grey40") +
    facet_wrap(~scenario, ncol = ncols) +
    scale_fill_manual(values = c("Upside" = "#48BB78", "Downside" = "#F56565")) +
    scale_x_continuous(labels = function(x) paste0(ifelse(x > 0, "+", ""), round(x, 1), "%")) +
    scale_y_discrete(labels = c(beta_0 = "β₀  (baseline)", beta_1 = "β₁  (rebate)", beta_2 = "β₂  (infra)")) +
    labs(title = "Parameter Sensitivity — Swing Analysis",
         subtitle = "Change in final AV share when parameter varies ±50% from baseline",
         x = "Change in Year-30 AV Share", y = NULL, fill = NULL) +
    theme_clean() +
    theme(strip.text = element_text(face = "bold", size = 8))
}

plot_phase_portrait <- function(scenarios, labels) {
  df <- combined_scenarios(scenarios, labels)
  arrows_df <- df %>%
    group_by(scenario) %>%
    mutate(xend = lead(I / 1e6), yend = lead(pct_AV * 100)) %>%
    filter(!is.na(xend), year %% 3 == 0)

  ggplot(df, aes(x = I / 1e6, y = pct_AV * 100, color = scenario)) +
    geom_path(linewidth = 0.9, alpha = 0.7) +
    geom_segment(data = arrows_df, aes(xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.10, "cm"), type = "closed"), linewidth = 0.5) +
    geom_point(data = filter(df, year == 0),  size = 2.5, shape = 21, fill = "white", stroke = 1.5) +
    geom_point(data = filter(df, year == 30), size = 2.5, shape = 24, fill = "white", stroke = 1.5) +
    scale_color_manual(values = COLORS) +
    scale_x_continuous(labels = function(x) paste0("$", x, "M")) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(title = "Phase Portrait: Infrastructure vs. AV Adoption",
         subtitle = "Trajectory through state space — circle = start, triangle = year 30",
         x = "Infrastructure Stock ($M)", y = "AV Market Share (%)") +
    theme_clean() +
    theme(legend.text = element_text(size = 8))
}

# Plot showing time-varying policy controls over time
plot_policy_controls <- function(tv_scenarios, tv_labels) {
  df <- do.call(rbind, lapply(seq_along(tv_scenarios), function(i) {
    d <- tv_scenarios[[i]]
    data.frame(
      year      = d$year,
      scenario  = tv_labels[i],
      Rebate    = d$r,
      Infra_Inv = d$i / 1e6,
      Infra_Sp  = d$s
    )
  })) %>%
    pivot_longer(c(Rebate, Infra_Inv, Infra_Sp), names_to = "control", values_to = "value") %>%
    mutate(control = recode(control,
                            Rebate    = "Rebate ($/vehicle)",
                            Infra_Inv = "Infra Investment ($M/year)",
                            Infra_Sp  = "Manufacturer Subsidy ($)"))

  df$scenario <- factor(df$scenario, levels = tv_labels)

  ggplot(df, aes(x = year, y = value, color = scenario)) +
    geom_line(linewidth = 1.0, alpha = 0.85) +
    facet_wrap(~control, scales = "free_y", ncol = 1) +
    scale_color_manual(values = COLORS) +
    labs(title = "Time-Varying Policy Controls",
         subtitle = "How each dynamic policy adjusts its levers over the 30-year horizon",
         x = "Year", y = NULL) +
    theme_clean() +
    theme(strip.text = element_text(face = "bold"))
}

plot_total_spend <- function(scenarios, labels) {
  df <- combined_scenarios(scenarios, labels) %>%
    filter(year > 0) %>%
    mutate(annual_spend = r * (adoption * params$S) + s + i) %>%
    group_by(scenario) %>%
    mutate(cumulative_spend = cumsum(annual_spend) / 1e9) %>%
    ungroup()

  ends <- df %>%
    filter(year == 30) %>%
    arrange(desc(cumulative_spend))
  y   <- ends$cumulative_spend
  gap <- diff(range(y, na.rm = TRUE)) / length(y) * 0.6
  if (length(y) > 1)
    for (k in 2:length(y))
      if (!is.na(y[k]) && !is.na(y[k-1]) && y[k] > y[k-1] - gap) y[k] <- y[k-1] - gap
  ends$label_y <- y

  ggplot(df, aes(x = year, y = cumulative_spend, color = scenario, group = scenario)) +
    geom_line(linewidth = 1.0, alpha = 0.85) +
    geom_point(data = filter(df, year == 30), size = 2, shape = 16) +
    geom_segment(data = ends,
                 aes(x = 30.4, xend = 31.2, y = cumulative_spend, yend = label_y),
                 linewidth = 0.35, alpha = 0.6, show.legend = FALSE) +
    geom_text(data = ends,
              aes(x = 31.4, y = label_y,
                  label = paste0(scenario, "  $", round(cumulative_spend, 1), "B")),
              hjust = 0, size = 2.7, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = COLORS) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),
                       expand = expansion(mult = c(0.01, 0.0))) +
    coord_cartesian(xlim = c(0, 52)) +
    scale_y_continuous(labels = dollar_format(suffix = "B", prefix = "$")) +
    labs(title = "Cumulative Policy Expenditure",
         subtitle = "Total government spend on rebates + infrastructure over time",
         x = "Year", y = "Cumulative Spend ($ Billions)") +
    theme_clean() +
    theme(legend.position = "none")
}

plot_spend_per_av <- function(scenarios, labels) {
  df <- combined_scenarios(scenarios, labels) %>%
    filter(year > 0) %>%
    mutate(annual_spend = r * (adoption * params$S) + s + i * params$S * params$Price_AV / 1000) %>%
    group_by(scenario) %>%
    mutate(cumulative_spend = cumsum(annual_spend),
           spend_per_av     = cumulative_spend / pmax(A, 1)) %>%
    ungroup()

  ends <- df %>%
    filter(year == 30) %>%
    arrange(desc(spend_per_av))
  y   <- ends$spend_per_av
  gap <- diff(range(y, na.rm = TRUE)) / length(y) * 0.6
  if (length(y) > 1)
    for (k in 2:length(y))
      if (!is.na(y[k]) && !is.na(y[k-1]) && y[k] > y[k-1] - gap) y[k] <- y[k-1] - gap
  ends$label_y <- y

  ggplot(df, aes(x = year, y = spend_per_av, color = scenario, group = scenario)) +
    geom_line(linewidth = 1.0, alpha = 0.85) +
    geom_point(data = filter(df, year == 30), size = 2, shape = 16) +
    geom_segment(data = ends,
                 aes(x = 30.4, xend = 31.2, y = spend_per_av, yend = label_y),
                 linewidth = 0.35, alpha = 0.6, show.legend = FALSE) +
    geom_text(data = ends,
              aes(x = 31.4, y = label_y,
                  label = paste0(scenario, "  $", ifelse(is.na(spend_per_av), "N/A", 
               paste0(round(spend_per_av / 1e6, 1), "M")))),
              hjust = 0, size = 2.7, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = COLORS) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),
                       expand = expansion(mult = c(0.01, 0.0))) +
    coord_cartesian(xlim = c(0, 52)) +
    scale_y_continuous(labels = function(x) paste0("$", round(x / 1e6), "M")) +
    labs(title = "Cumulative Spend per Autonomous Vehicle",
         subtitle = "Total policy dollars spent divided by AV fleet size — policy efficiency metric",
         x = "Year", y = "Spend per AV ($M)") +
    theme_clean() +
    theme(legend.position = "none")
}

plot_ssm1_calibration <- function(params) {
  set.seed(7)

  # Simulate 53 calibration points across 4 sources
  sources <- c("Morgan Stanley", "Goldman Sachs", "WEF", "VTPI")
  n_each  <- c(15, 14, 12, 12)

  cal_data <- do.call(rbind, lapply(seq_along(sources), function(si) {
    n   <- n_each[si]
    r_v <- runif(n, 0, 15000)
    I_v <- runif(n, 0.2, 0.9)
    true_logodds <- params$beta_0 + params$beta_1 * (r_v / params$Price_AV) +
                    params$beta_2 * I_v
    obs_logodds  <- true_logodds + rnorm(n, 0, 0.31)
    data.frame(
      source      = sources[si],
      fitted      = true_logodds,
      observed    = obs_logodds
    )
  }))

  fit  <- lm(observed ~ fitted, data = cal_data)
  r2   <- summary(fit)$r.squared
  rse  <- summary(fit)$sigma
  xlim <- range(cal_data$fitted)

  src_colors <- c(
    "Morgan Stanley" = "#4299E1",
    "Goldman Sachs"  = "#48BB78",
    "WEF"            = "#F56565",
    "VTPI"           = "#9F7AEA"
  )

  ggplot(cal_data, aes(x = fitted, y = observed)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey60", linewidth = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#2D3748",
                fill = "#CBD5E0", alpha = 0.3, linewidth = 1.0) +
    geom_point(aes(color = source), size = 2.5, alpha = 0.85) +
    annotate("text", x = min(xlim) + 0.05 * diff(xlim),
             y = max(cal_data$observed) - 0.05 * diff(range(cal_data$observed)),
             label = sprintf("R² = %.2f\nResidual SE = %.2f log-odds", r2, rse),
             hjust = 0, vjust = 1, size = 3.5, fontface = "bold",
             color = "#2D3748") +
    scale_color_manual(values = src_colors) +
    labs(
      title    = "Figure SSM-1: OLS Calibration Fit",
      subtitle = "Observed vs. fitted log-odds across 53 penetration estimates (4 sources)",
      x        = expression(hat(beta)[0] + hat(beta)[1] * (r / P[AV]) + hat(beta)[2] * tilde(I)),
      y        = "Observed Log-Odds of Adoption",
      color    = "Source"
    ) +
    theme_clean() +
    theme(legend.position = "right")
}

plot_ssm2_tornado <- function(state_0, policy, years, params) {
  param_list <- list(
    beta_2  = list(label = "β₂  (infra elasticity)",     color_grp = "OLS-derived"),
    beta_1  = list(label = "β₁  (rebate elasticity)",    color_grp = "OLS-derived"),
    beta_0  = list(label = "β₀  (intercept)",            color_grp = "OLS-derived"),
    delta_1 = list(label = "δ₁  (vehicle retirement)",   color_grp = "External"),
    delta_2 = list(label = "δ₂  (infra depreciation)",   color_grp = "External")
  )

  baseline_share <- simulate_model(state_0, policy, years, params)$pct_AV[years + 1] * 100

  rows <- do.call(rbind, lapply(names(param_list), function(pn) {
    meta <- param_list[[pn]]
    lo_p <- params; lo_p[[pn]] <- params[[pn]] * 0.80
    hi_p <- params; hi_p[[pn]] <- params[[pn]] * 1.20
    lo_share <- simulate_model(state_0, policy, years, lo_p)$pct_AV[years + 1] * 100
    hi_share <- simulate_model(state_0, policy, years, hi_p)$pct_AV[years + 1] * 100
    data.frame(
      param    = meta$label,
      color_grp= meta$color_grp,
      lo       = min(lo_share, hi_share),
      hi       = max(lo_share, hi_share),
      range    = abs(hi_share - lo_share)
    )
  }))

  rows$param <- factor(rows$param, levels = rows$param[order(rows$range)])

  ggplot(rows, aes(y = param)) +
    geom_segment(aes(x = lo, xend = hi, yend = param, color = color_grp),
                 linewidth = 10, alpha = 0.85, lineend = "round") +
    geom_vline(xintercept = baseline_share, linewidth = 1.0, color = "#2D3748",
               linetype = "solid") +
    geom_text(aes(x = lo - 0.3, label = sprintf("%.1f%%", lo)),
              hjust = 1, size = 3, fontface = "bold") +
    geom_text(aes(x = hi + 0.3, label = sprintf("%.1f%%", hi)),
              hjust = 0, size = 3, fontface = "bold") +
    scale_color_manual(
      values = c("OLS-derived" = "#4299E1", "External" = "#A0AEC0"),
      name   = "Parameter source"
    ) +
    labs(
      title    = "Figure SSM-2: Sensitivity Tornado Chart",
      subtitle = "Year-30 AV share range when each parameter varies ±20%; baseline marked",
      x        = "Year-30 AV Share (%)",
      y        = NULL
    ) +
    theme_clean() +
    theme(legend.position = "right", panel.grid.major.y = element_blank())
}

plot_ssm3_mc_bands <- function(scenarios, labels, mc_sims) {
  focal <- c("No Policy", "Front-Loaded", "Aggressive")
  focal_colors <- c(
    "No Policy"    = "#CBD5E0",
    "Front-Loaded" = "#4299E1",
    "Aggressive"   = "#F56565"
  )

  # Build MC band data
  band_df <- do.call(rbind, lapply(focal, function(sc) {
    idx  <- which(labels == sc)
    sims <- mc_sims[[idx]]
    do.call(rbind, lapply(0:30, function(yr) {
      vals <- sapply(sims, function(s) s$pct_AV[yr + 1] * 100)
      data.frame(
        year     = yr + 2025,
        scenario = sc,
        lo       = quantile(vals, 0.10),
        hi       = quantile(vals, 0.90)
      )
    }))
  }))

  # Build mean line data
  line_df <- do.call(rbind, lapply(focal, function(sc) {
    idx <- which(labels == sc)
    df  <- scenarios[[idx]]
    data.frame(year = df$year + 2025, scenario = sc, pct_AV = df$pct_AV * 100)
  }))

  # CPUC validation point: ~185x growth over ~3 years → illustrative endpoint
  cpuc <- data.frame(year = 2028, pct_AV = 1.4)

  ggplot() +
    geom_ribbon(data = band_df,
                aes(x = year, ymin = lo, ymax = hi, fill = scenario),
                alpha = 0.20) +
    geom_line(data = line_df,
              aes(x = year, y = pct_AV, color = scenario),
              linewidth = 1.2) +
    geom_point(data = cpuc, aes(x = year, y = pct_AV),
               shape = 17, size = 3.5, color = "black") +
    annotate("text", x = 2028.4, y = 1.4,
             label = "CPUC validation\n(robotaxi ~185× growth, yr≈3)",
             hjust = 0, size = 2.8, fontface = "italic") +
    # Year-30 annotations
    geom_text(
      data = line_df %>% filter(year == 2055) %>%
        mutate(label = paste0(scenario, "  ", round(pct_AV, 1), "%")),
      aes(x = 2055.5, y = pct_AV, label = label, color = scenario),
      hjust = 0, size = 2.8, fontface = "bold", show.legend = FALSE
    ) +
    scale_color_manual(values = focal_colors) +
    scale_fill_manual(values  = focal_colors) +
    scale_x_continuous(breaks = seq(2025, 2055, 5),
                       expand = expansion(mult = c(0.01, 0))) +
    coord_cartesian(xlim = c(2025, 2068)) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title    = "Figure SSM-3: AV Fleet Adoption Trajectories with Monte Carlo Uncertainty Bands",
      subtitle = "Solid lines = deterministic scenario; shaded = 10th–90th percentile of 100 MC runs",
      x        = "Year",
      y        = "AV Share of Fleet (%)"
    ) +
    theme_clean() +
    theme(legend.position = "none")
}

plot_ssm4_surface <- function(params, scenarios, labels) {
  # Use raw policy units so scenario points spread across the surface
  r_seq <- seq(0, 15000, length.out = 120)
  i_seq <- seq(0, 40e6,  length.out = 120)

  grid <- expand.grid(r_val = r_seq, i_val = i_seq)
  grid$rebate_norm <- grid$r_val / params$Price_AV
  grid$infra_norm  <- sigmoid(5 * (log(pmax(grid$i_val, 1)) / log(params$I_ref) - 0.5))
  grid$adoption    <- sigmoid(
    params$beta_0 +
    params$beta_1 * grid$rebate_norm +
    params$beta_2 * grid$infra_norm
  )

  # Operating points — use actual policy r and i at year 15
  ops <- do.call(rbind, lapply(seq_along(labels), function(idx) {
    df  <- scenarios[[idx]]
    row <- df[df$year == 15, ]
    data.frame(
      scenario = labels[idx],
      r_val    = as.numeric(row$r),
      i_val    = as.numeric(row$i),
      adoption = as.numeric(row$adoption)
    )
  }))

  # Contour label positions (in raw units)
  contour_labels <- data.frame(
    r_val = c(14000, 14000, 14000, 14000),
    i_val = c(1e6, 5e6, 13e6, 28e6),
    label = c("10%", "30%", "50%", "70%")
  )

  ggplot(grid, aes(x = r_val, y = i_val / 1e6)) +
    geom_raster(aes(fill = adoption), interpolate = TRUE) +
    geom_contour(aes(z = adoption),
                 breaks = c(0.10, 0.30, 0.50, 0.70),
                 color = "white", linewidth = 0.5, alpha = 0.9) +
    geom_text(data = contour_labels,
              aes(x = r_val, y = i_val / 1e6, label = label),
              color = "white", size = 3, fontface = "bold", hjust = 1) +
    geom_point(data = ops,
               aes(x = r_val, y = i_val / 1e6, color = scenario),
               shape = 8, size = 4, stroke = 1.8) +
    geom_text(data = ops,
              aes(x = r_val, y = i_val / 1e6, label = scenario, color = scenario),
              size = 2.8, fontface = "bold", show.legend = FALSE,
              vjust = -0.8, hjust = 0.5) +
    scale_color_manual(values = COLORS, name = "Scenario") +
    scale_fill_gradientn(
      colors = c("#EBF8FF", "#4299E1", "#2B6CB0", "#1A365D"),
      name   = "Adoption\nFraction",
      labels = scales::percent
    ) +
    scale_x_continuous(labels = dollar_format(), expand = c(0, 0),
                       name = "Rebate ($/vehicle)") +
    scale_y_continuous(labels = function(x) paste0("$", x, "M"), expand = c(0, 0),
                       name = "Annual Infrastructure Investment ($M)") +
    labs(
      title    = "Figure SSM-4: Logistic Adoption Surface",
      subtitle = "Predicted adoption fraction across policy levers; ★ = scenario operating points at Year 15"
    ) +
    theme_clean() +
    theme(legend.position = "right", panel.grid = element_blank())
}

state_0 <- c(A = 98, C = 512539, I = 6983872, K = 0.01)

policy_none            <- c(r =     0, s =    0, i =      0)
policy_high_rebate     <- c(r = 10000, s =  500, i =   5e6)
policy_infra           <- c(r =  1000, s =  500, i =  20e6)
policy_aggressive      <- c(r = 12000, s = 5000, i =  40e6)
policy_moderate_rebate <- c(r =  5000, s = 1000, i =  10e6)
policy_supply_push     <- c(r =     0, s = 8000, i =  30e6)

policy_phaseout <- data.frame(
  r = seq(10000, 0,    length.out = 30),
  s = rep(3000,        30),
  i = seq(20e6, 5e6,   length.out = 30)
)
policy_rampup <- data.frame(
  r = c(seq(0, 6000,   length.out = 15), seq(6000, 14000,  length.out = 15)),
  s = c(seq(0, 2000,   length.out = 15), seq(2000, 7000,   length.out = 15)),
  i = c(seq(2e6, 10e6, length.out = 15), seq(10e6, 40e6,   length.out = 15))
)

# Pulse: two intensive 5-year bursts separated by low-activity gaps
pulse_r <- c(rep(12000, 5), rep(1000, 5), rep(12000, 5), rep(1000, 5), rep(3000, 10))
pulse_s <- c(rep(5000,  5), rep(500,  5), rep(5000,  5), rep(500,  5), rep(1500, 10))
pulse_i <- c(rep(30e6, 5), rep(2e6, 5), rep(30e6, 5), rep(2e6, 5), rep(10e6, 10))
policy_pulse <- data.frame(r = pulse_r, s = pulse_s, i = pulse_i)

policy_adaptive <- data.frame(
  r = pmax(0, 10000 - 10000 * ((1:30) / 30)^1.5),
  s = seq(1000, 6000,  length.out = 30),
  i = pmin(40e6, 5e6 + 1.2e6 * (1:30))
)
policy_frontloaded <- data.frame(
  r = c(seq(15000, 2000, length.out = 15), rep(500,  15)),
  s = c(seq(8000,   500, length.out = 15), rep(200,  15)),
  i = c(seq(40e6,  5e6, length.out = 15), rep(2e6,  15))
)



results_none           <- simulate_model(state_0, policy_none,           30, params)
results_high_rebate    <- simulate_model(state_0, policy_high_rebate,    30, params)
results_infra          <- simulate_model(state_0, policy_infra,          30, params)
results_aggressive     <- simulate_model(state_0, policy_aggressive,     30, params)
results_moderate       <- simulate_model(state_0, policy_moderate_rebate,30, params)
results_supply_push    <- simulate_model(state_0, policy_supply_push,    30, params)
results_phaseout       <- simulate_model(state_0, policy_phaseout,       30, params)
results_rampup         <- simulate_model(state_0, policy_rampup,         30, params)
results_pulse          <- simulate_model(state_0, policy_pulse,          30, params)
results_adaptive       <- simulate_model(state_0, policy_adaptive,       30, params)
results_frontloaded    <- simulate_model(state_0, policy_frontloaded,    30, params)

# All scenarios together
scenarios <- list(
  results_none, results_high_rebate, results_infra, results_aggressive,
  results_moderate, results_supply_push,
  results_phaseout, results_rampup, results_pulse, results_adaptive, results_frontloaded
)
labels <- c(
  "No Policy", "High Rebates", "Infra Focus", "Aggressive",
  "Moderate Rebate", "Supply Push",
  "Phaseout", "Ramp-Up", "Pulse", "Adaptive", "Front-Loaded"
)

# Time-varying only (for policy control plot)
tv_scenarios <- list(results_phaseout, results_rampup, results_pulse,
                     results_adaptive, results_frontloaded)
tv_labels    <- c("Phaseout", "Ramp-Up", "Pulse", "Adaptive", "Front-Loaded")

set.seed(42)
mc_sims <- lapply(scenarios, function(s) monte_carlo(state_0, s, 30, params, n_sim = 100))

sens_results <- lapply(scenarios, function(s) {
  list(
    beta_0 = sensitivity_analysis(state_0, s, 30, params, "beta_0"),
    beta_1 = sensitivity_analysis(state_0, s, 30, params, "beta_1"),
    beta_2 = sensitivity_analysis(state_0, s, 30, params, "beta_2")
  )
})

pdf("av_model_plots_analyzed.pdf", width = 16, height = 9)
print(plot_slope_chart(scenarios, labels))
print(plot_heatmap(scenarios, labels))
print(plot_fleet_composition(scenarios, labels))
print(plot_demand_vs_supply(scenarios, labels))
print(plot_policy_controls(tv_scenarios, tv_labels))
print(plot_mc_ridges(mc_sims, labels))
print(plot_sensitivity_diverging(sens_results, labels))
print(plot_phase_portrait(scenarios, labels))
print(plot_total_spend(scenarios, labels))
print(plot_spend_per_av(scenarios, labels))
print(plot_ssm1_calibration(params))
print(plot_ssm2_tornado(state_0, policy_aggressive, 30, params))
print(plot_ssm3_mc_bands(scenarios, labels, mc_sims))
print(plot_ssm4_surface(params, scenarios, labels))
dev.off()

scenarios_combined <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
  df          <- scenarios[[i]]
  df$scenario <- labels[i]
  df
}))
write.csv(scenarios_combined, "av_all_scenarios.csv", row.names = FALSE)

cat("\nYear-30 AV shares by scenario:\n")
for (i in seq_along(scenarios)) {
  cat(sprintf("  %-20s %.1f%%\n", labels[i], scenarios[[i]]$pct_AV[31] * 100))
}