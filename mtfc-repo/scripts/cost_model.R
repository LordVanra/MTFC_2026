# =============================================================================
#  scripts/cost_model.R  —  DERIVED POLICY COST MODEL
# =============================================================================
#  Replaces (a) the hardcoded POLICY_COSTS table that lived inside
#  ConvectionDiffusion/plot_results.r, and (b) two mutually inconsistent and
#  dimensionally incorrect in-line expressions in StateModel/plot_model.r.
#
#  What was wrong
#  --------------
#  plot_model.r line 410 / 527:
#      annual_spend = r * (adoption * S) + s + i
#          `s` is a per-AV manufacturer subsidy in $/AV, but was added as a
#          bare scalar - i.e. the model charged $8,000 total per YEAR for
#          Supply Push instead of $8,000 per vehicle produced. This is exactly
#          the scenario the paper concludes is most cost-effective, so the
#          headline "5x more efficient" result was produced by under-charging
#          the winning scenario by a factor of ~(adoption * S).
#
#  plot_model.r line 450:
#      ... + i * S * Price_AV / 1000
#          `i` is already dollars per year (e.g. 40e6). Multiplying by
#          S * Price_AV / 1000 (~1.3e6) inflated infrastructure spend by six
#          orders of magnitude in the spend-per-AV figure only.
#
#  Correct accounting
#  ------------------
#      units_sold(t)  = adoption(t) * S_BOS          [AV/yr]
#      rebate_cost(t) = r(t) * units_sold(t)         [$/yr]
#      supply_cost(t) = s(t) * units_sold(t)         [$/yr]
#      infra_cost(t)  = i(t)                         [$/yr]
#      annual(t)      = rebate + supply + infra
# =============================================================================

annual_policy_cost <- function(r, s, i, adoption, S_BOS) {
  units_sold <- adoption * S_BOS
  rebate <- r * units_sold
  supply <- s * units_sold
  infra  <- i
  list(
    rebate = rebate,
    supply = supply,
    infra  = infra,
    total  = rebate + supply + infra,
    units_sold = units_sold
  )
}

# Add cost columns to a simulate_model() output data frame.
add_costs <- function(df, S_BOS, discount_rate = 0) {
  cp <- annual_policy_cost(df$r, df$s, df$i, df$adoption, S_BOS)
  disc <- 1 / (1 + discount_rate)^df$year

  df$units_sold     <- cp$units_sold
  df$cost_rebate    <- cp$rebate
  df$cost_supply    <- cp$supply
  df$cost_infra     <- cp$infra
  df$annual_spend   <- cp$total
  df$annual_spend_d <- cp$total * disc

  # year 0 carries no policy spend (state is initial condition)
  df$annual_spend[df$year == 0]   <- 0
  df$annual_spend_d[df$year == 0] <- 0

  df$cumulative_spend   <- cumsum(df$annual_spend)
  df$cumulative_spend_B <- df$cumulative_spend / 1e9
  df$spend_per_av       <- df$cumulative_spend / pmax(df$A, 1)
  df
}

# Total 30-year cost per scenario, in $B. Replaces the hardcoded table.
scenario_cost_table <- function(scenarios, labels, S_BOS, discount_rate = 0) {
  rows <- lapply(seq_along(scenarios), function(k) {
    d <- add_costs(scenarios[[k]], S_BOS, discount_rate)
    data.frame(
      scenario       = labels[k],
      total_cost_B   = sum(d$annual_spend)   / 1e9,
      total_cost_d_B = sum(d$annual_spend_d) / 1e9,
      rebate_share   = sum(d$cost_rebate) / max(sum(d$annual_spend), 1),
      supply_share   = sum(d$cost_supply) / max(sum(d$annual_spend), 1),
      infra_share    = sum(d$cost_infra)  / max(sum(d$annual_spend), 1),
      final_share_AV = d$pct_AV[nrow(d)],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
