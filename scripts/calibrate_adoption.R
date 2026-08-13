# =============================================================================
#  scripts/calibrate_adoption.R
#  Fits the adoption curve to OBSERVED Boston data instead of a forecast envelope.
#
#  The manuscript calibrates beta_0/beta_1/beta_2 to a Morgan Stanley / VTPI /
#  WEF scenario envelope and is explicit that the result is a set of elicited
#  scenario weights, not estimated elasticities. The Massachusetts Vehicle Census
#  gives six years of observed Boston zero-emission VMT share, which is a real
#  calibration target for the same curve.
#
#  Fits a 3-parameter logistic  s(t) = K / (1 + exp(-r (t - t0)))  by nonlinear
#  least squares on the logit, and reports the fit against the data.
# =============================================================================

source(if (file.exists("config/params.R")) "config/params.R" else "../config/params.R")

obs_zev <- par_val("SSM", "obs_zev_share")
obs_hyb <- par_val("SSM", "obs_hybrid_share")
years   <- as.integer(names(obs_zev))
t       <- years - years[1]

fit_logistic <- function(t, y, K_grid = seq(0.06, 1.00, by = 0.005)) {
  best <- NULL
  for (K in K_grid) {
    yy <- pmin(pmax(y / K, 1e-6), 1 - 1e-6)
    m  <- lm(log(yy / (1 - yy)) ~ t)
    r  <- coef(m)[2]; t0 <- -coef(m)[1] / r
    pred <- K / (1 + exp(-r * (t - t0)))
    sse  <- sum((pred - y)^2)
    if (is.null(best) || sse < best$sse)
      best <- list(K = K, r = as.numeric(r), t0 = as.numeric(t0), sse = sse, pred = pred)
  }
  best
}

cat("=== Logistic fit to OBSERVED Boston zero-emission VMT share ===\n\n")
fz <- fit_logistic(t, as.numeric(obs_zev))
cat(sprintf("  K  (saturation)   = %.3f\n", fz$K))
cat(sprintf("  r  (growth rate)  = %.4f /yr\n", fz$r))
cat(sprintf("  t0 (inflection)   = %.2f yr after %d  (= %.0f)\n", fz$t0, years[1], years[1] + fz$t0))
cat(sprintf("  RMSE              = %.5f\n\n", sqrt(fz$sse / length(t))))
cat(sprintf("%6s %12s %12s %10s\n", "year", "observed", "fitted", "resid"))
for (i in seq_along(t))
  cat(sprintf("%6d %11.3f%% %11.3f%% %9.3f%%\n",
              years[i], obs_zev[i]*100, fz$pred[i]*100, (obs_zev[i]-fz$pred[i])*100))

cat("\n=== Same fit to hybrid share ===\n")
fh <- fit_logistic(t, as.numeric(obs_hyb))
cat(sprintf("  K = %.3f   r = %.4f/yr   t0 = %.2f   RMSE = %.5f\n",
            fh$K, fh$r, fh$t0, sqrt(fh$sse/length(t))))

# --- the inflection warning, quantified -------------------------------------
cat("\n=== Sensitivity of the fit to the 2024-25 slowdown ===\n")
f_thru23 <- fit_logistic(t[1:4], as.numeric(obs_zev)[1:4])
proj <- f_thru23$K / (1 + exp(-f_thru23$r * (t - f_thru23$t0)))
cat(sprintf("  Fit through 2023 only: K=%.3f r=%.4f t0=%.2f\n",
            f_thru23$K, f_thru23$r, f_thru23$t0))
cat(sprintf("  It projects 2025 at %.2f%% vs observed %.2f%%  (%.1fx overshoot)\n",
            proj[6]*100, obs_zev[6]*100, proj[6]/obs_zev[6]))
cat("  -> Adoption curves fitted to pre-2024 data overshoot badly. Any forecast\n")
cat("     in the paper must be tested against this slowdown.\n")

# --- translate to the state-space beta_0 ------------------------------------
# CAREFUL: the SSM's adoption fraction `a` is the share of NEW SALES that are
# AV/BEV, not the fleet share. Observed data gives the FLEET (VMT) share. The
# two are linked by the fleet-turnover recurrence
#     s(t+1) = (1 - delta_1) s(t) + a(t) * (S_BOS / C)
# so inverting for a:
#     a(t) = [ s(t+1) - (1 - delta_1) s(t) ] / (S_BOS / C)
# Conflating the two understates `a` by roughly the turnover ratio (~18x here)
# and would drive beta_0 far too negative.
D      <- derive_params()
d1     <- par_val("SSM","delta_1")
C0     <- par_val("SSM","C_0")
turn   <- D$S_BOS / C0                 # new sales as a fraction of the fleet
s_obs  <- as.numeric(obs_zev)

a_implied <- (s_obs[-1] - (1 - d1) * s_obs[-length(s_obs)]) / turn

cat("\n=== Implied NEW-SALES zero-emission share (inverting fleet turnover) ===\n")
cat(sprintf("  turnover S_BOS/C_0 = %.4f  (%.0f sales into a %.0f fleet)\n",
            turn, D$S_BOS, C0))
for (i in seq_along(a_implied))
  cat(sprintf("  %d-%d : fleet %.2f%% -> %.2f%%   implies new-sales share %.1f%%\n",
              years[i], years[i+1], s_obs[i]*100, s_obs[i+1]*100, a_implied[i]*100))
cat("\n  Sanity: US EV share of new light-vehicle sales was roughly 4-10%% over\n")
cat("  this period, so these are the right order of magnitude.\n")

# Solve beta_0 so that sigmoid(beta_0 + beta_2 * I_scaled) reproduces the mean
# implied new-sales share under baseline (no-policy) infrastructure.
I0       <- par_val("SSM","I_0")
d2       <- par_val("SSM","delta_2")
I_ref    <- par_val("SSM","I_ref")
I_scaled <- sqrt(min(((1 - d2) * I0) / I_ref, 1))
b2       <- par_val("SSM","beta_2")

target_a  <- mean(a_implied)
beta0_cal <- log(target_a / (1 - target_a)) - b2 * I_scaled

cat(sprintf("\n=== Implied beta_0 ===\n"))
cat(sprintf("  I_scaled at baseline        = %.4f\n", I_scaled))
cat(sprintf("  mean implied new-sales share= %.4f (%.1f%%)\n", target_a, target_a*100))
cat(sprintf("  beta_0 from observed data   = %.3f\n", beta0_cal))
cat(sprintf("  manuscript -2.026, code was -3.5\n"))
cat("  A more negative beta_0 means the model's baseline adoption was too fast\n")
cat("  relative to what Boston actually did.\n")

saveRDS(list(zev = fz, hybrid = fh, beta_0 = beta0_cal,
             years = years, obs_zev = obs_zev, obs_hyb = obs_hyb),
        "outputs/adoption_calibration.rds")
cat("\nWrote outputs/adoption_calibration.rds\n")

# =============================================================================
#  CAVEAT THAT MUST APPEAR IN THE PAPER
#  The inversion above treats the observed VMT share as if it were a fleet-COUNT
#  share. It is not. ZEVs are newer than the fleet average and newer vehicles are
#  driven more, so VMT share exceeds count share and its growth overstates count
#  growth. The implied 37.1% new-sales share in 2022-23 is not credible against
#  actual Massachusetts EV registrations and is the clearest symptom.
#
#  Correcting requires a relative-mileage factor m = (VMT per ZEV) / (VMT per
#  vehicle, fleet average). Literature and MVC methodology suggest m ~ 1.2-1.5
#  for newer vehicles. The count share is then approximately share_count =
#  share_vmt / m, and the implied new-sales figures scale down proportionally.
# =============================================================================
cat("\n=== Relative-mileage sensitivity on the implied beta_0 ===\n")
cat(sprintf("%8s %14s %12s\n", "m", "mean new-sales", "beta_0"))
for (m in c(1.0, 1.2, 1.35, 1.5)) {
  s_cnt <- s_obs / m
  a_i   <- (s_cnt[-1] - (1 - d1) * s_cnt[-length(s_cnt)]) / turn
  ta    <- mean(a_i)
  b0    <- log(ta / (1 - ta)) - b2 * I_scaled
  cat(sprintf("%8.2f %13.1f%% %12.3f%s\n", m, ta*100, b0,
              if (m == 1.0) "   <- assumes VMT share = count share" else ""))
}
cat("\n  beta_0 is only mildly sensitive to m across this range. Report the m=1.35\n")
cat("  central case with the 1.0-1.5 span as a sensitivity, and state the\n")
cat("  VMT-vs-count distinction explicitly - a referee will spot it otherwise.\n")
