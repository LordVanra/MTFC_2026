# =============================================================================
#  scripts/convergence_test.R  —  STEADY-STATE VERIFICATION FOR THE CDM
# =============================================================================
#  The dispersion solver was previously run with T = 300 s over a 1000 m domain
#  at u = 3.5 m/s. Advective transit time across the domain is ~286 s, so the
#  field was sampled while the plume was still crossing it: every reported
#  concentration was a mid-transient value, and both `ground_max` and
#  `ground_mean` were sensitive to the arbitrary stopping time.
#
#  This script sweeps T and reports the relative change in the ground-level
#  field between successive values so that T_sim in config/params.R can be
#  justified rather than asserted. Run it whenever u, D_h, D_z or the domain
#  size changes.
# =============================================================================

.find_config <- function() {
  for (p in c("config/params.R", "../config/params.R")) if (file.exists(p)) return(p)
  stop("config/params.R not found - run from the repository root.")
}
source(.find_config())

# Minimal standalone copy of the solver (no plotting dependencies).
solve_2d_min <- function(u, w, D_h, D_z, k_dep, E_eff, road_length_m,
                         x_max, z_max, dx, dz, T, street_width) {
  nx <- as.integer(x_max / dx) + 1
  nz <- as.integer(z_max / dz) + 1
  x  <- seq(0, x_max, length.out = nx)
  dt <- min(dx / u * 0.9, 0.5 / (D_h/dx^2 + D_z/dz^2) * 0.9)
  nt <- as.integer(T / dt)
  C  <- matrix(0, nrow = nz, ncol = nx)
  S  <- matrix(0, nrow = nz, ncol = nx)
  mask <- x <= min(road_length_m, x_max)
  S[1, mask] <- E_eff / (street_width * dz)
  for (n in seq_len(nt)) {
    Cn     <- C
    adv_x  <- -u  * dt / (2*dx) * (cbind(Cn[, 2:nx], Cn[, nx]) - cbind(Cn[, 1], Cn[, 1:(nx-1)]))
    adv_z  <- -w  * dt / (2*dz) * (rbind(Cn[2:nz,], Cn[nz,]) - rbind(Cn[1,], Cn[1:(nz-1),]))
    diff_x <- D_h * dt / dx^2   * (cbind(Cn[, 2:nx], Cn[, nx]) - 2*Cn + cbind(Cn[, 1], Cn[, 1:(nx-1)]))
    diff_z <- D_z * dt / dz^2   * (rbind(Cn[2:nz,], Cn[nz,]) - 2*Cn + rbind(Cn[1,], Cn[1:(nz-1),]))
    C <- Cn + adv_x + adv_z + diff_x + diff_z + S * dt - k_dep * Cn * dt
    C[, 1] <- 0; C[, nx] <- C[, nx-1]
    C[1, !mask] <- C[2, !mask]; C[nz, ] <- 0; C[C < 0] <- 0
  }
  C[1, ]
}

run_convergence <- function(T_values = c(150, 300, 600, 1200, 2400, 4800)) {
  v <- function(n) par_val("CDM", n)

  # Representative corridor: 1 km, moderate flow, all-ICE fleet.
  ef_ice  <- par_val("EF","ice","exhaust") + par_val("EF","ice","brake") + par_val("EF","ice","tire")
  veh_per_s <- 1800 / 3600
  E_eff   <- ef_ice / 1000 * veh_per_s

  prev <- NULL
  cat(sprintf("Advective transit time across domain: %.0f s\n",
              v("x_domain_max") / v("u")))
  cat(sprintf("%8s %14s %14s %12s\n", "T (s)", "ground_max", "ground_mean", "rel.change"))
  out <- data.frame()
  for (T in T_values) {
    g <- solve_2d_min(v("u"), v("w"), v("D_h"), v("D_z"), v("k_dep"),
                      E_eff, 1000, v("x_domain_max"), v("z_domain_max"),
                      v("dx"), v("dz"), T, v("street_width"))
    rel <- if (is.null(prev)) NA_real_ else max(abs(g - prev)) / max(max(prev), 1e-30)
    cat(sprintf("%8.0f %14.4e %14.4e %12s\n", T, max(g), mean(g),
                if (is.na(rel)) "-" else sprintf("%.4f", rel)))
    out <- rbind(out, data.frame(T = T, ground_max = max(g),
                                 ground_mean = mean(g), rel_change = rel))
    prev <- g
  }
  tol <- v("converge_tol")
  ok  <- out$rel_change[out$T == v("T_sim")]
  cat(sprintf("\nConfigured T_sim = %.0f s; relative change at that T = %s (tol %.3f) -> %s\n",
              v("T_sim"),
              if (length(ok) && !is.na(ok)) sprintf("%.4f", ok) else "not in sweep",
              tol,
              if (length(ok) && !is.na(ok) && ok < tol) "CONVERGED" else "CHECK"))
  invisible(out)
}

if (sys.nframe() == 0) run_convergence()
