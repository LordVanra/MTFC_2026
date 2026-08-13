library(data.table)

emissions_data <- fread("emission_data.csv")

pollutants <- c("PM2.5", "NOx", "CO2")
for (p in pollutants) {
  col_km <- paste0(p, "_g_per_km")
  if (!col_km %in% names(emissions_data)) stop(paste("ERROR: Column", col_km, "not found"))
  emissions_data[[col_km]] <- as.numeric(emissions_data[[col_km]])
  emissions_data[[paste0(p, "_g_per_m")]] <- emissions_data[[col_km]] / 1000.0
}

vehicle     <- "Heavy_Duty_Diesel"
vehicle_row <- emissions_data[Vehicle_Type == vehicle]

E_PM25         <- vehicle_row$PM2.5_g_per_m
vehicles_per_s <- 1800 / 3600
E_eff          <- E_PM25 * vehicles_per_s

u     <- 3.5
w     <- 0.1
D_h   <- 50.0
D_z   <- 10.0
k_dep <- 3e-5
dx    <- 5.0
dz    <- 1.0

dt <- min(dx / u * 0.9, 0.5 / (D_h/dx^2 + D_z/dz^2) * 0.9)

T_sim           <- 600.0
x_domain        <- c(0, 1000)
z_domain        <- c(0, 100)
source_region_x <- c(0.0, 100.0)

solve_2d <- function(u, w, D_h, D_z, k_dep, source_region_x,
                     z_domain=c(0, 100), x_domain=c(0, 1000),
                     dx=5.0, dz=1.0, dt, T=600.0, E_eff) {

  nx <- as.integer((x_domain[2] - x_domain[1]) / dx) + 1
  nz <- as.integer((z_domain[2] - z_domain[1]) / dz) + 1
  nt <- as.integer(T / dt)
  x  <- seq(x_domain[1], x_domain[2], length.out = nx)
  z  <- seq(z_domain[1], z_domain[2], length.out = nz)
  C  <- matrix(0.0, nrow = nz, ncol = nx)

  S    <- matrix(0.0, nrow = nz, ncol = nx)
  mask <- (x >= source_region_x[1]) & (x <= source_region_x[2])
  S[1, mask] <- E_eff / dz

  for (n in seq_len(nt)) {
    Cn     <- C
    adv_x  <- -u  * dt / (2*dx) * (cbind(Cn[, 2:nx], Cn[, nx]) - cbind(Cn[, 1], Cn[, 1:(nx-1)]))
    adv_z  <- -w  * dt / (2*dz) * (rbind(Cn[2:nz, ], Cn[nz, ]) - rbind(Cn[1, ], Cn[1:(nz-1), ]))
    diff_x <- D_h * dt / dx^2   * (cbind(Cn[, 2:nx], Cn[, nx]) - 2*Cn + cbind(Cn[, 1], Cn[, 1:(nx-1)]))
    diff_z <- D_z * dt / dz^2   * (rbind(Cn[2:nz, ], Cn[nz, ]) - 2*Cn + rbind(Cn[1, ], Cn[1:(nz-1), ]))

    C <- Cn + adv_x + adv_z + diff_x + diff_z + S * dt - k_dep * Cn * dt

    C[, 1]      <- 0.0
    C[, nx]     <- C[, nx-1]
    C[1, !mask] <- C[2, !mask]
    C[nz, ]     <- 0.0
    C[C < 0]    <- 0.0
  }

  list(x=x, z=z, C=C)
}

result_2d <- solve_2d(u=u, w=w, D_h=D_h, D_z=D_z, k_dep=k_dep,
                      source_region_x=source_region_x,
                      z_domain=z_domain, x_domain=x_domain,
                      dx=dx, dz=dz, dt=dt, T=T_sim, E_eff=E_eff)
x2 <- result_2d$x
z2 <- result_2d$z
C2 <- result_2d$C

pdf("pm25_2d_results.pdf", width=8, height=6)

image(x2, z2, t(C2), col=hcl.colors(256, "plasma"),
      xlab="Downwind distance (m)", ylab="Height (m)",
      main="2D PM2.5 plume")
contour(x2, z2, t(C2), add=TRUE, col="white", nlevels=5)

ix_170  <- which.min(abs(x2 - 170))
ix_peak <- which.max(colSums(C2))
plot(C2[, ix_170], z2, type="l", col="steelblue",
     xlab="Concentration (g/m3)", ylab="Height (m)",
     main=paste("Vertical profile at x =", round(x2[ix_170]), "m"))
lines(C2[, ix_peak], z2, col="firebrick", lty=2)
legend("topright", legend=c(paste("x =", round(x2[ix_170]), "m"),
                             paste("peak x =", round(x2[ix_peak]), "m")),
       col=c("steelblue","firebrick"), lty=c(1,2))

plot(x2, C2[1,], type="l", col="darkgreen",
     xlab="Distance (m)", ylab="Ground concentration (g/m3)",
     main="Ground-level PM2.5 along downwind distance")
grid()

dev.off()

baseline     <- list(u=u, w=w, D_h=D_h, D_z=D_z, k_dep=k_dep)
param_sweeps <- list(
  u     = c(1.0, 3.5, 6.0),
  w     = c(0.0, 0.1, 0.3),
  D_h   = c(10.0, 50.0, 100.0),
  D_z   = c(5.0, 10.0, 30.0),
  k_dep = c(1e-5, 3e-5, 1e-4)
)

run_sweep <- function(param_name, values, baseline, source_region_x, E_eff,
                      z_domain, x_domain, dx, dz, T_sim) {
  peak_total  <- numeric(length(values))
  peak_ground <- numeric(length(values))
  for (i in seq_along(values)) {
    args <- baseline
    args[[param_name]] <- values[i]
    dt_i <- min(dx / args$u * 0.9, 0.5 / (args$D_h/dx^2 + args$D_z/dz^2) * 0.9)
    res  <- solve_2d(u=args$u, w=args$w, D_h=args$D_h, D_z=args$D_z,
                     k_dep=args$k_dep, source_region_x=source_region_x,
                     z_domain=z_domain, x_domain=x_domain,
                     dx=dx, dz=dz, dt=dt_i, T=T_sim, E_eff=E_eff)
    peak_total[i]  <- max(res$C)
    peak_ground[i] <- max(res$C[1, ])
  }
  list(peak_total=peak_total, peak_ground=peak_ground)
}

pdf("sensitivity_2d.pdf", width=8, height=6)
for (param_name in names(param_sweeps)) {
  values  <- param_sweeps[[param_name]]
  results <- run_sweep(param_name, values, baseline, source_region_x, E_eff,
                       z_domain, x_domain, dx, dz, T_sim)
  plot(values, results$peak_total, type="o", col="steelblue", pch=16,
       xlab=param_name, ylab="Peak concentration (g/m3)",
       main=paste("Sensitivity: full-domain peak vs", param_name))
  grid()
  plot(values, results$peak_ground, type="o", col="firebrick", pch=16,
       xlab=param_name, ylab="Ground peak concentration (g/m3)",
       main=paste("Sensitivity: ground-level peak vs", param_name))
  grid()
}
dev.off()