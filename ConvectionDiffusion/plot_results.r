library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(forcats)

policy_raw <- fread("av_all_scenarios.csv")
TARGET_YEARS <- c(0, 10, 20, 30)

load_emission_factors <- function(filepath) {
  data <- fread(filepath)
  for (p in c("PM2.5", "NOx", "CO2")) {
    col <- paste0(p, "_g_per_km")
    if (!col %in% names(data)) stop(paste("Missing column:", col))
    data[[col]] <- as.numeric(data[[col]])
  }
  return(data)
}

get_pm25_per_km <- function(emission_data, vehicle_type, is_av) {
  if (is_av) {
    row <- emission_data[Vehicle_Type == "Autonomous_EV_Lifecycle"]
  } else {
    row <- emission_data[Vehicle_Type == vehicle_type]
  }
  if (nrow(row) == 0) stop(paste("Vehicle type not found:", vehicle_type))
  return(row$PM2.5_g_per_km)
}

load_flow_matrix <- function(filepath) {
  raw <- fread(filepath, header = TRUE, check.names = FALSE)
  cell_ids <- as.character(raw[[1]])
  mat_data <- raw[, -1, with = FALSE]
  if (ncol(mat_data) != nrow(mat_data)) {
    raw <- fread(filepath, header = FALSE, skip = 1)
    cell_ids <- as.character(raw[[1]])
    mat_data <- raw[, -1, with = FALSE]
  }
  mat <- as.matrix(mat_data)
  nr <- nrow(mat); nc <- ncol(mat)
  rownames(mat) <- cell_ids[seq_len(nr)]
  colnames(mat) <- cell_ids[seq_len(nc)]
  mat[is.na(mat)] <- 0
  mat[mat < 0]    <- 0
  return(mat)
}

load_distance_matrix <- function(filepath) {
  raw <- fread(filepath)
  cell_ids <- as.character(raw[[1]])
  mat <- as.matrix(raw[, -1, with = FALSE])
  rownames(mat) <- cell_ids
  colnames(mat) <- cell_ids
  mat[is.na(mat)] <- 0
  mat[mat < 0]    <- 0
  return(mat)
}

solve_2d <- function(u = 3.5, w = 0.1, D_h = 50, D_z = 10, k_dep = 3e-5,
                     E_eff, road_length_m,
                     z_domain = c(0, 100), x_domain = c(0, 1000),
                     dx = 5, dz = 1, T = 600) {
  source_region_x <- c(0, min(road_length_m, x_domain[2]))
  nx <- as.integer((x_domain[2] - x_domain[1]) / dx) + 1
  nz <- as.integer((z_domain[2] - z_domain[1]) / dz) + 1
  x  <- seq(x_domain[1], x_domain[2], length.out = nx)
  z  <- seq(z_domain[1], z_domain[2], length.out = nz)
  dt <- min(dx / u * 0.9, 0.5 / (D_h/dx^2 + D_z/dz^2) * 0.9)
  nt <- as.integer(T / dt)
  C <- matrix(0, nrow = nz, ncol = nx)
  S <- matrix(0, nrow = nz, ncol = nx)
  mask <- (x >= source_region_x[1]) & (x <= source_region_x[2])
  S[1, mask] <- E_eff / dz
  for (n in seq_len(nt)) {
    Cn     <- C
    adv_x  <- -u  * dt / (2*dx) * (cbind(Cn[, 2:nx], Cn[, nx]) - cbind(Cn[, 1], Cn[, 1:(nx-1)]))
    adv_z  <- -w  * dt / (2*dz) * (rbind(Cn[2:nz,], Cn[nz,]) - rbind(Cn[1,], Cn[1:(nz-1),]))
    diff_x <- D_h * dt / dx^2   * (cbind(Cn[, 2:nx], Cn[, nx]) - 2*Cn + cbind(Cn[, 1], Cn[, 1:(nx-1)]))
    diff_z <- D_z * dt / dz^2   * (rbind(Cn[2:nz,], Cn[nz,]) - 2*Cn + rbind(Cn[1,], Cn[1:(nz-1),]))
    C <- Cn + adv_x + adv_z + diff_x + diff_z + S * dt - k_dep * Cn * dt
    C[, 1]      <- 0
    C[, nx]     <- C[, nx-1]
    C[1, !mask] <- C[2, !mask]
    C[nz, ]     <- 0
    C[C < 0]    <- 0
  }
  list(x = x, z = z, C = C,
       ground_max  = max(C[1, ]),
       domain_max  = max(C),
       ground_mean = mean(C[1, C[1,] > 0]))
}

run_scenario <- function(flow_matrix, dist_matrix, av_fraction,
                         emission_data, vehicle_type = "Passenger_Car_Gasoline",
                         scenario_label = "",
                         months_to_seconds = 30 * 24 * 3600) {
  pm25_normal <- get_pm25_per_km(emission_data, vehicle_type, is_av = FALSE)
  pm25_av     <- get_pm25_per_km(emission_data, vehicle_type, is_av = TRUE)
  cell_ids <- rownames(flow_matrix)
  results  <- list()
  for (i in seq_along(cell_ids)) {
    for (j in seq_along(cell_ids)) {
      total_veh <- flow_matrix[i, j]
      if (total_veh == 0) next
      dist_km <- dist_matrix[cell_ids[i], cell_ids[j]]
      if (is.na(dist_km) || dist_km == 0) next
      n_av     <- total_veh * av_fraction
      n_normal <- total_veh * (1 - av_fraction)
      veh_per_s_av     <- n_av     / months_to_seconds
      veh_per_s_normal <- n_normal / months_to_seconds
      E_eff <- (pm25_normal * veh_per_s_normal + pm25_av * veh_per_s_av) / 1000
      road_length_m <- dist_km * 1000
      sim <- solve_2d(E_eff = E_eff, road_length_m = road_length_m, dx = 20, dz = 5, T = 300)
      results[[length(results) + 1]] <- data.frame(
        origin        = cell_ids[i],
        destination   = cell_ids[j],
        dist_km       = dist_km,
        total_veh     = total_veh,
        n_av          = n_av,
        n_normal      = n_normal,
        av_fraction   = av_fraction,
        scenario      = scenario_label,
        E_eff         = E_eff,
        ground_max    = sim$ground_max,
        domain_max    = sim$domain_max,
        ground_mean   = sim$ground_mean
      )
    }
  }
  do.call(rbind, results)
}

theme_disp <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.background   = element_rect(fill = "#FFFFFF", color = NA),
      panel.background  = element_rect(fill = "#F8F9FA", color = NA),
      panel.grid.major  = element_line(color = "#E2E8F0", linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(face = "bold", size = 14, color = "#1A202C", margin = margin(b = 4)),
      plot.subtitle     = element_text(color = "#718096", size = 9, margin = margin(b = 10)),
      axis.text         = element_text(color = "#718096", size = 9),
      axis.title        = element_text(color = "#1A202C", size = 10),
      strip.background  = element_rect(fill = "#E2E8F0", color = NA),
      strip.text        = element_text(face = "bold", color = "#1A202C"),
      legend.background = element_rect(fill = "#F8F9FA", color = NA),
      legend.title      = element_text(color = "#1A202C"),
      legend.text       = element_text(color = "#718096"),
      plot.margin       = margin(14, 14, 14, 14)
    )
}

plot_ground_max_heatmap <- function(all_results) {
  df <- all_results %>%
    mutate(log_conc = log10(ground_max + 1e-15))
  ggplot(df, aes(x = destination, y = origin, fill = log_conc)) +
    geom_tile(color = "#E2E8F0", linewidth = 0.3) +
    scale_fill_gradientn(
      colors   = c("#F0F9FF", "#90CDF4", "#3182CE", "#C53030"),
      na.value = "#F8F9FA",
      name     = "Ground PM2.5\n(log10 g/m3)",
      guide    = guide_colorbar(barwidth = 0.8, barheight = 12)
    ) +
    facet_wrap(~scenario) +
    labs(title    = "Ground-Level Peak PM2.5 by Flow Pair",
         subtitle = "Each cell = one origin-destination corridor, color = peak ground concentration",
         x = "Destination", y = "Origin") +
    theme_disp() +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank())
}
plot_network_burden <- function(all_results) {
  df <- all_results %>%
    group_by(scenario, av_fraction) %>%
    summarise(
      total_ground_max = sum(ground_max),
      total_emission   = sum(E_eff),
      mean_ground      = mean(ground_max),
      .groups = "drop"
    )
  
  baseline <- df$total_emission[which.max(df$av_fraction == min(df$av_fraction))]
  df$pct_reduction <- (1 - df$total_emission / baseline) * 100

  # Sort scenarios by magnitude of reduction
  df <- df %>%
    mutate(scenario = fct_reorder(scenario, pct_reduction))

  ggplot(df, aes(x = scenario, y = pct_reduction, fill = pct_reduction)) +
    geom_col(width = 0.6, alpha = 0.9) +
    geom_text(aes(label = ifelse(pct_reduction == 0, "baseline",
                                 paste0("-", round(abs(pct_reduction), 1), "%"))),
              vjust = -0.5, fontface = "bold", color = "#1A202C", size = 3.5) +
    scale_fill_gradientn(
      colors = c("#C8DCF0", "#90B8D8", "#3182CE", "#1A4A7A"),
      name   = "% reduction\nvs baseline",
      guide  = guide_colorbar(barwidth = 1.2, barheight = 10,
                              ticks.colour = NA, frame.colour = NA)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15)),
      labels = function(x) paste0(x, "%")
    ) +
    labs(title    = "Total Network Emission Rate by Policy Scenario",
         subtitle = "Sum of E_eff across all active flow corridors; % reduction vs baseline scenario",
         x = "Policy Scenario", y = "% Reduction from Baseline") +
    theme_disp()
}

plot_corridor_av_benefit <- function(all_results) {
  baseline_frac <- min(all_results$av_fraction)
  baseline_scen <- all_results$scenario[which.min(all_results$av_fraction)]
  baseline_df <- all_results %>%
    filter(av_fraction == baseline_frac, scenario == baseline_scen) %>%
    select(origin, destination, ground_max_base = ground_max, dist_km, total_veh)
  df <- all_results %>%
    filter(!(av_fraction == baseline_frac & scenario == baseline_scen)) %>%
    left_join(baseline_df, by = c("origin", "destination", "dist_km", "total_veh")) %>%
    mutate(pct_reduction = (1 - ground_max / ground_max_base) * 100)
  ggplot(df, aes(x = dist_km, y = pct_reduction, color = log10(total_veh + 1), size = total_veh)) +
    geom_point(alpha = 0.75) +
    scale_color_gradientn(colors = c("#3182CE", "#D97706", "#C53030"), name = "Volume\n(log10 veh)") +
    scale_size_continuous(range = c(1, 6), guide = "none") +
    facet_wrap(~scenario) +
    labs(title    = "PM2.5 Reduction by Corridor vs Baseline",
         subtitle = "Each point = one flow pair; x = distance, y = % reduction in ground-level PM2.5",
         x = "Corridor Distance (km)", y = "% PM2.5 Reduction") +
    theme_disp()
}

plot_busiest_plume <- function(flow_matrix, dist_matrix, year_rows,
                               emission_data, vehicle_type = "Passenger_Car_Gasoline") {
  pm25_normal <- get_pm25_per_km(emission_data, vehicle_type, is_av = FALSE)
  pm25_av     <- get_pm25_per_km(emission_data, vehicle_type, is_av = TRUE)
  idx     <- which(flow_matrix == max(flow_matrix), arr.ind = TRUE)
  orig    <- rownames(flow_matrix)[idx[1,1]]
  dest    <- colnames(flow_matrix)[idx[1,2]]
  total   <- flow_matrix[idx[1,1], idx[1,2]]
  dist_km <- dist_matrix[orig, dest]
  plume_list <- list()
  for (k in seq_len(nrow(year_rows))) {
    av_frac  <- year_rows$pct_AV[k]
    scen_lbl <- year_rows$scenario[k]
    veh_per_s_av     <- (total * av_frac)     / (30 * 24 * 3600)
    veh_per_s_normal <- (total * (1-av_frac)) / (30 * 24 * 3600)
    E_eff <- (pm25_normal * veh_per_s_normal + pm25_av * veh_per_s_av) / 1000
    sim <- solve_2d(E_eff = E_eff, road_length_m = dist_km * 1000, dx = 20, dz = 5, T = 300)
    df_plume <- expand.grid(x = sim$x, z = sim$z)
    df_plume$conc  <- as.vector(t(sim$C))
    df_plume$label <- scen_lbl
    plume_list[[k]] <- df_plume
  }
  df_all <- do.call(rbind, plume_list) %>%
    mutate(label = factor(label, levels = unique(year_rows$scenario)))
  ggplot(df_all, aes(x = x, y = z, fill = conc)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(
      colors = c("#F0F9FF", "#90CDF4", "#3182CE", "#744210", "#C53030"),
      name   = "PM2.5\n(g/m3)",
      trans  = "sqrt",
      guide  = guide_colorbar(barwidth = 0.8, barheight = 12)
    ) +
    facet_wrap(~label, ncol = 1) +
    labs(title    = paste0("PM2.5 Plume - Busiest Corridor (", orig, " - ", dest, ")"),
         subtitle = paste0(formatC(total, format="f", digits=0, big.mark=","),
                           " vehicles/month, ", round(dist_km, 1), " km"),
         x = "Downwind Distance (m)", y = "Height (m)") +
    theme_disp() +
    theme(panel.grid = element_blank())
}

plot_ground_profile <- function(flow_matrix, dist_matrix, year_rows,
                                emission_data, vehicle_type = "Passenger_Car_Gasoline") {
  pm25_normal <- get_pm25_per_km(emission_data, vehicle_type, is_av = FALSE)
  pm25_av     <- get_pm25_per_km(emission_data, vehicle_type, is_av = TRUE)
  idx     <- which(flow_matrix == max(flow_matrix), arr.ind = TRUE)
  orig    <- rownames(flow_matrix)[idx[1,1]]
  dest    <- colnames(flow_matrix)[idx[1,2]]
  total   <- flow_matrix[idx[1,1], idx[1,2]]
  dist_km <- dist_matrix[orig, dest]
  n_scen  <- nrow(year_rows)
  COLORS  <- colorRampPalette(c("#3182CE", "#276749"))(n_scen)
  profile_list <- list()
  for (k in seq_len(n_scen)) {
    av_frac  <- year_rows$pct_AV[k]
    scen_lbl <- year_rows$scenario[k]
    veh_per_s_av     <- (total * av_frac)     / (30 * 24 * 3600)
    veh_per_s_normal <- (total * (1-av_frac)) / (30 * 24 * 3600)
    E_eff <- (pm25_normal * veh_per_s_normal + pm25_av * veh_per_s_av) / 1000
    sim <- solve_2d(E_eff = E_eff, road_length_m = dist_km * 1000, dx = 20, dz = 5, T = 300)
    profile_list[[k]] <- data.frame(x = sim$x, conc = sim$C[1, ], label = scen_lbl, color_idx = k)
  }
  df <- do.call(rbind, profile_list) %>%
    mutate(label = factor(label, levels = unique(year_rows$scenario)))
  ggplot(df, aes(x = x, y = conc * 1e6, color = label)) +
    geom_line(linewidth = 1.1, alpha = 0.9) +
    scale_color_manual(values = COLORS, name = "Scenario") +
    scale_y_continuous(labels = function(x) paste0(round(x, 2), " ug")) +
    labs(title    = paste0("Ground-Level PM2.5 Profile - ", orig, " - ", dest),
         subtitle = "Concentration along downwind axis at z = 0 (street level), converted to ug/m3",
         x = "Downwind Distance (m)", y = "PM2.5 Concentration (ug/m3)") +
    theme_disp()
}

plot_map_emission <- function(all_results, zip_coords) {
  node_burden <- all_results %>%
    group_by(zip = origin, scenario) %>%
    summarise(total_E = sum(E_eff), total_veh = sum(total_veh), .groups = "drop") %>%
    left_join(zip_coords, by = "zip") %>%
    filter(!is.na(lat))
  baseline_scen <- all_results$scenario[which.min(all_results$av_fraction)]
  flow_base <- all_results %>%
    filter(scenario == baseline_scen, total_veh > 0) %>%
    left_join(zip_coords %>% rename(orig_lat = lat, orig_lon = lon), by = c("origin" = "zip")) %>%
    left_join(zip_coords %>% rename(dest_lat = lat, dest_lon = lon), by = c("destination" = "zip")) %>%
    filter(!is.na(orig_lat), !is.na(dest_lat))
  ggplot() +
    geom_segment(data = flow_base,
                 aes(x = orig_lon, y = orig_lat, xend = dest_lon, yend = dest_lat,
                     alpha = log10(total_veh + 1)),
                 color = "#3182CE", linewidth = 0.4) +
    scale_alpha_continuous(range = c(0.05, 0.45), guide = "none") +
    ggnewscale::new_scale("alpha") +
    geom_point(data = node_burden,
               aes(x = lon, y = lat, size = total_veh, fill = total_E * 1e6, alpha = 0.9),
               shape = 21, color = "#1A202C", stroke = 0.4) +
    scale_fill_gradientn(colors = c("#EBF8FF", "#90CDF4", "#3182CE", "#C53030"),
                         name = "E_eff\n(ng/m/s)", trans = "sqrt") +
    scale_size_continuous(range = c(3, 14), guide = "none") +
    scale_alpha_continuous(range = c(0.85, 0.95), guide = "none") +
    geom_text(data = zip_coords, aes(x = lon, y = lat, label = zip),
              size = 2.2, color = "#2D3748", fontface = "bold", nudge_y = 0.004) +
    facet_wrap(~scenario) +
    coord_fixed(ratio = 1 / cos(mean(zip_coords$lat) * pi / 180)) +
    labs(title    = "Spatial Distribution of PM2.5 Emission Burden",
         subtitle = "Node size = outbound vehicle volume, color = total emission rate",
         x = "Longitude", y = "Latitude") +
    theme_disp() +
    theme(panel.grid.major = element_line(color = "#E2E8F0", linewidth = 0.3))
}

build_emission_heatmaps <- function(all_results, zip_coords) {
  node_burden <- all_results %>%
    group_by(zip = origin, scenario) %>%
    summarise(total_E = sum(E_eff), .groups = "drop") %>%
    left_join(zip_coords, by = "zip") %>%
    filter(!is.na(lat))
  lon_range <- range(zip_coords$lon)
  lat_range <- range(zip_coords$lat)
  res      <- 200
  sigma    <- 0.018
  lon_grid <- seq(lon_range[1] - 0.03, lon_range[2] + 0.03, length.out = res)
  lat_grid <- seq(lat_range[1] - 0.03, lat_range[2] + 0.03, length.out = res)
  grid     <- expand.grid(lon = lon_grid, lat = lat_grid)
  build_heatmap <- function(nodes_df) {
    z <- numeric(nrow(grid))
    for (k in seq_len(nrow(nodes_df))) {
      dx <- grid$lon - nodes_df$lon[k]
      dy <- grid$lat - nodes_df$lat[k]
      z  <- z + nodes_df$total_E[k] * exp(-(dx^2 + dy^2) / (2 * sigma^2))
    }
    grid$z <- z
    grid
  }
  heat_list <- lapply(unique(node_burden$scenario), function(scen) {
    sub           <- node_burden %>% filter(scenario == scen)
    heat          <- build_heatmap(sub)
    heat$scenario <- scen
    heat
  })
  list(heat_df = do.call(rbind, heat_list), lon_range = lon_range, lat_range = lat_range)
}

plot_map_emission_blended <- function(all_results, zip_coords) {
  hm        <- build_emission_heatmaps(all_results, zip_coords)
  heat_df   <- hm$heat_df
  lon_range <- hm$lon_range
  lat_range <- hm$lat_range
  ggplot() +
    geom_raster(data = heat_df, aes(x = lon, y = lat, fill = z), interpolate = TRUE) +
    scale_fill_gradientn(
      colors = c("#000000", "#1A0A00", "#5A1A08", "#8B3010", "#C87840", "#E8B878",
                "#F8EEE0", "#C8DCF0", "#90B8D8", "#D0ECF8", "#FFFFFF"),
      values = scales::rescale(c(0, 0.08, 0.18, 0.32, 0.50, 0.65,
                                0.75, 0.83, 0.90, 0.96, 1)),
      name   = "E_eff\n(g/m/s)",
      guide  = guide_colorbar(barwidth = 1.2, barheight = 14, ticks.colour = NA, frame.colour = NA)
    ) +
    facet_wrap(~scenario) +
    coord_fixed(ratio = 1 / cos(mean(zip_coords$lat) * pi / 180),
                xlim = lon_range + c(-0.03, 0.03),
                ylim = lat_range + c(-0.03, 0.03)) +
    labs(title    = "Spatial Distribution of PM2.5 Emission Burden (Absolute)",
         subtitle = "Gaussian emission surface per zone, warmer color = more pollution",
         x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 11) +
    theme(
      plot.background   = element_rect(fill = "#FFFFFF", color = NA),
      panel.background  = element_rect(fill = "#F0F4F8", color = NA),
      panel.grid.major  = element_line(color = "#D8E4EE", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(face = "bold", size = 14, color = "#1A202C", margin = margin(b = 4)),
      plot.subtitle     = element_text(color = "#718096", size = 9, margin = margin(b = 10)),
      axis.text         = element_text(color = "#718096", size = 9),
      axis.title        = element_text(color = "#1A202C", size = 10),
      strip.background  = element_rect(fill = "#E2ECF4", color = NA),
      strip.text        = element_text(face = "bold", color = "#1A202C"),
      legend.background = element_rect(fill = "#F0F4F8", color = NA),
      legend.title      = element_text(color = "#1A202C"),
      legend.text       = element_text(color = "#4A5568"),
      plot.margin       = margin(14, 14, 14, 14)
    )
}

plot_map_emission_blended_diff <- function(all_results, zip_coords) {
  hm        <- build_emission_heatmaps(all_results, zip_coords)
  heat_df   <- hm$heat_df
  lon_range <- hm$lon_range
  lat_range <- hm$lat_range
  cell_mean <- heat_df %>%
    group_by(lon, lat) %>%
    summarise(mean_z = mean(z), .groups = "drop")
  heat_df <- heat_df %>%
    left_join(cell_mean, by = c("lon", "lat")) %>%
    mutate(pct_diff = (z - mean_z) / (mean_z + 1e-30) * 100)
  abs_lim <- max(abs(quantile(heat_df$pct_diff, c(0.02, 0.98), na.rm = TRUE)), na.rm = TRUE)
  abs_lim <- ceiling(abs_lim / 5) * 5

  if (!is.finite(abs_lim) || abs_lim < 1) {
    message("plot_map_emission_blended_diff: all scenarios identical at this year - skipping diff plot.")
    return(ggplot() +
      labs(title    = "Spatial PM2.5 Emission Burden: Deviation from Cross-Scenario Mean",
           subtitle = "Not shown: all scenarios have identical AV fractions at this year - no variation to display.") +
      theme_void() +
      theme(plot.title    = element_text(face = "bold", size = 14, color = "#1A202C", margin = margin(b = 8)),
            plot.subtitle = element_text(color = "#718096", size = 11),
            plot.margin   = margin(40, 14, 40, 14))
    )
  }

  pivot <- scales::rescale(
    c(-abs_lim, -abs_lim*0.4, -abs_lim*0.1, 0, abs_lim*0.1, abs_lim*0.4, abs_lim),
    from = c(-abs_lim, abs_lim)
  )
  pivot <- pmin(pmax(pivot, 0), 1)
  pivot <- pivot + seq(0, 1e-9, length.out = length(pivot))

  ggplot() +
    geom_raster(data = heat_df, aes(x = lon, y = lat, fill = pct_diff), interpolate = TRUE) +
    scale_fill_gradientn(
      colors   = c("#1A6B99", "#5BAED6", "#C8DCF0", "#F5F5F5", "#F5C18A", "#D9601A", "#7A1F00"),
      values   = pivot,
      limits   = c(-abs_lim, abs_lim),
      breaks   = c(-abs_lim, -abs_lim/2, 0, abs_lim/2, abs_lim),
      labels   = function(x) paste0(round(x), "%"),
      name     = "vs. mean\nacross\nscenarios",
      guide    = guide_colorbar(barwidth = 1.4, barheight = 16,
                                ticks.colour = "#4A5568", frame.colour = "#4A5568",
                                title.position = "top")
    ) +
    facet_wrap(~scenario) +
    coord_fixed(ratio = 1 / cos(mean(zip_coords$lat) * pi / 180),
                xlim = lon_range + c(-0.03, 0.03),
                ylim = lat_range + c(-0.03, 0.03)) +
    labs(title    = "Spatial PM2.5 Emission Burden: Deviation from Cross-Scenario Mean",
         subtitle = "Orange/red = above average; blue = below average at each map cell",
         x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 11) +
    theme(
      plot.background   = element_rect(fill = "#FFFFFF", color = NA),
      panel.background  = element_rect(fill = "#F5F5F5", color = NA),
      panel.grid.major  = element_line(color = "#D8E4EE", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(face = "bold", size = 14, color = "#1A202C", margin = margin(b = 4)),
      plot.subtitle     = element_text(color = "#718096", size = 9, margin = margin(b = 10)),
      axis.text         = element_text(color = "#718096", size = 9),
      axis.title        = element_text(color = "#1A202C", size = 10),
      strip.background  = element_rect(fill = "#E2ECF4", color = NA),
      strip.text        = element_text(face = "bold", color = "#1A202C"),
      legend.background = element_rect(fill = "#F5F5F5", color = NA),
      legend.title      = element_text(color = "#1A202C", size = 9),
      legend.text       = element_text(color = "#4A5568", size = 8),
      plot.margin       = margin(14, 14, 14, 14)
    )
}


plot_emission_efficiency <- function(all_results) {
  df <- all_results %>%
    mutate(dist_bin = cut(dist_km, breaks = c(0, 2, 5, 10, 20, Inf),
                          labels = c("<2 km", "2-5 km", "5-10 km", "10-20 km", ">20 km"))) %>%
    group_by(scenario, dist_bin) %>%
    summarise(mean_E = mean(E_eff * 1e9), .groups = "drop") %>%
    filter(!is.na(dist_bin))
  n_scen <- length(unique(df$scenario))
  COLORS <- colorRampPalette(c("#C53030", "#3182CE", "#276749"))(n_scen)
  ggplot(df, aes(x = dist_bin, y = mean_E, fill = scenario, group = scenario)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.9) +
    scale_fill_manual(values = COLORS, name = "Scenario") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title    = "Mean Emission Rate by Corridor Distance Band",
         subtitle = "Average E_eff per corridor length category",
         x = "Corridor Distance", y = "Mean E_eff (ng/m/s)") +
    theme_disp()
}

plot_cdf_ground_max <- function(all_results) {
  df <- all_results %>%
    group_by(scenario) %>%
    arrange(ground_max) %>%
    mutate(ecdf_val = seq_along(ground_max) / n()) %>%
    ungroup()
  n_scenarios <- length(unique(df$scenario))
  COLORS <- colorRampPalette(c("#C53030", "#3182CE", "#276749"))(n_scenarios)
  ggplot(df, aes(x = ground_max * 1e6, y = ecdf_val, color = scenario)) +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    scale_color_manual(values = COLORS, name = "Scenario") +
    scale_x_log10(labels = function(x) paste0(signif(x, 2), " ug")) +
    scale_y_continuous(labels = scales::percent) +
    labs(title    = "Cumulative Distribution of Peak Ground-Level PM2.5",
         subtitle = "Each point = one corridor, rightward shift = worse air quality",
         x = "Peak Ground PM2.5 (ug/m3, log scale)", y = "Cumulative % of Corridors") +
    theme_disp()
}

plot_origin_ranking <- function(all_results) {
  df <- all_results %>%
    group_by(origin, scenario) %>%
    summarise(zone_ground_max = sum(ground_max * 1e6), .groups = "drop")
  baseline_scen <- df$scenario[which.min(df$zone_ground_max)]
  baseline_order <- df %>%
    filter(scenario == baseline_scen) %>%
    arrange(desc(zone_ground_max)) %>%
    pull(origin)
  df$origin <- factor(df$origin, levels = rev(unique(baseline_order)))
  n_scen <- length(unique(df$scenario))
  COLORS <- colorRampPalette(c("#C53030", "#3182CE", "#276749"))(n_scen)
  ggplot(df, aes(x = zone_ground_max, y = origin, fill = scenario)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.9) +
    scale_fill_manual(values = COLORS, name = "Scenario") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title    = "PM2.5 Burden by Origin Zone",
         subtitle = "Sum of peak ground concentrations across all outbound corridors per zone",
         x = "Total Ground PM2.5 Contribution (ug/m3)", y = "Origin Zone") +
    theme_disp()
}

plot_volume_vs_conc <- function(all_results) {
  ggplot(all_results, aes(x = total_veh, y = ground_max * 1e6, color = dist_km)) +
    geom_point(alpha = 0.6, size = 1.8) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = function(x) paste0(signif(x, 2), " ug")) +
    scale_color_gradientn(colors = c("#276749", "#D97706", "#C53030"),
                          name = "Distance\n(km)", trans = "sqrt") +
    facet_wrap(~scenario) +
    labs(title    = "Vehicle Volume vs. Peak Ground PM2.5 by Corridor",
         subtitle = "Both axes log scale, color = corridor distance",
         x = "Monthly Vehicle Volume", y = "Peak Ground PM2.5 (ug/m3)") +
    theme_disp()
}

emission_data <- load_emission_factors("emission_data.csv")
flow_matrix   <- load_flow_matrix("cumulative_flow_matrix.csv")
dist_matrix   <- load_distance_matrix("distance_matrix.csv")
zip_coords    <- fread("zip_coords.csv") %>% mutate(zip = as.character(zip))

dir.create("outputs", showWarnings = FALSE)

for (yr in TARGET_YEARS) {
  cat(sprintf("\n===== Year %d =====\n", yr))
  year_rows <- policy_raw[year == yr]
  if (nrow(year_rows) == 0) {
    warning(sprintf("No rows found for year %d - skipping.", yr))
    next
  }
  cat(sprintf("  Scenarios: %s\n", paste(year_rows$scenario, collapse = ", ")))
  all_results <- do.call(rbind, lapply(seq_len(nrow(year_rows)), function(k) {
    av_frac  <- year_rows$pct_AV[k]
    scen_lbl <- year_rows$scenario[k]
    cat(sprintf("    Running: %s (pct_AV = %.4f)\n", scen_lbl, av_frac))
    run_scenario(flow_matrix, dist_matrix,
                 av_fraction    = av_frac,
                 emission_data  = emission_data,
                 scenario_label = scen_lbl)
  }))
  fwrite(all_results, sprintf("outputs/dispersion_results_year%02d.csv", yr))
  pdf_path <- sprintf("outputs/pm25_av_scenarios_year%02d.pdf", yr)
  pdf(pdf_path, width = 14, height = 8)
  plot.new()
  title(main = sprintf("PM2.5 AV Scenarios - Year %d", yr),
        sub  = paste("Policies:", paste(year_rows$scenario, collapse = " | ")),
        cex.main = 2, cex.sub = 1.2)
  print(plot_network_burden(all_results))
  print(plot_ground_max_heatmap(all_results))
  print(plot_corridor_av_benefit(all_results))
  print(plot_busiest_plume(flow_matrix, dist_matrix, year_rows, emission_data))
  print(plot_ground_profile(flow_matrix, dist_matrix, year_rows, emission_data))
  print(plot_map_emission(all_results, zip_coords))
  print(plot_map_emission_blended(all_results, zip_coords))
  print(plot_map_emission_blended_diff(all_results, zip_coords))
  print(plot_emission_efficiency(all_results))
  print(plot_cdf_ground_max(all_results))
  print(plot_origin_ranking(all_results))
  print(plot_volume_vs_conc(all_results))
  dev.off()
  cat(sprintf("  Saved: %s\n", pdf_path))
}

cat("\nDone. PDFs saved to outputs/\n")