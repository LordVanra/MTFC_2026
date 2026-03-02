library(data.table)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(scales)

TRAFFIC_DARK   <- "#FFFFFF"
TRAFFIC_PANEL  <- "#F8F9FA"
TRAFFIC_GRID   <- "#E2E8F0"
TRAFFIC_TEXT   <- "#1A202C"
TRAFFIC_MUTED  <- "#718096"

COL_FREE    <- "#2B6CB0"
COL_CRIT    <- "#D97706"
COL_JAM     <- "#C53030"
COL_FLOW    <- "#276749"
COL_ACCENT  <- "#6B46C1"

CELL_COLORS <- c(
  "#2B6CB0","#2C7A7B","#276749","#744210","#702459",
  "#553C9A","#2D3748","#C53030","#C05621","#B7791F",
  "#4A5568","#1A365D","#1C4532","#3C366B","#742A2A",
  "#7B341E","#5F370E","#2D3748","#065666","#234E52",
  "#1A202C","#2D3748","#322659","#44337A","#63171B",
  "#1D4044","#14532D","#3B0764"
)

theme_traffic <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background    = element_rect(fill = "#FFFFFF",  color = NA),
      panel.background   = element_rect(fill = "#F8F9FA",  color = NA),
      panel.grid.major   = element_line(color = "#E2E8F0", linewidth = 0.4),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(color = "#1A202C", face = "bold",
                                        size = base_size + 3, margin = margin(b = 4)),
      plot.subtitle      = element_text(color = "#718096", size = base_size - 1,
                                        margin = margin(b = 12)),
      plot.caption       = element_text(color = "#718096", size = base_size - 3),
      axis.text          = element_text(color = "#718096", size = base_size - 2),
      axis.title         = element_text(color = "#1A202C", size = base_size - 1),
      axis.ticks         = element_line(color = "#E2E8F0"),
      legend.background  = element_rect(fill = "#F8F9FA",  color = NA),
      legend.text        = element_text(color = "#718096", size = base_size - 2),
      legend.title       = element_text(color = "#1A202C", size = base_size - 1),
      legend.key         = element_rect(fill = NA, color = NA),
      strip.background   = element_rect(fill = "#E2E8F0",  color = NA),
      strip.text         = element_text(color = "#1A202C", face = "bold",
                                        size = base_size - 1),
      plot.margin        = margin(16, 16, 16, 16)
    )
}

congestion_palette <- function() {
  scale_fill_gradientn(
    colors  = c(TRAFFIC_PANEL, COL_FREE, COL_CRIT, COL_JAM),
    values  = c(0, 0.25, 0.65, 1),
    name    = "Density\n(veh/km)",
    guide   = guide_colorbar(barwidth = 1, barheight = 10,
                             ticks.colour = TRAFFIC_GRID,
                             frame.colour = TRAFFIC_GRID)
  )
}

prep_density_df <- function(densityHistory, timeHistory, cellIDs) {
  df <- as.data.frame(densityHistory)
  colnames(df) <- as.character(cellIDs)
  df$time <- timeHistory
  df %>% pivot_longer(-time, names_to = "cell", values_to = "density") %>%
    mutate(cell = as.integer(cell))
}

# Cells on Y, time on X — reveals congestion waves propagating through network
plot_spacetime_heatmap <- function(densityHistory, timeHistory, cellIDs, jamDensity) {
  df <- prep_density_df(densityHistory, timeHistory, cellIDs) %>%
    mutate(pct_jam = density / jamDensity,
           cell    = factor(cell, levels = rev(sort(unique(cell)))))

  ggplot(df, aes(x = time, y = cell, fill = pct_jam)) +
    geom_tile() +
    scale_fill_gradientn(
      colors  = c(TRAFFIC_PANEL, "#0D2137", COL_FREE, COL_CRIT, COL_JAM),
      values  = c(0, 0.05, 0.3, 0.65, 1),
      limits  = c(0, 1),
      labels  = percent_format(),
      name    = "% Jam\nDensity",
      guide   = guide_colorbar(barwidth = 0.8, barheight = 12,
                               ticks.colour = TRAFFIC_GRID,
                               frame.colour = TRAFFIC_GRID)
    ) +
    scale_x_continuous(expand = c(0, 0), labels = function(x) paste0(x, "h")) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title    = "Space–Time Congestion Map",
      subtitle = "Each row = one cell · color = % of jam density · dark = free-flow · red = gridlock",
      x = "Simulation Time", y = "Cell ID"
    ) +
    theme_traffic() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 7))
}

# Distribution of density per cell across time — shows which cells are chronically congested
plot_density_ridges <- function(densityHistory, timeHistory, cellIDs, jamDensity) {
  library(ggridges)
  df <- prep_density_df(densityHistory, timeHistory, cellIDs) %>%
    mutate(pct_jam = density / jamDensity,
           cell    = factor(cell, levels = rev(sort(unique(cell)))))

  # Mean density per cell for coloring
  cell_means <- df %>% group_by(cell) %>% summarise(mean_pct = mean(pct_jam))
  df <- df %>% left_join(cell_means, by = "cell")

  ggplot(df, aes(x = pct_jam * 100, y = cell, fill = mean_pct)) +
    geom_density_ridges(alpha = 0.85, scale = 1.2, bandwidth = 1.5,
                        color = TRAFFIC_GRID, linewidth = 0.3) +
    scale_fill_gradientn(
      colors = c(COL_FREE, COL_CRIT, COL_JAM),
      name   = "Mean %\nJam Density",
      labels = percent_format(scale = 1)
    ) +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title    = "Density Distribution by Cell",
      subtitle = "Spread of occupancy states across simulation — color = chronic congestion level",
      x = "% of Jam Density", y = "Cell ID"
    ) +
    theme_traffic() +
    theme(axis.text.y = element_text(size = 7))
}

plot_network_pulse <- function(densityHistory, timeHistory, cellIDs, cellLengths,
                               flowHistory) {
  totalVeh <- rowSums(sweep(densityHistory, 2, cellLengths, "*"))

  # Total flow per timestep
  totalFlow <- sapply(flowHistory, function(fm) sum(fm))

  df <- data.frame(
    time      = timeHistory,
    occupancy = totalVeh,
    flow      = totalFlow
  ) %>%
    mutate(
      occ_norm  = occupancy / max(occupancy),
      flow_norm = flow      / max(flow)
    ) %>%
    pivot_longer(c(occ_norm, flow_norm), names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric,
                           occ_norm  = "Network Occupancy",
                           flow_norm = "Total Flow"))

  ggplot(df, aes(x = time, y = value, color = metric)) +
    geom_area(data = filter(df, metric == "Network Occupancy"),
              aes(fill = metric), alpha = 0.15, show.legend = FALSE) +
    geom_area(data = filter(df, metric == "Total Flow"),
              aes(fill = metric), alpha = 0.10, show.legend = FALSE) +
    geom_line(linewidth = 1.1, alpha = 0.9) +
    scale_color_manual(values = c("Network Occupancy" = COL_JAM,
                                  "Total Flow"        = COL_FLOW)) +
    scale_fill_manual(values  = c("Network Occupancy" = COL_JAM,
                                  "Total Flow"        = COL_FLOW)) +
    scale_x_continuous(labels = function(x) paste0(x, "h")) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      title    = "Network Pulse",
      subtitle = "Normalized occupancy vs. flow — divergence signals bottleneck formation",
      x = "Time", y = "Normalized Value (% of peak)"
    ) +
    theme_traffic() +
    theme(legend.position = "top")
}

plot_flow_matrix_styled <- function(cumulativeFlowMatrix, cellIDs) {
  df <- as.data.frame(cumulativeFlowMatrix)
  colnames(df) <- as.character(cellIDs)
  df$Origin <- as.character(cellIDs)

  df_melt <- melt(df, id.vars = "Origin", variable.name = "Destination", value.name = "Vehicles") %>%
    mutate(
      Origin      = factor(Origin,      levels = as.character(cellIDs)),
      Destination = factor(Destination, levels = as.character(cellIDs)),
      log_veh     = ifelse(Vehicles > 0, log10(Vehicles), NA),
      label       = ifelse(Vehicles >= 1000,
                           paste0(round(Vehicles/1000, 1), "K"),
                           ifelse(Vehicles > 0, as.character(round(Vehicles)), ""))
    )

  # Only label the top flows so it's not cluttered
  threshold <- quantile(df_melt$Vehicles[df_melt$Vehicles > 0], 0.75, na.rm = TRUE)
  df_melt$label_show <- ifelse(df_melt$Vehicles >= threshold, df_melt$label, "")

  ggplot(df_melt, aes(x = Destination, y = Origin, fill = log_veh)) +
    geom_tile(color = "#E2E8F0", linewidth = 0.4) +
    geom_text(aes(label = label_show), size = 2.2, color = "white", fontface = "bold") +
    scale_fill_gradientn(
      colors  = c("#EBF8FF", "#90CDF4", "#3182CE", "#C53030"),
      na.value = "#F8F9FA",
      name    = "Vehicles",
      labels  = function(x) {
        v <- round(10^x)
        ifelse(v >= 1000, paste0(round(v/1000, 0), "K"), as.character(v))
      },
      breaks  = log10(c(10, 100, 1000, 10000, 100000)),
      guide   = guide_colorbar(barwidth = 1, barheight = 14,
                               ticks.colour = "#E2E8F0",
                               frame.colour = "#E2E8F0")
    ) +
    labs(
      title    = "Cumulative Origin–Destination Flow Matrix",
      subtitle = "Total vehicles exchanged between cell pairs · empty cells = no direct connection · labels on top 25% flows",
      x = "Destination Cell", y = "Origin Cell"
    ) +
    theme_traffic() +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1, size = 7.5, color = "#1A202C"),
      axis.text.y = element_text(size = 7.5, color = "#1A202C"),
      panel.grid  = element_blank()
    ) +
    coord_equal()
}

plot_top_cells_timeseries <- function(densityHistory, timeHistory, cellIDs,
                                      jamDensity, n_top = 8) {
  avgDensity <- colMeans(densityHistory)
  topIdx     <- order(avgDensity, decreasing = TRUE)[1:n_top]
  topIDs     <- cellIDs[topIdx]

  df <- prep_density_df(densityHistory[, topIdx, drop = FALSE],
                        timeHistory, topIDs) %>%
    mutate(pct_jam = density / jamDensity,
           cell    = paste0("Cell ", cell))

  # Add congestion state bands
  bands <- data.frame(
    ymin  = c(0,    0.25, 0.65),
    ymax  = c(0.25, 0.65, 1.0),
    state = c("Free Flow", "Transitional", "Congested")
  )

  ggplot() +
    geom_rect(data = bands,
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = state),
              alpha = 0.06, inherit.aes = FALSE) +
    geom_hline(yintercept = c(0.25, 0.65), color = TRAFFIC_GRID,
               linewidth = 0.4, linetype = "dashed") +
    geom_line(data = df, aes(x = time, y = pct_jam, color = cell),
              linewidth = 0.9, alpha = 0.9) +
    scale_fill_manual(values = c("Free Flow"    = COL_FREE,
                                 "Transitional" = COL_CRIT,
                                 "Congested"    = COL_JAM),
                      guide  = "none") +
    scale_color_manual(values = CELL_COLORS[1:n_top], name = NULL) +
    scale_x_continuous(labels = function(x) paste0(x, "h")) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    facet_wrap(~cell, ncol = 4) +
    labs(
      title    = "Congestion Trajectory — Top Busiest Cells",
      subtitle = "% of jam density over time · shaded bands = free/transitional/congested regimes",
      x = "Time", y = "% of Jam Density"
    ) +
    theme_traffic() +
    theme(legend.position = "none")
}

# Shows which cells are bottlenecks (low receiving) vs generators (high sending)
plot_capacity_scatter <- function(densityHistory, cellIDs, cellLengths,
                                  maxVehicles, criticalVehicles, jamDensity, params) {
  nCells <- length(cellIDs)
  avg_send <- numeric(nCells)
  avg_recv <- numeric(nCells)

  for (j in 1:nCells) {
    sends <- sapply(1:nrow(densityHistory), function(t) {
      v <- densityHistory[t, j] * cellLengths[j]
      compute_sending_flow(v, cellLengths[j], criticalVehicles[j], params)
    })
    recvs <- sapply(1:nrow(densityHistory), function(t) {
      v <- densityHistory[t, j] * cellLengths[j]
      compute_receiving_flow(v, cellLengths[j], maxVehicles[j], criticalVehicles[j], params)
    })
    avg_send[j] <- mean(sends)
    avg_recv[j] <- mean(recvs)
  }

  avg_density <- colMeans(densityHistory)

  df <- data.frame(
    cell        = as.character(cellIDs),
    avg_send    = avg_send,
    avg_recv    = avg_recv,
    avg_density = avg_density,
    pct_jam     = avg_density / jamDensity
  )

  # Zoom to data range with 5% padding
  x_pad <- diff(range(df$avg_send)) * 0.05
  y_pad <- diff(range(df$avg_recv)) * 0.05
  xlims <- range(df$avg_send) + c(-x_pad, x_pad) * 3
  ylims <- range(df$avg_recv) + c(-y_pad, y_pad) * 3

  # Place diagonal label inside data range
  mid <- mean(c(max(xlims[1], ylims[1]), min(xlims[2], ylims[2])))

  ggplot(df, aes(x = avg_send, y = avg_recv)) +
    geom_abline(slope = 1, intercept = 0, color = "#CBD5E0",
                linewidth = 0.7, linetype = "dashed") +
    annotate("text", x = mid * 0.97, y = mid * 1.03,
             label = "Send = Receive", color = "#A0AEC0",
             size = 3, angle = 40) +
    geom_point(aes(color = pct_jam, size = avg_density), alpha = 0.85) +
    geom_text(aes(label = cell), color = "#2D3748", size = 2.4,
              vjust = -1.2, fontface = "bold") +
    scale_color_gradientn(
      colors = c(COL_FREE, COL_CRIT, COL_JAM),
      name   = "% Jam\nDensity",
      labels = percent_format()
    ) +
    scale_size_continuous(range = c(2, 8), guide = "none") +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    labs(
      title    = "Sending vs. Receiving Capacity",
      subtitle = "Points below diagonal = bottlenecks · color = congestion · zoomed to data range",
      x = "Avg Sending Capacity (veh/hr)", y = "Avg Receiving Capacity (veh/hr)"
    ) +
    theme_traffic()
}

plot_congestion_timeline <- function(densityHistory, timeHistory, cellIDs, jamDensity) {
  df <- prep_density_df(densityHistory, timeHistory, cellIDs) %>%
    mutate(
      state = case_when(
        density / jamDensity < 0.25 ~ "Free Flow",
        density / jamDensity < 0.65 ~ "Transitional",
        TRUE                         ~ "Congested"
      ),
      state = factor(state, levels = c("Free Flow", "Transitional", "Congested")),
      cell  = factor(cell, levels = sort(unique(cell)))
    )

  # Summarise fraction in each state per timestep
  df_pct <- df %>%
    group_by(time, state) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(time) %>%
    mutate(pct = n / sum(n))

  ggplot(df_pct, aes(x = time, y = pct, fill = state)) +
    geom_area(alpha = 0.9, color = NA) +
    scale_fill_manual(values = c(
      "Free Flow"    = COL_FREE,
      "Transitional" = COL_CRIT,
      "Congested"    = COL_JAM
    ), name = NULL) +
    scale_x_continuous(labels = function(x) paste0(x, "h")) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    labs(
      title    = "Network State Composition Over Time",
      subtitle = "Share of cells in each congestion regime at each timestep",
      x = "Time", y = "Share of Network"
    ) +
    theme_traffic() +
    theme(legend.position = "top")
}

# Throughput efficiency: how much flow does each cell generate per unit of stored density
plot_cell_efficiency <- function(densityHistory, timeHistory, cellIDs,
                                 cellLengths, flowHistory) {

  nCells <- length(cellIDs)

  # Total outflow per cell across all timesteps
  total_outflow <- numeric(nCells)
  for (step in seq_along(flowHistory)) {
    total_outflow <- total_outflow + rowSums(flowHistory[[step]])
  }

  avg_density  <- colMeans(densityHistory)
  avg_vehicles <- avg_density * cellLengths
  efficiency   <- total_outflow / pmax(avg_vehicles * length(flowHistory), 1)

  df <- data.frame(
    cell       = factor(cellIDs),
    efficiency = efficiency,
    avg_density = avg_density,
    total_outflow = total_outflow
  ) %>% arrange(desc(efficiency)) %>%
    mutate(cell = factor(cell, levels = cell),
           fill_col = efficiency / max(efficiency))

  ggplot(df, aes(x = cell, y = efficiency, fill = fill_col)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.2f", efficiency)),
              vjust = -0.5, color = TRAFFIC_TEXT, size = 2.8, fontface = "bold") +
    scale_fill_gradientn(
      colors = c(COL_JAM, COL_CRIT, COL_FREE),
      guide  = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title    = "Cell Throughput Efficiency",
      subtitle = "Total outflow / (avg vehicles × timesteps) — higher = more efficient throughput",
      x = "Cell ID", y = "Efficiency Index"
    ) +
    theme_traffic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}

plot_fundamental_diagram <- function(rho_c = 26.5, q_max = 2200, rho_j = 135,
                                     av_penetration = 0.70,
                                     cap_gain = 0.12, rho_c_gain = 0.10, rho_j_gain = 0.08) {
  v_f <- q_max / rho_c
  w   <- -q_max / (rho_j - rho_c)

  rho_c_av <- rho_c * (1 + rho_c_gain)
  q_max_av <- q_max * (1 + cap_gain)
  rho_j_av <- rho_j * (1 + rho_j_gain)

  base_df <- data.frame(
    rho = c(0, rho_c, rho_j),
    q   = c(0, q_max, 0),
    curve = "Baseline"
  )
  av_lbl <- sprintf("%.0f%% AV penetration", av_penetration * 100)
  av_df <- data.frame(
    rho = c(0, rho_c_av, rho_j_av),
    q   = c(0, q_max_av, 0),
    curve = av_lbl
  )
  df <- rbind(base_df, av_df)
  df$curve <- factor(df$curve, levels = c("Baseline", av_lbl))

  ann <- data.frame(
    rho   = c(rho_c,  rho_j,  rho_c_av, rho_j_av),
    q     = c(q_max,  0,      q_max_av, 0),
    label = c(sprintf("rho_c = %.1f\nq_max = %d", rho_c, as.integer(q_max)),
              sprintf("rho_j = %d", as.integer(rho_j)),
              sprintf("rho_c* = %.1f\nq_max* = %d", rho_c_av, as.integer(q_max_av)),
              sprintf("rho_j* = %d", as.integer(rho_j_av))),
    hjust = c(-0.1, 0.5, -0.1, 0.5),
    vjust = c(1.1, -0.5,  1.1, -0.5)
  )

  slope_lbl <- data.frame(
    rho   = c(10, 100),
    q     = c(10 * v_f * 0.55, q_max * 0.28),
    label = c(sprintf("v_f = %.0f km/h", v_f), sprintf("w = %.0f km/h", w)),
    angle = c(atan2(q_max, rho_c) * 180 / pi, atan2(-q_max, rho_j - rho_c) * 180 / pi)
  )

  ggplot(df, aes(x = rho, y = q, color = curve, linetype = curve)) +
    geom_line(linewidth = 1.3) +
    geom_point(data = rbind(base_df, av_df)[rbind(base_df, av_df)$rho %in%
                 c(0, rho_c, rho_j, rho_c_av, rho_j_av), ],
               aes(x = rho, y = q), size = 2.8, show.legend = FALSE) +
    geom_vline(xintercept = rho_c, linetype = "dotted", color = "#3182CE", linewidth = 0.5) +
    geom_text(data = ann,
              aes(x = rho, y = q, label = label, hjust = hjust, vjust = vjust),
              size = 3, inherit.aes = FALSE, color = "#2D3748", lineheight = 0.9) +
    geom_text(data = slope_lbl,
              aes(x = rho, y = q, label = label, angle = angle),
              size = 3, color = "#718096", inherit.aes = FALSE) +
    scale_color_manual(values = setNames(c("#3182CE", "#D97706"), c("Baseline", av_lbl))) +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c("Baseline", av_lbl))) +
    scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 25),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(limits = c(0, 2600), breaks = seq(0, 2500, 500),
                       expand = expansion(mult = c(0, 0.04))) +
    labs(title    = "Triangular Fundamental Diagram",
         subtitle = sprintf(
           "Baseline: v_f = %.0f km/h, w = %.0f km/h, rho_c = %.1f veh/km, q_max = %d veh/h, rho_j = %d veh/km\nDashed = AV-modified curve at p = %.0f%% penetration (shifted apex)",
           v_f, w, rho_c, as.integer(q_max), as.integer(rho_j), av_penetration * 100),
         x = "Density rho (veh/km)", y = "Flow q (veh/h)",
         color = NULL, linetype = NULL) +
    theme_traffic() +
    theme(legend.position = "top")
}

# ─────────────────────────────────────────────────────────────────────────────
# Figure CTM-2 : Congestion Heatmap of Boston ZIP Network
# ─────────────────────────────────────────────────────────────────────────────
plot_congestion_heatmap <- function(densityHistory, timeHistory, cellIDs, jamDensity) {
  avg_density   <- colMeans(densityHistory)
  rho_frac      <- avg_density / jamDensity
  bottleneck_id <- cellIDs[which.max(rho_frac)]

  df <- data.frame(
    cell          = factor(as.character(cellIDs), levels = as.character(sort(cellIDs))),
    rho_frac      = rho_frac,
    is_bottleneck = cellIDs == bottleneck_id
  )

  ggplot(df, aes(x = cell, y = 1, fill = rho_frac)) +
    geom_tile(color = "#E2E8F0", linewidth = 0.5) +
    geom_text(aes(label = ifelse(is_bottleneck,
                                 paste0(cell, "\n(bottleneck)"),
                                 as.character(cell))),
              size      = 2.6,
              color     = ifelse(df$rho_frac > 0.6, "white", "#1A202C"),
              fontface  = ifelse(df$is_bottleneck, "bold", "plain"),
              lineheight = 0.85) +
    scale_fill_gradientn(
      colors = c("#FEFCE8", "#FDE68A", "#FB923C", "#DC2626", "#7F1D1D"),
      values = scales::rescale(c(0, 0.3, 0.55, 0.80, 1.0)),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00 (jam)"),
      name   = "rho / rho_j",
      guide  = guide_colorbar(barwidth = 1, barheight = 10,
                              ticks.colour = "#4A5568", frame.colour = "#4A5568")
    ) +
    scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    labs(
      title    = "Congestion Heatmap - Boston ZIP Network",
      subtitle = sprintf(
        "Mean density (rho/rho_j) per cell across simulation · rho_j = %d veh/km · bottleneck = cell %s (%.2f x jam density)",
        as.integer(jamDensity), bottleneck_id, max(rho_frac)),
      x = "Cell ID", y = NULL
    ) +
    theme_traffic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          panel.grid  = element_blank())
}

params <- list(
  freeFlowSpeed       = 83,
  congestionWaveSpeed = 19.2,
  jamDensity          = 135,
  maxFlow             = 2200,
  timeStep            = 3600
)

load_distance_matrix <- function(filepath, cellIDs) {
  data    <- fread(filepath)
  rowIDs  <- as.character(data[[1]])
  data    <- data[, -1, with = FALSE]
  colnames(data) <- as.character(colnames(data))
  data    <- as.data.frame(data)
  rownames(data) <- rowIDs
  existingCells <- as.character(cellIDs[cellIDs %in% colnames(data)])
  data    <- data[existingCells, existingCells, drop = FALSE]
  mat     <- as.matrix(data)
  mat[is.na(mat)] <- 0
  mat[mat < 0]    <- 0
  return(mat)
}

compute_sending_flow <- function(vehicles, cellLength, criticalVehicles, params) {
  sf <- if (vehicles <= criticalVehicles)
    params$freeFlowSpeed * vehicles / cellLength
  else
    params$maxFlow
  sf * (params$timeStep / 3600)
}

compute_receiving_flow <- function(vehicles, cellLength, maxVehicles, criticalVehicles, params) {
  rf <- if (vehicles < criticalVehicles)
    params$maxFlow
  else
    params$congestionWaveSpeed * (maxVehicles - vehicles) / cellLength
  rf * (params$timeStep / 3600)
}

compute_flows <- function(currentVehicles, cellLengths, maxVehicles,
                          criticalVehicles, distanceMatrix, params) {
  nCells     <- length(currentVehicles)
  flowMatrix <- matrix(0, nrow = nCells, ncol = nCells)
  sendCaps   <- sapply(1:nCells, function(i)
    compute_sending_flow(currentVehicles[i], cellLengths[i], criticalVehicles[i], params))
  recvCaps   <- sapply(1:nCells, function(j)
    compute_receiving_flow(currentVehicles[j], cellLengths[j], maxVehicles[j], criticalVehicles[j], params))
  for (i in 1:nCells)
    for (j in 1:nCells)
      if (distanceMatrix[i,j] > 0)
        flowMatrix[i,j] <- min(sendCaps[i], recvCaps[j])
  for (i in 1:nCells) {
    tot <- sum(flowMatrix[i,])
    if (tot > sendCaps[i] && tot > 0) flowMatrix[i,] <- flowMatrix[i,] * sendCaps[i] / tot
  }
  for (j in 1:nCells) {
    tot <- sum(flowMatrix[,j])
    if (tot > recvCaps[j] && tot > 0) flowMatrix[,j] <- flowMatrix[,j] * recvCaps[j] / tot
  }
  flowMatrix
}

add_demand <- function(currentVehicles, maxVehicles, cellIDs, trafficDemand, currentTime, timeStep) {
  updatedVehicles <- currentVehicles
  if (!("time" %in% colnames(trafficDemand))) return(updatedVehicles)
  timeEnd        <- currentTime + timeStep / 3600
  relevantDemand <- trafficDemand[trafficDemand$time >= currentTime & trafficDemand$time < timeEnd,]
  if (nrow(relevantDemand) == 0) return(updatedVehicles)
  for (k in 1:nrow(relevantDemand)) {
    idx <- which(cellIDs == relevantDemand$cellID[k])
    if (length(idx) == 1) {
      space <- maxVehicles[idx] - updatedVehicles[idx]
      updatedVehicles[idx] <- updatedVehicles[idx] + min(relevantDemand$vehicles[k], space)
    }
  }
  updatedVehicles
}

update_state <- function(currentVehicles, flowMatrix, maxVehicles) {
  nCells          <- length(currentVehicles)
  updatedVehicles <- currentVehicles
  for (i in 1:nCells)
    updatedVehicles[i] <- pmax(0, pmin(
      currentVehicles[i] - sum(flowMatrix[i,]) + sum(flowMatrix[,i]),
      maxVehicles[i]))
  updatedVehicles
}

cellIDs <- c(2108, 2109, 2110, 2111, 2113, 2114, 2115, 2116, 2118, 2119,
             2120, 2121, 2122, 2124, 2125, 2126, 2127, 2128, 2129, 2130,
             2131, 2132, 2134, 2135, 2136, 2203, 2210, 2215)

distanceMatrix <- load_distance_matrix("distance_matrix.csv", cellIDs) * 1.60934
trafficDemand  <- fread("traffic_demand.csv")

duration  <- 600
nCells    <- length(cellIDs)
nSteps    <- floor(duration * 3600 / params$timeStep)

cellLengths      <- rep(mean(distanceMatrix[distanceMatrix > 0]), nCells)
maxVehicles      <- cellLengths * params$jamDensity
criticalVehicles <- cellLengths * params$maxFlow / params$freeFlowSpeed

currentVehicles <- 0.2 * params$jamDensity * cellLengths

densityHistory       <- matrix(0, nrow = nSteps, ncol = nCells)
timeHistory          <- numeric(nSteps)
flowHistory          <- vector("list", nSteps)
cumulativeFlowMatrix <- matrix(0, nrow = nCells, ncol = nCells,
                               dimnames = list(cellIDs, cellIDs))

for (step in 1:nSteps) {
  currentTime     <- (step - 1) * params$timeStep / 3600
  currentVehicles <- add_demand(currentVehicles, maxVehicles, cellIDs,
                                trafficDemand, currentTime, params$timeStep)
  flowMatrix      <- compute_flows(currentVehicles, cellLengths, maxVehicles,
                                   criticalVehicles, distanceMatrix, params)
  cumulativeFlowMatrix <- cumulativeFlowMatrix + flowMatrix
  currentVehicles <- update_state(currentVehicles, flowMatrix, maxVehicles)
  densityHistory[step,] <- currentVehicles / cellLengths
  flowHistory[[step]]   <- flowMatrix
  timeHistory[step]     <- currentTime
  if (step %% max(1, floor(nSteps/10)) == 0)
    cat(sprintf("Progress: %.1f%%\n", 100 * step / nSteps))
}
cat("Simulation complete!\n")

dir.create("outputs", showWarnings = FALSE)
write.csv(cumulativeFlowMatrix, "outputs/cumulative_flow_matrix.csv", row.names = TRUE)

pdf("outputs/ctm_analysis.pdf", width = 16, height = 9,
    bg = TRAFFIC_DARK)

print(plot_spacetime_heatmap(densityHistory, timeHistory, cellIDs, params$jamDensity))
print(plot_network_pulse(densityHistory, timeHistory, cellIDs, cellLengths, flowHistory))
print(plot_congestion_timeline(densityHistory, timeHistory, cellIDs, params$jamDensity))
print(plot_top_cells_timeseries(densityHistory, timeHistory, cellIDs, params$jamDensity))
print(plot_density_ridges(densityHistory, timeHistory, cellIDs, params$jamDensity))
print(plot_capacity_scatter(densityHistory, cellIDs, cellLengths, maxVehicles,
                            criticalVehicles, params$jamDensity, params))
print(plot_cell_efficiency(densityHistory, timeHistory, cellIDs, cellLengths, flowHistory))
print(plot_flow_matrix_styled(cumulativeFlowMatrix, cellIDs))
print(plot_fundamental_diagram())
print(plot_congestion_heatmap(densityHistory, timeHistory, cellIDs, params$jamDensity))

dev.off()
cat("Plots saved to outputs/ctm_analysis.pdf\n")