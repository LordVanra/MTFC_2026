library(data.table)
library(ggplot2)
library(reshape2)

params = list(
  freeFlowSpeed = 83,
  congestionWaveSpeed = 19.2,
  jamDensity = 135,
  maxFlow = 2200,
  timeStep = 3600
)

load_distance_matrix <- function(filepath, cellIDs) {
  data <- fread(filepath)
  
  rowIDs <- as.character(data[[1]])
  data <- data[, -1, with = FALSE]
  colnames(data) <- as.character(colnames(data))
  
  data <- as.data.frame(data)
  rownames(data) <- rowIDs
  
  existingCells <- as.character(cellIDs[cellIDs %in% colnames(data)])
  
  data <- data[existingCells, existingCells, drop = FALSE]
  
  mat <- as.matrix(data)
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  
  return(mat)
}

compute_sending_flow <- function(vehicles, cellLength, criticalVehicles, params) {
  if (vehicles <= criticalVehicles) {
    sendingFlow <- params$freeFlowSpeed * vehicles / cellLength
  } else {
    sendingFlow <- params$maxFlow
  }
  sendingFlow <- sendingFlow * (params$timeStep / 3600)
  return(sendingFlow)
}

compute_receiving_flow <- function(vehicles, cellLength, maxVehicles, criticalVehicles, params) {
  if (vehicles < criticalVehicles) {
    receivingFlow <- params$maxFlow
  } 
  else {
    receivingFlow <- params$congestionWaveSpeed * (maxVehicles - vehicles) / cellLength
  }
  receivingFlow <- receivingFlow * (params$timeStep / 3600)
  return(receivingFlow)
}

compute_flows <- function(currentVehicles, cellLengths, maxVehicles, 
                          criticalVehicles, distanceMatrix, params) {
  
  nCells <- length(currentVehicles)
  flowMatrix <- matrix(0, nrow = nCells, ncol = nCells)
  
  sendCaps <- numeric(nCells)
  recvCaps <- numeric(nCells)
  
  for (i in 1:nCells)
    sendCaps[i] <- compute_sending_flow(
      currentVehicles[i], cellLengths[i], criticalVehicles[i], params)
  
  for (j in 1:nCells)
    recvCaps[j] <- compute_receiving_flow(
      currentVehicles[j], cellLengths[j], 
      maxVehicles[j], criticalVehicles[j], params)
  
  for (i in 1:nCells) {
    for (j in 1:nCells) {
      if (distanceMatrix[i,j] > 0) {
        flowMatrix[i,j] <- min(sendCaps[i], recvCaps[j])
      }
    }
  }
  
  for (i in 1:nCells) {
    totalOutflow <- sum(flowMatrix[i, ])
    if (totalOutflow > sendCaps[i] && totalOutflow > 0) {
      flowMatrix[i, ] <- flowMatrix[i, ] * sendCaps[i] / totalOutflow
    }
  }
  
  for (j in 1:nCells) {
    totalInflow <- sum(flowMatrix[, j])
    if (totalInflow > recvCaps[j] && totalInflow > 0) {
      flowMatrix[, j] <- flowMatrix[, j] * recvCaps[j] / totalInflow
    }
  }
  
  return(flowMatrix)
}

add_demand <- function(currentVehicles, maxVehicles, cellIDs, trafficDemand, currentTime, timeStep) {
  
  updatedVehicles <- currentVehicles
  nCells <- length(currentVehicles)
  
  if (!("time" %in% colnames(trafficDemand)))
    return(updatedVehicles)
  
  timeEnd <- currentTime + timeStep / 3600
  
  relevantDemand <- trafficDemand[
    trafficDemand$time >= currentTime &
    trafficDemand$time < timeEnd, ]
  
  if (nrow(relevantDemand) == 0)
    return(updatedVehicles)
  
  for (k in 1:nrow(relevantDemand)) {
    
    originID <- relevantDemand$cellID[k]
    demand   <- relevantDemand$vehicles[k]
    
    cellIndex <- which(cellIDs == originID)
    
    if (length(cellIndex) == 1) {
      
      availableSpace <- maxVehicles[cellIndex] - updatedVehicles[cellIndex]
      
      vehiclesToAdd <- min(demand, availableSpace)
      
      updatedVehicles[cellIndex] <- 
        updatedVehicles[cellIndex] + vehiclesToAdd
    }
  }
  
  return(updatedVehicles)
}

update_state <- function(currentVehicles, flowMatrix, maxVehicles) {
  nCells <- length(currentVehicles)
  updatedVehicles <- currentVehicles
  
  for (i in 1:nCells) {
    outflow <- sum(flowMatrix[i, ])
    inflow <- sum(flowMatrix[, i])
    updatedVehicles[i] <- currentVehicles[i] - outflow + inflow
    updatedVehicles[i] <- pmax(0, pmin(updatedVehicles[i], maxVehicles[i]))
  }
  
  return(updatedVehicles)
}

save_plots <- function(densityHistory, timeHistory, cellIDs, cellLengths, output_dir = "outputs") {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE)
  
  # Prepare data
  df <- as.data.frame(densityHistory)
  colnames(df) <- paste0("Cell_", cellIDs)
  df$Time <- timeHistory
  df_melt <- melt(df, id.vars = "Time", variable.name = "Cell", value.name = "Density")
  
  # Plot 1: Density over time (simplified - show only subset of cells)
  # Select top 10 cells by average density for clarity
  avgDensity <- colMeans(densityHistory)
  topCells <- order(avgDensity, decreasing = TRUE)[1:min(10, length(cellIDs))]
  topCellNames <- paste0("Cell_", cellIDs[topCells])
  
  df_melt_subset <- df_melt[df_melt$Cell %in% topCellNames, ]
  
  png(file.path(output_dir, "density_over_time.png"), width = 1200, height = 600)
  p1 <- ggplot(df_melt_subset, aes(x = Time, y = Density, color = Cell)) +
    geom_line(linewidth = 0.8) + 
    theme_minimal() +
    labs(title = "Traffic Density over Time (Top 10 Busiest Cells)", 
         x = "Time (hours)", 
         y = "Density (vehicles/km)") +
    theme(legend.position = "right")
  print(p1)
  dev.off()
  
  # Plot 2: Total vehicles
  totalVehicles <- rowSums(densityHistory * 
                         matrix(cellLengths, 
                                nrow = nrow(densityHistory), 
                                ncol = ncol(densityHistory), 
                                byrow = TRUE))
  
  df_total <- data.frame(Time = timeHistory, TotalVehicles = totalVehicles)
  
  png(file.path(output_dir, "total_network_occupancy.png"), width = 1200, height = 600)
  p2 <- ggplot(df_total, aes(x = Time, y = TotalVehicles)) +
    geom_line(color = "blue", linewidth = 1) + 
    theme_minimal() +
    labs(title = "Total Network Occupancy", x = "Time (hours)", y = "Total Vehicles")
  print(p2)
  dev.off()
  
  # Plot 3: Improved Heatmap - aggregate time into bins
  nTimeBins <- 50  # Aggregate into 50 time bins for readability
  timeBinSize <- ceiling(nrow(densityHistory) / nTimeBins)
  
  # Aggregate density data
  df_heat_agg <- data.frame()
  for (i in seq(1, nrow(densityHistory), by = timeBinSize)) {
    endIdx <- min(i + timeBinSize - 1, nrow(densityHistory))
    avgDensity <- colMeans(densityHistory[i:endIdx, , drop = FALSE])
    avgTime <- mean(timeHistory[i:endIdx])
    
    df_temp <- data.frame(
      Time = avgTime,
      Cell = cellIDs,
      Density = avgDensity
    )
    df_heat_agg <- rbind(df_heat_agg, df_temp)
  }
  
  png(file.path(output_dir, "density_heatmap.png"), width = 1400, height = 800)
  p3 <- ggplot(df_heat_agg, aes(x = Time, y = factor(Cell), fill = Density)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma") +
    labs(title = "Density Heatmap (Time-Aggregated)", 
         x = "Time (hours)", 
         y = "Cell ID", 
         fill = "Density\n(veh/km)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  print(p3)
  dev.off()
  
  cat(sprintf("Plots saved to %s/\n", output_dir))
}

plot_flow_matrix_heatmap <- function(flowMatrix, cellIDs, output_path = "outputs/flow_matrix_heatmap.png") {
  # Convert matrix to data frame for ggplot
  df_flow <- as.data.frame(flowMatrix)
  colnames(df_flow) <- cellIDs
  df_flow$Origin <- cellIDs
  
  # Melt to long format
  df_flow_melt <- melt(df_flow, id.vars = "Origin", variable.name = "Destination", value.name = "Vehicles")
  df_flow_melt$Destination <- as.numeric(as.character(df_flow_melt$Destination))
  
  # Add log-scaled values for color mapping (add 1 to avoid log(0))
  df_flow_melt$VehiclesLog <- log10(df_flow_melt$Vehicles + 1)
  
  # Create the heatmap
  png(output_path, width = 1600, height = 1400, res = 120)
  
  p <- ggplot(df_flow_melt, aes(x = factor(Destination), y = factor(Origin), fill = VehiclesLog)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_viridis_c(
      option = "plasma",
      name = "Vehicles\n(log scale)",
      labels = function(x) {
        vals <- 10^x - 1
        ifelse(vals < 1000, 
               sprintf("%.0f", vals),
               sprintf("%.1fK", vals/1000))
      },
      breaks = log10(c(1, 10, 100, 1000, 10000, 100000, 1000000) + 1)
    ) +
    labs(
      title = "Cumulative Flow Matrix - Origin to Destination",
      subtitle = "Log-scaled color intensity showing vehicle flows between cells",
      x = "Destination Cell",
      y = "Origin Cell"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    coord_equal()
  
  print(p)
  dev.off()
  
  cat(sprintf("Flow matrix heatmap saved to %s\n", output_path))
  
  # Print summary statistics
  cat("\n=== FLOW MATRIX SUMMARY ===\n")
  nonZeroFlows <- df_flow_melt$Vehicles[df_flow_melt$Vehicles > 0]
  cat(sprintf("Total vehicles in network: %s\n", formatC(sum(df_flow_melt$Vehicles), format="f", big.mark=",", digits=0)))
  cat(sprintf("Number of active routes: %d / %d\n", length(nonZeroFlows), nrow(df_flow_melt)))
  cat(sprintf("Max flow on single route: %s vehicles\n", formatC(max(df_flow_melt$Vehicles), format="f", big.mark=",", digits=0)))
  cat(sprintf("Average flow (active routes): %s vehicles\n", formatC(mean(nonZeroFlows), format="f", big.mark=",", digits=0)))
  cat(sprintf("Median flow (active routes): %s vehicles\n", formatC(median(nonZeroFlows), format="f", big.mark=",", digits=0)))
  
  # Top 10 routes
  cat("\n=== TOP 10 BUSIEST ROUTES ===\n")
  df_top <- df_flow_melt[order(df_flow_melt$Vehicles, decreasing = TRUE), ]
  df_top <- df_top[df_top$Vehicles > 0, ]
  for (i in 1:min(10, nrow(df_top))) {
    cat(sprintf("%2d. Cell %d → %d: %s vehicles (%.1f veh/hr)\n", 
                i, 
                df_top$Origin[i], 
                df_top$Destination[i],
                formatC(df_top$Vehicles[i], format="f", big.mark=",", digits=0),
                df_top$Vehicles[i] / 600))
  }
  
  return(invisible(p))
}

# Function to load and plot flow matrix from CSV
plot_flow_matrix_from_csv <- function(csv_path = "outputs/cumulative_flow_matrix.csv", 
                                      output_path = "outputs/flow_matrix_heatmap.png") {
  # Load the CSV
  flowData <- read.csv(csv_path, row.names = 1, check.names = FALSE)
  flowMatrix <- as.matrix(flowData)
  
  # Get cell IDs from row/column names
  cellIDs <- as.numeric(rownames(flowMatrix))
  
  # Create the heatmap
  plot_flow_matrix_heatmap(flowMatrix, cellIDs, output_path)
}

cellIDs <- c(2108, 2109, 2110, 2111, 2113, 2114, 2115, 2116, 2118, 2119,
             2120, 2121, 2122, 2124, 2125, 2126, 2127, 2128, 2129, 2130,
             2131, 2132, 2134, 2135, 2136, 2203, 2210, 2215)

distanceMatrix <- load_distance_matrix("distance_matrix.csv", cellIDs) * 1.60934
trafficDemand <- fread("traffic_demand.csv")

duration <- 600
nCells <- length(cellIDs)
nSteps <- floor(duration * 3600 / params$timeStep)

cellLengths <- rep(mean(distanceMatrix[distanceMatrix > 0]), nCells)
maxVehicles <- cellLengths * params$jamDensity
criticalVehicles <- cellLengths * params$maxFlow / params$freeFlowSpeed

initialDensity <- 0.2 * params$jamDensity
currentVehicles <- initialDensity * cellLengths
currentDensity <- initialDensity

densityHistory <- matrix(0, nrow = nSteps, ncol = nCells)
timeHistory <- numeric(nSteps)
flowHistory <- vector("list", nSteps)

# Initialize cumulative flow matrix
cumulativeFlowMatrix <- matrix(0, nrow = nCells, ncol = nCells)
rownames(cumulativeFlowMatrix) <- cellIDs
colnames(cumulativeFlowMatrix) <- cellIDs

# Simulation loop
for (step in 1:nSteps) {
  currentTime <- (step - 1) * params$timeStep / 3600
  currentVehicles <- add_demand(currentVehicles, maxVehicles, cellIDs, trafficDemand, currentTime, params$timeStep)
  flowMatrix <- compute_flows(currentVehicles, cellLengths, maxVehicles, criticalVehicles, distanceMatrix, params)
  
  # Accumulate flows
  cumulativeFlowMatrix <- cumulativeFlowMatrix + flowMatrix
  
  currentVehicles <- update_state(currentVehicles, flowMatrix, maxVehicles)
  currentDensity <- currentVehicles / cellLengths
  
  densityHistory[step, ] <- currentDensity
  flowHistory[[step]] <- flowMatrix
  timeHistory[step] <- currentTime
  
  if (step %% max(1, floor(nSteps/10)) == 0) {
    cat(sprintf("Progress: %.1f%%\n", 100 * step / nSteps))
  }
}

cat("Simulation complete!\n")

# Save all plots to disk
save_plots(densityHistory, timeHistory, cellIDs, cellLengths)

# Save cumulative flow matrix
dir.create("outputs", showWarnings = FALSE)
write.csv(cumulativeFlowMatrix, "outputs/cumulative_flow_matrix.csv", row.names = TRUE)
cat("Cumulative flow matrix saved to outputs/cumulative_flow_matrix.csv\n")

# Create flow matrix heatmap
plot_flow_matrix_heatmap(cumulativeFlowMatrix, cellIDs)

# Summary statistics
totalVehicles <- densityHistory %*% diag(cellLengths)
maxVehicles_total <- max(rowSums(totalVehicles))
peakTime <- timeHistory[which.max(rowSums(totalVehicles))]

cat(sprintf("\nPeak network occupancy: %.0f vehicles at time %.2f hours\n", maxVehicles_total, peakTime))
cat("\nAverage density by cell:\n")
avgDensities <- colMeans(densityHistory)
for (i in 1:nCells) {
  cat(sprintf("Cell %d: %.2f vehicles/km\n", cellIDs[i], avgDensities[i]))
}