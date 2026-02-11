library(data.table)
library(ggplot2)
library(reshape2)

# -------------------------
# Load distance matrix
# -------------------------
load_distance_matrix <- function(filepath, cellIDs) {
  data <- fread(filepath)
  nCells <- length(cellIDs)
  distanceMatrix <- matrix(0, nrow = nCells, ncol = nCells)
  
  for (i in 1:nCells) {
    rowIdx <- which(data[[1]] == cellIDs[i])
    if (length(rowIdx) > 0) {
      for (j in 1:nCells) {
        if ((j+1) <= ncol(data)) {
          distanceMatrix[i,j] <- data[rowIdx, j+1, with = FALSE][[1]]
        }
      }
    }
  }
  return(distanceMatrix)
}

# -------------------------
# Calculate cell lengths
# -------------------------
calculate_cell_lengths <- function(distanceMatrix) {
  nCells <- nrow(distanceMatrix)
  cellLengths <- numeric(nCells)
  for (i in 1:nCells) {
    nonZeroDist <- distanceMatrix[i, distanceMatrix[i, ] > 0]
    if (length(nonZeroDist) > 0) {
      cellLengths[i] <- mean(nonZeroDist)
    } else {
      cellLengths[i] <- 1.0
    }
  }
  return(cellLengths)
}

# -------------------------
# Sending and receiving flows
# -------------------------
compute_sending_flow <- function(vehicles, cellLength, criticalVehicles, params) {
  if (vehicles <= criticalVehicles) {
    sendingFlow <- params$freeFlowSpeed * vehicles / cellLength
  } else {
    sendingFlow <- params$maxFlow
  }
  sendingFlow <- sendingFlow * (params$timeStep / 3600)  # vehicles per time step
  return(sendingFlow)
}

compute_receiving_flow <- function(vehicles, cellLength, maxVehicles, criticalVehicles, params) {
  if (vehicles < criticalVehicles) {
    receivingFlow <- params$maxFlow
  } else {
    receivingFlow <- params$congestionWaveSpeed * (maxVehicles - vehicles) / cellLength
  }
  receivingFlow <- receivingFlow * (params$timeStep / 3600)
  return(receivingFlow)
}

# -------------------------
# Compute flows
# -------------------------
compute_flows <- function(currentVehicles, cellLengths, maxVehicles, criticalVehicles, distanceMatrix, params) {
  nCells <- length(currentVehicles)
  flowMatrix <- matrix(0, nrow = nCells, ncol = nCells)
  
  for (i in 1:nCells) {
    for (j in 1:nCells) {
      if (distanceMatrix[i,j] > 0) {
        sendCap <- compute_sending_flow(currentVehicles[i], cellLengths[i], criticalVehicles[i], params)
        recvCap <- compute_receiving_flow(currentVehicles[j], cellLengths[j], maxVehicles[j], criticalVehicles[j], params)
        flowMatrix[i,j] <- min(sendCap, recvCap)
      }
    }
  }
  
  # Normalize if total outflow exceeds sending capacity
  for (i in 1:nCells) {
    totalOutflow <- sum(flowMatrix[i, ])
    sendCap <- compute_sending_flow(currentVehicles[i], cellLengths[i], criticalVehicles[i], params)
    if (totalOutflow > sendCap) {
      flowMatrix[i, ] <- flowMatrix[i, ] * sendCap / totalOutflow
    }
  }
  
  return(flowMatrix)
}

# -------------------------
# Add traffic demand
# -------------------------
add_demand <- function(currentVehicles, maxVehicles, cellIDs, trafficDemand, currentTime, timeStep) {
  updatedVehicles <- currentVehicles
  timeTolerance <- timeStep / 3600  # hours
  
  if (!("time" %in% colnames(trafficDemand))) return(updatedVehicles)
  
  relevantDemand <- trafficDemand[trafficDemand$time >= currentTime & trafficDemand$time < currentTime + timeTolerance, ]
  
  for (i in 1:nrow(relevantDemand)) {
    originCell <- relevantDemand$origin_cell[i]
    demand <- relevantDemand$demand[i]
    cellIdx <- which(cellIDs == originCell)
    if (length(cellIdx) > 0) {
      vehiclesToAdd <- demand * (timeStep / 3600)
      spaceAvailable <- maxVehicles[cellIdx] - updatedVehicles[cellIdx]
      vehiclesAdded <- min(vehiclesToAdd, spaceAvailable)
      updatedVehicles[cellIdx] <- updatedVehicles[cellIdx] + vehiclesAdded
    }
  }
  
  return(updatedVehicles)
}

# -------------------------
# Update vehicle state
# -------------------------
update_state <- function(currentVehicles, flowMatrix, maxVehicles) {
  nCells <- length(currentVehicles)
  updatedVehicles <- currentVehicles
  
  for (i in 1:nCells) {
    outflow <- sum(flowMatrix[i, ])
    inflow <- sum(flowMatrix[, i])
    updatedVehicles[i] <- currentVehicles[i] - outflow + inflow
    updatedVehicles[i] <- max(0, min(updatedVehicles[i], maxVehicles[i]))
  }
  
  return(updatedVehicles)
}

# -------------------------
# Plot results
# -------------------------
plot_results <- function(densityHistory, timeHistory, cellIDs, cellLengths) {
  # Convert to data frame for ggplot
  df <- as.data.frame(densityHistory)
  colnames(df) <- paste0("Cell_", cellIDs)
  df$Time <- timeHistory
  df_melt <- melt(df, id.vars = "Time", variable.name = "Cell", value.name = "Density")
  
  # Density over time
  p1 <- ggplot(df_melt, aes(x = Time, y = Density, color = Cell)) +
    geom_line() + theme_minimal() +
    labs(title = "Traffic Density over Time", x = "Time (hours)", y = "Density (vehicles/km)")
  print(p1)
  
  # Total vehicles
  totalVehicles <- densityHistory %*% diag(cellLengths)
  df_total <- data.frame(Time = timeHistory, TotalVehicles = rowSums(totalVehicles))
  
  p2 <- ggplot(df_total, aes(x = Time, y = TotalVehicles)) +
    geom_line(color = "blue") + theme_minimal() +
    labs(title = "Total Network Occupancy", x = "Time (hours)", y = "Total Vehicles")
  print(p2)
  
  # Heatmap
  df_heat <- melt(as.data.frame(densityHistory))
  df_heat$Time <- rep(timeHistory, ncol(densityHistory))
  df_heat$Cell <- rep(cellIDs, each = length(timeHistory))
  
  p3 <- ggplot(df_heat, aes(x = Time, y = factor(Cell), fill = value)) +
    geom_tile() + scale_fill_viridis_c() +
    labs(title = "Density Heatmap", x = "Time (hours)", y = "Cell", fill = "Density") +
    theme_minimal()
  print(p3)
}

# -------------------------
# Main simulation
# -------------------------
cellIDs <- c(2108, 2109, 2110, 2111, 2113, 2114, 2115, 2116, 2118, 2119,
             2120, 2121, 2122, 2124, 2125, 2126, 2127, 2128, 2129, 2130,
             2131, 2132, 2134, 2135, 2136, 2203, 2210, 2215)

distanceMatrix <- load_distance_matrix("distance_matrix.csv", cellIDs) * 1.60934
trafficDemand <- fread("traffic_demand.csv")

# Parameters
params <- list()
params$freeFlowSpeed <- 83
params$congestionWaveSpeed <- 19.2
params$jamDensity <- 135
params$maxFlow <- 2200
params$timeStep <- 86400
params$criticalDensity <- params$maxFlow / params$freeFlowSpeed

duration <- 600
nCells <- length(cellIDs)
nSteps <- floor(duration * 3600 / params$timeStep)

cellLengths <- calculate_cell_lengths(distanceMatrix)
maxVehicles <- cellLengths * params$jamDensity
criticalVehicles <- cellLengths * params$criticalDensity

initialDensity <- 0.2 * params$jamDensity
currentVehicles <- initialDensity * cellLengths
currentDensity <- initialDensity

densityHistory <- matrix(0, nrow = nSteps, ncol = nCells)
timeHistory <- numeric(nSteps)
flowHistory <- vector("list", nSteps)

# Simulation loop
for (step in 1:nSteps) {
  currentTime <- (step - 1) * params$timeStep / 3600
  currentVehicles <- add_demand(currentVehicles, maxVehicles, cellIDs, trafficDemand, currentTime, params$timeStep)
  flowMatrix <- compute_flows(currentVehicles, cellLengths, maxVehicles, criticalVehicles, distanceMatrix, params)
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

plot_results(densityHistory, timeHistory, cellIDs, cellLengths)

# Summary statistics
totalVehicles <- densityHistory %*% diag(cellLengths)
maxVehicles_total <- max(rowSums(totalVehicles))
peakTime <- timeHistory[which.max(rowSums(totalVehicles))]

cat(sprintf("Peak network occupancy: %.0f vehicles at time %.2f hours\n", maxVehicles_total, peakTime))
avgDensities <- colMeans(densityHistory)
for (i in 1:nCells) {
  cat(sprintf("Cell %d average density: %.2f vehicles/km\n", cellIDs[i], avgDensities[i]))
}
