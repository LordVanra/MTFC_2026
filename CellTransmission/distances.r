library(ggplot2)
library(reshape2)

plot_distance_matrix_heatmap <- function(distMatrix, cellIDs, output_path = "distance_matrix_heatmap.png") {
  # Convert matrix to data frame for ggplot
  df_dist <- as.data.frame(distMatrix)
  colnames(df_dist) <- cellIDs
  df_dist$Origin <- cellIDs
  
  # Melt to long format
  df_dist_melt <- melt(df_dist, id.vars = "Origin", variable.name = "Destination", value.name = "Distance")
  df_dist_melt$Destination <- as.numeric(as.character(df_dist_melt$Destination))
  
  # Create categories for coloring
  df_dist_melt$Category <- ifelse(df_dist_melt$Distance == 0, "Same Cell",
                                   ifelse(df_dist_melt$Distance < 0, "No Connection",
                                          "Connected"))
  
  # Create the heatmap
  png(output_path, width = 1600, height = 1400, res = 120)
  
  p <- ggplot(df_dist_melt, aes(x = factor(Destination), y = factor(Origin), fill = Distance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradientn(
      colors = c("#D32F2F", "#757575", "#FFF59D", "#4CAF50", "#2196F3"),
      values = scales::rescale(c(-1, -0.1, 0, 2, max(df_dist_melt$Distance))),
      name = "Distance\n(miles)",
      breaks = c(-1, 0, 1, 2, 3, 4, 5),
      labels = c("No Connection", "Same Cell", "1", "2", "3", "4", "5+"),
      na.value = "gray90"
    ) +
    labs(
      title = "Network Distance Matrix - Cell Connectivity",
      subtitle = "Red = No connection (-1), Yellow = Same cell (0), Green/Blue = Connected (distance in miles)",
      x = "Destination Cell",
      y = "Origin Cell"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    coord_equal()
  
  print(p)
  dev.off()
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)

# Default paths
input_csv <- ifelse(length(args) >= 1, args[1], "distance_matrix.csv")
output_png <- ifelse(length(args) >= 2, args[2], "outputs/distance_matrix_heatmap.png")

distData <- read.csv(input_csv, row.names = 1, check.names = FALSE)
distMatrix <- as.matrix(distData)

cellIDs <- as.numeric(rownames(distMatrix))

# Create output directory if needed
output_dir <- dirname(output_png)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Create the heatmap
plot_distance_matrix_heatmap(distMatrix, cellIDs, output_png)
