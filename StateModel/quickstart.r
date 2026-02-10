################################################################################
# QUICK START EXAMPLE - AV ADOPTION MODEL
#
# Minimal working example to test the model
# Run this first to verify everything works
################################################################################

# --- 1. Load the main model ---
source("av_adoption_model.R")

# --- 2. Run a simple 10-year simulation ---

cat("\n=== QUICK START: 10-YEAR PROJECTION ===\n\n")

# Define initial state
state_0 <- c(
  A = 5e6,       # 5 million AVs today
  C = 280e6,     # 280 million total vehicles
  I = 0.25       # 25% infrastructure readiness
)

# Define policy
policy <- c(
  r = 7500,      # $7,500 tax credit
  s = 3000,      # $3,000 manufacturer subsidy
  i = 0.08       # 8% infrastructure investment
)

# Run simulation
results <- simulate_model(state_0, policy, years = 10, params)

# Print results
cat("Year | AV Share | Adoption Rate | Infrastructure\n")
cat("-----|----------|---------------|---------------\n")
for (i in seq(1, 11, by = 2)) {
  cat(sprintf("%4d | %7.2f%% | %12.2f%% | %13.1f%%\n",
              results$year[i],
              results$pct_AV[i] * 100,
              results$adoption[i] * 100,
              results$I[i] * 100))
}

# Plot
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(results$year, results$pct_AV * 100,
     type = "l", lwd = 3, col = "#2E86AB",
     xlab = "Year", ylab = "AV Market Share (%)",
     main = "Projected AV Adoption")
grid(col = "gray80", lty = 2)

plot(results$year, results$adoption * 100,
     type = "l", lwd = 3, col = "#F18F01",
     xlab = "Year", ylab = "Annual Adoption Rate (%)",
     main = "New AV Purchases")
grid(col = "gray80", lty = 2)

par(mfrow = c(1, 1))

cat("\n✓ Quick start complete!\n")
cat("  Now try running the full av_adoption_model.R script\n")