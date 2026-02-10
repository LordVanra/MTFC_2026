# TS should be deleted eventually

################################################################################
# GUIDE: PREPARING THE MODEL FOR PUBLICATION
#
# Steps to turn this into a research-grade model for journal submission
################################################################################

# ==============================================================================
# 1. DATA CALIBRATION
# ==============================================================================

# Replace dummy parameters with real data:

# STEP 1: Collect historical data
# - EV/AV sales data (IHS Markit, NHTSA, automaker reports)
# - Prices and incentives (Kelley Blue Book, IRS tax credit data)
# - Infrastructure deployment (DOE Alternative Fuels Data Center)
# - Fleet composition (Bureau of Transportation Statistics)

# STEP 2: Estimate parameters using Maximum Likelihood
estimate_from_data <- function(historical_data) {
  # Example structure:
  # historical_data should have columns:
  #   year, av_sales, total_sales, avg_rebate, avg_subsidy, infrastructure
  
  # Estimate beta coefficients
  adoption_rate <- historical_data$av_sales / historical_data$total_sales
  price_advantage <- historical_data$avg_rebate / params$P_AV
  
  # Logistic regression
  logit_model <- glm(
    adoption_rate ~ price_advantage + infrastructure,
    data = historical_data,
    family = binomial(link = "logit"),
    weights = total_sales  # weight by market size
  )
  
  # Estimate decay rates using time-series
  # delta: vehicle survival models
  # rho: infrastructure maintenance literature
  
  # Return calibrated params
  list(
    beta_0 = coef(logit_model)[1],
    beta_r = coef(logit_model)[2],
    beta_i = coef(logit_model)[3],
    # ... estimate others from data
    fit = summary(logit_model)
  )
}

# STEP 3: Validate out-of-sample
# - Hold out last 2-3 years
# - Fit on training data
# - Check RMSE, MAE on test set
# - Report R² and coefficient standard errors


# ==============================================================================
# 2. MODEL VALIDATION
# ==============================================================================

# Goodness-of-fit tests to include in paper:

validation_suite <- function(model_results, actual_data) {
  # 1. In-sample fit
  rmse <- sqrt(mean((model_results$A - actual_data$A)^2))
  mae <- mean(abs(model_results$A - actual_data$A))
  mape <- mean(abs((model_results$A - actual_data$A) / actual_data$A)) * 100
  
  # 2. Residual diagnostics
  residuals <- model_results$A - actual_data$A
  
  par(mfrow = c(2, 2))
  
  # Time series of residuals
  plot(actual_data$year, residuals, type = "l",
       main = "Residuals Over Time", xlab = "Year", ylab = "Residual")
  abline(h = 0, col = "red", lty = 2)
  
  # QQ plot
  qqnorm(residuals, main = "Normal Q-Q Plot")
  qqline(residuals, col = "red")
  
  # ACF
  acf(residuals, main = "Autocorrelation")
  
  # Histogram
  hist(residuals, breaks = 20, main = "Residual Distribution", xlab = "Residual")
  
  par(mfrow = c(1, 1))
  
  # 3. Statistical tests
  # Durbin-Watson for autocorrelation
  # Jarque-Bera for normality
  # Breusch-Pagan for heteroskedasticity
  
  list(
    rmse = rmse,
    mae = mae,
    mape = mape,
    residuals = residuals
  )
}


# ==============================================================================
# 3. ROBUSTNESS CHECKS
# ==============================================================================

# Standard robustness tests for publication:

robustness_tests <- function() {
  # 1. Alternative specifications
  #    - Different functional forms (Gompertz, Bass diffusion)
  #    - Alternative lag structures
  #    - Non-linear infrastructure effects
  
  # 2. Subsample analysis
  #    - Urban vs rural
  #    - Early vs late adopters
  #    - Different regions/states
  
  # 3. Placebo tests
  #    - Randomize treatment assignment
  #    - Should show no effect
  
  # 4. Sensitivity to outliers
  #    - Jackknife: drop one year at a time
  #    - Robust regression (M-estimators)
  
  # 5. Parameter stability
  #    - Rolling window estimation
  #    - Chow test for structural breaks
}


# ==============================================================================
# 4. POLICY COUNTERFACTUALS
# ==============================================================================

# For policy analysis section of paper:

policy_experiments <- function() {
  # Define scenarios based on real policy proposals:
  
  scenarios <- list(
    # 1. Status quo baseline
    baseline = c(r = 7500, s = 0, i = 0.05),
    
    # 2. Eliminate tax credits (2025 proposal)
    no_credits = c(r = 0, s = 0, i = 0.05),
    
    # 3. Double infrastructure spending (Green New Deal)
    high_infra = c(r = 7500, s = 0, i = 0.15),
    
    # 4. European-style manufacturer mandates
    mandates = c(r = 5000, s = 8000, i = 0.08),
    
    # 5. Comprehensive policy package
    comprehensive = c(r = 10000, s = 5000, i = 0.20)
  )
  
  # Run all scenarios
  # Compare:
  #   - AV penetration at t=10, t=20, t=30
  #   - Time to 50% penetration
  #   - Total policy cost
  #   - Cost per AV adopted
  #   - Emissions impact (if you add that module)
  
  # Create comparison table for paper
}


# ==============================================================================
# 5. UNCERTAINTY QUANTIFICATION
# ==============================================================================

# Bootstrap confidence intervals for policy effects:

bootstrap_policy_effect <- function(data, policy, n_boot = 1000) {
  boot_results <- matrix(NA, nrow = n_boot, ncol = 31)
  
  for (b in 1:n_boot) {
    # Resample data
    boot_data <- data[sample(nrow(data), replace = TRUE), ]
    
    # Re-estimate parameters
    boot_params <- estimate_from_data(boot_data)
    
    # Simulate
    sim <- simulate_model(state_0, policy, 30, boot_params)
    boot_results[b, ] <- sim$pct_AV
  }
  
  # Compute confidence intervals
  ci_lower <- apply(boot_results, 2, quantile, probs = 0.025)
  ci_upper <- apply(boot_results, 2, quantile, probs = 0.975)
  
  return(list(lower = ci_lower, upper = ci_upper))
}


# ==============================================================================
# 6. TABLES FOR PUBLICATION
# ==============================================================================

# Generate LaTeX tables:

make_table_1_parameters <- function(params, std_errors) {
  # Table 1: Estimated Parameters
  cat("\\begin{table}[h]\n")
  cat("\\centering\n")
  cat("\\caption{Model Parameter Estimates}\n")
  cat("\\begin{tabular}{lcc}\n")
  cat("\\hline\n")
  cat("Parameter & Estimate & Std. Error \\\\\n")
  cat("\\hline\n")
  cat(sprintf("$\\beta_0$ (Intercept) & %.3f & (%.3f) \\\\\n", 
              params$beta_0, std_errors$beta_0))
  cat(sprintf("$\\beta_r$ (Rebate Effect) & %.3f & (%.3f) \\\\\n", 
              params$beta_r, std_errors$beta_r))
  cat(sprintf("$\\beta_i$ (Infrastructure Effect) & %.3f & (%.3f) \\\\\n", 
              params$beta_i, std_errors$beta_i))
  cat(sprintf("$\\delta$ (Retirement Rate) & %.3f & (%.3f) \\\\\n", 
              params$delta, std_errors$delta))
  cat(sprintf("$\\rho$ (Decay Rate) & %.3f & (%.3f) \\\\\n", 
              params$rho, std_errors$rho))
  cat("\\hline\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
}

make_table_2_scenarios <- function(scenario_results) {
  # Table 2: Policy Scenario Comparison
  # Shows AV share at t=10, 20, 30 and total cost
}

make_table_3_robustness <- function(robustness_results) {
  # Table 3: Robustness Checks
  # Alternative specs, subsamples, etc.
}


# ==============================================================================
# 7. FIGURES FOR PUBLICATION
# ==============================================================================

# Publication-quality figures (export as PDF):

make_figure_1_historical_fit <- function() {
  # Figure 1: Model Fit to Historical Data
  # Shows actual vs predicted with confidence bands
  
  pdf("figure1_model_fit.pdf", width = 8, height = 5)
  # [plotting code]
  dev.off()
}

make_figure_2_scenarios <- function() {
  # Figure 2: Policy Scenario Trajectories
  # Multiple lines showing different policies
  
  pdf("figure2_scenarios.pdf", width = 10, height = 6)
  # [plotting code]
  dev.off()
}

make_figure_3_sensitivity <- function() {
  # Figure 3: Sensitivity Analysis
  # Tornado diagram or parameter sweep
  
  pdf("figure3_sensitivity.pdf", width = 8, height = 6)
  # [plotting code]
  dev.off()
}


# ==============================================================================
# 8. MODEL EXTENSIONS TO CONSIDER
# ==============================================================================

# For a more complete model, consider adding:

extensions_to_add <- function() {
  # 1. Spatial heterogeneity
  #    - State-level variation in policies and infrastructure
  #    - Network effects across regions
  
  # 2. Consumer heterogeneity
  #    - Distribution of willingness to pay
  #    - Different consumer segments (urban/rural, income)
  
  # 3. Technology learning
  #    - Endogenous price reductions (learning curve)
  #    - Quality improvements over time
  
  # 4. Industry dynamics
  #    - Firm entry/exit
  #    - Economies of scale in production
  
  # 5. Emissions module
  #    - Track fleet emissions over time
  #    - Social cost of carbon calculations
  
  # 6. Network externalities
  #    - Value of AVs increases with prevalence
  #    - Infrastructure spillovers
  
  # 7. Regulatory constraints
  #    - Safety standards
  #    - CAFE standards
  
  # 8. Behavioral factors
  #    - Habit formation
  #    - Social influence/peer effects
}


# ==============================================================================
# 9. PAPER STRUCTURE TEMPLATE
# ==============================================================================

paper_outline <- function() {
  cat("
TITLE: State-Space Model of Autonomous Vehicle Adoption Under Policy Interventions

ABSTRACT (200 words)
  - Research question
  - Methodology (state-space model)
  - Key findings
  - Policy implications

1. INTRODUCTION
   - Motivation: AV adoption is policy-sensitive
   - Research gap: lack of integrated dynamic models
   - Contribution: tractable framework for policy analysis
   - Preview of results

2. LITERATURE REVIEW
   - Technology adoption models (Rogers, Bass)
   - Vehicle choice models (discrete choice)
   - Policy evaluation in transportation
   - State-space methods in economics

3. MODEL
   3.1 State variables and dynamics
   3.2 Policy instruments
   3.3 Behavioral foundations (sigmoid adoption)
   3.4 Supply-side constraints
   3.5 Identification and estimation strategy

4. DATA AND CALIBRATION
   4.1 Data sources
   4.2 Descriptive statistics
   4.3 Parameter estimation
   4.4 Model fit and validation

5. RESULTS
   5.1 Baseline projections
   5.2 Policy counterfactuals
   5.3 Cost-effectiveness analysis
   5.4 Sensitivity analysis

6. DISCUSSION
   6.1 Policy implications
   6.2 Limitations
   6.3 Future research

7. CONCLUSION

REFERENCES
APPENDIX A: Mathematical derivations
APPENDIX B: Robustness checks
APPENDIX C: Additional tables and figures
")
}


# ==============================================================================
# 10. CHECKLIST BEFORE SUBMISSION
# ==============================================================================

submission_checklist <- function() {
  cat("
PRE-SUBMISSION CHECKLIST:

□ Data
  □ All data sources documented and cited
  □ Replication data prepared
  □ Data cleaning code available
  
□ Estimation
  □ Standard errors calculated
  □ Confidence intervals reported
  □ Robustness checks completed
  
□ Results
  □ All tables formatted in LaTeX
  □ All figures high-resolution PDFs
  □ Numbers in text match tables
  
□ Code
  □ All code commented
  □ Replication script runs start-to-finish
  □ README file with software versions
  □ Code deposited (Dataverse, OSF, etc.)
  
□ Writing
  □ Abstract within word limit
  □ All references formatted
  □ Equations numbered
  □ Notation table included
  
□ Supplementary Materials
  □ Online appendix prepared
  □ Technical derivations complete
  □ Additional robustness checks
")
}


# ==============================================================================
# USAGE EXAMPLE
# ==============================================================================

cat("
TO PREPARE MODEL FOR PUBLICATION:

1. Collect real data and replace dummy parameters:
   real_data <- read.csv('av_historical_data.csv')
   calibrated_params <- estimate_from_data(real_data)

2. Run validation suite:
   validation <- validation_suite(model_results, real_data)

3. Generate all policy scenarios:
   scenarios <- policy_experiments()

4. Create publication tables:
   make_table_1_parameters(calibrated_params, std_errors)
   make_table_2_scenarios(scenarios)

5. Create publication figures:
   make_figure_1_historical_fit()
   make_figure_2_scenarios()
   make_figure_3_sensitivity()

6. Run robustness checks:
   robustness_tests()

7. Export replication package
")

cat("\n=== PUBLICATION GUIDE LOADED ===\n")