# This probably doesn't work

generate_synthetic_data <- function(n = 100, params) {
  set.seed(123)
  price_adv <- runif(n, 0, 0.3)
  infrastructure <- runif(n, 0, 1)
  z <- params$beta_0 + params$beta_1 * price_adv + params$beta_2 * infrastructure
  a_true <- sigmoid(z)
  adoption <- pmin(pmax(a_true + rnorm(n, 0, 0.05), 0), 1)
  data.frame(adoption = adoption, price_adv = price_adv, infrastructure = infrastructure)
}

estimate_parameters <- function(data) {
  model <- glm(adoption ~ price_adv + infrastructure, data = data,
               family = quasibinomial(link = "logit"))
  coefs <- coef(model)
  list(beta_0 = coefs[1], beta_1 = coefs[2], beta_2 = coefs[3], model = model)
}

cat("PARAMETER ESTIMATION (Synthetic Data)\n")
synthetic_data <- generate_synthetic_data(200, params)
est_params <- estimate_parameters(synthetic_data)
cat(sprintf("   beta_0: True = %.3f, Est = %.3f\n", params$beta_0, est_params$beta_0))
cat(sprintf("   beta_1: True = %.3f, Est = %.3f\n", params$beta_1, est_params$beta_1))
cat(sprintf("   beta_2: True = %.3f, Est = %.3f\n\n", params$beta_2, est_params$beta_2))
