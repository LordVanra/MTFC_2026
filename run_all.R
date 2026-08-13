# =============================================================================
#  run_all.R  —  SINGLE ENTRY POINT FOR THE FULL PIPELINE
# =============================================================================
#  Usage:  Rscript run_all.R            (from the repository root)
#
#  Order matters. The state-space model writes both the adoption trajectories
#  and the DERIVED policy costs that the dispersion model consumes. Running the
#  dispersion model first, or running any script in deprecated/, will either
#  fail loudly or use stale inputs.
#
#      config/params.R                       <- every parameter, one place
#             |
#      StateModel/plot_model.r               -> av_all_scenarios.csv
#             |                              -> outputs/policy_costs.csv
#      CellTransmission/plot_models.r        -> cumulative_flow_matrix.csv
#             |
#      ConvectionDiffusion/plot_results.r    -> outputs/dispersion_results_*.csv
#             |
#      outputs/paper_values.txt              <- every number the paper cites
# =============================================================================

t_start <- Sys.time()
options(warn = 1)
set.seed(42)   # global determinism; individual scripts re-seed where needed

REPO <- normalizePath(".")
stopifnot(file.exists(file.path(REPO, "config", "params.R")))
source(file.path(REPO, "config", "params.R"))
dir.create(file.path(REPO, "outputs"), showWarnings = FALSE)

banner <- function(msg) cat("\n", strrep("=", 70), "\n  ", msg, "\n",
                            strrep("=", 70), "\n", sep = "")

run_in <- function(dir, script) {
  banner(paste("RUN", file.path(dir, script)))
  old <- setwd(file.path(REPO, dir))
  on.exit(setwd(old), add = TRUE)
  source(script, echo = FALSE, local = new.env())
  invisible(NULL)
}

# -----------------------------------------------------------------------------
banner("PARAMETER PROVENANCE CHECK")
print_provenance()
D <- derive_params()
cat(sprintf("\nDerived: S_BOS = %.0f veh/yr  (S_US %.3g x C_0 %.0f / C_US %.3g)\n",
            D$S_BOS, par_val("SSM","S_US"), par_val("SSM","C_0"), par_val("SSM","C_US")))
cat(sprintf("Derived: EF_ICE = %.4f g/km,  EF_AV = %.4f g/km  (Eq. 15, explicit sum)\n",
            D$ef_ice, D$ef_av))

# -----------------------------------------------------------------------------
banner("CDM STEADY-STATE CONVERGENCE")
source(file.path(REPO, "scripts", "convergence_test.R"))
conv <- run_convergence()

# -----------------------------------------------------------------------------
run_in("StateModel", "plot_model.r")

# StateModel writes av_all_scenarios.csv next to itself; the dispersion model
# reads it from its own directory. Copy explicitly rather than relying on a
# stale duplicate being present.
file.copy(file.path(REPO, "StateModel", "av_all_scenarios.csv"),
          file.path(REPO, "ConvectionDiffusion", "av_all_scenarios.csv"),
          overwrite = TRUE)

# -----------------------------------------------------------------------------
run_in("CellTransmission", "plot_models.r")

file.copy(file.path(REPO, "CellTransmission", "outputs", "cumulative_flow_matrix.csv"),
          file.path(REPO, "ConvectionDiffusion", "cumulative_flow_matrix.csv"),
          overwrite = TRUE)

# -----------------------------------------------------------------------------
run_in("ConvectionDiffusion", "plot_results.r")

# -----------------------------------------------------------------------------
banner("DONE")
cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("Outputs in outputs/. Every cited number should be traceable to\n",
    "outputs/paper_values_ssm.txt and outputs/policy_costs.csv.\n", sep = "")
