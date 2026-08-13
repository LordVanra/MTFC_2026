# Deprecated scripts — DO NOT RUN

These are retained for provenance only. None of them produced any figure or
number in the manuscript, and running them corrupts the pipeline.

## `av_adoption_model.r`, `parameter_checking.r`
Superseded by `StateModel/plot_model.r`. Both are earlier development
versions that:

- use a **national** initial state (`A = 10e6, C = 250e6`) rather than Boston
  (`A = 98, C = 512539`);
- carry stale parameters (`delta_1 = 0.07`, `delta_2 = 0.05`, `S = 350000`,
  where `S` is annotated "15 million" but is neither);
- contain a dead-code defect in `step_model()`: `a_uncapped` is computed with a
  normalised infrastructure term and then immediately **overwritten** by a
  version using the raw index `I_i` (~7e6), which saturates the logistic to
  exactly 1.0 in every year of every scenario, so the rebate coefficient has no
  effect on output;
- add `sigmoid(k_i)` (never below 0.5) to the capacity recurrence, forcing
  `K -> 1` within two years even under zero policy;
- hard-override `n_sim <- 1` inside `monte_carlo()` while advertising 1000;
- and, critically, **write to `StateModel/av_all_scenarios.csv`** — the file the
  dispersion model consumes. Running either one silently replaces the
  pipeline's input with output from the defective model.

`plot_model.r` fixes each of these (its own comments record two of the fixes).

## `convection_diffusion.r`
A standalone demo that runs the entire dispersion routine on
`Heavy_Duty_Diesel` (0.05 g/km — an order of magnitude above the passenger
factor used in the pipeline) and writes `pm25_2d_results.pdf`. It is not part
of the pipeline and shares no code path with it. It also predates the
dimensional correction to the source term.
