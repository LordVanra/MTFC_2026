# =============================================================================
#  config/params.R  —  SINGLE SOURCE OF TRUTH FOR ALL MODEL PARAMETERS
# =============================================================================
#  Every model script sources this file. No parameter may be defined anywhere
#  else. Appendix tables in the manuscript are generated from this file via
#  scripts/export_param_tables.R, so code and paper cannot drift apart.
#
#  Each entry carries: value, units, and the source it is defensible from.
#  A parameter with source = "UNSOURCED" must not appear in the paper without
#  being flagged as a modelling assumption with a sensitivity range.
# =============================================================================

# ---- provenance helper -------------------------------------------------------
P <- function(value, units, source, note = NA_character_) {
  list(value = value, units = units, source = source, note = note)
}

# =============================================================================
# 1. STATE-SPACE MODEL (SSM)
# =============================================================================
SSM <- list(

  # --- fleet dynamics ---
  delta_1 = P(0.045, "1/yr",
              "S&P Global Mobility (2025); US avg vehicle age 12.8 yr",
              "Replacement-FLOW proxy, NOT 1/average-age. Revised from 0.07."),

  delta_2 = P(0.023, "1/yr",
              "BEA transportation capital stock; 40-60 yr service life",
              "CORRECTED from 0.05. Production code previously used 0.05 while
               the manuscript reported 0.023 - this was a live code/paper split."),

  S_US    = P(15.8e6, "vehicles/yr",
              "FRED TOTALSA (Jan 2026 SAAR ~15.4M); NADA 2024 (~15.9M)"),

  C_US    = P(283e6, "vehicles",
              "FHWA Highway Statistics - US light-duty registered vehicles",
              "Used only to scale national sales to Boston. VERIFY before submission."),

  # S_BOS is DERIVED, not asserted. See derive_params().
  # Previously hardcoded as 26970 with no shown derivation.

  Price_AV = P(48000, "USD",
               "Level-2 autonomy add-on $2.5-5k; full AV premium midpoint",
               "SOURCE IS WEAK (aggregator site). Flagged for replacement with a
                primary manufacturer/teardown source. Load-bearing: enters
                price_adv = r / Price_AV. Was 40000 in code vs 48000 in paper."),

  # --- adoption response coefficients ---
  # These are ELICITED SCENARIO WEIGHTS calibrated to a published adoption
  # envelope, NOT estimated causal elasticities. Do not report standard errors.
  beta_0 = P(-2.026, "log-odds",
             "Calibrated to Morgan Stanley / VTPI / WEF adoption envelope",
             "Code previously used -3.5 while paper reported -2.026."),

  beta_1 = P(3.250, "log-odds per unit price advantage",
             "Directional evidence: J.D. Power (2024) 87% federal EV credit uptake"),

  beta_2 = P(3.997, "log-odds per unit normalised infrastructure",
             "Directional evidence: ICCT (2025) charging build-out; VoxChina (2024)"),

  # --- infrastructure normalisation ---
  I_ref = P(66214286, "index units",
            "UNSOURCED - internal normalisation constant",
            "Was undocumented in the manuscript entirely. Must either be
             derived from a real infrastructure-stock figure or presented
             explicitly as a scaling choice with a sensitivity sweep."),

  infra_transform = P("sqrt", "-",
                      "Modelling choice",
                      "Concave diminishing-returns map I_scaled = sqrt(min(I/I_ref,1)).
                       NOTE: manuscript Eq.(2) and Appendix A.1 both print a
                       LOGISTIC map. Paper must be corrected to match this, or
                       this changed to match the paper. Options: 'sqrt' | 'linear' | 'logistic'."),

  # --- supply capacity ---
  cap_0 = P(0.20, "fraction of annual sales",
            "UNSOURCED - Year-1 adoption ceiling assumption"),

  gamma = P(0.00002, "capacity fraction per $/AV of supply subsidy",
            "UNSOURCED - calibration constant"),

  # --- initial state (Boston) ---
  A_0 = P(98,      "vehicles",    "Derived: AV registrations, Boston, base year"),
  C_0 = P(512539,  "vehicles",    "Boston registered light-duty vehicles"),
  I_0 = P(6983872, "index units", "Baseline transport capital stock proxy"),
  K_0 = P(0.01,    "fraction",    "Initial AV production capacity share"),

  horizon = P(30, "yr", "Modelling choice")
)

# =============================================================================
# 2. CELL TRANSMISSION MODEL (CTM)
# =============================================================================
CTM <- list(
  v_f     = P(50,   "km/h",    "MassDOT / Mass.gov urban arterial statutory limit"),
  rho_j   = P(120,  "veh/km",  "Urban arterial mid-range (115-155 veh/km/lane)",
                               "REPLACE 'Anonymous (2024)' citation before submission."),
  q_max   = P(1800, "veh/h",   "I-93 empirical 2000-2300 veh/h/lane, arterial discount applied"),
  n_cells = P(28,   "cells",   "Boston ZIP codes with non-zero ACS population"),

  # AV effect on the fundamental diagram
  alpha_v   = P(0.10, "-", "Rahman et al. (2021); Fagnant & Kockelman (2015)",
                           "Assumption with wide literature range - include in sensitivity sweep."),
  alpha_rho = P(0.08, "-", "Rahman et al. (2021)",
                           "Assumption - include in sensitivity sweep."),

  presence_factor = P(0.55, "-", "Urban transport-planning norm, analogical",
                                 "Tier-3 derived. Tested 0.45-0.65."),

  # O-D corridor selection. The manuscript contains a contradiction:
  # Sec 4.3.1 states 106 active corridors; Table 11 states top-50 of 378.
  od_pairs_total    = P(378, "pairs", "28x28 minus diagonal, adjacency-filtered"),
  od_pairs_retained = P(50,  "pairs", "Top-N by flow volume",
                                      "RESOLVE the 106-vs-50 contradiction and report
                                       the share of total modelled VMT retained."),

  time_step = P(1, "h", "Modelling choice")
)

# =============================================================================
# 3. CONVECTION-DIFFUSION MODEL (CDM)
# =============================================================================
CDM <- list(
  u       = P(3.5,   "m/s",   "NOAA local climatological data (New York proxy)",
                              "New York used as proxy for Boston - stated limitation."),
  w       = P(0.1,   "m/s",   "Modelling choice - weak vertical advection"),
  D_h     = P(50,    "m^2/s", "Pasquill & Smith (1983) urban boundary layer"),
  D_z     = P(10,    "m^2/s", "Turner (1970) mixed-layer estimate"),
  k_dep   = P(3e-5,  "1/s",   "EPA regulatory first-order dry deposition for PM2.5"),

  dx      = P(20,    "m",     "Grid resolution (downwind)"),
  dz      = P(5,     "m",     "Grid resolution (vertical)"),

  # --- DIMENSIONAL FIX -------------------------------------------------------
  # The source term was previously S = E_eff / dz, giving g/(m^2 s) injected
  # into a field carrying g/m^3. A line source of strength q [g/(m s)] spread
  # over a cell of cross-section (street_width x dz) has volumetric strength
  #     S = q / (street_width * dz)   [g/(m^3 s)]
  # street_width is now explicit and must be justified in the methods section.
  street_width = P(20, "m",
                   "Typical Boston arterial kerb-to-kerb width incl. parking lanes",
                   "NEW EXPLICIT PARAMETER introduced by the dimensional
                    correction. Was implicitly (and wrongly) 1 m."),

  # --- TEMPORAL BASIS FIX ----------------------------------------------------
  # Flows were divided by 30*24*3600 s, smearing a month of traffic across all
  # hours including 03:00, while the CTM runs on peak-hour capacity. The two
  # models were on different temporal bases. We now work in PEAK HOUR to match
  # the CTM, and convert to daily/annual means with an explicit diurnal factor.
  temporal_basis  = P("peak_hour", "-", "Consistency with CTM",
                                        "Options: 'peak_hour' | 'daily_mean'"),
  peak_hour_share = P(0.095, "-",
                      "Share of daily traffic in the peak hour; urban norm ~8-11%",
                      "VERIFY against MassDOT count data before submission."),
  diurnal_factor  = P(0.42, "-",
                      "Ratio of 24h-mean to peak-hour concentration",
                      "PLACEHOLDER - derive from MassDOT hourly counts."),

  # --- STEADY STATE FIX ------------------------------------------------------
  # T was 300 s over a 1000 m domain at u = 3.5 m/s (transit ~286 s), so the
  # field was sampled mid-transient. Convergence is now tested explicitly.
  T_sim          = P(2400, "s", "Set to reach steady state; verified by convergence test"),
  x_domain_max   = P(1000, "m", "Downwind extent"),
  z_domain_max   = P(100,  "m", "Vertical extent"),
  converge_tol   = P(0.01, "-", "Relative change in ground field between halvings"),

  # ground_mean was averaged only over cells > 0, so it depended on how far the
  # plume had travelled. Now averaged over a fixed window.
  ground_mean_window = P(c(0, 1000), "m", "Fixed averaging window for ground mean")
)

# =============================================================================
# 4. EMISSION FACTORS
# =============================================================================
#  Manuscript Eq.(15) states EF = exhaust + brake + tire, giving
#  ICE = 0.012 + 0.001 + 0.0006 = 0.0136 g/km and AV = 0.003 g/km.
#  The code instead read a SINGLE column from emission_data.csv
#  (ICE = 0.005, AV = 0.002) and never added the brake/tire rows that exist
#  in that same CSV. Eq.(15) as printed was not what ran.
#  These values now match the manuscript and are summed explicitly in code.
# =============================================================================
EF <- list(
  ice = list(
    exhaust = P(0.012,  "g/km", "EPA MOVES4 light-duty gasoline"),
    brake   = P(0.001,  "g/km", "EPA MOVES4 brake wear"),
    tire    = P(0.0006, "g/km", "EPA MOVES4 tire wear"),
    co2     = P(250,    "g/km", "EPA typical passenger vehicle")
  ),
  av = list(
    exhaust = P(0.0,    "g/km", "BEV tailpipe elimination"),
    brake   = P(0.0005, "g/km", "50% reduction vs ICE from regenerative braking",
                                "ASSUMPTION - regen share varies by drive cycle."),
    tire    = P(0.0006, "g/km", "NOT reduced: BEV mass offsets smoother driving"),
    residual_pm = P(0.0019, "g/km", "BEV-equivalent residual PM term",
                    "Chosen so total AV EF = 0.003 g/km per manuscript Eq.(15).
                     THIS IS A BACK-SOLVED PLUG. Either source it properly or
                     drop it and let the AV factor fall to 0.0011 g/km."),
    co2_bev      = P(100, "g/km", "BEV lifecycle, current grid mix"),
    co2_sensor   = P(8,   "g/km", "Sensor/computing overhead",
                      "CONTESTED in the literature - needs 2-3 independent sources.
                       Applied to lifecycle CO2 ONLY, never to PM2.5.")
  ),
  # road dust resuspension deliberately excluded - not available for all classes
  include_resuspension = P(FALSE, "-", "Documented exclusion")
)

# =============================================================================
# 5. VMT REBOUND
# =============================================================================
#  This drives the paper's headline counter-intuitive result and currently
#  rests on a single ride-hailing study. Treated as a first-class uncertain
#  parameter with an explicit sweep, not a fixed assumption.
# =============================================================================
REBOUND <- list(
  low      = P(0.30, "-", "Low-funding scenario rebound"),
  high     = P(0.40, "-", "High-funding scenario rebound (range 0.35-0.45)"),
  sweep    = P(seq(0.0, 0.6, by = 0.05), "-", "Sensitivity sweep range"),
  source   = P("Brasuell (2021) ride-hailing 97-118% VMT increase", "-",
               "EXTRAPOLATION from ride-hailing to private AV ownership.
                Needs induced-demand economics literature to support.")
)

# =============================================================================
# 6. POLICY COST MODEL
# =============================================================================
#  Costs were a HARDCODED table inside a plotting script, and the in-model
#  expression was dimensionally incoherent:
#      annual_spend = r*(adoption*S) + s + i*S*Price_AV/1000
#  - `s` is $/AV but was added as a bare scalar (never multiplied by units)
#  - `i` is already $/yr but was multiplied by S*Price_AV/1000 (~1.08e6),
#    inflating infrastructure spend by six orders of magnitude.
#  Costs are now DERIVED from the lever time paths. See scripts/cost_model.R.
# =============================================================================
COST <- list(
  discount_rate = P(0.00, "1/yr", "Undiscounted by default",
                    "Set >0 and report both if the target journal expects NPV."),
  units = P("USD", "-", "Reported in $B in tables")
)

# =============================================================================
# DERIVED QUANTITIES
# =============================================================================
derive_params <- function() {
  v <- function(x) x$value

  S_BOS <- v(SSM$S_US) * v(SSM$C_0) / v(SSM$C_US)

  ef_ice <- v(EF$ice$exhaust) + v(EF$ice$brake) + v(EF$ice$tire)
  ef_av  <- v(EF$av$exhaust)  + v(EF$av$brake)  + v(EF$av$tire) + v(EF$av$residual_pm)

  list(
    S_BOS  = S_BOS,
    ef_ice = ef_ice,
    ef_av  = ef_av,
    rho_c  = v(CTM$q_max) / v(CTM$v_f)
  )
}

# Flat accessor: par("SSM","delta_1") -> 0.045
par_val <- function(group, name, sub = NULL) {
  g <- get(group)
  x <- if (is.null(sub)) g[[name]] else g[[name]][[sub]]
  x$value
}

# Every parameter whose source is weak, absent, or back-solved.
unsourced_report <- function() {
  flags <- c()
  walk <- function(lst, path = "") {
    for (nm in names(lst)) {
      item <- lst[[nm]]
      p <- if (path == "") nm else paste0(path, "$", nm)
      if (is.list(item) && !is.null(item$source)) {
        bad <- grepl("UNSOURCED|WEAK|PLACEHOLDER|back-solved|plug|VERIFY|CONTESTED|REPLACE|RESOLVE",
                     paste(item$source, item$note), ignore.case = TRUE)
        if (bad) flags <<- c(flags, p)
      } else if (is.list(item)) walk(item, p)
    }
  }
  for (g in c("SSM", "CTM", "CDM", "EF", "REBOUND", "COST")) walk(get(g), g)
  flags
}
