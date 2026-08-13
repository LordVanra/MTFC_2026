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
              "S&P Global Mobility, 'U.S. Vehicle Age Rises Again to 12.8 Years in 2025',
               press release 21 May 2025. VERIFIED: this is the published scrappage
               rate (scrapped / vehicles-in-operation), not 1/average-age.
               https://press.spglobal.com/2025-05-21-U-S-Vehicle-Age-Rises-Again-to-12-8-Years-in-2025",
              "SOURCED. The paper's 'replacement-flow proxy, not inverse average age'
               caveat is now backed by the source's own definition. For a non-constant
               hazard see Greene & Leard (2023), Baker Center, Univ. of Tennessee:
               median lifetimes cars ~17 yr, SUVs ~20 yr, pickups ~25 yr."),

  delta_2 = P(0.0202, "1/yr",
              "Fraumeni, B.M. & Kornfeld, R. (2024), 'Constructing BEA Highways and
               Streets Net Wealth Stocks...', NBER Working Paper 32753. Geometric rate
               from 45-yr service life x 0.91 declining-balance (Hulten-Wykoff).
               https://www.nber.org/system/files/working_papers/w32753/w32753.pdf
               Stock: BEA Fixed Assets Table 7.1B, Structures > Highways and streets.",
              "SOURCED, and CORRECTED again: 2.02% is BEA's actual rate. The paper's
               2.3% was close but not the published figure; the code used 0.05."),

  S_US    = P(16.2e6, "vehicles/yr",
              "NADA Market Beat, 31 Dec 2025: 'New light-vehicle sales totaled 16.2
               million units in 2025, an increase of 2.4% compared to 2024.'
               https://www.nada.org/nada/nada-headlines/december-2025-market-beat-new-light-vehicle-sales-totaled-162-million-units",
              "CITATION CORRECTED: the paper cited FRED TOTALSA, which is TOTAL vehicle
               sales and includes heavy trucks. The light-duty series is ALTSALES
               (https://fred.stlouisfed.org/series/ALTSALES). Use ALTSALES or NADA."),

  C_US    = P(289e6, "vehicles",
              "S&P Global Mobility, 21 May 2025: 289 million light vehicles in
               operation in the U.S. as of 1 Jan 2025. Same release as delta_1.
               https://press.spglobal.com/2025-05-21-U-S-Vehicle-Age-Rises-Again-to-12-8-Years-in-2025",
              "SOURCED. Note FHWA Highway Statistics Table MV-1 (2024) gives 297.5M
               ALL motor vehicles but does NOT split light/heavy - MV-9 is currently
               suspended - so FHWA cannot supply a light-duty stock. S&P is the series
               BTS itself uses for NTS Table 1-26."),

  # S_BOS is DERIVED, not asserted. See derive_params().
  # Previously hardcoded as 26970 with no shown derivation.

  Price_AV = P(51377, "USD",
               "DERIVED: average new-vehicle transaction price $49,855 (Kelley Blue
                Book / Cox Automotive, July 2026 ATP report,
                https://www.coxautoinc.com/insights/july-2026-atp-report/)
                + Level-2 ADAS package premium midpoint $1,522 (range $750-$2,295;
                Goddard, McDonald, Wei & Batra (2022), 'Advanced Driver Assistance
                Systems in Top-Selling Vehicles in the United States', Findings,
                DOI 10.32866/001c.38291).",
               "REPLACES the aggregator citation. NOTE the previous $48,000 is
                suspiciously close to the plain new-vehicle ATP, suggesting the
                aggregator conflated an average car price with an AV price.
                UPPER BOUND for L4: ~$200,000/vehicle (Riggs & Richardson 2024,
                SSRN 4998828 - working paper, cite as such). Litman (2023, VTPI)
                gives the declining-premium trajectory but no single figure.
                Run the L4 anchor as a sensitivity case."),

  # --- adoption response coefficients ---
  # These are ELICITED SCENARIO WEIGHTS calibrated to a published adoption
  # envelope, NOT estimated causal elasticities. Do not report standard errors.
  beta_0 = P(-3.199, "log-odds",
             "CALIBRATED TO OBSERVED BOSTON DATA. Massachusetts Vehicle Census
              zero-emission VMT share 2020-2025, inverted through the fleet-turnover
              recurrence to a new-sales share, with a relative-mileage correction
              m = 1.35 (ZEVs are newer and driven more than the fleet average, so
              VMT share exceeds count share). See scripts/calibrate_adoption.R.",
             "REPLACES the Morgan Stanley / VTPI / WEF envelope calibration. This is
              the paper's single biggest methodological upgrade: the coefficient is
              now fitted to six years of the actual city rather than to other
              people's forecasts.
              Sensitivity across m = 1.0 to 1.5 gives beta_0 in [-3.32, -2.85];
              report m = 1.35 as central with that span.
              Manuscript reported -2.026; production code used -3.5."),

  beta_0_range = P(c(-3.319, -2.846), "log-odds",
                   "Span from the relative-mileage sensitivity, m = 1.5 to m = 1.0"),

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
  A_0 = P(98,      "vehicles",    "Derived: AV registrations, Boston, base year",
                                  "STILL UNSOURCED - no registry publishes L4 AV
                                   registrations for Boston. Consider setting to 0
                                   and starting adoption from the policy levers."),
  # --- OBSERVED ADOPTION, Massachusetts Vehicle Census -------------------------
  # Zero-emission share of Boston-garaged VMT, 2020-2025. This is REAL OBSERVED
  # Boston adoption data and is the calibration target the beta coefficients
  # never had. The manuscript calibrates to a national/industry forecast
  # envelope; this is the city itself.
  obs_zev_share = P(c("2020"=0.005548, "2021"=0.007764, "2022"=0.015331,
                      "2023"=0.035458, "2024"=0.046527, "2025"=0.049204), "-",
                    "Massachusetts Vehicle Census (MassDOT + MAPC), MVC Annual Zip
                     Code file, Boston-garaged vehicles, share of ANNUALVMT.
                     https://gis.data.mass.gov/datasets/dce8f7ce88f449f08150af4183b9a983
                     Derived by scripts/ingest_vehicle_census.R.",
                    "USE THIS TO CALIBRATE beta_0/beta_1/beta_2 instead of the
                     Morgan Stanley/VTPI envelope. Note it is a VMT share, not a
                     fleet-count share - ZEVs are newer and may be driven more or
                     less than average, so the two are not interchangeable.
                     Note also the 2024->2025 inflection: +1.11pp in 2023-24 but
                     only +0.27pp in 2024-25. A logistic fitted through 2023 would
                     badly overshoot."),
  obs_hybrid_share = P(c("2020"=0.03844, "2021"=0.04243, "2022"=0.05223,
                         "2023"=0.06283, "2024"=0.07830, "2025"=0.09611), "-",
                       "Same source. Hybrids are still growing faster than ZEVs in
                        absolute share terms.",
                       "The model treats the fleet as ICE-or-AV/BEV binary. Nearly
                        10% of Boston VMT is hybrid, which sits between the two
                        emission factors and is not represented."),
  obs_total_vmt = P(3.120e9, "miles/yr",
                    "MVC, Boston-garaged, 2025. Up 18.4% from 2020.",
                    "OBSERVED VMT GROWTH OF +18.4% IN FIVE YEARS with no AVs
                     involved. The model has no exogenous VMT growth term, so this
                     is a baseline trend the rebound analysis should be measured
                     against rather than attributed to automation."),

  C_0 = P(266899,  "vehicles",
          "DERIVED from the Massachusetts Vehicle Census 2025: 3.043e9 annual VMT
           across the 28 modelled ZIPs, divided by 11,400 mi/vehicle/yr (FHWA
           Highway Statistics VM-1 US average).
           https://gis.data.mass.gov/datasets/dce8f7ce88f449f08150af4183b9a983",
          "The MVC ZIP file publishes VMT but NOT vehicle counts, so this is a
           derivation, not a measurement. Urban Boston almost certainly drives less
           than the US average, which makes 266,899 an UNDERSTATEMENT of the fleet
           (fewer miles per vehicle would imply more vehicles for the same VMT).
           Cross-check: the 2014 MVC city count was 237,224 (Go Boston 2030,
           pp.54-55), so a 2025 figure of ~267k implies 1.2%/yr growth - plausible.
           Sanity: 266,899 / ~270,000 households = 0.99 veh/household.
           The manuscript's 512,539 remains ~1.9x too high.
           BETTER OPTION: the model is ultimately about VMT, not vehicle counts.
           Consider reframing the exposure calculation on observed ZIP-level VMT
           (outputs/mvc_zip_vmt.csv) and dropping the count entirely."),
  I_0 = P(6983872, "index units", "Baseline transport capital stock proxy"),
  K_0 = P(0.01,    "fraction",    "Initial AV production capacity share"),

  horizon = P(30, "yr", "Modelling choice")
)

# =============================================================================
# 2. CELL TRANSMISSION MODEL (CTM)
# =============================================================================
CTM <- list(
  v_f     = P(50,   "km/h",
              "Observed 85th-percentile operating speed on Boston streets = 31 mph
               (49.9 km/h), unchanged before/after the limit change: Hu, W. & Cicchino,
               J.B., 'Lowering the speed limit from 30 mph to 25 mph in Boston: effects
               on vehicle speeds', Injury Prevention (2019), reported by IIHS
               https://www.iihs.org/news/detail/city-drivers-slow-down-for-lower-speed-limit-in-boston",
              "VALUE UNCHANGED but now sourced to observed speeds rather than the
               statutory limit. Boston's posted default is 25 mph (40 km/h) since
               9 Jan 2017 (https://www.boston.gov/news/bostons-new-default-speed-limit-25-mph-effective-jan-9-2017),
               and does NOT cover state-owned roadways within the city. Using observed
               85th-percentile speed as free-flow speed is the defensible choice; cite
               the Injury Prevention paper directly, not the IIHS summary."),
  rho_j   = P(131,  "veh/km",
              "DERIVED from jam spacing: 1000 m / 7.62 m (25 ft) = 131 veh/km/lane.
               25 ft stopped spacing per Wu, X., Liu, H.X. & Geroliminis, N. (2011),
               'An empirical analysis on the arterial fundamental diagram',
               Transportation Research Part B 45(3):255-266, Sec.5
               (https://infoscience.epfl.ch/server/api/core/bitstreams/8eb2e456-94be-4df9-b17c-ac81e154b8e5/content);
               corroborated by FDOT Design Manual (2024) Sec.212.14.2 queue lengths.",
              "REPLACES the unusable 'Anonymous (2024)' citation. Present as an
               EXPLICIT DERIVATION in the text, not as a borrowed number - no source
               found stating an arterial jam density directly. Lloret-Batlle & Zheng
               (2023), TR-B 173:162-175, DOI 10.1016/j.trb.2023.02.007, is the natural
               citation AND makes exactly this criticism of the field (parameters
               'generally assumed with no estimation from data whatsoever').
               Obtain it via institutional access. Do NOT import freeway jam densities
               (100-110 veh/km/lane) - motorway values are lower than arterial."),
  q_max   = P(950,  "veh/h",
              "DERIVED per FHWA method: capacity = saturation flow x (g/C) x lanes.
               Base saturation flow 1,900 pc/h/lane (FHWA HPMS Field Manual App.N p.N-19,
               https://www.fhwa.dot.gov/OHIM/hpmsmanl/pdf/appn.pdf; confirmed HCM 2000
               p.16-10). Urban signalised default green ratio g/C = 0.50 (FHWA-PL-18-003
               Table 3 p.11, https://www.fhwa.dot.gov/policyinformation/pubs/pl18003/hpms_cap.pdf).
               1900 x 0.50 = 950 veh/h/lane.",
              "MAJOR CORRECTION: the previous 1,800 veh/h was a FREEWAY capacity with an
               undocumented 'arterial discount'. Signalised arterial through-capacity is
               about half that, because green time is shared. Show the g/C in the text."),
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
  street_width = P(16.5, "m",
                   "DERIVED from Boston Transportation Department, 'Boston Complete
                    Streets Design Guidelines' (2013), Minimum Widths for Roadway Lanes,
                    Arterial classification: travel lane 10 ft, parking lane 7 ft.
                    4 travel lanes + parking both sides = 4(10) + 2(7) = 54 ft = 16.5 m.
                    https://www.boston.gov/sites/default/files/file/2020/09/Minimum%20Widths%20for%20Roadway%20Lanes.pdf",
                   "NEW EXPLICIT PARAMETER from the dimensional correction (was
                    implicitly 1 m). BCS publishes MINIMA, not typicals, and no total
                    kerb-to-kerb dimension - the sum above is a derivation and must be
                    presented as one. Sensible sensitivity range 10.4 m (2 lanes +
                    parking) to 16.5 m. This parameter sets the absolute concentration
                    scale linearly, so it belongs in the sensitivity table."),

  # --- TEMPORAL BASIS FIX ----------------------------------------------------
  # Flows were divided by 30*24*3600 s, smearing a month of traffic across all
  # hours including 03:00, while the CTM runs on peak-hour capacity. The two
  # models were on different temporal bases. We now work in PEAK HOUR to match
  # the CTM, and convert to daily/annual means with an explicit diurnal factor.
  temporal_basis  = P("peak_hour", "-", "Consistency with CTM",
                                        "Options: 'peak_hour' | 'daily_mean'"),
  peak_hour_share = P(0.0724, "-",
                      "K-factor. EMPIRICAL, Boston area: MassDOT, 'Mount Auburn Street
                       Corridor Study, Appendix A: Traffic Counts' - Mount Auburn St EB
                       east of Homer Ave, Cambridge MA, Wed 27 Apr 2016. Peak hour
                       (07:00) 977 veh/h of 13,494 veh/day = 7.24%.
                       https://www.mass.gov/doc/mount-auburn-report-appendix-a-traffic-counts/download
                       Consistent with FHWA range 7-12% (Traffic Data Computation Method
                       Pocket Guide, FHWA-PL-18-027, 2018, pp.44-45) and with MassDOT's
                       urban 100th-highest-hour convention (PDDG Ch.3 Sec.3.4).",
                      "SOURCED, replacing the 0.095 placeholder. Caveat honestly: one
                       direction, one arterial, one weekday, single-day peak - not a
                       design-standard K. FHWA's HPMS default of 11% is an imputation
                       fallback for missing data, not an urban best estimate."),
  diurnal_factor  = P(0.576, "-",
                      "DERIVED, not independently sourced: mean-to-peak ratio =
                       1/(24 x K). From the same Mount Auburn count, 24-h mean 562.3
                       veh/h / peak 977 veh/h = 0.576, which equals 1/(24 x 0.0724)
                       by construction.",
                      "MUST stay algebraically tied to peak_hour_share - do not source
                       the two independently or the paper becomes internally
                       inconsistent. Computed in derive_params()."),

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
  # ===========================================================================
  #  ALL VALUES BELOW ARE FROM AN ACTUAL MOVES5.0.1 RUN.
  #  MOVES5.0.1, Default Scale, Emission Rates, Suffolk County MA (FIPS 25025),
  #  urban unrestricted access (roadTypeID 5), hour 07:00-07:59, weekday,
  #  average-speed bin 7 (~30 mph = 48 km/h, matching CTM free-flow 50 km/h),
  #  fleet-weighted 40% passenger car / 60% passenger truck.
  #  RunSpecs: moves/BOS_2025.mrs, moves/BOS_2050.mrs
  #  Raw output: moves/output/*_rateperdistance.txt
  #  Reproduce: Rscript moves/ingest_moves_rates.R
  # ===========================================================================
  ice = list(
    exhaust = P(0.0013719, "g/km",
                "MOVES5.0.1 run, gasoline, processID 1 (Running Exhaust), 2025, January.",
                "The manuscript's 0.012 g/km was ~9x too high. July is 0.0015427 -
                 the seasonal difference is only +12%, much smaller than expected.
                 By 2050 this falls to 0.0001838 g/km, a 87% decline from fleet
                 turnover and standards phase-in alone."),
    brake   = P(0.0031654, "g/km",
                "MOVES5.0.1 run, gasoline, processID 9 (Brakewear), 2025.",
                "Manuscript used 0.001 g/km - 3.2x too low. Note this is well above
                 the 0.00172 g/km implied by EPA-420-R-20-014's 2.77 mg/mile figure,
                 because that document reports a light-duty passenger car average
                 whereas this is a 40/60 car/truck fleet weighting at arterial speeds."),
    tire    = P(0.0010330, "g/km",
                "MOVES5.0.1 run, gasoline, processID 10 (Tirewear), 2025.",
                "Manuscript used 0.0006 g/km - 1.7x too low."),
    co2     = P(264, "g/km",
                "Argonne/DOE R&D GREET 2024 cradle-to-grave gasoline LDV:
                 425 g CO2e/mile = 264 g/km.
                 https://www.energy.gov/sites/default/files/2025-01/eere-greet-ldv-fact-sheet_january-2025.pdf",
                "NOT from the MOVES run - CO2 was not selected. Add pollutant 98
                 (CO2 Equivalent) to the RunSpec if a MOVES-consistent figure is wanted.")
  ),
  av = list(
    exhaust = P(0.0, "g/km",
                "MOVES5.0.1 run, electricity, processID 1: exactly zero.",
                "BEV tailpipe elimination, confirmed by the run rather than assumed."),
    brake   = P(0.0008783, "g/km",
                "MOVES5.0.1 run, electricity, processID 9 (Brakewear), 2025.",
                "This is a 72.3% reduction vs the gasoline brake rate - close to the
                 71% implied by MOVES5 Table 2.12 (0.115/0.400 g/hr) and well above
                 the 50% the manuscript assumed."),
    tire    = P(0.0010300, "g/km",
                "MOVES5.0.1 run, electricity, processID 10 (Tirewear), 2025.",
                "IMPORTANT CORRECTION to earlier guidance: MOVES5 does NOT raise tire
                 wear for BEVs. The electric rate (0.0010300) is essentially identical
                 to gasoline (0.0010330), a difference of 0.3%. MOVES5 differentiates
                 BRAKE wear by fuel type but not TIRE wear. The EMEP/EEA +7-10% mass
                 uplift is a European inventory convention that MOVES does not apply.
                 Do not claim a MOVES basis for a BEV tire-wear penalty."),
    residual_pm = P(0.0, "g/km",
                    "REMOVED - was a back-solved plug with no source.",
                    "The MOVES run confirms it was never needed: brake + tire alone
                     give 0.001908 g/km for a BEV."),
    co2_bev      = P(141, "g/km",
                     "Argonne/DOE R&D GREET 2024, BEV cradle-to-grave on the current
                      US grid: 227 g CO2e/mile = 141 g/km.",
                     "NOT from the MOVES run. New England grid (EPA eGRID2023, NEWE
                      246.4 vs US 349.7 g CO2e/kWh) would put this near 115-125 g/km."),
    co2_sensor   = P(8, "g/km",
                     "Gawron et al. (2018), Env. Sci. Technol. 52(5):3249-3256,
                      DOI 10.1021/acs.est.7b04576: L4 sensing/computing adds 3-20%.",
                     "CONTESTED - report the 4-28 g/km range, not a point value.
                      8 g/km implies ~1.6 kW continuous draw at 50 km/h, above the
                      1.2 kW ceiling in Sudhakar et al. (2022). Lifecycle CO2 only."),
    co2_sensor_range = P(c(4, 28), "g/km", "Gawron et al. (2018) 3-20% envelope")
  ),

  # ---- HYBRID: third fleet class ---------------------------------------------
  #  CONSTRUCTED, not measured. MOVES represents HEVs inside fuelTypeID 1
  #  (gasoline) and distinguishes them only by engine technology, which was not
  #  broken out in the executed RunSpec. A measured hybrid factor needs a third
  #  MOVES run with engine technology in the output detail.
  #  Construction, from the measured gasoline and electric rates:
  #    exhaust : 0.65 x gasoline  (HEV fuel-consumption reduction ~35%)
  #    brake   : gasoline less HALF the BEV regenerative credit (HEVs regen but
  #              with smaller batteries and more friction-brake reliance)
  #    tire    : same as gasoline (comparable mass)
  hybrid = list(
    exhaust = P(0.00089, "g/km",
                "CONSTRUCTED: 0.65 x measured gasoline exhaust (0.0013719).",
                "The 35% reduction is a fuel-consumption proxy, not a measured PM
                 rate. Sensitivity range 0.55-0.80 of the gasoline rate."),
    brake   = P(0.00202, "g/km",
                "CONSTRUCTED: gasoline (0.0031654) less 50% of the measured BEV
                 regenerative credit (0.0031654 - 0.0008783 = 0.0022871).",
                "The 50% regen share is an assumption. Sensitivity 25-75%."),
    tire    = P(0.0010330, "g/km",
                "Same as measured gasoline tire wear - comparable vehicle mass."),
    note    = P(NA, "-", "Total constructed hybrid PM2.5 = 0.00394 g/km",
                "Sits between gasoline 0.00557 and electric 0.00191, as expected.
                 REPLACE WITH A MEASURED RATE if a third MOVES run is done.")
  ),

  # ---- time-varying ICE factor: NOW POPULATED FROM THE 2050 RUN -------------
  # The ICE fleet cleans itself up substantially over the horizon. Holding the
  # baseline fixed would credit AV adoption with reductions that ordinary fleet
  # turnover delivers anyway.
  ice_total_2025 = P(0.005570, "g/km", "MOVES5.0.1, gasoline, 2025, January, fleet-weighted"),
  ice_total_2050 = P(0.004221, "g/km", "MOVES5.0.1, gasoline, 2050, January, fleet-weighted",
                     "24% lower than 2025. Running exhaust falls 87%; by 2050 it is
                      only 4% of total light-duty PM2.5 and the burden is almost
                      entirely non-exhaust."),
  av_total_2025  = P(0.001908, "g/km", "MOVES5.0.1, electricity, 2025, January"),
  av_total_2050  = P(0.001920, "g/km", "MOVES5.0.1, electricity, 2050, January",
                     "Essentially flat - BEV PM2.5 is all non-exhaust, which does not
                      improve with emissions standards."),
  ice_factor_interp = P(TRUE, "-",
                        "Modelling choice - linear interpolation 2025 to 2050",
                        "Now TRUE. The maximum achievable PM2.5 reduction from full
                         fleet electrification FALLS from 65.7% in 2025 to 54.5% in
                         2050, because the ICE fleet improves while the BEV floor
                         (non-exhaust) does not. A frozen baseline overstates
                         late-horizon AV benefits."),

  include_resuspension = P(FALSE, "-",
                           "MOVES treats road dust separately; not selected in the run.",
                           "The EU inventory (EMEP/EEA, Beddows & Harrison) reports
                            substantially higher total non-exhaust PM2.5. Report the
                            framework difference as a sensitivity.")
)

# =============================================================================
# 5. VMT REBOUND
# =============================================================================
#  This drives the paper's headline counter-intuitive result and currently
#  rests on a single ride-hailing study. Treated as a first-class uncertain
#  parameter with an explicit sweep, not a fixed assumption.
# =============================================================================
REBOUND <- list(
  # Central band retained, but it is now a DEFENDED RANGE with the disagreement
  # stated, rather than an extrapolation from a single ride-hailing news item.
  low      = P(0.30, "-",
               "Sun, R., Circella, G., Jaller, M., Qian, X. & Alemi, F. (2023),
                'Impacts of Connected and Automated Vehicles on Travel Demand and
                Emissions in California', Transportation Research Record 2678(4),
                DOI 10.1177/03611981231186984. Private CAV at 100% penetration:
                1,196-1,616M VMT vs 985M baseline = +21% to +64%.
                https://journals.sagepub.com/doi/10.1177/03611981231186984",
               "Statewide activity-based travel demand model + EMFAC, 7 scenarios to
                2050. Best single anchor: a full regional simulation, not an experiment."),
  high     = P(0.40, "-",
               "Same source (upper part of the +21% to +64% envelope). Empirical upper
                anchor: Harb, M., Malik, J., Circella, G. & Walker, J.L. (2022),
                'Simulating Life with Personally-Owned Autonomous Vehicles through a
                Naturalistic Experiment with Personal Drivers', UC ITS UC-ITS-2018-09,
                DOI 10.7922/G2WH2N96 - 43 Sacramento households given a chauffeur:
                +60% VMT, over half from zero-occupancy trips.
                https://rosap.ntl.bts.gov/view/dot/65477",
               "The assumed band sits BELOW the only revealed-behaviour evidence, so it
                can be framed as conservative."),
  sweep    = P(seq(0.0, 0.6, by = 0.05), "-",
               "Sensitivity sweep spanning the low-case literature through the
                empirical upper anchor."),

  # Mechanism, for deriving rather than asserting the value
  price_elasticity = P(-0.4, "-",
    "Taiebat, M., Stolper, S. & Xu, M. (2019), 'Forecasting the Impact of Connected and
     Automated Vehicles on Energy Use: A Microeconomic Study of Induced Travel and Energy
     Rebound', Applied Energy 247:297-308. https://arxiv.org/abs/1902.00382
     Combined price elasticity of VMT demand = -0.4; projects +2% to +47% travel demand
     under full CAV adoption; households more sensitive to TIME cost than fuel cost.",
    "This is the bridge from 'AVs cut the time price of travel' to a VMT elasticity.
     Its +47% upper bound brackets the assumed high case."),

  # Decomposition components, so the band can be built bottom-up
  empty_repositioning = P(0.08, "-",
    "Fagnant, D.J. & Kockelman, K.M. (2015), 'Operations of a Shared Autonomous Vehicle
     Fleet for the Austin, Texas Market', Transportation Research Record 2536:98-106,
     DOI 10.3141/2536-12: +8% VMT from unoccupied repositioning alone.
     http://caee.webhost.utexas.edu/prof/Kockelman/public_html/TRB15SAVsinAustin.pdf"),

  # Induced-demand prior from the closely analogous capacity-expansion case
  capacity_elasticity = P(1.03, "-",
    "Duranton, G. & Turner, M.A. (2011), 'The Fundamental Law of Road Congestion:
     Evidence from US Cities', American Economic Review 101(6):2616-2652,
     DOI 10.1257/aer.101.6.2616. Elasticity of VKT to interstate lane-km = 1.03 (IV).
     https://matthewturner.org/papers/published/Duranton_Turner_AER_2011.pdf
     Corroborated: Hsu & Zhang (2014) J. Urban Econ. 81:65-76 find 1.24-1.34 in Japan;
     Hymel (2019) Transport Policy 76:57-66 finds unit elasticity; Volker & Handy (2021,
     UC-ITS-2021-04, DOI 10.7922/G22F7KQH) adopt 1.0 interstates / 0.75 arterials for
     California regulatory use."),

  # THE CHALLENGE - must be confronted in the text, not omitted
  contrary_evidence = P(0.0595, "-",
    "Naz, F. & Mattingly, S.P. (2024), 'Assessing Automated Vehicle Induced VMT: Meta
     Analysis of Current Research', SSRN DOI 10.2139/ssrn.5030045 (journal version in
     Case Studies on Transport Policy). Pooled mean +5.95% VMT across 26 articles and
     195 effect sizes. https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5030045",
    "THIS IS THE MOST SERIOUS THREAT TO THE ASSUMPTION - the pooled literature central
     estimate is roughly one-sixth of the assumed low case. The defence is that the pool
     is dominated by SHARED-AV and low-penetration simulations with conservative
     value-of-time reductions, whereas this paper models private ownership at high
     penetration. That defence must appear in the text. Also confront:
     Childress et al. (2015) TRR 2493:99-106, Puget Sound ABM, only +3-5% VMT under a
     -35% VOT assumption; and ITF/OECD (2016) Lisbon simulation, -23% VKT over the full
     day under FULL POOLING - which shows the SIGN FLIPS with fleet structure, so the
     assumed band is conditional on private ownership, not a property of automation.
     https://www.itf-oecd.org/sites/default/files/docs/shared-mobility-liveable-cities.pdf"),

  # Ride-hailing evidence, properly cited
  tnc_induced = P(1.57, "-",
    "Schaller, B. (2021), 'Can sharing a ride make for less traffic? Evidence from Uber
     and Lyft and implications for cities', Transport Policy 102:1-10,
     DOI 10.1016/j.tranpol.2020.12.015: VMT increase where personal auto trips shift to
     ride-hail - BOSTON +157%, SF +134%, NYC +114%, Chicago +97%. Deadhead 40-48%.
     http://www.schallerconsult.com/rideservices/sharingride.pdf",
    "REPLACES the Brasuell (2021) Planetizen news item - this is the peer-reviewed
     original behind the '97-118%' figures the paper was citing second-hand, and it
     has a Boston-specific number. Corroborating: Henao & Marshall (2019) Transportation
     46(6):2173-2194 (DOI 10.1007/s11116-018-9923-2), Denver, +83.5% VMT, deadheading
     >=40.8%; Erhardt et al. (2019) Science Advances 5(5):eaau2670
     (DOI 10.1126/sciadv.aau2670), San Francisco delay +62% vs +22% counterfactual.
     CAUTION: TNC evidence bounds the mechanism but is NOT directly transferable to
     private AV ownership - say so.")
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

  # Time-varying factors from the two MOVES years. ef_ice_t(t) and ef_av_t(t)
  # linearly interpolate 2025 -> 2050 and hold flat beyond.
  horizon <- v(SSM$horizon)
  frac    <- pmin(seq(0, horizon) / 25, 1)     # 2025 + 25 yr = 2050
  ef_ice_t <- v(EF$ice_total_2025) + frac * (v(EF$ice_total_2050) - v(EF$ice_total_2025))
  ef_av_t  <- v(EF$av_total_2025)  + frac * (v(EF$av_total_2050)  - v(EF$av_total_2025))
  if (!v(EF$ice_factor_interp)) { ef_ice_t[] <- v(EF$ice_total_2025); ef_av_t[] <- v(EF$av_total_2025) }

  list(
    S_BOS  = S_BOS,
    ef_ice = ef_ice,
    ef_av  = ef_av,
    rho_c  = v(CTM$q_max) / v(CTM$v_f),
    # Algebraically tied to K so the two can never drift apart in the paper.
    diurnal_factor = 1 / (24 * v(CDM$peak_hour_share)),
    ef_reduction_potential = 1 - ef_av / ef_ice,
    ef_ice_t = ef_ice_t,
    ef_av_t  = ef_av_t,
    ef_reduction_potential_t = 1 - ef_av_t / ef_ice_t
  )
}

# Flat accessor: par("SSM","delta_1") -> 0.045
par_val <- function(group, name, sub = NULL) {
  g <- get(group)
  x <- if (is.null(sub)) g[[name]] else g[[name]][[sub]]
  x$value
}

# Parameters that still need work before submission. Matches explicit ALL-CAPS
# status tokens only, so ordinary prose ("replacement-flow", "verified") does not
# trip the flag.
TOKENS <- c("UNSOURCED", "STILL UNSOURCED", "PLACEHOLDER", "REFRESH NEEDED",
            "MUST BE REPLACED", "VERIFY BEFORE", "CONTESTED", "SOURCE IS WEAK")

unsourced_report <- function() {
  flags <- list()
  walk <- function(lst, path = "") {
    for (nm in names(lst)) {
      item <- lst[[nm]]
      pth <- if (path == "") nm else paste0(path, "$", nm)
      if (is.list(item) && !is.null(item$source)) {
        txt <- paste(item$source, item$note)
        hit <- TOKENS[vapply(TOKENS, function(t) grepl(t, txt, fixed = TRUE), logical(1))]
        if (length(hit)) flags[[pth]] <<- paste(hit, collapse = "; ")
      } else if (is.list(item)) walk(item, pth)
    }
  }
  for (g in c("SSM", "CTM", "CDM", "EF", "REBOUND", "COST")) walk(get(g), g)
  flags
}

print_provenance <- function() {
  f <- unsourced_report()
  if (!length(f)) { cat("All parameters sourced.\n"); return(invisible(NULL)) }
  cat("Parameters still needing source work (", length(f), "):\n", sep = "")
  for (nm in names(f)) cat(sprintf("  %-28s %s\n", nm, f[[nm]]))
  invisible(f)
}
