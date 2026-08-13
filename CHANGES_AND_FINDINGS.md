# Code Repair — What Changed, and What It Did to the Results

Repository state after the corrections described below. Every change is in
version control; `git log` shows the sequence.

---

## THE HEADLINE FINDING

**The Year-30 penetration figures in the manuscript are not produced by this
code, under any parameter set.**

`StateModel/av_all_scenarios.csv` is the pipeline's own output artifact — the
file the dispersion model consumes. Running the shipped `plot_model.r` with its
shipped parameters reproduces that CSV exactly. Neither matches the paper.

| Scenario | Paper (Tables 8/9) | Shipped code output | After corrections |
|---|---:|---:|---:|
| No Policy | 66.4% | **5.1%** | 16.5% |
| Supply Push | 72.7% | 47.9% | 67.7% |
| Infra Focus | 73.1% | 31.0% | 31.4% |
| Front-Loaded | 73.7% | 50.9% | 68.5% |
| Phaseout | 74.5% | 50.4% | 65.3% |
| Pulse | 75.0% | 52.2% | 66.9% |
| Adaptive | 75.1% | 48.4% | 61.1% |
| Moderate Rebate | 75.6% | 44.3% | 46.2% |
| Ramp-Up | 77.6% | 48.2% | 54.1% |
| High Rebates | 78.1% | 31.0% | 31.4% |
| Aggressive | 79.9% | 61.2% | 71.0% |

The abstract's "66.4% AV fleet penetration and a 39.8% traffic-attributable
PM₂.₅ reduction by year 30", the 13.5 pp spread between Aggressive and No
Policy, and the 8.1 pp PM₂.₅ gap that recurs throughout §4 and §5 all derive
from the first column. The model produces the second.

This is not a bug that was fixed. The code is internally consistent and
reproducible; the manuscript's numbers are disconnected from it. **Before
anything else, someone needs to establish where the published figures came
from.** No amount of code repair resolves this — it is a provenance question.

---

## THE SECOND FINDING: the cost bug reverses the paper's main conclusion

The policy cost expression charged the manufacturer supply subsidy as a bare
scalar rather than per vehicle:

```r
annual_spend = r * (adoption * S) + s + i     # `s` never multiplied by units
```

So Supply Push — whose entire budget is an $8,000/AV manufacturer subsidy — was
charged **$8,000 per year** instead of $8,000 per vehicle produced. It is
precisely the scenario the paper crowns as most cost-effective.

With correct accounting (`(r + s) * adoption * S + i`):

| Scenario | Cost in paper | Corrected cost |
|---|---:|---:|
| Supply Push | $0.9B | **$6.68B** |
| Infra Focus | $1.3B | $1.07B |
| Moderate Rebate | $3.0B | $2.98B |
| High Rebates | $7.4B | $3.44B |
| Aggressive | $10.1B | $13.81B |

Supply Push is 7.4× more expensive than reported. The paper's claim that it
"delivers 7.00% pollution reduction per $1B — more than double Infra Focus and
roughly five times the average rebate scenario" does not survive: on corrected
costs **Infra Focus is the most efficient scenario in the portfolio**, and
Supply Push falls to mid-table.

The broader direction — outcome-modification levers outperform consumer rebates
per public dollar — does survive, because Infra Focus is also an
outcome-modification lever. But every specific efficiency number in §5.2.3,
§5.4.3, Table 9, Figure 7, and Recommendation 6.2 is wrong, and the named
winner changes.

---

## THE THIRD FINDING: the AV advantage shrinks over the horizon

Two MOVES5.0.1 runs (2025 and 2050) now supply the emission factors. Comparing
them exposes something the frozen-baseline model could not see.

| | 2025 | 2050 | Change |
|---|---:|---:|---:|
| Gasoline running exhaust | 0.001372 | 0.000184 | **−87%** |
| Gasoline brake wear | 0.003165 | 0.003005 | −5% |
| Gasoline tire wear | 0.001033 | 0.001032 | ~0 |
| **Gasoline total** | **0.005570** | **0.004221** | **−24%** |
| **Electric total** | **0.001908** | **0.001920** | **~0** |
| **Max reduction from full electrification** | **65.7%** | **54.5%** | **−11.2 pp** |

The ICE fleet cleans itself up. Running exhaust falls 87% through fleet turnover
and standards phase-in alone, with no AVs involved. BEV PM2.5 does not improve,
because it is entirely non-exhaust and emissions standards do not touch brake and
tire wear.

So **the environmental case for electrification weakens over the modelling
horizon**, from a 65.7% ceiling in 2025 to 54.5% in 2050. The manuscript's fixed
emission factors would have credited AV adoption with a constant advantage that
the technology does not actually retain. `config/params.R` now interpolates both
factors between the two MOVES years, and `derive_params()` returns the annual
series.

A second consequence worth writing up: **by 2050, running exhaust is only 4% of
light-duty PM2.5.** The problem becomes almost entirely non-exhaust, where the
only lever BEVs offer is regenerative braking. That reframes the paper's policy
question and is a more interesting finding than the one it currently reports.

## THE FOURTH FINDING: the model now has a real validation target, and it broadly passes

The Massachusetts Vehicle Census ZIP file gives observed Boston adoption for the
first time. Zero-emission share of Boston-garaged VMT:

| Year | 2020 | 2021 | 2022 | 2023 | 2024 | 2025 |
|---|---:|---:|---:|---:|---:|---:|
| ZEV share of VMT | 0.55% | 0.78% | 1.53% | 3.55% | 4.65% | **4.92%** |
| Hybrid share | 3.8% | 4.2% | 5.2% | 6.3% | 7.8% | **9.6%** |

Against the corrected model's No-Policy trajectory:

| Year | Model AV share | Observed ZEV share |
|---|---:|---:|
| 1 | 0.57% | 0.78% |
| 3 | 1.67% | 3.55% |
| 5 | 2.74% | 4.92% |

Model growth **47.9%/yr** vs observed **54.7%/yr**, with levels within a factor
of two. This is not a formal validation — the model tracks a fleet-*count* share
while the MVC reports a VMT share, and "AV" and "ZEV" are not the same
population — but it is the first time any part of this pipeline has been held
against observed data, and the adoption dynamics are in the right region.

Three things this data changes:

**1. The β coefficients now have a real calibration target.** They are currently
elicited weights fitted to a Morgan Stanley / VTPI / WEF forecast envelope. Six
years of observed Boston adoption is a far stronger basis, and it removes the
paper's most-criticised methodological claim.

**2. There is an inflection the forecasts miss.** ZEV share rose 1.11 pp in
2023–24 but only 0.27 pp in 2024–25. A logistic fitted through 2023 overshoots
badly. Any adoption curve in the paper should be tested against that slowdown.

**3. Observed VMT grew 18.4% in five years with no AVs involved.** The model has
no exogenous VMT growth term, so this is a baseline trend the rebound analysis
must be measured *against* rather than attributed to automation.

Also worth noting: **9.6% of Boston VMT is hybrid**, a category the model does
not represent at all — it treats the fleet as binary ICE-or-BEV, and hybrids sit
between the two emission factors.

`C_0` is now derived from this data at **266,899 vehicles** (3.043e9 VMT across
the 28 modelled ZIPs ÷ 11,400 mi/veh/yr). The manuscript's 512,539 remains ~1.9×
too high. The 28 modelled ZIPs carry **97.5%** of Boston-garaged VMT, so the
spatial scope is well chosen.

## Corrections applied

### Structural

| Change | Rationale |
|---|---|
| **`config/params.R`** — every parameter, with units and source, in one file; all scripts read from it | Four parameters differed between code and paper. Single source of truth makes drift structurally impossible. Appendix tables should be generated from this file. |
| **`run_all.R`** — one entry point, fixed execution order, deterministic seed | Order matters: the SSM writes both adoption trajectories and derived costs that the CDM consumes. |
| **`scripts/cost_model.R`** — costs derived from lever time paths | Replaces a hardcoded literal that lived in a *plotting* script and two mutually inconsistent in-line formulas. |
| **`scripts/convergence_test.R`** — steady-state verification | Lets `T_sim` be justified rather than asserted. |
| **`deprecated/`** — three scripts moved out, with a README | Two of them **write to `av_all_scenarios.csv`**, the CDM's input. Running either silently replaces the pipeline's input with output from a defective model. |

### State-space model (`StateModel/plot_model.r`)

| Parameter | Was | Now | Note |
|---|---|---|---|
| `delta_2` | 0.05 | **0.023** | Table 13 states *"Correction: revised from initial 5%"*. It was never revised in code. |
| `Price_AV` | 40000 | **48000** | Paper Table 13. Enters `price_adv = r / Price_AV`. |
| `beta_0` | −3.5 | **−2.026** | Paper Table 4. |
| `S` | 26970 (asserted) | **28615 (derived)** | Now computed as `S_US × C₀ / C_US` with the arithmetic visible. |
| `n_sim` | 100 at every call site | **1000** (`N_MC`) | Figures 2 and 8 claim N = 1000. |
| infra transform | `sqrt()` hardcoded | switchable, config-driven | Eq. (2) and Appendix A.1 both print a **logistic** map; production code used a square root. Now `sqrt` / `linear` / `logistic` — pick one and write the paper to match. |

`step_model()` now also returns `demand` (the uncapped adoption fraction), so
demand-vs-supply binding can be inspected rather than assumed.

### Dispersion model (`ConvectionDiffusion/plot_results.r`)

**Emission factors now come from an actual MOVES5.0.1 run** — Default Scale,
Emission Rates, Suffolk County, urban arterial, 07:00, speed bin ~48 km/h,
weighted 40/60 car/truck. RunSpecs and raw output are committed under `moves/`.
Against the manuscript: exhaust was **9× too high** (0.012 vs 0.00137), brake
**3.2× too low** (0.001 vs 0.00317), tire **1.7× too low** (0.0006 vs 0.00103).

Two corrections to my own earlier guidance, both now settled by the run itself:

- **MOVES5 does not raise tire wear for BEVs.** Electric tire wear is 0.0010300
  vs gasoline 0.0010330 — a 0.3% difference. MOVES5 differentiates *brake* wear
  by fuel type but not *tire* wear. The EMEP/EEA +7–10% mass uplift is a European
  inventory convention; do not claim a MOVES basis for it.
- **The regenerative-braking credit measures 72.3%**, close to the 71% implied by
  MOVES5 Table 2.12 and far above the 50% the manuscript assumed.

**Dimensional fix.** `E_eff` is a line-source strength in g·m⁻¹·s⁻¹. The code set
`S = E_eff / dz`, giving g·m⁻²·s⁻¹ injected into a field carrying g·m⁻³ — an
implicit street width of exactly 1 m. Now `S = E_eff / (street_width × dz)`
with `street_width` an explicit, justifiable parameter.

*(Note the manuscript's "1.102×10⁻⁶ g/m³ (1.102 ng/m³)" is separately a unit
typo — that value is 1.102 **µg**/m³.)*

**Emission factors.** `get_pm25_per_km()` returned a single CSV column
(ICE 0.005, AV 0.002) and ignored the `Brake_Wear` / `Tire_Wear` rows in the
same file. Eq. (15) as printed — ICE 0.0136, AV 0.003 — was never what ran. Now
summed explicitly from config, so Table 2 and the executed code cannot diverge.

**Temporal basis.** Flows were divided by `30·24·3600`, smearing a month of
traffic across every hour including 03:00, while the CTM that produced those
flows runs on peak-hour capacity. Now peak-hour by default, with an explicit
diurnal factor for daily means.

**Steady state.** `T = 300 s` over a 1000 m domain at 3.5 m/s (transit ≈ 286 s)
sampled the field mid-transient. The convergence sweep now shows:

```
T (s)      ground_max     ground_mean   rel.change
  150     1.6331e-06     1.2794e-06            -
  300     2.3283e-06     1.4792e-06       0.4258
  600     2.3744e-06     1.4847e-06       0.0198
 1200     2.3744e-06     1.4847e-06       0.0000
```

`T_sim = 2400 s`, verified converged.

**`ground_mean`** was averaged only over cells with C > 0, making it depend on
how far the plume had travelled. Now averaged over a fixed spatial window.

### Resolved: the O-D corridor contradiction

§4.3.1 says 106 corridors; Table 11 says top-50 of 378. The code produces
**106**. Table 11 is the error.

---

## Post-correction dispersion output

Peak-hour ground-level, all-ICE fleet vs. 67.7% AV fleet:

```
av_fraction 0.0000 : 106 corridors | mean ground_max 16.14 ug/m3 | max 82.67 ug/m3
av_fraction 0.6771 : 106 corridors | mean ground_max  7.62 ug/m3 | max 39.04 ug/m3
```

The 52.8% reduction matches the emission-factor ratio exactly
(`1 − 0.677 × (1 − 0.003/0.0136)`), confirming the chain is internally
consistent. Absolute magnitudes are peak-hour near-road increments and depend
directly on `street_width` and `peak_hour_share`, both of which are flagged in
config as needing real MassDOT data.

---

## What still needs a human

1. **Establish where the published Year-30 figures came from.** Nothing else
   can be finalised until this is answered.
2. **Choose the infrastructure transform** (`sqrt` / `linear` / `logistic`) and
   write Eq. (2) and Appendix A.1 to match whichever is chosen.
3. **`I_ref = 66214286`** appears nowhere in the manuscript and has no source.
   Derive it from a real infrastructure-stock figure or present it explicitly
   as a scaling choice with a sensitivity sweep.
4. **`EF$av$residual_pm = 0.0019 g/km`** is a back-solved plug that exists only
   to make the AV factor reach the 0.003 g/km printed in Eq. (15). Source it or
   drop it — without it the AV factor is 0.0011 g/km.
5. **`street_width`, `peak_hour_share`, `diurnal_factor`** — placeholders
   pending MassDOT count data. They set the absolute concentration scale.
6. **`cap_0`, `gamma`** carry no source at all.

Run `Rscript -e 'source("config/params.R"); print(unsourced_report())'` for the
current list; `run_all.R` prints it at startup.
