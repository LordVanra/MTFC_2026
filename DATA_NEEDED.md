# What Still Needs Gathering

Ordered by what unblocks the paper, not by effort.

---

## TIER 1 — blocks publication

### 1. Where did the published Year-30 figures come from?

Not a dataset — a question for whoever ran the original model. It remains the
single largest open item. The manuscript reports No Policy at 66.4%; the code
produces 5.1%, and reproduces its own output CSV exactly. Until someone can say
where 66.4% came from, it is unclear whether this is a correction pass or a new
set of results, and that determines how the paper is written.

### 2. Massachusetts Vehicle Census — ZIP-code file

**One download. Fixes four separate problems.**

<https://gis.data.mass.gov/datasets/dce8f7ce88f449f08150af4183b9a983>
→ Download → CSV (~4.3 MB). I cannot retrieve it — the ArcGIS content endpoint
returns binary to my fetcher and the hub download API rejects the request.

Fields: `ZipCode`, `Municipality`, `Year`, `FuelClass`, `VehicleType`,
`AdvancedVehicleType`, `ModelYear`, `GVWRCategory`, **`AnnualVMT`**.

What it resolves:

| Currently | With this file |
|---|---|
| `C_0` = 237,224 from **2014** | Current-year Boston count |
| Car/truck VMT split **assumed 40/60** — this weights the MOVES emission factors | Actual split from `VehicleType` |
| Baseline ICE/BEV mix **assumed** | Actual from `FuelClass` |
| Vehicle stock is a single city-wide number | Per-ZIP distribution across your 28 modelled ZIPs |
| **Nothing in the pipeline is validated against observed data** | `AnnualVMT` per ZIP is a ground-truth check on the CTM |

That last row is the one that matters. It is the difference between a model that
is internally consistent and a model that has been tested against reality.

---

## TIER 2 — converts the paper's biggest weakness into a strength

### 3. Observed PM₂.₅ from Boston monitors

The dispersion model has never been compared to a measurement. Even a rough
order-of-magnitude and spatial-pattern check would materially change how a
referee reads the paper.

- **EPA AQS** pre-generated data files: <https://aqs.epa.gov/aqsweb/airdata/download_files.html> — hourly PM2.5, by year, filter to Massachusetts / Suffolk County
- **MassDEP** air quality data: current and historical Boston monitors

You are modelling the *traffic-attributable increment*, not total ambient PM₂.₅,
so do not expect the levels to match. What you can test: does the modelled
spatial ordering of corridors match the ordering across monitors, and is the
modelled increment a plausible fraction of measured total?

### 4. More MassDOT traffic counts

`peak_hour_share` currently rests on **one arterial, one direction, one weekday
in 2016**. That is thin for a parameter that scales every concentration.

MassDOT publishes counts through its Transportation Data Management System.
Pull 5–10 Boston-area arterial count stations and take a distribution rather
than a point value. The same data cross-checks the CTM's corridor volumes.

---

## TIER 3 — improves defensibility, not blocking

### 5. Observed corridor speeds
For the speed-coupling fix, the MOVES rate table is now indexed by speed bin
(`moves/output/ef_by_speed.csv`), so what remains is mapping CTM speeds onto
those bins. That is internal work. Real observed speeds would additionally let
you check the CTM's speeds are plausible.

### 6. BEA Fixed Assets Table 7.1B, current year
<https://www.bea.gov/itable/fixed-assets> — Structures → Highways and streets.
Scaled to Boston, this lets `I_ref` be *derived* rather than declared. Optional
but removes an unsourced parameter.

### 7. A sample of Boston arterial widths
`street_width` = 16.5 m is built from Boston Complete Streets *minima*. Measuring
a handful of the actual modelled corridors would replace a construction with an
observation. It scales concentration linearly, so it is worth an hour.

---

## NOT FINDABLE — stop looking, present as assumptions

`cap_0`, `gamma`, `A_0`, and the AV sensor overhead have no external source and
searching harder will not produce one. They need the assumptions-table treatment
described in `TODO_OPEN_PARAMETERS.md`: state the value, state that it is a
modelling choice, give the tested range, report the effect on the result.

`gamma` deserves its own sensitivity figure rather than a table row — it is the
only channel through which Supply Push acts, so it is load-bearing for the
central comparison.

---

## Suggested split

You chase Tier 1 item 2 (one download, biggest payoff). In parallel I can:

- build the assumptions table and sensitivity sweeps for the five unsourceable
  parameters
- wire the speed-binned MOVES table into the CTM so emission factors vary with
  modelled corridor speed
- regenerate all figures from the corrected pipeline
- draft the methods and limitations sections around what the model now does

None of that needs new data.
