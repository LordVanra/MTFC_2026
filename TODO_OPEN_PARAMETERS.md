# Instructions for the Seven Open Parameters

Two require you to run software or use a browser. Five cannot be sourced at all
and need to be *presented* correctly instead. Do them in this order — Task 1 is
the one a reviewer is most likely to challenge.

---

# TASK 1 — Run MOVES and replace the exhaust emission factor

**Time: half a day, mostly waiting.** This is the single most important open item.

## First, an important change: it is MOVES5 now, not MOVES4

EPA released **MOVES5.0.1 in April 2026**, and the Federal Register notice makes
MOVES5 **required** for SIP and transportation-conformity purposes in states
other than California. MOVES4 is superseded.

This affects your manuscript beyond the one number. The paper cites "EPA MOVES4"
throughout — Table 2, Table 10, §3.3.2, Eq. (15). Using a superseded version is
a legitimate reviewer objection, and it also *solves a problem for you*:

> MOVES4 applies identical brake and tire wear rates to every fuel type. That
> means MOVES4 cannot support your regenerative-braking credit **at all**.
> MOVES5 differentiates by fuel type (gasoline car 0.400 g/hr vs electric
> 0.115 g/hr). Moving to MOVES5 turns an unsupportable assumption into a
> sourced one.

So: run MOVES5, and change every "MOVES4" reference in the paper to MOVES5.

## Steps

**1. Install.** Download and run:
`https://www.epa.gov/system/files/other-files/2026-04/moves5.0.1-setup-20260402.exe`

It is Windows-only, which suits you. The installer bundles MariaDB — let it
install its own database server rather than pointing it at an existing one. Give
it ~20 GB of disk. Do not install it under a path with spaces or OneDrive sync.

**2. Choose the scale — and you probably do not want County Scale.**

The Scale panel says it outright: *"This scale requires user-supplied local data
for most activity and fleet inputs."* That is the whole reason the County Data
Manager exists and why so many of its tabs show red X marks. County Scale is
built for SIP and conformity submissions, where a regulator must supply local
fleet counts, local VMT, and a local speed distribution.

**You are not doing a conformity analysis. You need emission *rates*, not a
Suffolk County inventory.** Fighting the CDM to produce a number you will then
divide by distance anyway is the long way round.

### Recommended: Default Scale + Emission Rates

| Panel | Setting |
|---|---|
| Model | Onroad |
| Domain/Scale | **Default Scale** |
| Calculation Type | **Emission Rates** (not Inventory) |

This combination:

- **requires no user-supplied local data** — no County Data Manager at all, and
  none of the red X tabs
- returns **g per unit distance directly**, indexed by average-speed bin, which
  is exactly the lookup table you need for the speed-coupling fix. Instead of one
  emission factor applied to every corridor, you get a rate that varies with
  speed — so AV-induced flow smoothing actually shows up in the emissions rather
  than being asserted
- still lets you select Massachusetts → Suffolk County in Geographic Bounds; it
  simply uses EPA's default allocation factors instead of local data

The trade-off is real and must be disclosed: the fleet is a nationally-derived
default allocated to Suffolk County, not a locally-surveyed Boston fleet. For a
scenario-comparison paper that is defensible, because your result depends on the
*ratio* between ICE and BEV rates far more than on the absolute level. Write in
the methods: *"Emission rates were generated with MOVES5.0.1 at Default Scale
with EPA default allocation to Suffolk County, MA. Local fleet data were not
substituted; results are therefore representative rather than
regulatory-grade."*

Note the panel's own warning — Default Scale is unsuitable for SIP or conformity
work. That restriction does not bind you, but say so explicitly rather than
letting a referee wonder whether you noticed.

### Upgrade path

If you later obtain Massachusetts fleet data — MassDEP maintains it, and the
Vehicle Census from Task 2 gives you fleet composition — switch to County Scale
and localise Age Distribution, Source Type Population, and Vehicle Type VMT.
That is a genuine improvement, not a prerequisite. Do not let it block the run.

### The rest of the RunSpec

| Panel | Setting |
|---|---|
| Time Spans | Analysis year; **Months** = January **and** July; **Hours** = all 24; Day = weekday |
| Geographic Bounds | Massachusetts → **Suffolk County** (FIPS 25025) |
| Vehicles/Equipment | **Passenger Car** and **Passenger Truck**; fuels **Gasoline** and **Electricity** |
| Road Type | **Urban Restricted Access** and **Urban Unrestricted Access** (the latter is your arterials) |
| Pollutants and Processes | PM2.5 as **three** pollutants — see below — with Running Exhaust, Brakewear, Tirewear, and Start Exhaust |

## Which hour? — select all 24, then filter

**Do not pick a single hour in the RunSpec.** Selecting all 24 hours costs
almost nothing extra in County-scale runtime (MOVES computes the rate table
once and applies it across the selected time spans), and it buys you four
things you otherwise cannot get without re-running:

1. **The diurnal factor becomes measured, not assumed.** `CDM$diurnal_factor` is
   currently derived algebraically as 1/(24K) from a single traffic count. With
   all 24 hours you can compute the emissions-weighted 24-hour mean directly.
2. **Hour choice becomes a sensitivity, not a hidden assumption.** You can report
   how much the headline result moves between peak-hour and daily-mean framing.
3. **A referee asking for a 24-hour or annual mean** — likely, since the PM₂.₅
   NAAQS are 24-hour and annual — does not force a new run.
4. **You can match the dispersion model's stability classes.** The CDM uses
   Class B (summer day), D (neutral baseline), and F (winter night). Having all
   hours lets each meteorological case be paired with its own emission rate
   instead of reusing one number across all three.

**If you are forced to choose one hour:** the weekday morning peak, **07:00–07:59**.
That matches the Mount Auburn count that sets `CDM$peak_hour_share = 0.0724`
(peak at 07:00), and it matches the CTM, which runs on peak-hour capacity. In
the MOVES GUI the Time Spans panel lists hours with explicit clock labels, so
select the one labelled 7:00–7:59 AM rather than working from an hour ID.

## Which month? — January and July, and January if only one

Month matters more than hour for PM₂.₅, through two channels: exhaust PM rises
at low ambient temperature, and MOVES applies **seasonal fuel formulations**
(summer vs winter RVP) that change exhaust composition. A January rate and a
July rate differ more than 7 AM differs from 8 AM.

**Run January and July.** They bracket both channels, and they map onto the
meteorological cases the dispersion model already carries:

| CDM case | Stability class | Mixing height | MOVES month |
|---|---|---|---|
| Winter night | F | 200 m | **January** |
| Winter day | D (baseline) | 800 m | **January** |
| Summer day | B | 1,600 m | **July** |

**If only one: January.** The manuscript's headline peak concentration (~7.5
µg/m³ on the 02127→02210 corridor) is explicitly reported "under winter
stability conditions," and Class F winter night is your worst-case dispersion.
Pairing a worst-case dispersion regime with a summer emission rate would be
internally inconsistent — and inconsistent in the direction that flatters the
result, which is the kind of thing referees look for.

**Rule: the month in MOVES must match the meteorological case in the CDM.**
Never apply a July emission factor to a Class F winter-night run.

### A bigger issue the month question exposes: which *year*?

While setting Time Spans you also pick the analysis year, and this deserves more
thought than the month does.

Your emission factors are currently **constant across all 30 years**. Real ICE
fleet emission factors are not. As older vehicles retire and Tier 3 / LEV III
phase in, the average gasoline vehicle's PM₂.₅ rate falls on its own, with no AVs
involved. MOVES models exactly this: run it for 2026 and again for 2050 and the
ICE fleet rate will differ materially.

This matters for your central claim. If the ICE baseline is frozen at today's
rate while the AV share grows, **the model attributes to AV adoption a reduction
that would have happened anyway through ordinary fleet turnover.** The No-Policy
scenario should already show declining traffic-attributable PM₂.₅ from fleet
modernisation alone — and every scenario's benefit is measured against that
baseline, so freezing it inflates all of them.

**Recommendation:** run MOVES for your base year and at least one future year
(e.g. 2026 and 2050), then interpolate the ICE factor across the horizon rather
than holding it fixed. That is a modest amount of extra work and it removes a
serious confound. If you choose not to, it must appear in the limitations,
stated in the direction it biases — because it biases toward your conclusion.

### The trap to avoid: speed coupling

This is the part most likely to draw a methodological objection.

MOVES running-exhaust rates are a function of **average speed** — emissions per
kilometre rise substantially in stop-and-go conditions. MOVES will use EPA's
default `avgSpeedDistribution` for Suffolk County unless you replace it.

But your CTM *computes its own speeds*, and those speeds change across scenarios
as AV penetration alters the fundamental diagram (`alpha_v`, `alpha_rho`). If you
take a MOVES factor generated at EPA's default speed distribution and apply it
to CTM-derived flows, you have emission rates from one speed regime multiplied by
traffic from another — and, worse, the mechanism by which AVs are supposed to
reduce emissions through smoother flow never actually enters the emission factor.

Two ways to handle it, in order of preference:

- **Feed CTM speeds into MOVES.** Export your CTM's speed distribution by road
  type, bin it into MOVES' 16 average-speed bins, and load it in the County Data
  Manager's Average Speed Distribution tab. Then the emission factor and the
  traffic model describe the same road. If you run MOVES at Emission Rates scale
  instead of Inventory, you get a speed-indexed lookup table and can apply the
  right rate per cell per scenario — this is the cleanest option and turns a
  weakness into a genuine strength of the paper.
- **At minimum, check and disclose.** Note the average speed MOVES assumed for
  urban unrestricted access, compare it against your CTM's modelled speeds, and
  state the discrepancy in the limitations. If they differ by more than a few
  km/h, say what direction that biases the result.

**3. County Data Manager.**

### Where it actually is

It is **not** a top menu item, and there is **no "Generate" button** — both of
those were wrong in an earlier version of this document.

In MOVES5 the County Data Manager is reached from a panel in the **left
navigation list**, below the output panels, called **"Create Input Database"**.
Click that panel, then the **"Enter/Edit Data"** button on it.

### Why you may not see it yet

The Create Input Database panel stays disabled, showing a **red X**, until
*every other panel in the RunSpec has a green checkmark*. That is almost
certainly what is happening.

Scan the left navigation list for red X marks and finish those panels first:

- Description
- Scale — must be **County** (Default scale has no CDM at all)
- Time Spans
- Geographic Bounds — must resolve to a **single county**
- Onroad Vehicles / Vehicles and Equipment
- Road Type
- Pollutants and Processes
- General Output
- Output Emissions Detail

Two common blockers: Scale left on Default, and a Vehicles/Equipment selection
containing an invalid fuel/source-type combination. Green across the board and
the CDM unlocks.

### The actual per-tab workflow

There is no one-click generate. For each tab:

1. **Export Default Data** → save the Excel file (into `moves/cdm_exports/`)
2. Open it and check the values are sane
3. Back in the CDM, **Browse** to that file
4. Select the worksheet name from the dropdown
5. **Import**
6. Confirm a green checkmark appears with no error message

Accepting EPA defaults still requires that export-then-import round trip — the
data has to physically land in your county input database. The upside is you end
up with a spreadsheet of every input, which is exactly the audit trail the paper
needs. Tabs where no default exists offer **Create Template** instead, giving a
blank correctly-formatted sheet to fill in.

Start on the **Database** tab and create the input database first, e.g.
`bos_2026_jan_in` — distinct from the *output* database `bos_2026_jan_out`.

### The tabs, and which matter for you

| Tab | Priority |
|---|---|
| Database | First — create the input DB |
| RunSpec Summary | Read-only check |
| **Average Speed Distribution** | **Highest.** This is the speed-coupling problem. Defaults here are precisely what decouple your emission factors from your CTM speeds. |
| Age Distribution | High — drives fleet turnover, and it is what makes the 2026-vs-2050 comparison meaningful. MassDEP can supply Massachusetts data. |
| Source Type Population | High — cross-check against the Vehicle Census figure from Task 2 |
| Vehicle Type VMT | High |
| Fuel | Medium — carries the seasonal RVP that makes January differ from July |
| Meteorology | Medium — the temperature channel behind the month choice |
| Road Type Distribution | Medium |
| I/M Programs | Massachusetts runs an active I/M program; the default should be right, but confirm |
| Starts | Needed, since you selected the start-exhaust process |
| Hotelling, Idle | Heavy-duty concepts; defaults are fine for a light-duty study |
| Retrofit Data | Skip |

Export every tab to `moves/cdm_exports/` and commit them. That is what lets you
write "EPA defaults were used for X, Y and Z; A and B were localised" and have
the claim be checkable.

### General Output panel

| Field | Set to | Why |
|---|---|---|
| Server | `localhost` | Blank by default. Type it, then click **Refresh**. |
| Database | **Create Database…** → `bos_2026_jan_out` | You will do up to four runs (Jan/Jul × base/horizon year). Name each one so you can tell them apart later. |
| Mass Units | **Grams** | Already correct. |
| Energy Units | Joules | Irrelevant unless you report energy. Leave it. |
| Distance Units | **Kilometers** | Change this from Miles. Your entire codebase is g/km. Taking kilometres directly removes a unit conversion — and a botched unit conversion is exactly what produced the 0.012 g/km error in the first place. |

**Activity checkboxes — tick these:**

- ☑ **Distance Traveled** — *mandatory*. Without it an Inventory run gives you total grams with no denominator, and you cannot form g/km at all.
- ☑ **Population** — lets you cross-check MOVES' assumed Suffolk County fleet size against the Massachusetts Vehicle Census figure from Task 2. If they disagree badly, your County Data Manager inputs are wrong.
- ☑ **Source Hours Operating** — gives you g/hr, which is the unit EPA's MOVES5 brake-wear table (Table 2.12) is published in. This is how you verify the 71% regenerative-braking figure against your own run rather than taking my word for it.
- ☑ **Starts** — you are excluding start emissions (a line-source road model only carries running emissions). Ticking this quantifies what you left out, so the exclusion can be stated with a number in the limitations instead of just asserted.

Leave Hotelling Hours and Source Hours Parked unticked — both are heavy-duty
and parked-vehicle concepts, irrelevant to a light-duty running-emissions study.

### Output Emissions Detail panel

| Setting | Set to | Why |
|---|---|---|
| Time | **Hour** | Correct as shown. You selected all 24 hours, so you need hourly resolution to recover the diurnal profile. |
| Geographic | **COUNTY** | Correct as shown. |

The greyed-out boxes — Fuel Type, Emission Process, Road Type, Source Use Type —
are forced on at this scale. That is fortunate, because you need every one of
them:

- **Fuel Type** separates gasoline from electricity — this is what gives you the ICE and AV factors from a single run
- **Emission Process** separates Running Exhaust / Brakewear / Tirewear, which is Eq. (15)'s three terms
- **Road Type** separates urban unrestricted (your arterials) from restricted
- **Source Use Type** separates passenger car from passenger truck

**Untick these three** — they are all currently on and none earns its place:

- ☐ **Model Year** — multiplies your output rows by ~40 model years. Fleet-turnover effects are better captured by running two *calendar* years (the year issue below) than by disaggregating model year within one.
- ☐ **SCC** — source classification codes are for regulatory inventory reporting, not analysis. Pure row bloat.
- ☐ **Regulatory Class** — you are not doing a standards-compliance analysis.

Leave Fuel Subtype unticked. Leave all Nonroad boxes unticked — you are modelling
onroad vehicles only.

Unticking those three will cut the output table size by roughly two orders of
magnitude and makes the post-processing query much simpler.

**4. Execute.** Expect 20–60 minutes.

**5. Extract the number.** In the output database, the table you want is
`rateperdistance` (if you ran Emission Rates scale) or `movesoutput` joined to
`activitytype` for distance (if Inventory). You want grams per vehicle-kilometre
by `processID`:

- `processID = 1` → Running Exhaust
- `processID = 9` → Brakewear
- `processID = 10` → Tirewear
- `pollutantID = 110` → PM2.5 (total, elemental + organic)

Sum the three for your ICE factor. The electric-fuel rows give your AV factor
directly, which is better than my derivation.

**6. Put it in the code.** Open `config/params.R` and replace:

```r
EF$ice$exhaust    # currently 0.0015 g/km  <- placeholder, replace
EF$ice$brake      # currently 0.00172 g/km <- confirm against your run
EF$ice$tire       # currently 0.00080 g/km <- confirm against your run
EF$av$brake       # currently 0.000499     <- use the electric rows
EF$av$tire        # currently 0.000868     <- use the electric rows
```

In the `source` field, write: *"MOVES5.0.1, County Scale run, Suffolk County MA,
[year], [month], [hour], urban unrestricted access; rateperdistance output.
RunSpec archived at `moves/[filename].mrs`."*

**7. Archive the RunSpec.** Commit the `.mrs` file and your exported CountyDataManager
tables to `moves/` in the repo. This is what makes the number reproducible, and
it is what a referee will ask for.

## What to do if you cannot run MOVES

Cite a published MOVES-based study for a comparable urban fleet, and say plainly
in the limitations that emission factors are taken from the literature rather
than a project-specific run. That is honest and survives review. What does not
survive is citing "EPA MOVES4" for a number that is not a MOVES output.

---

# TASK 2 — Get the current Boston vehicle count (and a validation dataset)

**Time: 15 minutes.** And it turns out to be worth much more than the one number.

## The better dataset

While checking, I found the MVC is published at **ZIP-code level**, not just
municipality:

**MVC Annual Zip Code** — https://gis.data.mass.gov/datasets/dce8f7ce88f449f08150af4183b9a983

Fields: `ZipCode`, `Municipality`, `MPO`, `Year`, `FuelClass`,
`AdvancedVehicleType`, `VehicleType`, `VehicleUse`, `ModelYear`,
`GVWRCategory`, **`AnnualVMT`**, `GarageCode`.

Your model is built on 28 Boston ZIP codes. This dataset is at exactly that
resolution, and it carries **AnnualVMT per ZIP**.

That second point is the important one. Your paper's biggest structural weakness
is that **nothing in the pipeline is validated against observed data**. This
dataset gives you:

- `C_0` per ZIP — replacing a single city-wide number with a real spatial distribution
- fleet composition by `FuelClass` — a real baseline for the ICE/BEV mix instead of an assumption
- `AnnualVMT` per ZIP — **a ground-truth check on your cell-transmission model's traffic distribution**

A model-vs-observed comparison on ZIP-level VMT would close the single biggest
gap in the paper. I would treat this as a headline opportunity, not a chore.

## Steps

1. Open https://gis.data.mass.gov/datasets/dce8f7ce88f449f08150af4183b9a983 in Chrome.
2. Click **Download** → **CSV** (~4.3 MB). I could not retrieve it programmatically — the ArcGIS content endpoint returns binary through my fetcher.
3. Filter to `Municipality = "Boston"`, take the most recent `Year`.
4. Sum vehicle counts by `ZipCode` for your 28 modelled ZIPs; sum across all ZIPs for the city total.
5. Put the city total in `config/params.R` as `SSM$C_0`, with source *"Massachusetts Vehicle Census, MVC Annual Zip Code, [year], MassDOT/MAPC, accessed [date]"*.
6. Save the CSV to `data/` and commit it.

Also grab the municipality-level file if you want a cross-check:
https://gis.data.mass.gov/datasets/0e7873d506d547d88b555545e3627ca8

**Sanity check before you accept the number:** divide by ~270,000 Boston
households. You should get roughly 0.8–1.0 vehicles per household. If you get
1.9, you have the wrong column — that is exactly how the 512,539 error happened.

---

# TASKS 3–7 — The five that cannot be sourced

`I_ref`, `cap_0`, `gamma`, `A_0`, and the AV sensor overhead have no external
source, and searching harder will not produce one. Three of them are internal
calibration constants of your own model; one is a quantity no registry publishes;
one is genuinely contested in the literature.

**This is fine.** Models are allowed to contain assumptions. What is not allowed
is assumptions that look like sourced facts. The fix is presentational, and it
is the same fix in each case:

> State the value. State that it is a modelling choice. Give the range you
> tested. Report how much the result moves across that range.

A parameter handled that way is unobjectionable. The same parameter with a
citation-shaped footnote pointing at nothing is what gets a paper rejected.

## 3. `I_ref` = 66,214,286 — infrastructure normalisation constant

Worst of the five, because it appears **nowhere in the manuscript at all** while
materially affecting adoption through `β₂ · sqrt(I/I_ref)`.

Two options:

- **Derive it.** Take BEA Fixed Assets Table 7.1B (Structures → Highways and streets), scale to Boston by population or lane-miles, and show the arithmetic. Then it becomes a derived quantity with a real anchor.
- **Declare it.** State that infrastructure is normalised to a reference stock chosen so that the baseline index maps to ~0.1 of saturation, that the value is a scaling convention, and that adoption trajectories are reported against it.

Either way it must appear in the paper, and it must be in the sensitivity table.
Sweep it ±50% — if Year-30 share moves more than a couple of points, say so.

## 4. `cap_0` = 0.20 — Year-1 adoption ceiling

A supply constraint: at most 20% of new sales can be AVs in year 1. Present as a
stated assumption. Justify qualitatively from manufacturing ramp-up
constraints — you can cite Litman (2023, VTPI) for the general deployment
timeline without pretending he gives this number. Sweep 0.10–0.40.

## 5. `gamma` = 0.00002 — capacity response to subsidy

Pure calibration constant: how much production capacity expands per dollar per
vehicle of manufacturer subsidy. There is no empirical basis and you should not
invent one. Say: *"γ is a calibration constant governing the responsiveness of
production capacity to supply-side subsidy; it is set so that a $8,000/AV
subsidy relaxes the capacity constraint within roughly N years. Results are
reported across γ ∈ [range]."*

Note this one matters more than it looks — it is the only channel through which
Supply Push does anything, so its sensitivity sweep is load-bearing for your
main comparison. Give it a full figure, not a table row.

## 6. `A_0` = 98 — initial Boston AV fleet

No registry publishes Level-4 AV registrations for Boston, because there
essentially are none. **My recommendation: set it to 0.**

A 30-year logistic that starts from 98 versus 0 will be visually identical after
year 3, and 0 is defensible without any source. If you keep 98, you must explain
where it came from, and I do not think you can.

If you want a non-zero start, use a *stated* seed: "adoption is seeded at 0.01%
of the fleet to avoid a degenerate zero state," which is a numerical convention,
not an empirical claim.

## 7. AV sensor/computing overhead = 8 g CO₂/km — contested

Different from the other four: there *is* literature, it just disagrees.

- Gawron et al. (2018): L4 sensing/computing raises energy and GHG **3–20%** → 4–28 g CO₂e/km
- Sudhakar, Sze & Karaman (2022): a 0.84 kW compute budget → ~4–6 g/km
- Onat et al. (2023): **+8% lifecycle** once rebound and +40% manufacturing are counted
- Gawron also finds net **reductions** are possible once operational benefits are included

Report **4–28 g/km as a range**, with 8 as the central case, and say the
literature disagrees on sign once operational effects are included.

One concrete thing that will impress a referee: **state the implied power draw.**
8 g CO₂/km on the New England grid is ~32 Wh/km, which at 50 km/h is **~1.6 kW
continuous** — above the 1.2 kW ceiling Sudhakar et al. identify. That is worth a
sentence. It makes the assumption auditable and shows you understand what the
number physically means.

---

# After you finish

```bash
Rscript -e 'source("config/params.R"); print_provenance()'
```

Every remaining flag should be one you have consciously decided to keep. Then:

```bash
Rscript run_all.R
```

and check `outputs/policy_costs.csv` against whatever the paper says.

## Build the assumptions table from the code

For each of the five judgement parameters, the paper needs one row: value,
status ("modelling choice"), range tested, and effect on the headline result.
Generate it from `config/params.R` rather than typing it — that is the whole
point of the config file, and it guarantees the appendix cannot drift from what
executed.

## Sequence note

Do Task 2 before Task 1 if you want to be efficient. The MVC fleet-composition
data tells you the actual gasoline/hybrid/BEV mix in Boston, which is an input
you can feed into the MOVES County Data Manager to localise the run rather than
accepting EPA defaults.
