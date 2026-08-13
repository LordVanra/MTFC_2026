# MOVES5 runs

Everything needed to reproduce the emission factors in `config/params.R`.

## Contents

| File | What it is |
|---|---|
| `extract_emission_factors.sql` | Post-processing for **Inventory** runs. Divides emissions by distance to get g/km. |
| `extract_rates_emissionrates.sql` | Post-processing for **Emission Rates** runs. Reads `rateperdistance` directly and produces the speed-binned table. |
| `*.mrs` | RunSpec files — **commit these**. They are what make the numbers reproducible. |
| `cdm_exports/` | County Data Manager tables exported to Excel, so it is on record which inputs were EPA defaults and which were localised. |
| `output/*.csv` | Exported `ef_rates` tables, one per run. |

## Corrected RunSpecs

`BOS_2025.mrs` and `BOS_2050.mrs` are ready to load (File → Open). They fix the
issues found in the first draft RunSpec — see "RunSpec review" below.

**Two runs, not four.** MOVES accepts multiple months in a single run, and
`monthID` is carried in the output, so January and July are separable
afterwards. Only the *year* needs a separate run, because that is what captures
ICE fleet turnover.

## RunSpec review — what was wrong in the first draft

| # | Issue | Severity |
|---|---|---|
| 1 | `scenarioid` contained all four labels concatenated and truncated at 40 chars: `BOS_2026_JAN BOS_2026_JUL BOS_2050_JAN B`. MOVESScenarioID is **one label for one run**. | **Blocking** |
| 2 | `truncateoutput`, `truncateactivity`, `truncatebaserates` all `true`. Each run would have wiped the previous run's output. Combined with issue 1, running all four into one database would have left only the last. | **Blocking** |
| 3 | `outputvmtdata=false` and `outputsho=false` — Distance Traveled and Source Hours Operating were not ticked. Without them there is no denominator for a cross-check and no way to derive the average speed MOVES assumed. | High |
| 4 | `outputshidling=true` — Hotelling Hours is a heavy-duty concept, irrelevant here. | Cosmetic |
| 5 | Year set to **2025** while the output database was named `bos_2026_out`. | Naming |
| 6 | All five road types selected, including Rural Restricted and Rural Unrestricted. Suffolk County has essentially no rural roads; the rows are noise. Off-Network is retained because start and idle emissions live there. | Low |
| 7 | `scaleinputdatabase` still pointed at `bos_2026_in`. Not needed at Default Scale. | Low |

Fixed in the supplied RunSpecs. Two decisions left for you:

**Which fuels count as "ICE"?** The RunSpec selects Gasoline (1), Diesel (2),
Ethanol E-85 (5) and Electricity (9). `config/params.R` currently treats ICE as
gasoline only. Either restrict the ICE factor to fuel type 1 and say so, or
compute a fleet-weighted ICE factor across 1, 2 and 5 — the latter is more
realistic for Boston and is what a referee would expect from a fleet model.

**Start Exhaust (process 2) is not selected.** Without it you cannot quantify
what the running-emissions-only assumption leaves out. Add it if you want that
number for the limitations section; skip it if you are content to state the
exclusion qualitatively.

## Run matrix and MOVESScenarioID

`MOVESScenarioID` is a free-text label, 40 characters or fewer, stamped into
every row of the rate tables. It is required for Emission Rates runs. Because it
is carried in the output, **all four runs can share a single output database**
and be separated by scenario ID afterwards — simpler than juggling four
databases.

| MOVESScenarioID | Year | Months | Output database | Purpose |
|---|---|---|---|---|
| `BOS_2025` | 2025 | Jan + Jul | `bos_2025_out` | Base year |
| `BOS_2050` | 2050 | Jan + Jul | `bos_2050_out` | Fleet-turnover endpoint for the ICE baseline |

Separate output databases per run, because the truncate flags are easy to leave
on by accident. Keep scenario IDs to letters, digits and underscores — the ID
ends up in table keys and joins.

The two 2050 runs are what stop the model from crediting AV adoption with
reductions that ordinary fleet turnover would have delivered anyway.

## Correction to earlier guidance — pollutant selection

In the Pollutants and Processes panel, PM2.5 is **not one entry**. Selecting
only "Primary Exhaust PM2.5 - Total" gives you exhaust and nothing else — no
brake wear, no tire wear, i.e. two of Eq. (15)'s three terms silently missing.

Select all three:

| Pollutant | ID | Pair with process |
|---|---|---|
| Primary Exhaust PM2.5 - Total | 110 | Running Exhaust (1) |
| Primary PM2.5 - Brakewear Particulate | 116 | Brakewear (9) |
| Primary PM2.5 - Tirewear Particulate | 117 | Tirewear (10) |

MOVES enforces prerequisite chains — selecting PM2.5 Total may require Total
Energy Consumption and the carbon/sulfate components as well. Let it add them;
they cost nothing and you can ignore the extra rows.

Add **Start Exhaust (process 2)** to the selection even though your line-source
model excludes it. Query 7(c) then tells you how large the exclusion is, so it
can be stated with a number rather than asserted.

## Getting to the County Data Manager

Not a menu item. It is the **"Create Input Database"** panel in the left
navigation list, then the **"Enter/Edit Data"** button.

It stays greyed out with a red X until every other RunSpec panel shows a green
checkmark. If you cannot see it, that is why — look for red X marks in the left
list, usually Scale (must be County, not Default) or an invalid combination in
Vehicles/Equipment.

Per tab the workflow is **Export Default Data → Browse → select worksheet →
Import → green check**. Even accepting EPA defaults requires that round trip;
there is no one-click generate. Save the exports to `cdm_exports/`.

## Troubleshooting

### "ERROR: Unable to validate County Data Status"

Check the database name first. If it contains `out` — e.g. `bos_2026_out_bet` —
you are pointing the **input** database field at your **output** database. MOVES
finds output tables where it expects county input tables and cannot validate.

Input and output are two separate databases:

| Where | Name |
|---|---|
| General Output panel | `bos_2026_jan_out` |
| CDM Database tab | `bos_2026_jan_in` |

Type a new name in the CDM Database field and click **Create Database**.

### Do you even need County Scale?

Probably not. The Scale panel states that County Scale "requires user-supplied
local data for most activity and fleet inputs" — that is why the CDM tabs show
red X marks. County Scale exists for SIP and conformity submissions.

For emission *rates*, use **Default Scale + Emission Rates** instead. No CDM, no
local data, and the output is g per unit distance indexed by speed bin — which
is what the speed-coupling fix needs anyway. See TODO_OPEN_PARAMETERS.md Task 1.

### "Unable to use this entry as a County Domain database. The database does not have the required county."

Answer **No** (discard). Keeping a mismatched input database only defers the
failure to run time.

The message means MOVES checked the database you named against the county in
Geographic Bounds and did not find that county's data in it. Causes, in the
order they usually turn out to be true:

1. **You picked an existing database from the dropdown instead of creating a new
   one.** The dropdown lists every database on the server — MOVES' own default
   database, output databases, leftovers from other work. Almost none of them are
   valid county input databases for your county. Type a name that does not exist
   yet and click **Create Database**.
2. **Geographic Bounds has no county actually selected.** Highlighting Suffolk
   County in the list is not the same as adding it — you have to click the arrow
   button to move it into the selected pane on the right. If nothing was added,
   MOVES has no "required county" to match against. Go back to Geographic Bounds
   and confirm Suffolk County appears in the selected list, and that the panel
   shows a green checkmark.
3. **You are pointing at the output database.** `bos_2026_jan_out` is where
   results go; the input database is a separate one, `bos_2026_jan_in`.
4. **The name is left over from an earlier attempt with a different county
   selected.** Drop it, or use a fresh suffix (`bos_2026_jan_in_v2`).

Also confirm Server reads `localhost`, and that County scale has **exactly one**
county selected — not zero, not several.

## Running it

**1. Confirm which scale you are on.** This decides whether the County Data
Manager matters at all.

- **Default Scale** — the CDM is optional. Any red X tabs inside it are
  irrelevant; you do not need to fill them. Close the CDM and move on.
- **County Scale** — every red X tab must be filled (Export Default Data →
  Browse → Import) before MOVES will run. If you are still here and did not
  intend to be, switch to Default Scale and save the effort.

**2. Save the RunSpec before running.** File → Save As → `moves/bos_2026_jan.mrs`.

Commit it. This single file is what lets anyone reproduce the emission factors,
and it is the first thing a referee asking about your inputs will want. Save one
per run in the matrix.

**3. Check the left navigation list is all green.** MOVES will refuse to execute
with any red X outstanding. The Create Input Database panel should be green or
grey, not red.

**4. Action → Execute** (or the green arrow on the toolbar).

Expect 20–60 minutes. Progress appears in the log window at the bottom, not in a
dialog. **Read the log even on a successful run** — MOVES reports missing-data
warnings there rather than stopping, and a run that "succeeded" while silently
substituting defaults for something important looks identical to a clean one.

**5. Repeat for the other three scenario IDs.** Change only the year and month
in Time Spans and the MOVESScenarioID; leave everything else alone. Point all
four at the same output database.

## After the run

1. **Post Processing → Run MySQL Script on Output Database**, choose
   `extract_rates_emissionrates.sql` if you ran **Emission Rates**, or
   `extract_emission_factors.sql` if you ran **Inventory**. Picking the wrong
   one returns nothing, because the two calculation types write to different
   tables.
2. Run the reference queries at the top first and confirm the ID codes against
   your own database. Do not trust the comments — MOVES has renumbered things
   between versions before.
3. Export `ef_rates` to `output/<dbname>.csv`.
4. Read the four sanity checks in section 7 before using anything.
5. Update `config/params.R`:

```r
EF$ice$exhaust   # gasoline, processID 1
EF$ice$brake     # gasoline, processID 9
EF$ice$tire      # gasoline, processID 10
EF$av$brake      # electricity, processID 9
EF$av$tire       # electricity, processID 10
                 # EF$av$exhaust stays 0 - BEVs have no tailpipe
```

In each `source` field record: MOVES version, database name, county, year,
month, hour, road type, and the RunSpec filename. That single line is what
makes the number auditable.

6. `Rscript run_all.R` and compare against `outputs/policy_costs.csv`.

## Watch for

**If electric-vehicle distance is ~zero in the base year** (sanity check 7a),
you cannot derive an AV factor from that run. Two options: use the 2050 run for
the AV factor, or fall back to the MOVES5 brake-wear ratio (0.115/0.400 g/hr)
applied to the gasoline rate, which is what `config/params.R` currently does.

**If MOVES' fleet population is smaller than the City of Boston MVC count**
(7b), something is wrong — MOVES covers all of Suffolk County, which is
strictly larger than Boston.

**If MOVES' average speed differs materially from your CTM speeds** (7d), the
emission factors and the traffic model are describing different roads. See the
speed-coupling section in `TODO_OPEN_PARAMETERS.md`.
