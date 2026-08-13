# Parameter Sourcing — What Was Found, What Changed, What Is Still Open

Every URL below was actually retrieved. Nothing here comes from an aggregator,
SEO "statistics" site, or news write-up of a study. Where a value is a
*derivation* rather than a quoted figure, it says so — those must be presented
as derivations in the paper, not as borrowed numbers.

---

## Corrections that change results materially

| Parameter | Was | Now | Magnitude |
|---|---|---|---|
| **Boston vehicle stock `C_0`** | 512,539 | **237,224** | **2.2× too high** |
| **ICE exhaust PM₂.₅** | 0.012 g/km | **0.0015 g/km** | **~8× too high** |
| **Arterial capacity `q_max`** | 1,800 veh/h | **950 veh/h** | **~1.9× too high** |
| BEV lifecycle CO₂ | 100 g/km | 141 g/km | 30% optimistic |
| AV price | $48,000 (SEO site) | $51,377 (derived) | citation replaced |
| Infra depreciation `δ₂` | 0.023 (paper) / 0.05 (code) | **0.0202** | BEA's actual rate |
| ICE brake wear | 0.001 g/km | 0.00172 g/km | 42% low |
| ICE tire wear | 0.0006 g/km | 0.00080 g/km | 25% low |
| BEV brake reduction | 50% assumed | **71%** (MOVES5) | understated |
| BEV tire wear | "not reduced" | **+8.5%** (heavier) | wrong direction |
| BEV residual PM | 0.0019 g/km | **0 — deleted** | was a back-solved plug |
| Peak-hour share | 0.095 placeholder | **0.0724** (measured) | now empirical |
| Street width | 20 m assumed | **16.5 m** (derived) | now sourced |

**Net effect on the emission chain:** EF_ICE falls from 0.0136 to **0.00402
g/km**, EF_AV from 0.003 to **0.00137 g/km**. The maximum achievable PM₂.₅
reduction from full fleet conversion drops from 78% to **66%**.

---

## 1. Boston vehicle stock — the 512,539 figure is wrong

**237,224 registered passenger vehicles, Boston, Q4 2014.**
Massachusetts Vehicle Census (MassDOT + MAPC), as published in City of Boston,
*Go Boston 2030: Vision and Action Plan* (2017), pp. 54–55.
<https://www.boston.gov/sites/default/files/document-file-03-2017/go_boston_2030_-_3_boston_today_spreads.pdf>

512,539 implies ~1.9 vehicles per Boston household in a city where ACS shows
~35% of households are car-free. 237,224 gives ~0.88 veh/household, which is
consistent with ACS. Note that the same Go Boston chart plots 639,594 = Boston
*population* on a second axis; the bad figure may be a garbled derivative.

**Still open:** 2014 is the newest figure retrievable without a browser. The
current-year Boston row is available from the MassDOT Vehicle Census dashboard
— <https://geodot-massdot.hub.arcgis.com/pages/vehicle-census> — and is a
two-minute lookup for someone with Chrome open. Methodology:
<https://www.mapc.org/wp-content/uploads/2021/07/MA_VehicleCensus_v3_Documentation.pdf>

## 2. National fleet and sales

- **289 million light vehicles in operation** (US, 1 Jan 2025) — S&P Global Mobility, 21 May 2025. <https://press.spglobal.com/2025-05-21-U-S-Vehicle-Age-Rises-Again-to-12-8-Years-in-2025>
  FHWA Highway Statistics Table MV-1 (2024) gives 297.5M *all* motor vehicles but **cannot** supply a light-duty split — MV-9 is currently suspended. S&P is the series BTS itself uses for NTS Table 1-26.
- **16.2 million new light-vehicle sales, CY2025** — NADA Market Beat, 31 Dec 2025. <https://www.nada.org/nada/nada-headlines/december-2025-market-beat-new-light-vehicle-sales-totaled-162-million-units>
  **Citation error found:** the paper cites FRED **TOTALSA**, which is *total* vehicle sales including heavy trucks. The light-duty series is **ALTSALES**.
- **Scrappage 4.5%/yr** — S&P Global Mobility, same release. This is almost certainly where the paper's 4.5% came from, and it *is* a replacement-flow measure (scrapped ÷ VIO), so the paper's "not 1/average-age" caveat is now defensible rather than asserted. Average age 12.8 yr (cars 14.5, light trucks 11.9).
  For a non-constant hazard: Greene & Leard (2023), Baker Center, Univ. of Tennessee — median lifetimes cars ~17 yr, SUVs ~20 yr, pickups ~25 yr. <https://baker.utk.edu/wp-content/uploads/2023/03/Formatted-Vehicle-Scrappage-and-Survival_Single-Column_7Mar23.pdf>

## 3. AV price — the $48,000 was probably an average car price

Kelley Blue Book / Cox Automotive: average new-vehicle **transaction price
$49,855** (July 2026). <https://www.coxautoinc.com/insights/july-2026-atp-report/>

The paper's $48,000 is within 4% of the plain ATP, which strongly suggests the
aggregator conflated an average new-car price with an AV price.

Replacement is built explicitly: ATP + **Level-2 ADAS premium $750–$2,295**
(Goddard, McDonald, Wei & Batra 2022, *Findings*, DOI 10.32866/001c.38291 —
peer-reviewed, US, MSRP-based). Midpoint gives $51,377.
**L4 upper bound for sensitivity:** ~$200,000/vehicle (Riggs & Richardson 2024,
SSRN 4998828 — *working paper*, cite as such).

## 4. Emission factors — the exhaust figure is an order of magnitude out

- **ICE brake 2.77 mg/mi = 0.00172 g/km**, **tire 1.28 mg/mi = 0.00080 g/km** — EPA-420-R-20-014 (2020). <https://www.epa.gov/sites/default/files/2020-11/documents/420r20014.pdf> EPA's *Overview of MOVES4* confirms brake/tire rates are unchanged from MOVES3, so these **are** the MOVES4 values.
- **ICE exhaust:** the paper's 0.012 g/km = 19.3 mg/mile. EPA's own light-duty benchmarks are ~1–3 mg/mile (Tier 3 standard 3 mg/mi, LEV III 1 mg/mi from 2025). <https://19january2021snapshot.epa.gov/sites/static/files/2017-04/documents/04_-_light-duty_pm_emission_rates_update.pdf> The paper's value looks like a pre-Tier-2 or non-DPF diesel rate.
- **Regenerative braking = 71% reduction, not 50%** — EPA MOVES5 (EPA-420-R-24-012, Nov 2024) Table 2.12: gasoline car 0.400 g/hr vs electric 0.115 g/hr. **Critical:** MOVES4 applies the *same* brake/tire rates to all fuel types, so a regen credit **cannot** be cited to MOVES4 — it must be cited to MOVES5.
- **BEV tire wear is +7–10% higher**, not unchanged — EMEP/EEA Guidebook 2023 ch.1.A.3.b.vi–vii. Corroborated by Beddows & Harrison (2021, DOI 10.1016/j.atmosenv.2020.117886) and OECD (2020).
- **The 0.0019 g/km "residual PM" is deleted.** No primary source contains an unallocated residual PM category for BEVs. It existed only to make the AV factor reach the 0.003 printed in Eq. (15).

**Framework warning:** the EU inventory (EMEP/EEA, Beddows & Harrison) gives
total non-exhaust PM₂.₅ roughly **3× higher** than the MOVES-consistent
construction. Pick one framework, state which, and report the other as a
sensitivity — do not average them silently.

## 5. Traffic parameters

- **Capacity 950 veh/h/lane** — derived per FHWA's own method: saturation flow **1,900 pc/h/lane** (HPMS Field Manual App. N p. N-19, confirmed HCM 2000 p. 16-10) × green ratio **g/C = 0.50** (FHWA-PL-18-003 Table 3). The previous 1,800 was a *freeway* capacity with an undocumented "arterial discount"; signalised arterial through-capacity is about half, because green time is shared.
- **Free-flow speed 50 km/h** — value unchanged but now sourced to *observed* 85th-percentile speed of 31 mph on Boston streets (Hu & Cicchino, *Injury Prevention* 2019) rather than the statutory limit. Boston's posted default is 25 mph since Jan 2017 and does **not** cover state-owned roadways in the city.
- **Peak-hour share K = 7.24%** — empirical, MassDOT Mount Auburn Street Corridor Study, Appendix A: 977 veh/h peak of 13,494 veh/day, Cambridge MA, Wed 27 Apr 2016. <https://www.mass.gov/doc/mount-auburn-report-appendix-a-traffic-counts/download> Consistent with FHWA's 7–12% range (FHWA-PL-18-027 pp. 44–45). Caveat honestly: one direction, one arterial, one weekday.
- **Diurnal factor 0.576** — derived as 1/(24K), algebraically tied to K in `derive_params()` so the two cannot drift apart.
- **Street width 16.5 m** — derived from Boston Complete Streets Design Guidelines arterial minima (4 travel lanes @ 10 ft + parking both sides @ 7 ft = 54 ft). <https://www.boston.gov/sites/default/files/file/2020/09/Minimum%20Widths%20for%20Roadway%20Lanes.pdf> These are **minima**, and BCS publishes no total kerb-to-kerb dimension, so this is a construction. Sensitivity range 10.4–16.5 m. It scales concentration linearly.
- **Jam density 131 veh/km/lane** — derived from 25 ft stopped spacing (Wu, Liu & Geroliminis 2011, *TR-B* 45(3):255–266, §5). **No source states an arterial jam density directly.** Lloret-Batlle & Zheng (2023, *TR-B* 173:162–175, DOI 10.1016/j.trb.2023.02.007) is the natural citation and makes exactly this criticism — that the field assumes this parameter "with no estimation from data whatsoever." Worth obtaining via institutional access. Do **not** import freeway jam densities (100–110 veh/km/lane); motorway values are lower than arterial.

## 6. VMT rebound — now a defended range instead of one news article

The single citation (Brasuell 2021, a Planetizen news item) is replaced.

**Supporting the +30–45% band:**
- **Sun et al. (2023)**, *TRR* 2678(4), DOI 10.1177/03611981231186984 — California statewide activity-based model + EMFAC to 2050. Private CAV at 100% penetration: **+21% to +64%** VMT. The assumed band sits squarely inside. **Best single anchor** — a full regional simulation.
- **Harb et al. (2022)**, UC-ITS-2018-09, DOI 10.7922/G2WH2N96 — 43 Sacramento households given a chauffeur: **+60% VMT, over half zero-occupancy**. <https://rosap.ntl.bts.gov/view/dot/65477> The assumed band sits *below* the only revealed-behaviour evidence, so it can be framed as conservative. (Precursor: Harb et al. 2018, *Transportation* 45(6), N=13, +83%.)
- **Taiebat, Stolper & Xu (2019)**, *Applied Energy* 247:297–308 — VMT price elasticity **−0.4**, households more sensitive to *time* than fuel cost, +2% to +47% range. <https://arxiv.org/abs/1902.00382> This is the mechanism that lets you *derive* the number rather than assert it.
- **Duranton & Turner (2011)**, *AER* 101(6):2616–2652 — capacity elasticity **1.03**. Hsu & Zhang (2014) find 1.24–1.34; Hymel (2019) unit elasticity; Volker & Handy (2021) adopt 1.0/0.75 for California regulatory use. Establishes that a large rebound is the *prior*, not the exception.
- **Fagnant & Kockelman (2015)**, *TRR* 2536:98–106 — **+8%** from empty repositioning alone, for bottom-up decomposition.

**Ride-hailing, properly cited:**
- **Schaller (2021)**, *Transport Policy* 102:1–10 — the peer-reviewed original behind the "97–118%" the paper cited second-hand, and it has a **Boston-specific +157%**. Deadhead 40–48%.
- Henao & Marshall (2019), *Transportation* 46(6) — Denver, +83.5%, deadheading ≥40.8%.
- Erhardt et al. (2019), *Science Advances* 5(5):eaau2670 — SF delay +62% vs +22% counterfactual.
- **Caveat that must appear in the text:** TNC evidence bounds the mechanism but is not directly transferable to private AV ownership.

**What genuinely threatens the assumption — must be confronted, not omitted:**
- **Naz & Mattingly (2024)**, SSRN 5030045 — meta-analysis, **pooled mean +5.95%** across 26 articles / 195 effect sizes. This is roughly **one-sixth** of the assumed low case. The defence is that the pool is dominated by shared-AV, low-penetration simulations with conservative VOT reductions, whereas this paper models private ownership at high penetration. That defence has to be written down.
- **Childress et al. (2015)**, *TRR* 2493:99–106 — Puget Sound ABM, only **+3–5%** under a −35% VOT assumption.
- **ITF/OECD (2016)**, Lisbon simulation — **−23% VKT** under full pooling. The sign flips with fleet structure, so the band is conditional on private ownership, not a property of automation.

## 7. Infrastructure depreciation

**2.02%/yr**, geometric, from a 45-year service life × 0.91 declining balance
(Hulten–Wykoff). Fraumeni & Kornfeld (2024), NBER WP 32753 — Kornfeld is BEA, so
this is the authoritative description of BEA's own method.
<https://www.nber.org/system/files/working_papers/w32753/w32753.pdf>
Stock: **BEA Fixed Assets Table 7.1B**, Structures → Highways and streets.
The paper's 2.3% was close but is not the published rate.

---

## Still open — seven items

Run `Rscript -e 'source("config/params.R"); print_provenance()'`.

| Parameter | Status | What it needs |
|---|---|---|
| `EF$ice$exhaust` | **MUST BE REPLACED** | Run MOVES4 for the target year/county and cite the `rateperdistance` output. The 0.0015 placeholder is defensible but a reviewer will expect a real run. |
| `SSM$C_0` | **REFRESH NEEDED** | Current-year Boston row from the MassDOT Vehicle Census dashboard. Browser lookup, two minutes. |
| `SSM$I_ref` | **UNSOURCED** | Appears nowhere in the manuscript. Derive from a real infrastructure-stock figure (BEA 7.1B scaled to Boston) or present as an explicit scaling choice with a sweep. |
| `SSM$cap_0` | **UNSOURCED** | Year-1 adoption ceiling. No source found. |
| `SSM$gamma` | **UNSOURCED** | Capacity response to subsidy. Pure calibration constant. |
| `SSM$A_0` | **UNSOURCED** | No registry publishes L4 AV registrations for Boston. Consider setting to 0. |
| `EF$av$co2_sensor` | **CONTESTED** | Report the 4–28 g/km range, not 8 as a point value. State the implied power draw (8 g/km ≈ 1.6 kW at 50 km/h, above Sudhakar's 1.2 kW ceiling). |

**Could not verify, do not cite without checking:** Harper et al. (2016) latent-demand
figure; Mercedes DRIVE PILOT US price; Tesla FSD outright purchase price; BTS Table 1-26
numeric values 2022–25; current BEA Table 7.1B net stock; HCM 6th ed. jam density.
ScienceDirect, SAGE, Science.org and ResearchGate all block automated retrieval — several
citations above were verified through RePEc or institutional repository mirrors, so DOIs
should be confirmed against the publisher before the reference list is finalised.
