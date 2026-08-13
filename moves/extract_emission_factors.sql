-- ============================================================================
--  extract_emission_factors.sql
--  Pulls PM2.5 g/km by fuel type and process from a MOVES5 County-scale run.
--
--  HOW TO RUN
--    MOVES GUI -> Post Processing -> Run MySQL Script on Output Database
--    Select this file, choose your output database, run.
--    Results land in the tables created at the bottom; export them to CSV.
--
--    Or from a shell:
--      mysql -u moves -pmoves <your_output_db> < extract_emission_factors.sql
--
--  ASSUMPTIONS
--    - Distance Units were set to KILOMETERS in the General Output panel.
--      If you left them as Miles, every g_per_km below is actually g_per_mile.
--    - "Distance Traveled" was ticked in the Activity panel.
--    - Road type 5 = Urban Unrestricted Access (arterials). Change to 4 for
--      Urban Restricted Access (freeways) if you want that instead.
-- ============================================================================

-- ---------------------------------------------------------------------------
-- Reference: confirm the codes in YOUR output database rather than trusting
-- these comments. Run these three first and read the results.
-- ---------------------------------------------------------------------------
-- SELECT * FROM activitytype;          -- which activityTypeID is Distance Traveled
-- SELECT * FROM pollutant;             -- pollutantID names
-- SELECT * FROM emissionprocess;       -- processID names
-- SELECT * FROM fueltype;              -- fuelTypeID names
-- SELECT * FROM roadtype;              -- roadTypeID names

-- Expected (verify above):
--   pollutantID 110 = Primary Exhaust PM2.5 - Total
--   pollutantID 116 = Primary PM2.5 - Brakewear Particulate
--   pollutantID 117 = Primary PM2.5 - Tirewear Particulate
--   processID     1 = Running Exhaust
--   processID     9 = Brakewear
--   processID    10 = Tirewear
--   fuelTypeID    1 = Gasoline
--   fuelTypeID    9 = Electricity
--   roadTypeID    5 = Urban Unrestricted Access
--   activityTypeID 1 = Distance Traveled

SET @ROAD := 5;          -- urban arterials
SET @DIST_ACT := 1;      -- Distance Traveled

-- ---------------------------------------------------------------------------
-- 1. Emissions, grams, by year / month / hour / fuel / pollutant / process
-- ---------------------------------------------------------------------------
DROP TABLE IF EXISTS ef_emissions;
CREATE TABLE ef_emissions AS
SELECT
    yearID,
    monthID,
    hourID,
    fuelTypeID,
    pollutantID,
    processID,
    SUM(emissionQuant) AS grams
FROM movesoutput
WHERE roadTypeID = @ROAD
  AND pollutantID IN (110, 116, 117)
GROUP BY yearID, monthID, hourID, fuelTypeID, pollutantID, processID;

-- ---------------------------------------------------------------------------
-- 2. Activity, distance, by year / month / hour / fuel
-- ---------------------------------------------------------------------------
DROP TABLE IF EXISTS ef_distance;
CREATE TABLE ef_distance AS
SELECT
    yearID,
    monthID,
    hourID,
    fuelTypeID,
    SUM(activity) AS distance
FROM movesactivityoutput
WHERE roadTypeID = @ROAD
  AND activityTypeID = @DIST_ACT
GROUP BY yearID, monthID, hourID, fuelTypeID;

-- ---------------------------------------------------------------------------
-- 3. Rates. THIS IS THE TABLE YOU EXPORT.
--    g_per_km is emissions divided by the distance travelled by the SAME
--    fuel type - not by total fleet distance. Getting that wrong is the
--    single easiest way to produce a badly wrong emission factor.
-- ---------------------------------------------------------------------------
DROP TABLE IF EXISTS ef_rates;
CREATE TABLE ef_rates AS
SELECT
    e.yearID,
    e.monthID,
    e.hourID,
    e.fuelTypeID,
    CASE e.fuelTypeID WHEN 1 THEN 'Gasoline'
                      WHEN 2 THEN 'Diesel'
                      WHEN 9 THEN 'Electricity'
                      ELSE CONCAT('fuel_', e.fuelTypeID) END AS fuel_name,
    e.pollutantID,
    e.processID,
    CASE e.processID WHEN 1  THEN 'RunningExhaust'
                     WHEN 9  THEN 'Brakewear'
                     WHEN 10 THEN 'Tirewear'
                     ELSE CONCAT('process_', e.processID) END AS process_name,
    e.grams,
    d.distance,
    CASE WHEN d.distance > 0 THEN e.grams / d.distance ELSE NULL END AS g_per_km
FROM ef_emissions e
JOIN ef_distance d
  ON  e.yearID     = d.yearID
  AND e.monthID    = d.monthID
  AND e.hourID     = d.hourID
  AND e.fuelTypeID = d.fuelTypeID
ORDER BY e.yearID, e.monthID, e.hourID, e.fuelTypeID, e.processID;

-- ---------------------------------------------------------------------------
-- 4. The three numbers the paper needs, for the 07:00 peak hour.
--    Change hourID if you are using a different peak. Note hour labels in the
--    GUI are authoritative - confirm which hourID corresponds to 7:00-7:59.
-- ---------------------------------------------------------------------------
SELECT
    yearID,
    monthID,
    fuel_name,
    process_name,
    ROUND(g_per_km, 6) AS g_per_km
FROM ef_rates
WHERE hourID = 8            -- verify this is 7:00-7:59 in your run
ORDER BY yearID, monthID, fuelTypeID, processID;

-- ---------------------------------------------------------------------------
-- 5. Totals: one ICE factor and one AV factor, summed across the three
--    processes. These map straight onto config/params.R.
-- ---------------------------------------------------------------------------
SELECT
    yearID,
    monthID,
    hourID,
    fuel_name,
    ROUND(SUM(g_per_km), 6) AS total_pm25_g_per_km
FROM ef_rates
WHERE fuelTypeID IN (1, 9)
GROUP BY yearID, monthID, hourID, fuelTypeID
ORDER BY yearID, monthID, hourID, fuelTypeID;

-- ---------------------------------------------------------------------------
-- 6. Diurnal profile: emissions-weighted, for CDM$diurnal_factor.
--    Replaces the value currently derived algebraically as 1/(24K).
-- ---------------------------------------------------------------------------
SELECT
    yearID,
    monthID,
    hourID,
    ROUND(SUM(grams), 2) AS total_grams,
    ROUND(SUM(grams) / MAX(SUM(grams)) OVER (PARTITION BY yearID, monthID), 4) AS share_of_peak_hour
FROM ef_emissions
WHERE pollutantID IN (110, 116, 117)
GROUP BY yearID, monthID, hourID
ORDER BY yearID, monthID, hourID;

-- ---------------------------------------------------------------------------
-- 7. Sanity checks. Read these before trusting anything above.
-- ---------------------------------------------------------------------------

-- (a) Is there ANY electric-vehicle distance? If this is zero or tiny, the
--     electricity rows in ef_rates are unreliable and you cannot derive an AV
--     factor from this run - use a future analysis year instead.
SELECT yearID, fuel_name, ROUND(SUM(distance), 1) AS total_distance
FROM ef_rates GROUP BY yearID, fuelTypeID ORDER BY yearID, fuelTypeID;

-- (b) Fleet population, to cross-check against the Massachusetts Vehicle
--     Census figure. MOVES models Suffolk County; the MVC figure is the City
--     of Boston, so MOVES should be somewhat LARGER. If it is smaller, or
--     larger by more than ~50%, your County Data Manager inputs are wrong.
SELECT yearID, SUM(activity) AS population
FROM movesactivityoutput
WHERE activityTypeID = 6      -- verify against SELECT * FROM activitytype
GROUP BY yearID;

-- (c) What you are excluding by modelling running emissions only.
--     Start exhaust (process 2) is real and is not in your line-source model.
SELECT processID, ROUND(SUM(emissionQuant), 2) AS grams
FROM movesoutput
WHERE pollutantID IN (110, 116, 117)
GROUP BY processID ORDER BY processID;

-- (d) Average speed actually used, for the speed-coupling check. Compare this
--     against the speeds your CTM produces. If they differ materially, the
--     emission factors and the traffic model describe different roads.
SELECT
    o.yearID, o.hourID,
    ROUND(SUM(o.activity) / NULLIF(SUM(h.activity), 0), 2) AS avg_speed_km_per_h
FROM movesactivityoutput o
JOIN movesactivityoutput h
  ON  o.yearID = h.yearID AND o.monthID = h.monthID AND o.hourID = h.hourID
  AND o.roadTypeID = h.roadTypeID AND o.fuelTypeID = h.fuelTypeID
WHERE o.activityTypeID = @DIST_ACT
  AND h.activityTypeID = 4      -- Source Hours Operating; verify the ID
  AND o.roadTypeID = @ROAD
GROUP BY o.yearID, o.hourID
ORDER BY o.yearID, o.hourID;
