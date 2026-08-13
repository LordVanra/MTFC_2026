-- ============================================================================
--  extract_rates_emissionrates.sql
--  For runs made with Calculation Type = EMISSION RATES.
--  (Use extract_emission_factors.sql instead if you ran Inventory.)
--
--  Emission Rates runs write to rateperdistance / rateperhour / ratepervehicle
--  rather than movesoutput. The rate is given directly - there is no dividing
--  by activity - but the rows are DISAGGREGATED BY MODEL YEAR AND SPEED BIN,
--  so a fleet-average factor requires weighting. See section 3.
--
--  MOVESScenarioID lets several runs share one output database. Suggested IDs:
--    BOS_2026_JAN   BOS_2026_JUL   BOS_2050_JAN   BOS_2050_JUL
-- ============================================================================

SET @ROAD := 5;     -- Urban Unrestricted Access (arterials)
SET @HOUR := 8;     -- verify this is 07:00-07:59 in your run

-- ---------------------------------------------------------------------------
-- 0. What is in here?
-- ---------------------------------------------------------------------------
SELECT MOVESScenarioID, yearID, monthID, COUNT(*) AS rows_
FROM rateperdistance
GROUP BY MOVESScenarioID, yearID, monthID;

-- Speed bin reference - this is what makes Emission Rates worth the extra work
SELECT * FROM avgspeedbin;

-- ---------------------------------------------------------------------------
-- 1. Rates by speed bin. THE SPEED-COUPLING TABLE.
--    Export this. Join it to your CTM's per-corridor speeds so each corridor
--    gets the rate matching its own modelled speed, instead of one flat factor.
-- ---------------------------------------------------------------------------
DROP TABLE IF EXISTS ef_by_speed;
CREATE TABLE ef_by_speed AS
SELECT
    r.MOVESScenarioID,
    r.yearID, r.monthID, r.hourID,
    r.fuelTypeID,
    CASE r.fuelTypeID WHEN 1 THEN 'Gasoline'
                      WHEN 9 THEN 'Electricity'
                      ELSE CONCAT('fuel_', r.fuelTypeID) END AS fuel_name,
    r.processID,
    CASE r.processID WHEN 1  THEN 'RunningExhaust'
                     WHEN 9  THEN 'Brakewear'
                     WHEN 10 THEN 'Tirewear'
                     ELSE CONCAT('process_', r.processID) END AS process_name,
    r.avgSpeedBinID,
    b.avgBinSpeed,
    SUM(r.ratePerDistance) AS rate_per_distance
FROM rateperdistance r
LEFT JOIN avgspeedbin b ON b.avgSpeedBinID = r.avgSpeedBinID
WHERE r.roadTypeID = @ROAD
  AND r.pollutantID IN (110, 116, 117)
GROUP BY r.MOVESScenarioID, r.yearID, r.monthID, r.hourID,
         r.fuelTypeID, r.processID, r.avgSpeedBinID, b.avgBinSpeed;

-- ---------------------------------------------------------------------------
-- 2. Single fleet-average factor per fuel type, for config/params.R.
--    NOTE: this is an UNWEIGHTED mean across model years. It is fine as a first
--    pass but is not the fleet-weighted number - see section 3 before quoting
--    it in the paper.
-- ---------------------------------------------------------------------------
SELECT
    MOVESScenarioID, yearID, monthID, fuel_name, process_name,
    ROUND(AVG(rate_per_distance), 6) AS rate_unweighted
FROM ef_by_speed
WHERE hourID = @HOUR
GROUP BY MOVESScenarioID, yearID, monthID, fuelTypeID, processID
ORDER BY MOVESScenarioID, fuelTypeID, processID;

-- ---------------------------------------------------------------------------
-- 3. Fleet-weighted factor. THIS is the number for the paper.
--    Weights each model year by its share of travel, using the run's own
--    age distribution and relative mileage. If sourcetypeagedistribution is
--    absent (Default Scale keeps it in the default database, not the output
--    database), read it from the MOVES default DB and adjust the schema
--    prefix below.
-- ---------------------------------------------------------------------------
-- SELECT
--     r.fuelTypeID, r.processID,
--     SUM(r.ratePerDistance * a.ageFraction) / SUM(a.ageFraction) AS rate_weighted
-- FROM rateperdistance r
-- JOIN <movesdb>.sourcetypeagedistribution a
--   ON  a.sourceTypeID = r.sourceTypeID
--   AND a.yearID       = r.yearID
--   AND a.ageID        = r.yearID - r.modelYearID
-- WHERE r.roadTypeID = @ROAD AND r.pollutantID IN (110,116,117)
-- GROUP BY r.fuelTypeID, r.processID;

-- ---------------------------------------------------------------------------
-- 4. Sanity check: is there any electric-vehicle rate at all?
-- ---------------------------------------------------------------------------
SELECT MOVESScenarioID, yearID, fuel_name, COUNT(*) AS n_rows,
       ROUND(AVG(rate_per_distance), 6) AS mean_rate
FROM ef_by_speed
GROUP BY MOVESScenarioID, yearID, fuelTypeID
ORDER BY MOVESScenarioID, yearID, fuelTypeID;
