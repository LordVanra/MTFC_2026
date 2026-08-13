# =============================================================================
#  moves/ingest_moves_rates.R
#  Turns exported MOVES5 rateperdistance tables into the emission factors used
#  by config/params.R, and into a speed-indexed table for the CTM coupling.
#
#  Inputs : moves/output/bos_*_01_rateperdistance.txt  (tab separated)
#  Outputs: moves/output/ef_by_speed.csv     - rate by fuel/process/speed bin
#           moves/output/ef_summary.csv      - headline factors
# =============================================================================

ROAD      <- 5      # Urban Unrestricted Access (arterials)
HOUR      <- 8      # 07:00-07:59
SPEED_BIN <- 7      # ~30 mph = 48 km/h, matches CTM free-flow 50 km/h

# Light-duty VMT split. US light trucks have exceeded cars for years; 60/40 is
# the conventional split. THIS IS AN ASSUMPTION because the exported
# rateperdistance extract does not include VMT by source type for weighting.
W_CAR   <- 0.40
W_TRUCK <- 0.60

files <- list.files("moves/output", pattern = "rateperdistance.txt$", full.names = TRUE)
stopifnot(length(files) > 0)

raw <- do.call(rbind, lapply(files, function(f) {
  d <- read.delim(f, stringsAsFactors = FALSE)
  d$srcfile <- basename(f)
  d
}))

names(raw) <- sub("^ratePerDistance$", "rate", names(raw))
d <- raw[raw$roadTypeID == ROAD, ]

fuel_name <- c("1" = "Gasoline", "2" = "Diesel", "5" = "E85", "9" = "Electric")
proc_name <- c("1" = "RunningExhaust", "9" = "Brakewear", "10" = "Tirewear")

d$fuel    <- fuel_name[as.character(d$fuelTypeID)]
d$process <- proc_name[as.character(d$processID)]

# ---- speed-indexed table (for the CTM coupling) ----------------------------
agg <- aggregate(rate ~ yearID + monthID + hourID + fuelTypeID + fuel +
                        processID + process + avgSpeedBinID + sourceTypeID,
                 data = d, FUN = sum)

# MOVES average-speed-bin midpoints
bin_mph <- c(2.5,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75)
agg$bin_mph  <- bin_mph[agg$avgSpeedBinID]
agg$bin_kmh  <- round(agg$bin_mph * 1.609344, 1)

write.csv(agg[order(agg$yearID, agg$monthID, agg$fuelTypeID,
                    agg$processID, agg$avgSpeedBinID), ],
          "moves/output/ef_by_speed.csv", row.names = FALSE)

# ---- headline factors -------------------------------------------------------
sel <- agg[agg$hourID == HOUR & agg$avgSpeedBinID == SPEED_BIN, ]

fleet_rate <- function(year, month, fuelid, procid) {
  s <- sel[sel$yearID == year & sel$monthID == month &
           sel$fuelTypeID == fuelid & sel$processID == procid, ]
  car   <- sum(s$rate[s$sourceTypeID == 21])
  truck <- sum(s$rate[s$sourceTypeID == 31])
  W_CAR * car + W_TRUCK * truck
}

out <- expand.grid(year = sort(unique(sel$yearID)), month = c(1, 7),
                   fuelid = c(1, 9), stringsAsFactors = FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(out)), function(i) {
  y <- out$year[i]; m <- out$month[i]; f <- out$fuelid[i]
  ex <- fleet_rate(y, m, f, 1); br <- fleet_rate(y, m, f, 9); ti <- fleet_rate(y, m, f, 10)
  data.frame(year = y, month = m, fuel = fuel_name[as.character(f)],
             exhaust = ex, brake = br, tire = ti, total = ex + br + ti)
}))

write.csv(res, "moves/output/ef_summary.csv", row.names = FALSE)

cat("\n=== Fleet-weighted PM2.5 (g/km), urban arterial, 07:00, ~48 km/h ===\n")
cat(sprintf("    car/truck VMT weights %.2f / %.2f (assumption, not from MOVES)\n\n", W_CAR, W_TRUCK))
print(format(res, digits = 4), row.names = FALSE)

for (y in sort(unique(res$year))) {
  ice <- res$total[res$year == y & res$month == 1 & res$fuel == "Gasoline"]
  ev  <- res$total[res$year == y & res$month == 1 & res$fuel == "Electric"]
  cat(sprintf("\n%d: EF_ICE = %.5f, EF_AV = %.5f -> max reduction %.1f%%\n",
              y, ice, ev, (1 - ev/ice) * 100))
}
