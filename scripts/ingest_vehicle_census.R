# =============================================================================
#  scripts/ingest_vehicle_census.R
#  Massachusetts Vehicle Census (MassDOT + MAPC), ZIP-code annual file.
#  Source: https://gis.data.mass.gov/datasets/dce8f7ce88f449f08150af4183b9a983
#
#  Produces three things the model did not previously have:
#    1. outputs/mvc_zip_vmt.csv       - observed VMT per ZIP, for CTM validation
#    2. outputs/mvc_adoption.csv      - observed ZEV/hybrid VMT share 2020-2025,
#                                       a real calibration target for the SSM
#    3. outputs/mvc_summary.txt       - headline figures
#
#  NOTE the file carries ANNUALVMT only - no vehicle counts. C_0 must still be
#  sourced elsewhere, or derived by dividing VMT by an assumed per-vehicle
#  annual mileage (done below and clearly flagged as derived).
# =============================================================================

f <- "data/MVC_Annual_Zip_Code.csv"
stopifnot(file.exists(f))
d <- read.csv(f, stringsAsFactors = FALSE, colClasses = c(ZIPCODE = "character"))
names(d) <- toupper(names(d))
d$ANNUALVMT <- as.numeric(d$ANNUALVMT)
d$MUNICIPALITY <- toupper(trimws(d$MUNICIPALITY))

# The 28 ZIPs the manuscript models: 02108-02136 plus 02199, 02203, 02210, 02215.
# NOTE 02203/02210/02215 are 022xx, not 021xx - an easy filter mistake.
MODELLED <- c(sprintf("021%02d", 8:36), "02199", "02203", "02210", "02215")

bos <- d[d$MUNICIPALITY == "BOSTON", ]

# ---- 1. adoption trajectory (all Boston-garaged vehicles) -------------------
ad <- aggregate(ANNUALVMT ~ YEAR + FUELCLASS, data = bos, FUN = sum)
ad <- reshape(ad, idvar = "YEAR", timevar = "FUELCLASS", direction = "wide")
names(ad) <- sub("^ANNUALVMT\\.", "", names(ad))
ad$TOTAL    <- rowSums(ad[, setdiff(names(ad), "YEAR")], na.rm = TRUE)
ad$zev_share <- ad$`Zero-Emission` / ad$TOTAL
ad$hyb_share <- ad$Hybrid / ad$TOTAL
ad <- ad[order(ad$YEAR), ]
write.csv(ad, "outputs/mvc_adoption.csv", row.names = FALSE)

# ---- 2. per-ZIP VMT, most recent year ---------------------------------------
latest <- max(bos$YEAR)
z <- aggregate(ANNUALVMT ~ ZIPCODE + FUELCLASS, data = bos[bos$YEAR == latest, ], FUN = sum)
zw <- reshape(z, idvar = "ZIPCODE", timevar = "FUELCLASS", direction = "wide")
names(zw) <- sub("^ANNUALVMT\\.", "", names(zw))
zw[is.na(zw)] <- 0
zw$TOTAL <- rowSums(zw[, setdiff(names(zw), "ZIPCODE")])
zw$modelled <- zw$ZIPCODE %in% MODELLED
zw <- zw[order(-zw$TOTAL), ]
write.csv(zw, "outputs/mvc_zip_vmt.csv", row.names = FALSE)

# ---- 3. summary --------------------------------------------------------------
tot_all <- sum(zw$TOTAL)
tot_mod <- sum(zw$TOTAL[zw$modelled])

# Derived vehicle count. US average annual VMT per light vehicle ~11,400 mi
# (FHWA Highway Statistics VM-1). Urban Boston drives less than average, so
# this is an UPPER bound on the implied fleet.
VMT_PER_VEH <- 11400
c0_derived  <- tot_mod / VMT_PER_VEH

sink("outputs/mvc_summary.txt")
cat("Massachusetts Vehicle Census - Boston, ", latest, "\n\n", sep = "")
cat(sprintf("Total Boston-garaged VMT        : %15.0f mi/yr\n", tot_all))
cat(sprintf("VMT in the 28 modelled ZIPs     : %15.0f mi/yr  (%.1f%% of city)\n",
            tot_mod, 100 * tot_mod / tot_all))
cat(sprintf("Derived fleet at %d mi/veh/yr : %15.0f vehicles\n\n", VMT_PER_VEH, c0_derived))
cat("Observed fuel-class VMT share, Boston:\n")
print(ad[, c("YEAR", "zev_share", "hyb_share", "TOTAL")], row.names = FALSE, digits = 4)
cat("\nZEV VMT share 2020 -> 2025: ",
    sprintf("%.2f%% -> %.2f%%", 100*ad$zev_share[1], 100*ad$zev_share[nrow(ad)]), "\n")
cat("Total VMT change 2020 -> 2025: ",
    sprintf("%+.1f%%", 100*(ad$TOTAL[nrow(ad)]/ad$TOTAL[1] - 1)), "\n")
sink()

cat(readLines("outputs/mvc_summary.txt"), sep = "\n")
cat("\n\nTop 10 modelled ZIPs by VMT:\n")
print(head(zw[zw$modelled, c("ZIPCODE","Not Advanced","Hybrid","Zero-Emission","TOTAL")], 10),
      row.names = FALSE)
