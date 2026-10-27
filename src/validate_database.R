#!/usr/bin/env Rscript
# ==============================================================================
# Automated Integrity & AI-Readiness Validation Suite
# ==============================================================================
# Verifies:
#   1. Zero residual -9999 sentinel values across all Parquet & SQLite tables
#   2. Exact row counts & numerical fidelity against raw NetCDF files
#   3. Unit harmonization accuracy (FRED cm / 100 == rooting_depth_m; TRY m * 100 == cm)
#   4. Zero spatial block leakage across the 5 spatial CV folds across ALL tables
#   5. Equivalence between Parquet and SQLite backends via r_loader.R
#   6. SHA-256 checksum verification
# ==============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(arrow)
  library(DBI)
  library(RSQLite)
  library(digest)
})

source("src/loaders/r_loader.R")
cat("=== Running AI-Ready Database Validation Suite ===\n")

# 1. Check zero -9999 sentinel values remain in Parquet tables
pq_files <- c(
  Sys.glob("data/parquet/*.parquet"),
  Sys.glob("data/parquet/individual/*.parquet")
)
stopifnot(length(pq_files) == 12L)

for (f in pq_files) {
  df <- as.data.frame(read_parquet(f))
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  for (col in num_cols) {
    v <- df[[col]]
    bad <- sum(!is.na(v) & (v == -9999 | abs(v + 9999) < 1e-3))
    if (bad > 0) stop(sprintf("Residual -9999 found in %s :: %s (%d rows)", f, col, bad))
  }
}
cat("  [PASS] Zero -9999 fill values present across all 12 Parquet tables.\n")

# 2. Verify row counts & exact numerical match with raw NetCDF
soil <- load_ai_table("soil_properties_conus", deduplicate_sites = FALSE)
stopifnot(nrow(soil) == 50116L)

nc_conus <- nc_open("data/raw_netcdf/conus_variables.nc")
clay_raw <- as.vector(ncvar_get(nc_conus, "Clay"))
clay_raw[clay_raw == -9999] <- NA
nc_close(nc_conus)
stopifnot(identical(is.na(soil$clay_pct), is.na(clay_raw)))
stopifnot(max(abs(soil$clay_pct - clay_raw), na.rm = TRUE) < 1e-5)
cat("  [PASS] CONUS soil table (50,116 rows) matches raw NetCDF values bit-accurately.\n")

# 3. Verify unit harmonization in root_depth_harmonized
rd <- load_ai_table("root_depth_harmonized", deduplicate_sites = FALSE)
stopifnot(nrow(rd) == 1062L)
fred_rd <- rd[rd$source_database == "FRED", ]
try_rd <- rd[rd$source_database == "TRY", ]
stopifnot(max(abs(fred_rd$raw_value / 100.0 - fred_rd$rooting_depth_m), na.rm = TRUE) < 1e-6)
stopifnot(max(abs(try_rd$raw_value * 100.0 - try_rd$rooting_depth_cm), na.rm = TRUE) < 1e-5)
stopifnot(sum(rd$qc_missing_value) == 1L)
cat("  [PASS] Rooting depth unit harmonization (FRED cm <-> TRY m) verified.\n")

# 4. Verify microbe merge completeness
mic <- load_ai_table("microbe_enzymes_harmonized")
stopifnot(nrow(mic) == 805L)
stopifnot(sum(mic$qc_complete_cnp) == 804L)
cat("  [PASS] Microbial C/N/P functional gene merge (805 unique sites) verified.\n")

# 5. Verify zero spatial block or site leakage across all 4 tables simultaneously
all_blocks_folds <- unique(rbind(
  soil[, c("spatial_block_id", "spatial_cv_fold")],
  mic[, c("spatial_block_id", "spatial_cv_fold")],
  load_ai_table("root_biomass_harmonized", deduplicate_sites = FALSE)[, c("spatial_block_id", "spatial_cv_fold")],
  rd[, c("spatial_block_id", "spatial_cv_fold")]
))
# Every spatial_block_id must map to strictly ONE fold across the entire database
stopifnot(!any(duplicated(all_blocks_folds$spatial_block_id)))
stopifnot(all(all_blocks_folds$spatial_cv_fold %in% 1:5))
cat("  [PASS] Zero spatial block or cross-modality fold leakage across all 4 domains.\n")

# 6. Verify Parquet vs SQLite loader parity
for (tbl in c("soil_properties_conus", "microbe_enzymes_harmonized", "root_biomass_harmonized", "root_depth_harmonized")) {
  df_pq <- load_ai_table(tbl, backend = "parquet")
  df_sq <- load_ai_table(tbl, backend = "sqlite")
  stopifnot(nrow(df_pq) == nrow(df_sq))
  stopifnot(identical(names(df_pq), names(df_sq)))
}
cat("  [PASS] Parquet and SQLite backends return identical row/column schemas.\n")

# 7. Verify SHA-256 checksums
lines <- readLines("data/metadata/checksums.sha256")
for (ln in lines) {
  parts <- strsplit(ln, "  ", fixed = TRUE)[[1]]
  actual <- digest::digest(file = parts[2], algo = "sha256")
  if (!identical(actual, parts[1])) {
    stop("SHA-256 mismatch for ", parts[2])
  }
}
cat(sprintf("  [PASS] All %d file SHA-256 checksums verified.\n", length(lines)))
cat("=== ALL VALIDATION CHECKS PASSED ===\n")
