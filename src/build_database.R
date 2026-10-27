#!/usr/bin/env Rscript
# ==============================================================================
# Build AI-Ready CONUS Soil, Plant Root & Microbial Functional Database
# ==============================================================================
# Converts raw CF-1.8 NetCDF point datasets in data/raw_netcdf/ into:
#   1. Harmonized, QC-flagged, spatially-folded Apache Parquet tables (data/parquet/)
#   2. Direct 1-to-1 cleaned Parquet mirrors of each NetCDF file (data/parquet/individual/)
#   3. Indexed single-file SQLite relational database (data/sqlite/conus_soil_root_microbe.sqlite)
#   4. Machine-readable variable dictionary, summary stats & SHA-256 checksums (data/metadata/)
#   5. Spatial coverage & completeness overview figure (figures/spatial_coverage_overview.png)
# ==============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(arrow)
  library(DBI)
  library(RSQLite)
  library(jsonlite)
  library(digest)
  library(maps)
  library(ggplot2)
  library(cowplot)
})

cat("=== Starting AI-Ready Database Build Pipeline ===\n")

# Ensure output directories exist
dirs <- c(
  "data/parquet",
  "data/parquet/individual",
  "data/sqlite",
  "data/metadata",
  "figures"
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

read_netcdf_table <- function(nc_path) {
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc))

  cols <- list()
  # Include dimension 'obs' if not already a variable
  if ("obs" %in% names(nc$dim)) {
    cols[["obs"]] <- as.integer(nc$dim$obs$vals)
  }

  for (vname in names(nc$var)) {
    raw_vals <- as.vector(ncvar_get(nc, vname))
    prec <- nc$var[[vname]]$prec
    # Replace -9999 fill values with NA
    is_fill <- !is.na(raw_vals) & (raw_vals == -9999 | abs(raw_vals - (-9999)) < 1e-3)
    raw_vals[is_fill] <- NA
    if (prec %in% c("int", "short", "byte")) {
      cols[[vname]] <- as.integer(round(raw_vals))
    } else {
      cols[[vname]] <- as.numeric(raw_vals)
    }
  }
  df <- as.data.frame(cols, stringsAsFactors = FALSE)
  return(df)
}

make_site_id <- function(lat, lon) {
  coord_str <- sprintf("%.6f,%.6f", round(lat, 6), round(lon, 6))
  hashes <- vapply(coord_str, function(s) {
    substr(digest::digest(s, algo = "sha1", serialize = FALSE), 1, 10)
  }, character(1), USE.NAMES = FALSE)
  paste0("SITE_", hashes)
}

make_spatial_block <- function(lat, lon, cell_deg = 2.0) {
  b_lat <- as.integer(floor(lat / cell_deg) * cell_deg)
  b_lon <- as.integer(floor(lon / cell_deg) * cell_deg)
  sprintf("B_%+03d_%+04d", b_lat, b_lon)
}

check_outside_conus_border <- function(lat, lon) {
  w <- maps::map.where("world", lon, lat)
  # Flag points that explicitly fall onto non-USA land polygons (e.g., Canada, Mexico)
  !is.na(w) & !grepl("^USA", w)
}

# ------------------------------------------------------------------------------
# 1. Read all 8 raw NetCDF files & export 1-to-1 cleaned Parquet mirrors
# ------------------------------------------------------------------------------
cat("[1/6] Reading raw NetCDF files and writing individual Parquet mirrors...\n")

raw_files <- c(
  "conus_variables.nc",
  "fred_belowground_biomass_points.nc",
  "fred_rooting_depth_points.nc",
  "try_belowground_biomass_points.nc",
  "try_rooting_depth_points.nc",
  "microbe_C_degrading_points.nc",
  "microbe_N_degrading_points.nc",
  "microbe_P_degrading_points.nc"
)

raw_dfs <- list()
for (fname in raw_files) {
  fpath <- file.path("data/raw_netcdf", fname)
  if (!file.exists(fpath)) stop("Missing raw NetCDF file: ", fpath)
  df <- read_netcdf_table(fpath)
  raw_dfs[[fname]] <- df

  out_pq <- file.path("data/parquet/individual", sub("\\.nc$", ".parquet", fname))
  arrow::write_parquet(df, out_pq, compression = "snappy")
}

# ------------------------------------------------------------------------------
# 2. Construct Unified Spatial Block Grid & Balanced 5-Fold Spatial CV Assignment
# ------------------------------------------------------------------------------
cat("[2/6] Computing leak-free 2x2 degree spatial blocks and 5-fold CV assignments...\n")

# Collect all coordinates with primary weights from CONUS soil + secondary weights from root/microbe
all_blocks_soil <- make_spatial_block(
  raw_dfs[["conus_variables.nc"]]$latitude,
  raw_dfs[["conus_variables.nc"]]$longitude
)
all_blocks_other <- c(
  make_spatial_block(raw_dfs[["microbe_N_degrading_points.nc"]]$latitude, raw_dfs[["microbe_N_degrading_points.nc"]]$longitude),
  make_spatial_block(raw_dfs[["fred_belowground_biomass_points.nc"]]$latitude, raw_dfs[["fred_belowground_biomass_points.nc"]]$longitude),
  make_spatial_block(raw_dfs[["try_belowground_biomass_points.nc"]]$latitude, raw_dfs[["try_belowground_biomass_points.nc"]]$longitude),
  make_spatial_block(raw_dfs[["fred_rooting_depth_points.nc"]]$latitude, raw_dfs[["fred_rooting_depth_points.nc"]]$longitude),
  make_spatial_block(raw_dfs[["try_rooting_depth_points.nc"]]$latitude, raw_dfs[["try_rooting_depth_points.nc"]]$longitude)
)

unique_blocks <- sort(unique(c(all_blocks_soil, all_blocks_other)))
soil_counts <- table(factor(all_blocks_soil, levels = unique_blocks))
other_counts <- table(factor(all_blocks_other, levels = unique_blocks))

block_df <- data.frame(
  spatial_block_id = unique_blocks,
  soil_n = as.integer(soil_counts),
  other_n = as.integer(other_counts),
  stringsAsFactors = FALSE
)
# Deterministic sort: largest soil count first, then other_n, then block_id hash
block_df$hash_order <- vapply(block_df$spatial_block_id, function(b) {
  digest::digest(b, algo = "crc32", serialize = FALSE)
}, character(1))
block_df <- block_df[order(-block_df$soil_n, -block_df$other_n, block_df$hash_order), ]

# Greedy balanced assignment across k = 5 folds
k_folds <- 5L
fold_soil_totals <- rep(0L, k_folds)
fold_other_totals <- rep(0L, k_folds)
block_df$spatial_cv_fold <- integer(nrow(block_df))

for (i in seq_len(nrow(block_df))) {
  # Choose fold with minimum soil_totals, tie-breaking by minimum other_totals, then fold index
  ord <- order(fold_soil_totals, fold_other_totals, seq_len(k_folds))
  best_f <- ord[1]
  block_df$spatial_cv_fold[i] <- best_f
  fold_soil_totals[best_f] <- fold_soil_totals[best_f] + block_df$soil_n[i]
  fold_other_totals[best_f] <- fold_other_totals[best_f] + block_df$other_n[i]
}

fold_lookup <- setNames(block_df$spatial_cv_fold, block_df$spatial_block_id)

# ------------------------------------------------------------------------------
# 3. Build Harmonized AI-Ready Tables
# ------------------------------------------------------------------------------
cat("[3/6] Building harmonized, unit-standardized, QC-flagged tables...\n")

# --- 3A. CONUS Soil Properties (50,116 rows) ---
c_raw <- raw_dfs[["conus_variables.nc"]]
c_site <- make_site_id(c_raw$latitude, c_raw$longitude)
c_block <- make_spatial_block(c_raw$latitude, c_raw$longitude)
c_fold <- unname(fold_lookup[c_block])

soil_prop_cols <- c(
  "Clay", "Sand", "Silt", "Aggregate_stability",
  "Water_retension_3bar", "Water_retension_15_bar",
  "BD_Third_bar", "BD_oven_dry", "TC", "OC", "Estimated_OC",
  "T_nitrogen", "CEC", "CN", "pH"
)
n_valid_props <- as.integer(rowSums(!is.na(c_raw[, soil_prop_cols])))

coord_exact <- sprintf("%.6f_%.6f", c_raw$latitude, c_raw$longitude)
dup_coord_tab <- table(coord_exact)
qc_dup_coords <- as.logical(dup_coord_tab[coord_exact] > 1L)

qc_pct_gt_100 <- (!is.na(c_raw$Clay) & c_raw$Clay > 100) |
  (!is.na(c_raw$Sand) & c_raw$Sand > 100) |
  (!is.na(c_raw$Silt) & c_raw$Silt > 100) |
  (!is.na(c_raw$Aggregate_stability) & c_raw$Aggregate_stability > 100)

tex_sum <- c_raw$Clay + c_raw$Sand + c_raw$Silt
qc_tex_sum <- !is.na(tex_sum) & (tex_sum < 95.0 | tex_sum > 105.0)

qc_wr_inv <- !is.na(c_raw$Water_retension_3bar) &
  !is.na(c_raw$Water_retension_15_bar) &
  (c_raw$Water_retension_15_bar > c_raw$Water_retension_3bar)

qc_oc_gt_tc <- !is.na(c_raw$OC) & !is.na(c_raw$TC) & ((c_raw$OC - c_raw$TC) > 0.5)

qc_any_anomaly <- qc_pct_gt_100 | qc_tex_sum | qc_wr_inv | qc_oc_gt_tc

soil_df <- data.frame(
  obs_id = as.integer(c_raw$obs),
  site_id = c_site,
  latitude = c_raw$latitude,
  longitude = c_raw$longitude,
  depth_top_cm = 0L,
  depth_bottom_cm = 30L,
  dataset_source = as.integer(c_raw$dataset_source),
  spatial_block_id = c_block,
  spatial_cv_fold = as.integer(c_fold),
  clay_pct = c_raw$Clay,
  sand_pct = c_raw$Sand,
  silt_pct = c_raw$Silt,
  aggregate_stability_pct = c_raw$Aggregate_stability,
  water_retention_33kpa_pct = c_raw$Water_retension_3bar,
  water_retention_1500kpa_pct = c_raw$Water_retension_15_bar,
  bulk_density_third_bar_g_cm3 = c_raw$BD_Third_bar,
  bulk_density_oven_dry_g_cm3 = c_raw$BD_oven_dry,
  total_carbon_pct = c_raw$TC,
  organic_carbon_pct = c_raw$OC,
  estimated_organic_carbon_pct = c_raw$Estimated_OC,
  total_nitrogen_pct = c_raw$T_nitrogen,
  cec_cmolc_kg = c_raw$CEC,
  cn_ratio = c_raw$CN,
  ph = c_raw$pH,
  n_valid_properties = n_valid_props,
  qc_duplicate_coords = qc_dup_coords,
  qc_pct_gt_100 = qc_pct_gt_100,
  qc_texture_sum_out_of_range = qc_tex_sum,
  qc_wr_inverted = qc_wr_inv,
  qc_oc_gt_tc = qc_oc_gt_tc,
  qc_any_anomaly = qc_any_anomaly,
  stringsAsFactors = FALSE
)

# --- 3B. Microbial Functional Enzyme Genes Harmonized (805 unique sites) ---
mc <- raw_dfs[["microbe_C_degrading_points.nc"]]
mn <- raw_dfs[["microbe_N_degrading_points.nc"]]
mp <- raw_dfs[["microbe_P_degrading_points.nc"]]

mc$coord_key <- sprintf("%.6f_%.6f", mc$latitude, mc$longitude)
mn$coord_key <- sprintf("%.6f_%.6f", mn$latitude, mn$longitude)
mp$coord_key <- sprintf("%.6f_%.6f", mp$latitude, mp$longitude)

# Merge onto the full union of coordinates (mn has all 805)
all_m_keys <- unique(c(mn$coord_key, mc$coord_key, mp$coord_key))
idx_c <- match(all_m_keys, mc$coord_key)
idx_n <- match(all_m_keys, mn$coord_key)
idx_p <- match(all_m_keys, mp$coord_key)

m_lat <- ifelse(!is.na(idx_n), mn$latitude[idx_n],
                ifelse(!is.na(idx_c), mc$latitude[idx_c], mp$latitude[idx_p]))
m_lon <- ifelse(!is.na(idx_n), mn$longitude[idx_n],
                ifelse(!is.na(idx_c), mc$longitude[idx_c], mp$longitude[idx_p]))
m_block <- make_spatial_block(m_lat, m_lon)

microbe_df <- data.frame(
  microbe_site_idx = seq_along(all_m_keys),
  site_id = make_site_id(m_lat, m_lon),
  latitude = m_lat,
  longitude = m_lon,
  spatial_block_id = m_block,
  spatial_cv_fold = as.integer(unname(fold_lookup[m_block])),
  c_degrading_rel_abundance = mc$avg_freq_scaled[idx_c],
  n_degrading_rel_abundance = mn$avg_freq_scaled[idx_n],
  p_degrading_rel_abundance = mp$avg_freq_scaled[idx_p],
  obs_id_c = as.integer(mc$obs[idx_c]),
  obs_id_n = as.integer(mn$obs[idx_n]),
  obs_id_p = as.integer(mp$obs[idx_p]),
  qc_complete_cnp = !is.na(mc$avg_freq_scaled[idx_c]) &
    !is.na(mn$avg_freq_scaled[idx_n]) &
    !is.na(mp$avg_freq_scaled[idx_p]),
  qc_outside_conus_border = check_outside_conus_border(m_lat, m_lon),
  stringsAsFactors = FALSE
)

# --- 3C. Plant Root Belowground Biomass Harmonized (FRED + TRY = 256 rows) ---
fb <- raw_dfs[["fred_belowground_biomass_points.nc"]]
tb <- raw_dfs[["try_belowground_biomass_points.nc"]]

fb_site <- make_site_id(fb$latitude, fb$longitude)
tb_site <- make_site_id(tb$latitude, tb$longitude)
shared_bio_sites <- intersect(fb_site, tb_site)

fb_dup_within <- duplicated(fb_site) | duplicated(fb_site, fromLast = TRUE)
tb_dup_within <- duplicated(tb_site) | duplicated(tb_site, fromLast = TRUE)

bio_lat <- c(fb$latitude, tb$latitude)
bio_lon <- c(fb$longitude, tb$longitude)
bio_site <- c(fb_site, tb_site)
bio_block <- make_spatial_block(bio_lat, bio_lon)
bio_val <- c(fb$mean_belowground_biomass, tb$mean_belowground_biomass)

# Primary deduplicated record per unique site_id (first valid observation, FRED ordered first)
is_prim_bio <- !is.na(bio_val) & !duplicated(ifelse(is.na(bio_val), NA_character_, bio_site))

root_biomass_df <- data.frame(
  record_id = seq_along(bio_lat),
  site_id = bio_site,
  source_database = c(rep("FRED", nrow(fb)), rep("TRY", nrow(tb))),
  source_obs_id = c(as.integer(fb$obs), as.integer(tb$obs)),
  latitude = bio_lat,
  longitude = bio_lon,
  spatial_block_id = bio_block,
  spatial_cv_fold = as.integer(unname(fold_lookup[bio_block])),
  belowground_biomass_g_m2 = bio_val,
  raw_value = bio_val,
  raw_units = c(rep("unknown", nrow(fb)), rep("g m-2; NA", nrow(tb))),
  n_obs_reported = c(rep(NA_integer_, nrow(fb)), as.integer(tb$n_obs)),
  qc_fred_inferred_unit = c(rep(TRUE, nrow(fb)), rep(FALSE, nrow(tb))),
  qc_within_db_duplicate_coord = c(fb_dup_within, tb_dup_within),
  qc_shared_fred_try_coord = bio_site %in% shared_bio_sites,
  qc_outside_conus_border = check_outside_conus_border(bio_lat, bio_lon),
  is_primary_site_record = is_prim_bio,
  stringsAsFactors = FALSE
)

# --- 3D. Plant Rooting Depth Harmonized (FRED [cm] + TRY [m] = 1,062 rows) ---
fr <- raw_dfs[["fred_rooting_depth_points.nc"]]
tr <- raw_dfs[["try_rooting_depth_points.nc"]]

fr_site <- make_site_id(fr$latitude, fr$longitude)
tr_site <- make_site_id(tr$latitude, tr$longitude)
shared_rd_sites <- intersect(fr_site, tr_site)

fr_dup_within <- duplicated(fr_site) | duplicated(fr_site, fromLast = TRUE)
tr_dup_within <- duplicated(tr_site) | duplicated(tr_site, fromLast = TRUE)

rd_lat <- c(fr$latitude, tr$latitude)
rd_lon <- c(fr$longitude, tr$longitude)
rd_site <- c(fr_site, tr_site)
rd_block <- make_spatial_block(rd_lat, rd_lon)

# Harmonize units: FRED raw is cm -> divide by 100 for m; TRY raw is m -> multiply by 100 for cm
rd_m <- c(fr$mean_rooting_depth / 100.0, tr$mean_rooting_depth)
rd_cm <- c(fr$mean_rooting_depth, tr$mean_rooting_depth * 100.0)
rd_raw <- c(fr$mean_rooting_depth, tr$mean_rooting_depth)

is_prim_rd <- !is.na(rd_m) & !duplicated(ifelse(is.na(rd_m), NA_character_, rd_site))

root_depth_df <- data.frame(
  record_id = seq_along(rd_lat),
  site_id = rd_site,
  source_database = c(rep("FRED", nrow(fr)), rep("TRY", nrow(tr))),
  source_obs_id = c(as.integer(fr$obs), as.integer(tr$obs)),
  latitude = rd_lat,
  longitude = rd_lon,
  spatial_block_id = rd_block,
  spatial_cv_fold = as.integer(unname(fold_lookup[rd_block])),
  rooting_depth_m = rd_m,
  rooting_depth_cm = rd_cm,
  raw_value = rd_raw,
  raw_units = c(rep("cm", nrow(fr)), rep("m", nrow(tr))),
  n_obs_reported = c(rep(NA_integer_, nrow(fr)), as.integer(tr$n_obs)),
  qc_missing_value = is.na(rd_raw),
  qc_within_db_duplicate_coord = c(fr_dup_within, tr_dup_within),
  qc_shared_fred_try_coord = rd_site %in% shared_rd_sites,
  qc_outside_conus_border = check_outside_conus_border(rd_lat, rd_lon),
  is_primary_site_record = is_prim_rd,
  stringsAsFactors = FALSE
)

# Write harmonized Parquet files
arrow::write_parquet(soil_df, "data/parquet/soil_properties_conus.parquet", compression = "snappy")
arrow::write_parquet(microbe_df, "data/parquet/microbe_enzymes_harmonized.parquet", compression = "snappy")
arrow::write_parquet(root_biomass_df, "data/parquet/root_biomass_harmonized.parquet", compression = "snappy")
arrow::write_parquet(root_depth_df, "data/parquet/root_depth_harmonized.parquet", compression = "snappy")

# ------------------------------------------------------------------------------
# 4. Generate Variable Dictionary & Summary Metadata
# ------------------------------------------------------------------------------
cat("[4/6] Generating variable dictionary and dataset metadata...\n")

describe_table <- function(df, table_name, descriptions, units_map, orig_var_map) {
  rows <- lapply(names(df), function(col) {
    v <- df[[col]]
    n_tot <- length(v)
    n_val <- sum(!is.na(v))
    miss_pct <- round(100 * (n_tot - n_val) / n_tot, 2)
    dtype <- class(v)[1]
    if (is.numeric(v) && !is.logical(v) && n_val > 0) {
      vmin <- round(min(v, na.rm = TRUE), 6)
      vmed <- round(median(v, na.rm = TRUE), 6)
      vmax <- round(max(v, na.rm = TRUE), 6)
    } else if (is.logical(v) && n_val > 0) {
      vmin <- 0
      vmed <- round(mean(v, na.rm = TRUE), 4)
      vmax <- 1
    } else {
      vmin <- NA_real_
      vmed <- NA_real_
      vmax <- NA_real_
    }
    data.frame(
      table_name = table_name,
      column_name = col,
      original_nc_variable = ifelse(col %in% names(orig_var_map), orig_var_map[[col]], ""),
      data_type = dtype,
      units = ifelse(col %in% names(units_map), units_map[[col]], "1"),
      total_rows = n_tot,
      valid_count = n_val,
      missing_pct = miss_pct,
      min_or_false = vmin,
      median_or_rate = vmed,
      max_or_true = vmax,
      description = ifelse(col %in% names(descriptions), descriptions[[col]], col),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

soil_desc <- c(
  obs_id = "Observation index (1 to 50,116) matching raw conus_variables.nc",
  site_id = "Deterministic SHA-1 site identifier derived from 6-decimal (latitude, longitude)",
  latitude = "WGS84 Latitude",
  longitude = "WGS84 Longitude",
  depth_top_cm = "Top of soil sampling horizon (0 cm)",
  depth_bottom_cm = "Bottom of soil sampling horizon (30 cm)",
  dataset_source = "Upstream master compilation source identifier (1 or 2)",
  spatial_block_id = "2x2 degree geographic grid cell ID for spatial cross-validation",
  spatial_cv_fold = "Balanced 5-fold spatial cross-validation split (1 to 5, shared across all tables)",
  clay_pct = "Soil clay mass fraction (<0.002 mm)",
  sand_pct = "Soil sand mass fraction (0.05-2.0 mm)",
  silt_pct = "Soil silt mass fraction (0.002-0.05 mm)",
  aggregate_stability_pct = "Soil wet aggregate stability percentage",
  water_retention_33kpa_pct = "Water retention at 1/3 bar (33 kPa, field capacity)",
  water_retention_1500kpa_pct = "Water retention at 15 bar (1500 kPa, permanent wilting point)",
  bulk_density_third_bar_g_cm3 = "Soil bulk density measured at 1/3 bar moisture",
  bulk_density_oven_dry_g_cm3 = "Soil bulk density (oven dry at 105 C)",
  total_carbon_pct = "Total soil carbon concentration",
  organic_carbon_pct = "Measured soil organic carbon concentration",
  estimated_organic_carbon_pct = "Estimated soil organic carbon concentration",
  total_nitrogen_pct = "Total soil nitrogen concentration",
  cec_cmolc_kg = "Cation exchange capacity (numerically equal to meq/100g)",
  cn_ratio = "Soil Carbon-to-Nitrogen ratio (C:N)",
  ph = "Soil pH",
  n_valid_properties = "Count of non-null soil properties available at this observation (0 to 15)",
  qc_duplicate_coords = "TRUE if exact (latitude, longitude) appears more than once in soil_properties_conus",
  qc_pct_gt_100 = "TRUE if Clay, Sand, Silt, or Aggregate_stability exceeds 100%",
  qc_texture_sum_out_of_range = "TRUE if Clay + Sand + Silt is outside [95%, 105%] when all three are present",
  qc_wr_inverted = "TRUE if 15-bar wilting point water retention exceeds 1/3-bar field capacity",
  qc_oc_gt_tc = "TRUE if Organic Carbon exceeds Total Carbon by more than 0.5%",
  qc_any_anomaly = "TRUE if any physical value QC flag (pct_gt_100, texture_sum, wr_inverted, oc_gt_tc) is TRUE"
)
soil_units <- c(
  obs_id = "id", site_id = "id", latitude = "degrees_north", longitude = "degrees_east",
  depth_top_cm = "cm", depth_bottom_cm = "cm", dataset_source = "category",
  spatial_block_id = "id", spatial_cv_fold = "fold (1-5)",
  clay_pct = "%", sand_pct = "%", silt_pct = "%", aggregate_stability_pct = "%",
  water_retention_33kpa_pct = "%", water_retention_1500kpa_pct = "%",
  bulk_density_third_bar_g_cm3 = "g cm-3", bulk_density_oven_dry_g_cm3 = "g cm-3",
  total_carbon_pct = "%", organic_carbon_pct = "%", estimated_organic_carbon_pct = "%",
  total_nitrogen_pct = "%", cec_cmolc_kg = "cmol(+) kg-1", cn_ratio = "1", ph = "1",
  n_valid_properties = "count", qc_duplicate_coords = "boolean", qc_pct_gt_100 = "boolean",
  qc_texture_sum_out_of_range = "boolean", qc_wr_inverted = "boolean",
  qc_oc_gt_tc = "boolean", qc_any_anomaly = "boolean"
)
soil_orig <- c(
  obs_id = "obs", latitude = "latitude", longitude = "longitude", dataset_source = "dataset_source",
  clay_pct = "Clay", sand_pct = "Sand", silt_pct = "Silt",
  aggregate_stability_pct = "Aggregate_stability",
  water_retention_33kpa_pct = "Water_retension_3bar",
  water_retention_1500kpa_pct = "Water_retension_15_bar",
  bulk_density_third_bar_g_cm3 = "BD_Third_bar",
  bulk_density_oven_dry_g_cm3 = "BD_oven_dry",
  total_carbon_pct = "TC", organic_carbon_pct = "OC",
  estimated_organic_carbon_pct = "Estimated_OC",
  total_nitrogen_pct = "T_nitrogen", cec_cmolc_kg = "CEC", cn_ratio = "CN", ph = "pH"
)

microbe_desc <- c(
  microbe_site_idx = "Unique site row index (1 to 805)",
  site_id = "Deterministic SHA-1 site identifier derived from 6-decimal (latitude, longitude)",
  latitude = "WGS84 Latitude",
  longitude = "WGS84 Longitude",
  spatial_block_id = "2x2 degree geographic grid cell ID for spatial cross-validation",
  spatial_cv_fold = "Balanced 5-fold spatial cross-validation split (1 to 5)",
  c_degrading_rel_abundance = "Scaled relative abundance (0-1) of C-degrading microbial enzyme encoding genes",
  n_degrading_rel_abundance = "Scaled relative abundance (0-1) of N-degrading microbial enzyme encoding genes",
  p_degrading_rel_abundance = "Scaled relative abundance (0-1) of P-degrading microbial enzyme encoding genes",
  obs_id_c = "Original obs index in microbe_C_degrading_points.nc",
  obs_id_n = "Original obs index in microbe_N_degrading_points.nc",
  obs_id_p = "Original obs index in microbe_P_degrading_points.nc",
  qc_complete_cnp = "TRUE when all three C, N, and P enzyme abundances are non-null (804 of 805 sites)",
  qc_outside_conus_border = "TRUE if point falls on non-US land outside CONUS borders"
)
microbe_units <- c(
  microbe_site_idx = "id", site_id = "id", latitude = "degrees_north", longitude = "degrees_east",
  spatial_block_id = "id", spatial_cv_fold = "fold (1-5)",
  c_degrading_rel_abundance = "relative abundance (0-1)",
  n_degrading_rel_abundance = "relative abundance (0-1)",
  p_degrading_rel_abundance = "relative abundance (0-1)",
  obs_id_c = "id", obs_id_n = "id", obs_id_p = "id",
  qc_complete_cnp = "boolean", qc_outside_conus_border = "boolean"
)
microbe_orig <- c(
  latitude = "latitude", longitude = "longitude",
  c_degrading_rel_abundance = "avg_freq_scaled (C)",
  n_degrading_rel_abundance = "avg_freq_scaled (N)",
  p_degrading_rel_abundance = "avg_freq_scaled (P)"
)

bio_desc <- c(
  record_id = "Harmonized row identifier (1 to 256)",
  site_id = "Deterministic SHA-1 site identifier derived from 6-decimal (latitude, longitude)",
  source_database = "Trait database source ('FRED' or 'TRY')",
  source_obs_id = "Original obs index in source NetCDF file",
  latitude = "WGS84 Latitude",
  longitude = "WGS84 Longitude",
  spatial_block_id = "2x2 degree geographic grid cell ID for spatial cross-validation",
  spatial_cv_fold = "Balanced 5-fold spatial cross-validation split (1 to 5)",
  belowground_biomass_g_m2 = "Harmonized mean belowground plant biomass in g m-2",
  raw_value = "Unmodified numeric value from source NetCDF",
  raw_units = "Original unit attribute string recorded in source NetCDF",
  n_obs_reported = "Reported observation count per location (provided by TRY; NA for FRED)",
  qc_fred_inferred_unit = "TRUE for FRED records where g m-2 unit was verified via exact cross-database site matches with TRY",
  qc_within_db_duplicate_coord = "TRUE if coordinate appears multiple times within the same source database",
  qc_shared_fred_try_coord = "TRUE if site_id exists in both FRED and TRY belowground biomass files",
  qc_outside_conus_border = "TRUE if coordinate falls in Canada or Mexico outside strict US borders",
  is_primary_site_record = "TRUE for a single canonical non-null record per unique site_id (recommended filter to prevent cross-DB duplicate leakage)"
)
bio_units <- c(
  record_id = "id", site_id = "id", source_database = "category", source_obs_id = "id",
  latitude = "degrees_north", longitude = "degrees_east",
  spatial_block_id = "id", spatial_cv_fold = "fold (1-5)",
  belowground_biomass_g_m2 = "g m-2", raw_value = "g m-2", raw_units = "text",
  n_obs_reported = "count", qc_fred_inferred_unit = "boolean",
  qc_within_db_duplicate_coord = "boolean", qc_shared_fred_try_coord = "boolean",
  qc_outside_conus_border = "boolean", is_primary_site_record = "boolean"
)
bio_orig <- c(
  source_obs_id = "obs", latitude = "latitude", longitude = "longitude",
  belowground_biomass_g_m2 = "mean_belowground_biomass", raw_value = "mean_belowground_biomass",
  n_obs_reported = "n_obs"
)

rd_desc <- c(
  record_id = "Harmonized row identifier (1 to 1,062)",
  site_id = "Deterministic SHA-1 site identifier derived from 6-decimal (latitude, longitude)",
  source_database = "Trait database source ('FRED' or 'TRY')",
  source_obs_id = "Original obs index in source NetCDF file",
  latitude = "WGS84 Latitude",
  longitude = "WGS84 Longitude",
  spatial_block_id = "2x2 degree geographic grid cell ID for spatial cross-validation",
  spatial_cv_fold = "Balanced 5-fold spatial cross-validation split (1 to 5)",
  rooting_depth_m = "Harmonized mean plant rooting depth in meters (FRED cm divided by 100; TRY m as-is)",
  rooting_depth_cm = "Harmonized mean plant rooting depth in centimeters (FRED cm as-is; TRY m multiplied by 100)",
  raw_value = "Unmodified numeric value from source NetCDF in raw_units",
  raw_units = "Original unit attribute string recorded in source NetCDF ('cm' for FRED, 'm' for TRY)",
  n_obs_reported = "Reported observation count per location (provided by TRY; NA for FRED)",
  qc_missing_value = "TRUE if mean_rooting_depth is missing/NA (1 record in TRY with n_obs = 0)",
  qc_within_db_duplicate_coord = "TRUE if coordinate appears multiple times within the same source database",
  qc_shared_fred_try_coord = "TRUE if site_id exists in both FRED and TRY rooting depth files",
  qc_outside_conus_border = "TRUE if coordinate falls in Canada or Mexico outside strict US borders",
  is_primary_site_record = "TRUE for a single canonical non-null record per unique site_id (recommended filter to prevent cross-DB duplicate leakage)"
)
rd_units <- c(
  record_id = "id", site_id = "id", source_database = "category", source_obs_id = "id",
  latitude = "degrees_north", longitude = "degrees_east",
  spatial_block_id = "id", spatial_cv_fold = "fold (1-5)",
  rooting_depth_m = "m", rooting_depth_cm = "cm", raw_value = "raw_units", raw_units = "text",
  n_obs_reported = "count", qc_missing_value = "boolean",
  qc_within_db_duplicate_coord = "boolean", qc_shared_fred_try_coord = "boolean",
  qc_outside_conus_border = "boolean", is_primary_site_record = "boolean"
)
rd_orig <- c(
  source_obs_id = "obs", latitude = "latitude", longitude = "longitude",
  rooting_depth_m = "mean_rooting_depth", rooting_depth_cm = "mean_rooting_depth",
  raw_value = "mean_rooting_depth", n_obs_reported = "n_obs"
)

var_dict <- rbind(
  describe_table(soil_df, "soil_properties_conus", soil_desc, soil_units, soil_orig),
  describe_table(microbe_df, "microbe_enzymes_harmonized", microbe_desc, microbe_units, microbe_orig),
  describe_table(root_biomass_df, "root_biomass_harmonized", bio_desc, bio_units, bio_orig),
  describe_table(root_depth_df, "root_depth_harmonized", rd_desc, rd_units, rd_orig)
)

write.csv(var_dict, "data/metadata/variable_dictionary.csv", row.names = FALSE)

dataset_catalog <- data.frame(
  table_name = c(
    "soil_properties_conus",
    "microbe_enzymes_harmonized",
    "root_biomass_harmonized",
    "root_depth_harmonized"
  ),
  domain = c("Soil Physics & Chemistry", "Microbial Functional Genomics", "Plant Fine/Coarse Root Traits", "Plant Rooting Depth Traits"),
  row_count = c(nrow(soil_df), nrow(microbe_df), nrow(root_biomass_df), nrow(root_depth_df)),
  unique_sites = c(
    length(unique(soil_df$site_id)),
    length(unique(microbe_df$site_id)),
    sum(root_biomass_df$is_primary_site_record),
    sum(root_depth_df$is_primary_site_record)
  ),
  crs = "EPSG:4326 (WGS84)",
  depth_representation = c("0-30 cm topsoil", "Surface soil metagenomes", "Whole-profile belowground biomass", "Maximum / mean rooting depth"),
  parquet_path = c(
    "data/parquet/soil_properties_conus.parquet",
    "data/parquet/microbe_enzymes_harmonized.parquet",
    "data/parquet/root_biomass_harmonized.parquet",
    "data/parquet/root_depth_harmonized.parquet"
  ),
  source_netcdf_files = c(
    "conus_variables.nc",
    "microbe_C_degrading_points.nc; microbe_N_degrading_points.nc; microbe_P_degrading_points.nc",
    "fred_belowground_biomass_points.nc; try_belowground_biomass_points.nc",
    "fred_rooting_depth_points.nc; try_rooting_depth_points.nc"
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 5. Build Indexed Relational SQLite Database
# ------------------------------------------------------------------------------
cat("[5/6] Building indexed SQLite database...\n")

sqlite_path <- "data/sqlite/conus_soil_root_microbe.sqlite"
if (file.exists(sqlite_path)) unlink(sqlite_path)
con <- dbConnect(RSQLite::SQLite(), sqlite_path)

dbWriteTable(con, "dataset_catalog", dataset_catalog, overwrite = TRUE)
dbWriteTable(con, "variable_dictionary", var_dict, overwrite = TRUE)
dbWriteTable(con, "soil_properties_conus", soil_df, overwrite = TRUE)
dbWriteTable(con, "microbe_enzymes_harmonized", microbe_df, overwrite = TRUE)
dbWriteTable(con, "root_biomass_harmonized", root_biomass_df, overwrite = TRUE)
dbWriteTable(con, "root_depth_harmonized", root_depth_df, overwrite = TRUE)

# Create SQL indexes for fast ML fold / coordinate / site lookups
index_stmts <- c(
  "CREATE INDEX idx_soil_fold ON soil_properties_conus(spatial_cv_fold);",
  "CREATE INDEX idx_soil_site ON soil_properties_conus(site_id);",
  "CREATE INDEX idx_soil_coords ON soil_properties_conus(latitude, longitude);",
  "CREATE INDEX idx_microbe_fold ON microbe_enzymes_harmonized(spatial_cv_fold);",
  "CREATE INDEX idx_microbe_site ON microbe_enzymes_harmonized(site_id);",
  "CREATE INDEX idx_bio_fold ON root_biomass_harmonized(spatial_cv_fold);",
  "CREATE INDEX idx_bio_prim ON root_biomass_harmonized(is_primary_site_record);",
  "CREATE INDEX idx_rd_fold ON root_depth_harmonized(spatial_cv_fold);",
  "CREATE INDEX idx_rd_prim ON root_depth_harmonized(is_primary_site_record);",
  # Legacy view matching raw conus_variables.nc column names
  "CREATE VIEW conus_variables_legacy AS
   SELECT
     obs_id AS obs,
     latitude,
     longitude,
     dataset_source,
     clay_pct AS Clay,
     sand_pct AS Sand,
     silt_pct AS Silt,
     aggregate_stability_pct AS Aggregate_stability,
     water_retention_33kpa_pct AS Water_retension_3bar,
     water_retention_1500kpa_pct AS Water_retension_15_bar,
     bulk_density_third_bar_g_cm3 AS BD_Third_bar,
     bulk_density_oven_dry_g_cm3 AS BD_oven_dry,
     total_carbon_pct AS TC,
     organic_carbon_pct AS OC,
     estimated_organic_carbon_pct AS Estimated_OC,
     total_nitrogen_pct AS T_nitrogen,
     cec_cmolc_kg AS CEC,
     cn_ratio AS CN,
     ph AS pH,
     spatial_cv_fold,
     qc_any_anomaly
   FROM soil_properties_conus;"
)
for (sql in index_stmts) dbExecute(con, sql)
dbExecute(con, "ANALYZE;")
dbDisconnect(con)

# ------------------------------------------------------------------------------
# 6. Generate Publication-Quality Overview Figure & Checksums
# ------------------------------------------------------------------------------
cat("[6/6] Rendering spatial coverage overview figure and SHA-256 checksums...\n")

us_map <- map_data("state")

# Subsample CONUS soil slightly only if needed for crisp plotting, or plot all 50,116 with small alpha!
p_soil <- ggplot() +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
               fill = "#f8f9fa", color = "#ced4da", linewidth = 0.25) +
  geom_point(data = soil_df, aes(x = longitude, y = latitude, color = factor(spatial_cv_fold)),
             size = 0.22, alpha = 0.45) +
  coord_quickmap(xlim = c(-125, -66.5), ylim = c(24.5, 49.5)) +
  scale_color_brewer(palette = "Set2", name = "Spatial Fold") +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(
    title = "A. CONUS Soil Properties (n = 50,116; 0-30 cm)",
    subtitle = "Colored by leak-free 2x2\u00b0 5-fold spatial cross-validation split",
    x = NULL, y = "Latitude (\u00b0N)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, color = "#495057"),
    panel.grid.minor = element_blank()
  )

bio_prim <- root_biomass_df[root_biomass_df$is_primary_site_record, ]
rd_prim <- root_depth_df[root_depth_df$is_primary_site_record, ]

lbl_microbe <- sprintf("Microbial C/N/P Genes (n = %s)", format(nrow(microbe_df), big.mark = ","))
lbl_rd <- sprintf("Rooting Depth Sites (n = %s)", format(nrow(rd_prim), big.mark = ","))
lbl_bio <- sprintf("Belowground Biomass Sites (n = %s)", format(nrow(bio_prim), big.mark = ","))

overlay_pts <- rbind(
  data.frame(longitude = microbe_df$longitude, latitude = microbe_df$latitude,
             layer = lbl_microbe),
  data.frame(longitude = rd_prim$longitude, latitude = rd_prim$latitude,
             layer = lbl_rd),
  data.frame(longitude = bio_prim$longitude, latitude = bio_prim$latitude,
             layer = lbl_bio)
)
overlay_pts$layer <- factor(overlay_pts$layer, levels = c(lbl_microbe, lbl_rd, lbl_bio))

p_bio <- ggplot() +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
               fill = "#f8f9fa", color = "#ced4da", linewidth = 0.25) +
  geom_point(data = overlay_pts, aes(x = longitude, y = latitude, color = layer, shape = layer),
             size = 1.35, alpha = 0.85) +
  coord_quickmap(xlim = c(-125, -66.5), ylim = c(24.5, 49.5)) +
  scale_color_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3"),
    name = NULL
  ) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2.5, alpha = 1))) +
  labs(
    title = "B. Microbial Functional Genes & Plant Root Traits",
    subtitle = "Harmonized multi-domain point observations (deduplicated unique sites)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, color = "#495057"),
    panel.grid.minor = element_blank()
  )

# Panel C: Variable completeness bar chart across the 15 soil properties
soil_vars_dict <- var_dict[var_dict$table_name == "soil_properties_conus" &
  var_dict$column_name %in% c(
    "clay_pct", "sand_pct", "silt_pct", "ph", "cec_cmolc_kg",
    "organic_carbon_pct", "water_retention_1500kpa_pct", "total_nitrogen_pct",
    "cn_ratio", "water_retention_33kpa_pct", "bulk_density_third_bar_g_cm3",
    "bulk_density_oven_dry_g_cm3", "total_carbon_pct",
    "estimated_organic_carbon_pct", "aggregate_stability_pct"
  ), ]
soil_vars_dict$completeness_pct <- 100 - soil_vars_dict$missing_pct
soil_vars_dict$label_clean <- c(
  clay_pct = "Clay (%)",
  sand_pct = "Sand (%)",
  silt_pct = "Silt (%)",
  aggregate_stability_pct = "Aggregate Stability (%)",
  water_retention_33kpa_pct = "Water Retention 33 kPa (%)",
  water_retention_1500kpa_pct = "Water Retention 1500 kPa (%)",
  bulk_density_third_bar_g_cm3 = "Bulk Density 1/3 bar (g/cm\u00b3)",
  bulk_density_oven_dry_g_cm3 = "Bulk Density Oven Dry (g/cm\u00b3)",
  total_carbon_pct = "Total Carbon (%)",
  organic_carbon_pct = "Organic Carbon (%)",
  estimated_organic_carbon_pct = "Estimated OC (%)",
  total_nitrogen_pct = "Total Nitrogen (%)",
  cec_cmolc_kg = "CEC (cmol+/kg)",
  cn_ratio = "C:N Ratio",
  ph = "Soil pH"
)[soil_vars_dict$column_name]
soil_vars_dict <- soil_vars_dict[order(soil_vars_dict$completeness_pct), ]
soil_vars_dict$label_clean <- factor(soil_vars_dict$label_clean, levels = soil_vars_dict$label_clean)

p_comp <- ggplot(soil_vars_dict, aes(x = completeness_pct, y = label_clean)) +
  geom_col(fill = "#2b6cb0", width = 0.72) +
  geom_text(
    aes(label = sprintf("%s (%0.1f%%)", format(valid_count, big.mark = ","), completeness_pct)),
    hjust = -0.05, size = 2.8, color = "#2d3748"
  ) +
  scale_x_continuous(limits = c(0, 115), expand = c(0, 0)) +
  labs(
    title = "C. CONUS Soil Property Observation Completeness (Out of 50,116 Sites)",
    x = "Non-Missing Coverage (%)", y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

top_row <- plot_grid(p_soil, p_bio, ncol = 2, rel_widths = c(1, 1))
final_fig <- plot_grid(top_row, p_comp, ncol = 1, rel_heights = c(1.15, 0.85))

ggsave(
  filename = "figures/spatial_coverage_overview.png",
  plot = final_fig,
  width = 12,
  height = 8.5,
  dpi = 220,
  bg = "white"
)

# Compute SHA-256 checksums across all data artifacts
all_data_files <- sort(c(
  Sys.glob("data/parquet/*.parquet"),
  Sys.glob("data/parquet/individual/*.parquet"),
  Sys.glob("data/sqlite/*.sqlite"),
  Sys.glob("data/raw_netcdf/*.nc"),
  "data/metadata/variable_dictionary.csv"
))
sha_vals <- vapply(all_data_files, function(f) {
  digest::digest(file = f, algo = "sha256")
}, character(1))

checksum_lines <- sprintf("%s  %s", sha_vals, all_data_files)
writeLines(checksum_lines, "data/metadata/checksums.sha256")

summary_json <- list(
  name = "CONUS Soil, Plant Root, and Microbial Functional AI-Ready Database",
  version = "1.0.0",
  generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  crs = "EPSG:4326",
  spatial_cv = list(
    scheme = "2.0 x 2.0 degree geographic block cross-validation",
    n_folds = 5L,
    total_unique_blocks = nrow(block_df),
    soil_obs_per_fold = as.list(setNames(as.integer(table(soil_df$spatial_cv_fold)), paste0("fold_", 1:5)))
  ),
  tables = list(
    soil_properties_conus = list(
      rows = nrow(soil_df),
      unique_sites = length(unique(soil_df$site_id)),
      qc_duplicate_coords_rows = sum(soil_df$qc_duplicate_coords),
      qc_any_anomaly_rows = sum(soil_df$qc_any_anomaly),
      qc_pct_gt_100_rows = sum(soil_df$qc_pct_gt_100),
      qc_texture_sum_out_of_range_rows = sum(soil_df$qc_texture_sum_out_of_range),
      qc_wr_inverted_rows = sum(soil_df$qc_wr_inverted),
      qc_oc_gt_tc_rows = sum(soil_df$qc_oc_gt_tc)
    ),
    microbe_enzymes_harmonized = list(
      rows = nrow(microbe_df),
      unique_sites = length(unique(microbe_df$site_id)),
      complete_cnp_sites = sum(microbe_df$qc_complete_cnp),
      qc_outside_conus_border_rows = sum(microbe_df$qc_outside_conus_border)
    ),
    root_biomass_harmonized = list(
      rows = nrow(root_biomass_df),
      fred_rows = sum(root_biomass_df$source_database == "FRED"),
      try_rows = sum(root_biomass_df$source_database == "TRY"),
      unique_primary_sites = sum(root_biomass_df$is_primary_site_record),
      shared_fred_try_rows = sum(root_biomass_df$qc_shared_fred_try_coord),
      qc_outside_conus_border_rows = sum(root_biomass_df$qc_outside_conus_border)
    ),
    root_depth_harmonized = list(
      rows = nrow(root_depth_df),
      fred_rows = sum(root_depth_df$source_database == "FRED"),
      try_rows = sum(root_depth_df$source_database == "TRY"),
      unique_primary_sites = sum(root_depth_df$is_primary_site_record),
      missing_value_rows = sum(root_depth_df$qc_missing_value),
      shared_fred_try_rows = sum(root_depth_df$qc_shared_fred_try_coord),
      qc_outside_conus_border_rows = sum(root_depth_df$qc_outside_conus_border)
    )
  )
)
writeLines(
  jsonlite::toJSON(summary_json, auto_unbox = TRUE, pretty = TRUE),
  "data/metadata/dataset_summary.json"
)

cat("=== Build Complete! ===\n")
