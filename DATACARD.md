# Datasheet / AI Data Card: CONUS Soil, Plant Root & Microbial Functional Database

This document follows the *Datasheets for Datasets* framework (Gebru et al., 2021) to document provenance, harmonization decisions, quality control flags, and machine learning best practices.

---

## 1. Motivation & Intended AI/ML Use Cases

### Purpose
Terrestrial ecosystem and carbon-cycle modeling across the Contiguous United States (CONUS) frequently suffers from fragmented observational data: soil physics/chemistry, plant root traits, and microbial functional genomics are distributed in separate files with inconsistent missing-value sentinels (`-9999`), mismatched physical units (`cm` vs. `m`), unflagged cross-database site duplicates, and no standardized spatial validation splits.

This database unifies **52,239 raw point records** across three coupled belowground domains into clean, self-describing **Apache Parquet** and **SQLite** tables while preserving the original **CF-1.8 NetCDF** files.

### Recommended ML Tasks
1. **Multi-Task Digital Soil Mapping (DSM) & Foundation Models**: Jointly predicting 15 partially observed topsoil (`0–30 cm`) physical and chemical properties using masked multi-task losses (`target_mask` in `src/loaders/python_loader.py`).
2. **Pedotransfer Function (PTF) Learning**: Learning nonlinear relationships among texture (`clay_pct`, `sand_pct`, `silt_pct`), organic carbon (`organic_carbon_pct`), bulk density, cation exchange capacity (`cec_cmolc_kg`), and water retention (`33 kPa` and `1500 kPa`).
3. **Coupled Soil–Root–Microbe Biogeochemical Modeling**: Linking topsoil edaphic properties with plant rooting depth, belowground biomass, and C-, N-, and P-degrading microbial enzyme gene relative abundances.
4. **Spatial Extrapolation & Uncertainty Quantification**: Benchmarking Gaussian Processes, Random Forests, Graph Neural Networks, and Implicit Neural Representations under strict spatial block cross-validation (`spatial_cv_fold`).

### Out-of-Scope Uses
- Subsoil (>30 cm) profile mass-stock accounting directly from `soil_properties_conus` without accounting for depth extrapolation (`soil_properties_conus` represents harmonized `0–30 cm` topsoil).
- Site-level regulatory, engineering, or parcel-scale agronomic decisions without local soil sampling.

---

## 2. Dataset Composition

All coordinates use **WGS84 (`EPSG:4326`)**. Every site is assigned a deterministic 6-decimal coordinate SHA-1 identifier (`site_id`, e.g., `SITE_6b97a66128`) and a 2° × 2° geographic block (`spatial_block_id`).

| Harmonized Table | Domain | Total Rows | Unique Valid Sites | Target Variables | Upstream NetCDF Sources |
| :--- | :--- | ---: | ---: | :--- | :--- |
| `soil_properties_conus` | Topsoil Physics & Chemistry (`0–30 cm`) | 50,116 | 46,764 | 15 edaphic properties | `conus_variables.nc` |
| `microbe_enzymes_harmonized` | Soil Microbial Functional Genomics | 805 | 805 | C, N, P enzyme gene abundance (`0–1`) | `microbe_{C,N,P}_degrading_points.nc` |
| `root_biomass_harmonized` | Plant Belowground Biomass | 256 | 173 | Belowground biomass (`g m-2`) | `fred_belowground_biomass_points.nc`, `try_belowground_biomass_points.nc` |
| `root_depth_harmonized` | Plant Rooting Depth | 1,062 | 596 | Rooting depth (`m` and `cm`) | `fred_rooting_depth_points.nc`, `try_rooting_depth_points.nc` |

---

## 3. Harmonization & Preprocessing Decisions

All preprocessing is deterministic and scripted in [`src/build_database.R`](src/build_database.R). Raw numeric values are **never destructively clipped or deleted**.

1. **Sentinel Fill Value Decoding**:
   - All `-9999` / `-9999.f` NetCDF `_FillValue` entries are converted to native `NULL` / `NA` / `NaN` in Parquet and SQLite. Automated tests (`src/validate_database.R`) verify zero `-9999` values remain across all 12 Parquet tables.
2. **Unit Standardization Across FRED and TRY**:
   - **Rooting Depth**: Raw `fred_rooting_depth_points.nc` stores values in **centimeters (`cm`)** (median `145.5 cm`, max `5300 cm`), whereas `try_rooting_depth_points.nc` stores values in **meters (`m`)** (median `1.24 m`, max `53.0 m`). `root_depth_harmonized` provides both `rooting_depth_m` and `rooting_depth_cm` alongside unmodified `raw_value` and `raw_units`.
   - **Belowground Biomass**: Raw `fred_belowground_biomass_points.nc` recorded `units = "unknown"`. Cross-database coordinate matching against `try_belowground_biomass_points.nc` (`g m-2`) revealed 84 shared site pairs, of which 62 have identical numeric values to `< 0.01` precision (e.g., max `10420.0 g m-2` in both), confirming FRED values are in **`g m-2`**. These rows carry `qc_fred_inferred_unit = TRUE`.
3. **Microbial C/N/P Functional Gene Alignment**:
   - `microbe_C_degrading_points.nc` (`n = 804`) and `microbe_P_degrading_points.nc` (`n = 804`) share identical coordinates; `microbe_N_degrading_points.nc` (`n = 805`) includes those exact 804 locations plus one additional Colorado site (`40.983582, -106.630974`). Merging on exact 6-decimal coordinates produces an 805-row multi-target table with provenance indices `obs_id_c`, `obs_id_n`, and `obs_id_p`.
4. **Column Naming Standardization**:
   - Harmonized tables use explicit `snake_case` with embedded units (e.g., `water_retention_33kpa_pct` for `Water_retension_3bar`, `cec_cmolc_kg` for `CEC` where `1 meq/100g == 1 cmol(+)/kg`). Exact legacy NetCDF column names remain available in `data/parquet/individual/` and via the SQLite view `conus_variables_legacy`.

---

## 4. Preventing Data Leakage in Machine Learning

Geospatial soil and trait compilations exhibit three major leakage risks if split randomly:

1. **Co-located Soil Pedons**: In `soil_properties_conus`, 4,466 rows share exact coordinates across 1,114 multi-observation locations (`qc_duplicate_coords == TRUE`).
2. **Cross-Database Root Trait Overlap (FRED vs. TRY)**: Because FRED and TRY compile overlapping ecological primary literature, **165 of 256 rows (64.5%)** in `root_biomass_harmonized` and **930 of 1,062 rows (87.6%)** in `root_depth_harmonized` share coordinates across both databases (`qc_shared_fred_try_coord == TRUE`). Concatenating FRED + TRY without deduplicating will leak identical sites across splits.
   - **Solution**: Filter `is_primary_site_record == TRUE` (enabled **by default** in `python_loader.py` and `r_loader.R`), which selects one canonical valid record per unique `site_id` (`173` biomass sites; `596` rooting depth sites).
3. **Spatial Autocorrelation Across Nearby Sites & Modalities**:
   - **Solution**: All four tables share a unified **2.0° × 2.0° geographic block grid** (`265` populated blocks) partitioned into 5 balanced folds (`spatial_cv_fold` $\in \{1, 2, 3, 4, 5\}$, with `10,023–10,024` soil rows per fold). Any soil, root, or microbial point falling in the same 2° × 2° block is guaranteed to belong to the **exact same fold**.

---

## 5. Non-Destructive Quality Control (QC) Flags

Rather than silently dropping observations, potential physical edge cases are annotated with boolean flags:

| Table | Flag Column | Flagged Rows | Definition & Guidance |
| :--- | :--- | ---: | :--- |
| `soil_properties_conus` | `qc_duplicate_coords` | 4,466 (8.91%) | Multiple observations at identical `(latitude, longitude)`. Group by `site_id` or `spatial_cv_fold` during validation. |
| `soil_properties_conus` | `qc_pct_gt_100` | 8 (0.02%) | `sand_pct > 100` (5 rows, max `105.6%`) or `aggregate_stability_pct > 100` (3 rows, max `103.8%`), typical of depth-spline or lab rounding artifacts. |
| `soil_properties_conus` | `qc_texture_sum_out_of_range` | 96 (0.19%) | `clay_pct + sand_pct + silt_pct` outside `[95%, 105%]` when all three fractions are reported. |
| `soil_properties_conus` | `qc_wr_inverted` | 240 (0.48%) | 15-bar wilting point water retention exceeds 1/3-bar field capacity (`water_retention_1500kpa_pct > water_retention_33kpa_pct`). |
| `soil_properties_conus` | `qc_oc_gt_tc` | 66 (0.13%) | Measured `organic_carbon_pct` exceeds `total_carbon_pct` by `> 0.5%` (often due to separate analytical methods, e.g., Walkley–Black vs. dry combustion). |
| `soil_properties_conus` | `qc_any_anomaly` | 398 (0.79%) | Union of the 4 physical anomaly flags above (`exclude_qc_anomalies=True` in loaders drops these 398 rows). |
| `root_depth_harmonized` | `qc_missing_value` | 1 (0.09%) | TRY `obs = 559` (`46.00°N, -77.55°W`) has `n_obs = 0` and missing `mean_rooting_depth`. |
| `root_*` / `microbe_*` | `qc_outside_conus_border` | 10 (biomass) / 73 (depth) | Upstream FRED/TRY bounding-box clip included near-border points in southern Canada or northern Mexico. Set `conus_only=True` when joining with CONUS-only rasters (e.g., PRISM, gSSURGO, NLCD, POLARIS). |

> **Note on Water Retention > 100%**: `water_retention_33kpa_pct` reaches up to `837.3%` (`142` rows `> 100%`) and `water_retention_1500kpa_pct` up to `230.9%` (`44` rows `> 100%`). In histosols/peat and high-organic-matter horizons, gravimetric water content (% of oven-dry soil mass) physically exceeds 100% and is therefore **not** flagged as an error unless inverted (`1500 kPa > 33 kPa`).

---

## 6. Maintenance, Reproducibility & Checksums

- **Rebuild Pipeline**: `Rscript src/build_database.R` regenerates all Parquet files, the SQLite database, `variable_dictionary.csv`, `dataset_summary.json`, `checksums.sha256`, and `figures/spatial_coverage_overview.png` from `data/raw_netcdf/`.
- **Automated Verification**: `Rscript src/validate_database.R` and `python3 src/loaders/python_loader.py --self-test` run end-to-end schema, unit, spatial-leakage, and SHA-256 integrity checks.
