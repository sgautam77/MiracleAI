# ==============================================================================
# AI-Ready R Data Loader for CONUS Soil, Plant Root & Microbial Database
# ==============================================================================
# Supports loading via Apache Parquet (`arrow`) or SQLite (`DBI` + `RSQLite`),
# spatial cross-validation splits (`spatial_cv_fold`), cross-database site
# deduplication, QC filtering, and optional conversion to `sf` spatial objects.
# ==============================================================================

load_ai_table <- function(
  table_name = "soil_properties_conus",
  repo_root = ".",
  deduplicate_sites = TRUE,
  exclude_qc_anomalies = FALSE,
  conus_only = FALSE,
  as_sf = FALSE,
  backend = c("parquet", "sqlite")
) {
  valid_tables <- c(
    "soil_properties_conus",
    "microbe_enzymes_harmonized",
    "root_biomass_harmonized",
    "root_depth_harmonized"
  )
  if (!(table_name %in% valid_tables)) {
    stop("Unknown table_name: ", table_name, ". Choose one of: ", paste(valid_tables, collapse = ", "))
  }
  backend <- match.arg(backend)

  if (backend == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Package 'arrow' is required for backend = 'parquet'. Use backend = 'sqlite' instead.")
    }
    pq_path <- file.path(repo_root, "data", "parquet", paste0(table_name, ".parquet"))
    df <- as.data.frame(arrow::read_parquet(pq_path))
  } else {
    if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RSQLite", quietly = TRUE)) {
      stop("Packages 'DBI' and 'RSQLite' are required for backend = 'sqlite'.")
    }
    db_path <- file.path(repo_root, "data", "sqlite", "conus_soil_root_microbe.sqlite")
    con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    df <- DBI::dbReadTable(con, table_name)
  }

  if (isTRUE(deduplicate_sites) && "is_primary_site_record" %in% names(df)) {
    df <- df[as.logical(df$is_primary_site_record), , drop = FALSE]
  }
  if (isTRUE(exclude_qc_anomalies)) {
    if ("qc_any_anomaly" %in% names(df)) {
      df <- df[!as.logical(df$qc_any_anomaly), , drop = FALSE]
    }
    if ("qc_missing_value" %in% names(df)) {
      df <- df[!as.logical(df$qc_missing_value), , drop = FALSE]
    }
  }
  if (isTRUE(conus_only) && "qc_outside_conus_border" %in% names(df)) {
    df <- df[!as.logical(df$qc_outside_conus_border), , drop = FALSE]
  }

  rownames(df) <- NULL

  if (isTRUE(as_sf)) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required when as_sf = TRUE.")
    }
    df <- sf::st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  }
  return(df)
}

get_spatial_cv_split <- function(
  table_name = "soil_properties_conus",
  test_fold = 1L,
  val_fold = NULL,
  repo_root = ".",
  ...
) {
  if (!(test_fold %in% 1:5)) stop("test_fold must be an integer in 1:5.")
  df <- load_ai_table(table_name = table_name, repo_root = repo_root, ...)

  if (is.null(val_fold)) {
    return(list(
      train = df[df$spatial_cv_fold != test_fold, , drop = FALSE],
      test = df[df$spatial_cv_fold == test_fold, , drop = FALSE]
    ))
  }
  if (!(val_fold %in% 1:5) || val_fold == test_fold) {
    stop("val_fold must be in 1:5 and distinct from test_fold.")
  }
  list(
    train = df[!(df$spatial_cv_fold %in% c(test_fold, val_fold)), , drop = FALSE],
    val = df[df$spatial_cv_fold == val_fold, , drop = FALSE],
    test = df[df$spatial_cv_fold == test_fold, , drop = FALSE]
  )
}

get_variable_dictionary <- function(repo_root = ".") {
  read.csv(file.path(repo_root, "data", "metadata", "variable_dictionary.csv"), stringsAsFactors = FALSE)
}
