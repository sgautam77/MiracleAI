#!/usr/bin/env python3
"""
AI-Ready Python & PyTorch Data Loader for the CONUS Soil, Root & Microbial Database.

Supports:
  - Fast columnar loading from Apache Parquet (via pandas / pyarrow)
  - Zero-dependency fallback via Python's built-in `sqlite3` module
  - Leak-free 5-fold spatial cross-validation splits (`spatial_cv_fold`)
  - Automatic cross-database site deduplication (`is_primary_site_record`)
  - Optional QC filtering (`qc_any_anomaly`, `qc_outside_conus_border`)
  - Multi-task PyTorch `Dataset` with missing-target boolean masks (`target_mask`)
"""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

SOIL_TARGET_COLUMNS: Tuple[str, ...] = (
    "clay_pct",
    "sand_pct",
    "silt_pct",
    "aggregate_stability_pct",
    "water_retention_33kpa_pct",
    "water_retention_1500kpa_pct",
    "bulk_density_third_bar_g_cm3",
    "bulk_density_oven_dry_g_cm3",
    "total_carbon_pct",
    "organic_carbon_pct",
    "estimated_organic_carbon_pct",
    "total_nitrogen_pct",
    "cec_cmolc_kg",
    "cn_ratio",
    "ph",
)

MICROBE_TARGET_COLUMNS: Tuple[str, ...] = (
    "c_degrading_rel_abundance",
    "n_degrading_rel_abundance",
    "p_degrading_rel_abundance",
)

VALID_TABLES: Tuple[str, ...] = (
    "soil_properties_conus",
    "microbe_enzymes_harmonized",
    "root_biomass_harmonized",
    "root_depth_harmonized",
)


class SoilBiogeochemDatabase:
    """Unified interface to the CONUS Soil, Plant Root & Microbial AI-Ready Database."""

    def __init__(self, repo_root: Union[str, Path, None] = None) -> None:
        if repo_root is None:
            # Default to two levels above src/loaders/python_loader.py
            self.repo_root = Path(__file__).resolve().parents[2]
        else:
            self.repo_root = Path(repo_root).resolve()

        self.parquet_dir = self.repo_root / "data" / "parquet"
        self.sqlite_path = self.repo_root / "data" / "sqlite" / "conus_soil_root_microbe.sqlite"
        self.summary_path = self.repo_root / "data" / "metadata" / "dataset_summary.json"

        if not self.sqlite_path.exists() and not self.parquet_dir.exists():
            raise FileNotFoundError(
                f"Could not find database files under {self.repo_root / 'data'}. "
                "Run `Rscript src/build_database.R` or check `repo_root`."
            )

    def summary(self) -> Dict[str, Any]:
        """Return dataset-level metadata and summary counts."""
        with open(self.summary_path, "r", encoding="utf-8") as f:
            return json.load(f)

    def load_table(
        self,
        table_name: str = "soil_properties_conus",
        *,
        deduplicate_sites: bool = True,
        exclude_qc_anomalies: bool = False,
        conus_only: bool = False,
        backend: str = "auto",
    ) -> Any:
        """
        Load a harmonized table as a pandas DataFrame (if pandas is installed)
        or a list of row dicts (stdlib sqlite3 fallback).

        Parameters
        ----------
        table_name : str
            One of 'soil_properties_conus', 'microbe_enzymes_harmonized',
            'root_biomass_harmonized', or 'root_depth_harmonized'.
        deduplicate_sites : bool, default True
            For root trait tables combining FRED and TRY, filter to
            `is_primary_site_record == True` to prevent cross-DB duplicate leakage.
        exclude_qc_anomalies : bool, default False
            If True, drop rows flagged by `qc_any_anomaly` (soil) or `qc_missing_value` (root depth).
        conus_only : bool, default False
            If True, filter out near-border Canada/Mexico records (`qc_outside_conus_border == False`).
        backend : {'auto', 'parquet', 'sqlite'}, default 'auto'
            Storage backend to read from.
        """
        if table_name not in VALID_TABLES:
            raise ValueError(f"Unknown table '{table_name}'. Expected one of {VALID_TABLES}.")

        has_pandas = False
        try:
            import pandas as pd  # type: ignore

            has_pandas = True
        except ImportError:
            pd = None  # type: ignore

        if backend == "auto":
            if has_pandas:
                try:
                    import pyarrow  # noqa: F401 # type: ignore

                    backend = "parquet"
                except ImportError:
                    backend = "sqlite"
            else:
                backend = "sqlite"

        if backend == "parquet":
            if not has_pandas:
                raise ImportError("pandas and pyarrow are required for backend='parquet'.")
            pq_path = self.parquet_dir / f"{table_name}.parquet"
            df = pd.read_parquet(pq_path)
            return self._filter_dataframe(
                df,
                deduplicate_sites=deduplicate_sites,
                exclude_qc_anomalies=exclude_qc_anomalies,
                conus_only=conus_only,
            )

        # SQLite path (works with pandas or pure stdlib)
        where_clauses: List[str] = []
        if deduplicate_sites and table_name in ("root_biomass_harmonized", "root_depth_harmonized"):
            where_clauses.append("is_primary_site_record = 1")
        if exclude_qc_anomalies:
            if table_name == "soil_properties_conus":
                where_clauses.append("qc_any_anomaly = 0")
            elif table_name == "root_depth_harmonized":
                where_clauses.append("qc_missing_value = 0")
        if conus_only and table_name != "soil_properties_conus":
            where_clauses.append("qc_outside_conus_border = 0")

        where_sql = f" WHERE {' AND '.join(where_clauses)}" if where_clauses else ""
        query = f"SELECT * FROM {table_name}{where_sql};"

        with sqlite3.connect(self.sqlite_path) as conn:
            if has_pandas:
                return pd.read_sql_query(query, conn)
            conn.row_factory = sqlite3.Row
            cur = conn.execute(query)
            return [dict(row) for row in cur.fetchall()]

    @staticmethod
    def _filter_dataframe(
        df: Any,
        *,
        deduplicate_sites: bool,
        exclude_qc_anomalies: bool,
        conus_only: bool,
    ) -> Any:
        mask = pd_true_mask(len(df), df.index)
        if deduplicate_sites and "is_primary_site_record" in df.columns:
            mask &= df["is_primary_site_record"].astype(bool)
        if exclude_qc_anomalies:
            if "qc_any_anomaly" in df.columns:
                mask &= ~df["qc_any_anomaly"].astype(bool)
            if "qc_missing_value" in df.columns:
                mask &= ~df["qc_missing_value"].astype(bool)
        if conus_only and "qc_outside_conus_border" in df.columns:
            mask &= ~df["qc_outside_conus_border"].astype(bool)
        return df.loc[mask].reset_index(drop=True)

    def get_spatial_split(
        self,
        table_name: str = "soil_properties_conus",
        *,
        test_fold: int = 1,
        val_fold: Optional[int] = None,
        deduplicate_sites: bool = True,
        exclude_qc_anomalies: bool = False,
        conus_only: bool = False,
    ) -> Union[Tuple[Any, Any], Tuple[Any, Any, Any]]:
        """
        Split a table using the precomputed 2x2-degree `spatial_cv_fold` (1..5)
        to prevent spatial autocorrelation leakage between train, val, and test sets.
        """
        if test_fold not in (1, 2, 3, 4, 5):
            raise ValueError("test_fold must be in {1, 2, 3, 4, 5}.")
        if val_fold is not None and (val_fold not in (1, 2, 3, 4, 5) or val_fold == test_fold):
            raise ValueError("val_fold must be in {1..5} and distinct from test_fold.")

        data = self.load_table(
            table_name,
            deduplicate_sites=deduplicate_sites,
            exclude_qc_anomalies=exclude_qc_anomalies,
            conus_only=conus_only,
        )

        # Pandas DataFrame branch
        if hasattr(data, "loc"):
            test_df = data.loc[data["spatial_cv_fold"] == test_fold].reset_index(drop=True)
            if val_fold is None:
                train_df = data.loc[data["spatial_cv_fold"] != test_fold].reset_index(drop=True)
                return train_df, test_df
            val_df = data.loc[data["spatial_cv_fold"] == val_fold].reset_index(drop=True)
            train_df = data.loc[
                ~data["spatial_cv_fold"].isin([test_fold, val_fold])
            ].reset_index(drop=True)
            return train_df, val_df, test_df

        # Pure Python list-of-dicts branch
        test_rows = [r for r in data if r["spatial_cv_fold"] == test_fold]
        if val_fold is None:
            train_rows = [r for r in data if r["spatial_cv_fold"] != test_fold]
            return train_rows, test_rows
        val_rows = [r for r in data if r["spatial_cv_fold"] == val_fold]
        train_rows = [r for r in data if r["spatial_cv_fold"] not in (test_fold, val_fold)]
        return train_rows, val_rows, test_rows


def pd_true_mask(length: int, index: Any) -> Any:
    import pandas as pd  # type: ignore

    return pd.Series([True] * length, index=index, dtype=bool)


def build_pytorch_dataset(
    df: Any,
    target_cols: Sequence[str] = SOIL_TARGET_COLUMNS,
    feature_cols: Sequence[str] = ("latitude", "longitude"),
) -> Any:
    """
    Wrap a loaded DataFrame into a PyTorch Dataset for single- or multi-task learning.
    Returns dictionaries with:
      - 'features': float32 tensor of shape (F,)
      - 'targets': float32 tensor of shape (T,) with NaNs replaced by 0.0
      - 'target_mask': bool tensor of shape (T,) True where ground-truth is observed
    """
    import numpy as np  # type: ignore
    import torch  # type: ignore
    from torch.utils.data import Dataset  # type: ignore

    class MultiTaskGeospatialDataset(Dataset):
        def __init__(self, frame: Any) -> None:
            x_np = frame[list(feature_cols)].to_numpy(dtype=np.float32, copy=True)
            y_raw = frame[list(target_cols)].to_numpy(dtype=np.float32, copy=True)
            mask_np = ~np.isnan(y_raw)
            y_clean = np.where(mask_np, y_raw, 0.0)

            self.features = torch.from_numpy(x_np)
            self.targets = torch.from_numpy(y_clean)
            self.target_mask = torch.from_numpy(mask_np)
            self.folds = torch.from_numpy(
                frame["spatial_cv_fold"].to_numpy(dtype=np.int64, copy=True)
            )

        def __len__(self) -> int:
            return int(self.features.shape[0])

        def __getitem__(self, idx: int) -> Dict[str, Any]:
            return {
                "features": self.features[idx],
                "targets": self.targets[idx],
                "target_mask": self.target_mask[idx],
                "spatial_cv_fold": self.folds[idx],
            }

    return MultiTaskGeospatialDataset(df)


def _run_cli() -> None:
    parser = argparse.ArgumentParser(description="Inspect and test the AI-Ready Soil Database.")
    parser.add_argument("--info", action="store_true", help="Print dataset summary JSON.")
    parser.add_argument("--self-test", action="store_true", help="Run loader integrity checks.")
    args = parser.parse_args()

    db = SoilBiogeochemDatabase()
    if args.info or not args.self_test:
        print(json.dumps(db.summary(), indent=2))

    if args.self_test:
        print("\nRunning Python loader self-test...")
        for tbl in VALID_TABLES:
            train, val, test = db.get_spatial_split(tbl, test_fold=1, val_fold=2)
            n_tr, n_va, n_te = len(train), len(val), len(test)
            assert n_tr > 0 and n_va > 0 and n_te > 0, f"Empty split in {tbl}"
            print(f"  [OK] {tbl:<28} -> train={n_tr:5d}, val={n_va:5d}, test={n_te:5d}")
        print("All Python loader checks passed!")


if __name__ == "__main__":
    _run_cli()
