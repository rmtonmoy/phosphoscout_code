#!/usr/bin/env python
"""
MCP server exposing CELLxGENE Census differential expression tools.

Tools:

census_list_obs_columns
  List relevant column names in the Census obs (cell metadata) schema.
  Returns a curated list of commonly used columns: cell_type, disease,
  is_primary_data, and tissue_general. Use this to discover correct column names
  for filter strings.

census_get_filter_options
  Return possible values for high-level filters (tissue_general, cell_type, disease)
  from the Census obs metadata. Use this to discover valid filter values
  before constructing group1_filter and group2_filter for census_diffexp.

census_get_disease_for_tissue_and_cell_type
  Given tissue_general and cell_type filters, return which disease values have
  non-zero cell counts. Use this to discover valid disease options for a specific
  tissue and cell type combination.

census_get_cell_type_for_tissue_and_disease
  Given tissue_general and disease filters, return which cell_type values have
  non-zero cell counts. Use this to discover valid cell type options for a specific
  tissue and disease combination.

census_get_tissue_for_cell_type_and_disease
  Given cell_type and disease filters, return which tissue_general values have
  non-zero cell counts. Use this to discover valid tissue options for a specific
  cell type and disease combination.

census_diffexp
  Run differential expression analysis between two cell groups using Welch's t-test
  and Cohen's d on ln(CPTT+1) normalized expression. Filters low-coverage cells and
  requires explicit gene lists for efficient querying. Use census_get_filter_options
  to discover valid values for group1_filter and group2_filter parameters.

Supports:
- census_version (e.g. "latest", "stable", "2025-01-30", "2024-07-01")
- optional local `uri` for local SOMA copies
- SOMA value_filter syntax for cell selection

Example filter: "tissue_general == 'lung' and cell_type == 'B cell' and is_primary_data == True"

Note: If you encounter errors about columns not existing, use census_list_obs_columns
to see the relevant column names. Common mistakes include using 'donor_sex' instead
of 'sex', or other similar naming variations.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Set

import numpy as np
import pandas as pd
from fastmcp import FastMCP
from pydantic import BaseModel, model_validator
from scipy import sparse, stats
import cellxgene_census

mcp = FastMCP("census_diffexp_mcp")


# --------------------------
# Helpers
# --------------------------

ORG_ALIASES = {
    "homo sapiens": "Homo sapiens",
    "Homo sapiens": "Homo sapiens",
    "human": "Homo sapiens",
    "mus musculus": "Mus musculus",
    "Mus musculus": "Mus musculus",
    "mouse": "Mus musculus",
}


def _normalize_organism(organism: str) -> str:
    return ORG_ALIASES.get(organism, organism)


def _organism_to_summary_format(organism: str) -> str:
    """Convert normalized organism name to the format used in summary_cell_counts table.
    
    The summary table uses lowercase with underscores (e.g., 'homo_sapiens'),
    while normalized names use title case with spaces (e.g., 'Homo sapiens').
    """
    org = _normalize_organism(organism)
    return org.lower().replace(" ", "_")


def _open_census(
    census_version: str = "latest",
    uri: Optional[str] = None,
):
    """
    Open Census with either a specific census_version or a direct URI.

    - If `uri` is provided, it takes precedence (local / S3 / custom).
    - Otherwise uses `census_version` (e.g. "latest", "2025-01-30").
    """
    if uri:
        return cellxgene_census.open_soma(uri=uri)
    else:
        return cellxgene_census.open_soma(census_version=census_version)


def _add_is_primary_clause(value_filter: Optional[str], only_primary: bool) -> Optional[str]:
    """
    Ensure we always restrict to is_primary_data == True when requested,
    unless the user explicitly mentions is_primary_data in their filter.
    """
    if not only_primary:
        return value_filter

    if value_filter and "is_primary_data" in value_filter:
        # user is explicitly handling it
        return value_filter

    if value_filter and value_filter.strip():
        return f"({value_filter}) and is_primary_data == True"
    else:
        return "is_primary_data == True"


def _bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """
    Benjamini–Hochberg FDR correction.

    Parameters
    ----------
    p_values : array-like of p-values

    Returns
    -------
    q_values : np.ndarray of adjusted p-values
    """
    p = np.asarray(p_values, dtype=float)
    n = p.size
    if n == 0:
        return p

    order = np.argsort(p)
    ranked = p[order]
    # q_i = p_i * n / i
    q = ranked * n / (np.arange(1, n + 1))
    # enforce monotone non-increasing from tail
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    result = np.empty_like(q)
    result[order] = q
    return result

def _filter_low_coverage_cells(adata, min_genes: int = 500):
    """
    Filter out cells with fewer than `min_genes` expressed genes.

    Uses obs["nnz"] when available (Census schema: number of non-zero genes).
    Falls back to computing (X > 0).sum(axis=1) if needed.
    """
    import numpy as _np

    if "nnz" in adata.obs.columns:
        nnz = adata.obs["nnz"].to_numpy()
    elif "n_measured_vars" in adata.obs.columns:
        nnz = adata.obs["n_measured_vars"].to_numpy()
    else:
        X = adata.X
        if sparse.issparse(X):
            nnz = _np.asarray((X > 0).sum(axis=1)).ravel()
        else:
            nnz = (X > 0).sum(axis=1)

    nnz = nnz.astype(int)
    mask = nnz >= min_genes

    return adata[mask].copy()


def _count_cells_with_min_genes(
    census,
    organism: str,
    obs_value_filter: Optional[str],
    min_genes: int = 500,
) -> Optional[int]:
    column_names = ["nnz", "n_measured_vars"]
    df = cellxgene_census.get_obs(
        census,
        organism,
        value_filter=obs_value_filter,
        column_names=column_names,
    )
    if "nnz" in df.columns:
        nnz = df["nnz"].to_numpy()
    elif "n_measured_vars" in df.columns:
        nnz = df["n_measured_vars"].to_numpy()
    else:
        return None
    return int((nnz.astype(int) >= min_genes).sum())


def _prepare_expression_matrix(adata) -> sparse.spmatrix | np.ndarray:
    """
    Take AnnData.X and return ln(CPTT + 1) per cell:

      - CPTT = counts per 10,000, using per-cell library size from obs["raw_sum"]
        when available (Census schema), otherwise fall back to sum of X.
      - Then apply natural log: ln(CPTT + 1).

    For sparse X, all operations are performed in sparse form:
      - Row scaling via multiply(factors[:, None])
      - log1p applied only to non-zero entries (X.data)
    """
    X = adata.X

    # Normalize type and format
    if sparse.issparse(X):
        X = X.tocsr().astype(np.float32)
    else:
        X = np.asarray(X, dtype=np.float32)

    # 1) Per-cell library sizes
    if "raw_sum" in adata.obs.columns:
        libsizes = adata.obs["raw_sum"].to_numpy().astype(np.float32)
    else:
        if sparse.issparse(X):
            libsizes = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
        else:
            libsizes = X.sum(axis=1).astype(np.float32)

    # 2) Convert counts to CPTT (counts per 10,000)
    scale = 1e4
    libsizes = np.clip(libsizes, 1e-8, None)
    factors = (scale / libsizes).astype(np.float32)

    if sparse.issparse(X):
        # Row-wise scaling while staying sparse
        X = X.multiply(factors[:, None])
        # 3) ln(CPTT + 1) on non-zero entries only
        X.data = np.log1p(X.data)
        return X  # CSR sparse matrix
    else:
        # Dense path (typically small slices)
        X *= factors[:, None]
        X = np.log1p(X)
        return X



def _get_gene_names(adata, gene_field: str = "feature_name") -> np.ndarray:
    """
    Return an array of gene names for the var axis.

    Prefer adata.var[gene_field] if present, otherwise fall back to var.index.
    """
    var_df = adata.var
    if gene_field in var_df.columns:
        return var_df[gene_field].astype(str).values
    else:
        return var_df.index.astype(str).values


def _build_var_value_filter(genes: List[str], gene_field: str = "feature_name") -> str:
    """
    Build a var_value_filter string for SOMA queries to filter genes at fetch time.
    
    Parameters
    ----------
    genes : List[str]
        List of gene names to filter.
    gene_field : str
        Column name in var to match against (default: "feature_name").
    
    Returns
    -------
    str
        SOMA value_filter string, e.g. "feature_name in ['AQP5', 'TUBB4B']".
    """
    if not genes:
        raise ValueError("genes list cannot be empty")
    
    genes_quoted = [f"'{g}'" for g in genes]
    return f"{gene_field} in [{', '.join(genes_quoted)}]"


def _extract_column_names_from_filter(value_filter: Optional[str]) -> Set[str]:
    """
    Extract column names from a SOMA value_filter string.
    
    This is a simple parser that looks for patterns like "column_name == ..."
    or "column_name in ...". It's not a full parser but handles common cases.
    
    Parameters
    ----------
    value_filter : Optional[str]
        SOMA value_filter string, e.g. "tissue_general == 'lung' and cell_type == 'B cell'".
    
    Returns
    -------
    Set[str]
        Set of column names mentioned in the filter.
    """
    if not value_filter:
        return set()
    
    # Pattern to match column names before ==, !=, in, not in, etc.
    # Matches: word characters, underscores, dots before comparison operators
    pattern = r'([a-zA-Z_][a-zA-Z0-9_]*(?:\.[a-zA-Z_][a-zA-Z0-9_]*)?)\s*(?:==|!=|in|not in|>|<|>=|<=)'
    matches = re.findall(pattern, value_filter)
    
    # Also handle cases where column names might be at the start of the filter
    # or after "and"/"or"
    pattern2 = r'(?:^|\s+and\s+|\s+or\s+)([a-zA-Z_][a-zA-Z0-9_]*(?:\.[a-zA-Z_][a-zA-Z0-9_]*)?)\s*(?:==|!=|in|not in|>|<|>=|<=)'
    matches2 = re.findall(pattern2, value_filter)
    
    all_matches = set(matches) | set(matches2)
    return all_matches


def _get_relevant_obs_columns() -> List[str]:
    """
    Return a constant list of relevant obs column names.
    
    Returns
    -------
    List[str]
        List of relevant column names: cell_type, disease, is_primary_data, tissue_general.
    """
    return ["cell_type", "disease", "is_primary_data", "tissue_general"]


# --------------------------
# Tools
# --------------------------

def _census_list_obs_values(
    organism: str,
    column: str,
    value_filter: Optional[str] = None,
    census_version: str = "latest",
    uri: Optional[str] = None,
    max_values: int = 100,
) -> Dict[str, Any]:
    """
    List distinct values (with counts) for a given obs column in the Census.

    Use this to discover valid values for filters like tissue_general, cell_type, disease, etc.

    Parameters
    ----------
    organism : str
        "Homo sapiens" or "Mus musculus" (case-insensitive allowed for 'human'/'mouse').
    column : str
        obs column name, e.g. "tissue_general", "cell_type", "disease", "sex".
    value_filter : Optional[str]
        SOMA value_filter for obs, e.g. "tissue_general == 'lung' and sex == 'female'".
        If None, all cells are considered.
    census_version : str
        Census version, e.g. "latest", "2025-01-30", "2024-07-01".
    uri : Optional[str]
        Optional SOMA URI for a local / S3 copy instead of a public release.
    max_values : int
        Maximum number of distinct values to return, sorted by decreasing count.

    Returns
    -------
    dict with:
      - column: str
      - organism: str
      - census_version: str
      - n_cells_scanned: int
      - values: list[{value: str | null, count: int}]
    """
    org = _normalize_organism(organism)

    with _open_census(census_version=census_version, uri=uri) as census:
        df = cellxgene_census.get_obs(
            census,
            org,
            value_filter=value_filter,
            column_names=[column],
        )

    n_cells = int(df.shape[0])
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in obs for organism '{org}'.")

    vc = df[column].value_counts(dropna=False).sort_values(ascending=False)

    values = []
    for val, count in vc.head(max_values).items():
        if pd.isna(val):
            display_val = None
        else:
            display_val = str(val)
        values.append({"value": display_val, "count": int(count)})

    return {
        "organism": org,
        "census_version": census_version,
        "column": column,
        "n_cells_scanned": n_cells,
        "values": values,
    }


def _census_count_cells(
    organism: str,
    value_filter: Optional[str] = None,
    census_version: str = "latest",
    uri: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Count how many cells match an obs filter, without pulling expression data.

    Use this before running DE to avoid accidentally selecting huge groups.

    Parameters
    ----------
    organism : str
        "Homo sapiens" or "Mus musculus".
    value_filter : Optional[str
        SOMA value_filter for obs, e.g.
        "tissue_general == 'lung' and cell_type == 'B cell' and is_primary_data == True".
    census_version : str
        Census version ("latest", "2025-01-30", etc.).
    uri : Optional[str]
        Optional SOMA URI for a local / S3 copy.

    Returns
    -------
    dict with:
      - organism
      - census_version
      - value_filter
      - n_cells
    """
    org = _normalize_organism(organism)

    with _open_census(census_version=census_version, uri=uri) as census:
        df = cellxgene_census.get_obs(
            census,
            org,
            value_filter=value_filter,
            column_names=["soma_joinid"],
        )

    return {
        "organism": org,
        "census_version": census_version,
        "value_filter": value_filter,
        "n_cells": int(df.shape[0]),
    }


def _census_diffexp_impl(
    organism: str,
    group1_filter: str,
    group2_filter: str,
    genes: List[str],
    census_version: str = "2025-11-08",
    uri: Optional[str] = None,
    gene_field: str = "feature_name",
    only_primary_data: bool = True,
    lfc_threshold: float = 0.0,
    min_effect_size: float = 0.0,
    top_n: int = 50,
) -> Dict[str, Any]:
    """
    Run differential expression for two cell groups in CELLxGENE Census.

    What this function DOES
    -----------------------
    - Cell selection:
      * Applies the user-provided SOMA obs filters for group1 and group2
        (e.g. "tissue_general == 'blood' and cell_type == 'B cell'").
      * If `only_primary_data=True`, automatically appends
        "and is_primary_data == True" to each group filter unless the user
        has already mentioned `is_primary_data`.
      * After fetching each AnnData slice, drops low-coverage cells with
        fewer than a minimum number of expressed genes (e.g. < 500), so we
        roughly match the low-coverage filtering described in the
        CELLxGENE Gene Expression documentation.

    - Gene selection:
      * Requires an explicit `genes` list and constructs a `var_value_filter`
        so only those genes are fetched from Census. This avoids downloading
        the full gene matrix and keeps memory usage manageable.

    - Expression transform:
      * Starts from raw counts.
      * Performs simple library-size normalization to counts-per-10,000
        (CPTT) per cell, using per-cell library size (e.g. obs["raw_sum"]
        when available).
      * Applies a log1p transform to obtain ln(CPTT + 1).
      * All downstream statistics (means, variances, effect sizes, tests)
        are computed on this ln(CPTT + 1) scale.
      * For sparse matrices, all operations stay in sparse space and only
        non-zero entries are transformed.

    - Differential expression statistics:
      * Uses Welch's t-test (unequal variance) on ln(CPTT + 1) values,
        implemented from summary statistics (means/variances) to avoid
        densifying X.
      * Effect size is Cohen's d on the same scale.
      * `log_fold_change` is the difference in group means on the
        ln(CPTT + 1) scale, so positive values mean higher expression
        in group1 and negative values mean higher expression in group2.
      * P-values are Benjamini–Hochberg corrected to produce q-values.
      * Results can be further filtered by `lfc_threshold`,
        `min_effect_size`, and `top_n`.

    What this function does NOT do (and why)
    ----------------------------------------
    Compared to the full CELLxGENE Gene Expression processing pipeline,
    this implementation is intentionally simpler:

    - It does NOT:
      * Restrict to a curated set of sequencing assays or perform
        assay-specific handling; all assays that pass the obs filters
        are treated uniformly.
      * Apply gene-length pre-normalization for full-length assays
        (e.g. Smart-seq) before normalization.
      * Mask ultra-low normalized expression values (e.g. setting
        ln(CPTT + 1) ≤ threshold to missing), or compute means only over
        “expressing” cells; group means are taken over all cells kept
        after the low-coverage filter.

    These additional steps are omitted to keep the MCP tool fast,
    lightweight, and easy to understand. As a result, the direction and
    qualitative ranking of differential expression should generally
    agree with the official CELLxGENE Gene Expression UI, but absolute
    effect sizes and fold changes will not be numerically identical.

    Parameters
    ----------
    organism : str
        "Homo sapiens" or "Mus musculus" (case-insensitive 'human'/'mouse'
        also ok).
    group1_filter : str
        SOMA value_filter for group 1 obs, e.g.
        "tissue_general == 'lung' and cell_type == 'B cell'".
    group2_filter : str
        SOMA value_filter for group 2 obs.
    genes : List[str]
        Required list of gene names to analyze (matched against `gene_field`,
        typically "feature_name", i.e. gene symbols). Genes are filtered at
        fetch time to avoid downloading unnecessary data.
    census_version : str
        Census version, e.g. "latest", "2025-01-30", "2024-07-01".
    uri : Optional[str]
        Optional SOMA URI for local / S3 Census (overrides census_version).
    gene_field : str
        Column in adata.var to use as gene name (default: "feature_name").
        If not present, falls back to var.index.
    only_primary_data : bool
        If True (default), automatically append "and is_primary_data == True"
        to both group filters unless they already mention is_primary_data.
    lfc_threshold : float
        Minimum absolute log fold change (on the ln(CPTT + 1) scale) to keep
        in output (after t-test).
    min_effect_size : float
        Minimum absolute Cohen's d to keep in output.
    top_n : int
        Return only the top N genes by |effect_size| after filtering.

    Returns
    -------
    dict with:
      - organism
      - census_version
      - group1: {filter, n_cells}
      - group2: {filter, n_cells}
      - parameters: {...}
      - results: list of:
          {
            "gene": str,
            "log_fold_change": float,
            "effect_size": float,
            "t_stat": float,
            "p_value": float,
            "q_value": float,
          }
    """
    org = _normalize_organism(organism)

    if not genes:
        raise ValueError("genes parameter is required and cannot be empty")

    g1_filter_full = _add_is_primary_clause(group1_filter, only_primary_data)
    g2_filter_full = _add_is_primary_clause(group2_filter, only_primary_data)

    var_filter = _build_var_value_filter(genes, gene_field=gene_field)

    with _open_census(census_version=census_version, uri=uri) as census:
        pre_n1 = _count_cells_with_min_genes(
            census,
            org,
            g1_filter_full,
            min_genes=500,
        )
        pre_n2 = _count_cells_with_min_genes(
            census,
            org,
            g2_filter_full,
            min_genes=500,
        )
        if pre_n1 == 0 or pre_n2 == 0:
            return {
                "organism": org,
                "census_version": census_version,
                "group1": {"filter": g1_filter_full, "n_cells": int(pre_n1 or 0)},
                "group2": {"filter": g2_filter_full, "n_cells": int(pre_n2 or 0)},
                "parameters": {
                    "genes": genes,
                    "gene_field": gene_field,
                    "only_primary_data": only_primary_data,
                    "lfc_threshold": lfc_threshold,
                    "min_effect_size": min_effect_size,
                    "top_n": top_n,
                },
                "results": [],
                "no_data_found": True,
                "message": (
                    "Subset too small after min_genes precheck. "
                    "Adjust filters."
                ),
            }

        adata1 = cellxgene_census.get_anndata(
            census=census,
            organism=org,
            obs_value_filter=g1_filter_full,
            var_value_filter=var_filter,
        )
        adata1 = _filter_low_coverage_cells(adata1, min_genes=500)

        adata2 = cellxgene_census.get_anndata(
            census=census,
            organism=org,
            obs_value_filter=g2_filter_full,
            var_value_filter=var_filter,
        )
        adata2 = _filter_low_coverage_cells(adata2, min_genes=500)

    n1 = int(adata1.n_obs)
    n2 = int(adata2.n_obs)

    if n1 == 0 or n2 == 0:
        return {
            "organism": org,
            "census_version": census_version,
            "group1": {"filter": g1_filter_full, "n_cells": n1},
            "group2": {"filter": g2_filter_full, "n_cells": n2},
            "parameters": {
                "genes": genes,
                "gene_field": gene_field,
                "only_primary_data": only_primary_data,
                "lfc_threshold": lfc_threshold,
                "min_effect_size": min_effect_size,
                "top_n": top_n,
            },
            "results": [],
            "no_data_found": True,
            "message": "Subset too small after low-coverage filtering.",
        }

    gene_names = _get_gene_names(adata1, gene_field=gene_field)

    # Prepare ln(CPTT+1) matrices (sparse or dense)
    X1 = _prepare_expression_matrix(adata1)
    del adata1
    X2 = _prepare_expression_matrix(adata2)
    del adata2

    if len(gene_names) == 0:
        return {
            "organism": org,
            "census_version": census_version,
            "group1": {
                "filter": g1_filter_full,
                "n_cells": n1,
            },
            "group2": {
                "filter": g2_filter_full,
                "n_cells": n2,
            },
            "parameters": {
                "genes": genes,
                "gene_field": gene_field,
                "only_primary_data": only_primary_data,
                "lfc_threshold": lfc_threshold,
                "min_effect_size": min_effect_size,
                "top_n": top_n,
            },
            "results": [],
            "no_data_found": True,
            "message": "No genes matched the provided gene list.",
        }

    n1_f = float(n1)
    n2_f = float(n2)

    def _mean_and_var_from_matrix(X, n_obs: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute mean and unbiased variance per gene from ln(CPTT+1) matrix X.

        Works for both sparse and dense matrices, treating missing entries
        as zeros (i.e., mean over all cells).
        """
        if sparse.issparse(X):
            # Mean over all cells (including zeros)
            mean = np.asarray(X.mean(axis=0)).ravel().astype(np.float32)

            # E[x^2] via element-wise square
            X_sq = X.multiply(X)
            mean_sq = np.asarray(X_sq.mean(axis=0)).ravel().astype(np.float32)

            # Unbiased variance: n/(n-1) * (E[x^2] - (E[x])^2)
            n = float(n_obs)
            var = (n / max(n - 1.0, 1.0)) * (mean_sq - mean ** 2)
            # Numerical guard: no negative variances
            var = np.clip(var, 0.0, None)
            return mean, var
        else:
            X = np.asarray(X, dtype=np.float32)
            mean = X.mean(axis=0)
            var = X.var(axis=0, ddof=1)
            return mean, var

    mean1, var1 = _mean_and_var_from_matrix(X1, n1)
    mean2, var2 = _mean_and_var_from_matrix(X2, n2)

    # Log fold-change on ln(CPTT+1) scale
    log_fc = mean1 - mean2

    # Cohen's d
    pooled_var = ((n1_f - 1.0) * var1 + (n2_f - 1.0) * var2) / max(
        n1_f + n2_f - 2.0, 1.0
    )
    pooled_var = np.clip(pooled_var, 1e-8, None)
    cohen_d = log_fc / np.sqrt(pooled_var)

    # Welch t-statistics
    se = np.sqrt(var1 / n1_f + var2 / n2_f + 1e-8)
    t_stat = log_fc / se

    # Welch–Satterthwaite degrees of freedom
    df_num = (var1 / n1_f + var2 / n2_f) ** 2
    denom1 = (var1 ** 2) / (n1_f ** 2 * max(n1_f - 1.0, 1.0))
    denom2 = (var2 ** 2) / (n2_f ** 2 * max(n2_f - 1.0, 1.0))
    df_den = denom1 + denom2
    df = df_num / (df_den + 1e-8)
    df = np.clip(df, 1.0, 1e9)

    # Two-sided p-values from t distribution
    p_val = 2.0 * stats.t.sf(np.abs(t_stat), df)

    # BH FDR
    q_val = _bh_fdr(p_val)

    df_res = pd.DataFrame(
        {
            "gene": gene_names,
            "log_fold_change": log_fc,
            "effect_size": cohen_d,
            "t_stat": t_stat,
            "p_value": p_val,
            "q_value": q_val,
        }
    )

    df_res.replace({np.inf: np.nan, -np.inf: np.nan}, inplace=True)
    df_res = df_res.dropna(
        subset=["log_fold_change", "effect_size", "p_value", "q_value"]
    )

    if lfc_threshold > 0:
        df_res = df_res[df_res["log_fold_change"].abs() >= float(lfc_threshold)]
    if min_effect_size > 0:
        df_res = df_res[df_res["effect_size"].abs() >= float(min_effect_size)]

    df_res = df_res.reindex(df_res["effect_size"].abs().sort_values(ascending=False).index)

    if top_n > 0:
        df_res = df_res.head(int(top_n))

    results = []
    for _, row in df_res.iterrows():
        results.append(
            {
                "gene": str(row["gene"]),
                "log_fold_change": float(row["log_fold_change"]),
                "effect_size": float(row["effect_size"]),
                "t_stat": float(row["t_stat"]),
                "p_value": float(row["p_value"]),
                "q_value": float(row["q_value"]),
            }
        )

    return {
        "organism": org,
        "census_version": census_version,
        "group1": {
            "filter": g1_filter_full,
            "n_cells": n1,
        },
        "group2": {
            "filter": g2_filter_full,
            "n_cells": n2,
        },
        "parameters": {
            "genes": genes,
            "gene_field": gene_field,
            "only_primary_data": only_primary_data,
            "lfc_threshold": lfc_threshold,
            "min_effect_size": min_effect_size,
            "top_n": top_n,
        },
        "results": results,
    }
    
def _census_get_filter_options_impl(
    organism: str,
    census_version: str = "2025-11-08",
    uri: Optional[str] = None,
) -> Dict[str, Any]:
    """Return unique values for tissue_general, cell_type, and disease from Census summary data."""
    org = _normalize_organism(organism)
    org_summary_format = _organism_to_summary_format(organism)

    with _open_census(census_version=census_version, uri=uri) as census:
        # Use the small aggregated table instead of scanning all obs rows
        summary = census["census_info"]["summary_cell_counts"]
        pdf = summary.read().concat().to_pandas()

    # Restrict to this organism if the column is present
    # The summary table uses lowercase with underscores (e.g., 'homo_sapiens')
    if "organism" in pdf.columns:
        pdf = pdf[pdf["organism"] == org_summary_format]

    def _labels_for_category(category: str) -> List[str]:
        """Get sorted unique labels for a given summary 'category'."""
        if "category" not in pdf.columns or "label" not in pdf.columns:
            return []
        sub = pdf[pdf["category"] == category]
        if sub.empty:
            return []
        vals = (
            sub["label"]
            .dropna()
            .astype(str)
            .drop_duplicates()
            .sort_values()
            .tolist()
        )
        return vals

    tissues = _labels_for_category("tissue_general")
    cell_types = _labels_for_category("cell_type")
    diseases = _labels_for_category("disease")

    return {
        "organism": org,
        "census_version": census_version,
        "filters": {
            "tissue_general": tissues,
            "cell_type": cell_types,
            "disease": diseases,
        },
    }


def _census_get_third_dimension_impl(
    organism: str,
    census_version: str = "2025-11-08",
    uri: Optional[str] = None,
    tissue_general: Optional[str] = None,
    cell_type: Optional[str] = None,
    disease: Optional[str] = None,
    only_primary_data: bool = True,
    candidates_of_interest: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Given two of the three filters (tissue_general, cell_type, disease),
    return which values of the third dimension have non-zero cell counts.
    
    Parameters
    ----------
    organism : str
        "Homo sapiens" or "Mus musculus" (case-insensitive aliases ok).
    census_version : str
        Census version, e.g. "latest", "2025-01-30".
    uri : Optional[str]
        Optional SOMA URI for local/S3 Census (overrides census_version).
    tissue_general : Optional[str]
        Tissue filter value (provide exactly 2 of the 3 filters).
    cell_type : Optional[str]
        Cell type filter value (provide exactly 2 of the 3 filters).
    disease : Optional[str]
        Disease filter value (provide exactly 2 of the 3 filters).
    only_primary_data : bool
        If True (default), restrict to is_primary_data == True.
    candidates_of_interest : Optional[List[str]]
        Optional list of candidate values for the third dimension. If provided,
        only return counts for these candidates. Entries with 0 count are not
        returned even if they are in this list.
    
    Returns
    -------
    dict with:
      - organism
      - census_version
      - provided_filters: dict with the two provided filters
      - missing_dimension: str (one of "tissue_general", "cell_type", "disease")
      - n_cells_total: int (total cells matching the two provided filters)
      - values: list of {value: str | null, count: int} for the missing dimension,
        sorted by count descending. Only includes entries with count > 0.
        If candidates_of_interest is provided, only includes candidates from that list.
        Note: If a candidate has no entry (count is 0 or not found), it is not
        returned in the results, even if it was provided in candidates_of_interest.
    """
    org = _normalize_organism(organism)
    
    provided = []
    if tissue_general is not None:
        provided.append(("tissue_general", tissue_general))
    if cell_type is not None:
        provided.append(("cell_type", cell_type))
    if disease is not None:
        provided.append(("disease", disease))
    
    if len(provided) != 2:
        raise ValueError(
            f"Must provide exactly 2 of the 3 filters (tissue_general, cell_type, disease). "
            f"Provided: {[p[0] for p in provided]}"
        )
    
    provided_dict = {k: v for k, v in provided}
    all_dims = {"tissue_general", "cell_type", "disease"}
    missing_dim = (all_dims - set(provided_dict.keys())).pop()
    
    filter_parts = []
    for dim, value in provided:
        filter_parts.append(f"{dim} == '{value}'")
    
    value_filter = " and ".join(filter_parts)
    if only_primary_data:
        value_filter = f"({value_filter}) and is_primary_data == True"
    
    with _open_census(census_version=census_version, uri=uri) as census:
        df = cellxgene_census.get_obs(
            census,
            org,
            value_filter=value_filter,
            column_names=[missing_dim],
        )
    
    n_cells_total = int(df.shape[0])
    
    if missing_dim not in df.columns:
        raise ValueError(f"Column '{missing_dim}' not found in obs for organism '{org}'.")
    
    vc = df[missing_dim].value_counts(dropna=False).sort_values(ascending=False)
    
    values = []
    for val, count in vc.items():
        if int(count) == 0:
            continue
        
        if pd.isna(val):
            display_val = None
        else:
            display_val = str(val)
        
        if candidates_of_interest is not None:
            if display_val not in candidates_of_interest:
                continue
        
        values.append({"value": display_val, "count": int(count)})
    
    return {
        "organism": org,
        "census_version": census_version,
        "provided_filters": provided_dict,
        "missing_dimension": missing_dim,
        "n_cells_total": n_cells_total,
        "values": values,
    }


@mcp.tool
def census_get_filter_options(
    organism: str,
    census_version: str = "2025-11-08",
    uri: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Return possible values for the high-level filters we support:
      - tissue_general
      - cell_type
      - disease

    This queries the obs metadata from Census to get unique values for these columns.

    Use this tool to discover valid filter values before constructing
    group1_filter and group2_filter for census_diffexp.

    Parameters
    ----------
    organism : str
        "Homo sapiens" or "Mus musculus" (case-insensitive aliases ok).
    census_version : str
        Census version, e.g. "stable", "latest", "2025-01-30".
    uri : Optional[str]
        Optional SOMA URI for local/S3 Census (overrides census_version).

    Returns
    -------
    dict with:
      {
        "organism": <normalized>,
        "census_version": ...,
        "filters": {
          "tissue_general": [str, ...],
          "cell_type": [str, ...],
          "disease": [str, ...],
        }
      }
    """
    return _census_get_filter_options_impl(
        organism=organism,
        census_version=census_version,
        uri=uri,
    )


@mcp.tool
def census_diffexp(
    organism: str,
    group1_filter: str,
    group2_filter: str,
    genes: List[str],
    census_version: str = "2025-11-08",
    uri: Optional[str] = None,
    gene_field: str = "feature_name",
    only_primary_data: bool = True,
    lfc_threshold: float = 0.0,
    min_effect_size: float = 0.0,
    top_n: int = 50,
) -> Dict[str, Any]:
    """
    Run differential expression for two cell groups in CELLxGENE Census.

    What this tool does
    -------------------
    - Uses Welch's t-test (unequal variance) on ln(CPTT + 1), where:
        * Raw counts are library-size–normalized to counts-per-10,000 (CPTT).
        * Then a log1p transform is applied: ln(CPTT + 1).
    - Effect size is Cohen's d on the ln(CPTT + 1) scale.
    - `log_fold_change` is mean(group1) - mean(group2) on this scale:
        * Positive => higher expression in group1.
        * Negative => higher expression in group2.
    - Genes are filtered at fetch time using a `var_value_filter` for efficiency.

    Filter expression language (group1_filter / group2_filter)
    ----------------------------------------------------------
    `group1_filter` and `group2_filter` are SOMA `obs_value_filter` strings.
    They support a small expression language, including:

    - Comparison operators:
        * `==`, `!=`, `<`, `<=`, `>`, `>=`
    - Boolean logic:
        * `and`, `or`, `not`
    - Membership (compound OR on a single field):
        * `in [...]` with a Python-like list of values

      For example, to select multiple cell types:

          "tissue_general == 'lung' and cell_type in ['B cell', 'T cell', 'NK cell']"

      Or combining disease states:

          "tissue_general == 'brain' and disease in ['Alzheimer disease', 'Parkinson disease']"

    - Parentheses to control precedence when mixing `and` / `or`:

          "(cell_type in ['B cell', 'T cell']) and (disease == 'healthy')"

    NOTE: String values must be wrapped in single quotes inside the filter.

    Primary-data handling
    ---------------------
    If `only_primary_data=True` (default), the tool will automatically append
    `and is_primary_data == True` to each of `group1_filter` and `group2_filter`
    unless the filter already mentions `is_primary_data`.

    Use `census_get_filter_options` to discover valid values for `tissue_general`,
    `cell_type`, and `disease` when constructing `group1_filter` and `group2_filter`.

    Parameters
    ----------
    organism : str
        "Homo sapiens" or "Mus musculus" (case-insensitive "human"/"mouse" also ok).
    group1_filter : str
        SOMA `obs_value_filter` for group 1, e.g.:
            "tissue_general == 'lung' and cell_type in ['B cell', 'T cell']"
        Use `census_get_filter_options` to discover valid filter values.
    group2_filter : str
        SOMA `obs_value_filter` for group 2.
        Use `census_get_filter_options` to discover valid filter values.
    genes : List[str]
        Required list of gene names to analyze (matched against `gene_field`,
        typically "feature_name", i.e. gene symbols). Genes are filtered at
        fetch time to avoid downloading unnecessary data.
    census_version : str
        Census version, e.g. "latest", "stable", "2025-01-30", "2024-07-01".
    uri : Optional[str]
        Optional SOMA URI for local / S3 Census (overrides `census_version`).
    gene_field : str
        Column in `adata.var` to use as gene name (default: "feature_name").
        If not present, falls back to `var.index`.
    only_primary_data : bool
        If True (default), automatically append "and is_primary_data == True"
        to both group filters unless they already mention `is_primary_data`.
    lfc_threshold : float
        Minimum absolute log fold change (on ln(CPTT + 1) scale) to keep
        in output (applied after t-test).
    min_effect_size : float
        Minimum absolute Cohen's d to keep in output.
    top_n : int
        Return only the top N genes by |effect_size| after filtering.

    Returns
    -------
    dict with:
      - organism
      - census_version
      - group1: {filter, n_cells}
      - group2: {filter, n_cells}
      - parameters: {...}
      - results: list of:
          {
            "gene": str,
            "log_fold_change": float,
            "effect_size": float,
            "t_stat": float,
            "p_value": float,
            "q_value": float,
          }

    Interpreting the returned statistics
    ------------------------------------
    The top-level result contains:

      - `group1`: { "filter": <effective filter string>, "n_cells": int }
      - `group2`: { "filter": <effective filter string>, "n_cells": int }
      - `results`: list of per-gene dictionaries, one per gene.

    Each entry in `results` has:

      - `gene` : str
          Gene name as given by `gene_field` (default: "feature_name").

      - `log_fold_change` : float
          Difference in group means on the ln(CPTT + 1) scale:

              log_fold_change = mean_group1 - mean_group2

          Interpretation:
            * `log_fold_change > 0`  → gene is higher on average in group1.
            * `log_fold_change < 0`  → gene is higher on average in group2.

          Approximate fold-change on the CPTT scale can be obtained by:

              fold_change ≈ exp(log_fold_change)

          Example:
            * log_fold_change ≈ +0.69 → ~2× higher in group1.
            * log_fold_change ≈ -0.69 → ~2× higher in group2.

      - `effect_size` : float
          Cohen's d on the ln(CPTT + 1) scale:

              effect_size = (mean_group1 - mean_group2) / pooled_sd

          where `pooled_sd` is the pooled standard deviation across groups.

          Interpretation (rough guidelines):
            * |d| ≈ 0.2 → small effect
            * |d| ≈ 0.5 → medium effect
            * |d| ≥ 0.8 → large effect

          This is often the best field for ranking "strength" of the
          differential signal, especially when sample sizes differ.

      - `t_stat` : float
          Welch's t-statistic for the difference in means on the ln(CPTT + 1)
          scale. The sign matches `log_fold_change`:

              t_stat > 0  → higher in group1
              t_stat < 0  → higher in group2

          You usually do not need this directly; it is used to compute p-values.

      - `p_value` : float
          Two-sided p-value from Welch's t-test for the null hypothesis of
          equal means between group1 and group2. Smaller values indicate that
          the observed difference in means is unlikely under the null.

      - `q_value` : float
          Benjamini–Hochberg FDR-corrected p-value across all tested genes.
          Interpretation is standard FDR:

            * q_value <= 0.05 → at most ~5% expected false discoveries among
              genes called "significant" at this threshold.
    """
    return _census_diffexp_impl(
        organism=organism,
        group1_filter=group1_filter,
        group2_filter=group2_filter,
        genes=genes,
        census_version=census_version,
        uri=uri,
        gene_field=gene_field,
        only_primary_data=only_primary_data,
        lfc_threshold=lfc_threshold,
        min_effect_size=min_effect_size,
        top_n=top_n,
    )


@mcp.tool
def census_get_disease_for_tissue_and_cell_type(
    organism: str,
    tissue_general: str,
    cell_type: str,
    census_version: str = "2025-11-08",
    uri: Optional[str] = None,
    only_primary_data: bool = True,
    candidates_of_interest: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Given tissue_general and cell_type filters, return which disease values have non-zero cell counts.

    This tool helps discover valid disease options for a specific tissue and cell type combination.
    Both tissue_general and cell_type are REQUIRED parameters (not optional).

    Parameters
    ----------
    organism : str
        Organism name: "Homo sapiens" or "Mus musculus" (case-insensitive aliases ok).
    tissue_general : str
        REQUIRED. Tissue filter value (e.g., "lung", "blood", "brain").
        Use census_get_filter_options to discover valid values.
    cell_type : str
        REQUIRED. Cell type filter value (e.g., "T cell", "B cell", "neuron").
        Use census_get_filter_options to discover valid values.
    census_version : str, default "2025-11-08"
        Census version identifier. Examples: "latest", "stable", "2025-01-30", "2024-07-01".
        Ignored if uri is provided.
    uri : Optional[str], default None
        Optional SOMA URI for local or S3-based Census. If provided, overrides census_version.
    only_primary_data : bool, default True
        If True, restrict results to cells where is_primary_data == True.
        If False, include all cells regardless of primary data status.
    candidates_of_interest : Optional[List[str]], default None
        Optional list of disease candidate values to check. If provided,
        only return counts for these specific candidates. Entries with 0 count are not
        returned even if they are in this list. If None, returns all non-zero disease values.

    Returns
    -------
    dict with:
      - organism : str
      - census_version : str
      - provided_filters : dict with tissue_general and cell_type
      - missing_dimension : str (always "disease")
      - n_cells_total : int
      - values : list[dict]  # each: {"value": str|None, "count": int}, count > 0 only
    """
    return _census_get_third_dimension_impl(
        organism=organism,
        census_version=census_version,
        uri=uri,
        tissue_general=tissue_general,
        cell_type=cell_type,
        disease=None,
        only_primary_data=only_primary_data,
        candidates_of_interest=candidates_of_interest,
    )


@mcp.tool
def census_get_cell_type_for_tissue_and_disease(
    organism: str,
    tissue_general: str,
    disease: str,
    census_version: str = "2025-11-08",
    uri: Optional[str] = None,
    only_primary_data: bool = True,
    candidates_of_interest: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Given tissue_general and disease filters, return which cell_type values have non-zero cell counts.

    This tool helps discover valid cell type options for a specific tissue and disease combination.
    Both tissue_general and disease are REQUIRED parameters (not optional).

    Parameters
    ----------
    organism : str
        Organism name: "Homo sapiens" or "Mus musculus" (case-insensitive aliases ok).
    tissue_general : str
        REQUIRED. Tissue filter value (e.g., "lung", "blood", "brain").
        Use census_get_filter_options to discover valid values.
    disease : str
        REQUIRED. Disease filter value (e.g., "normal", "Alzheimer disease", "cancer").
        Use census_get_filter_options to discover valid values.
    census_version : str, default "2025-11-08"
        Census version identifier. Examples: "latest", "stable", "2025-01-30", "2024-07-01".
        Ignored if uri is provided.
    uri : Optional[str], default None
        Optional SOMA URI for local or S3-based Census. If provided, overrides census_version.
    only_primary_data : bool, default True
        If True, restrict results to cells where is_primary_data == True.
        If False, include all cells regardless of primary data status.
    candidates_of_interest : Optional[List[str]], default None
        Optional list of cell_type candidate values to check. If provided,
        only return counts for these specific candidates. Entries with 0 count are not
        returned even if they are in this list. If None, returns all non-zero cell_type values.

    Returns
    -------
    dict with:
      - organism : str
      - census_version : str
      - provided_filters : dict with tissue_general and disease
      - missing_dimension : str (always "cell_type")
      - n_cells_total : int
      - values : list[dict]  # each: {"value": str|None, "count": int}, count > 0 only
    """
    return _census_get_third_dimension_impl(
        organism=organism,
        census_version=census_version,
        uri=uri,
        tissue_general=tissue_general,
        cell_type=None,
        disease=disease,
        only_primary_data=only_primary_data,
        candidates_of_interest=candidates_of_interest,
    )


@mcp.tool
def census_get_tissue_for_cell_type_and_disease(
    organism: str,
    cell_type: str,
    disease: str,
    census_version: str = "2025-11-08",
    uri: Optional[str] = None,
    only_primary_data: bool = True,
    candidates_of_interest: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Given cell_type and disease filters, return which tissue_general values have non-zero cell counts.

    This tool helps discover valid tissue options for a specific cell type and disease combination.
    Both cell_type and disease are REQUIRED parameters (not optional).

    Parameters
    ----------
    organism : str
        Organism name: "Homo sapiens" or "Mus musculus" (case-insensitive aliases ok).
    cell_type : str
        REQUIRED. Cell type filter value (e.g., "T cell", "B cell", "neuron").
        Use census_get_filter_options to discover valid values.
    disease : str
        REQUIRED. Disease filter value (e.g., "normal", "Alzheimer disease", "cancer").
        Use census_get_filter_options to discover valid values.
    census_version : str, default "2025-11-08"
        Census version identifier. Examples: "latest", "stable", "2025-01-30", "2024-07-01".
        Ignored if uri is provided.
    uri : Optional[str], default None
        Optional SOMA URI for local or S3-based Census. If provided, overrides census_version.
    only_primary_data : bool, default True
        If True, restrict results to cells where is_primary_data == True.
        If False, include all cells regardless of primary data status.
    candidates_of_interest : Optional[List[str]], default None
        Optional list of tissue_general candidate values to check. If provided,
        only return counts for these specific candidates. Entries with 0 count are not
        returned even if they are in this list. If None, returns all non-zero tissue_general values.

    Returns
    -------
    dict with:
      - organism : str
      - census_version : str
      - provided_filters : dict with cell_type and disease
      - missing_dimension : str (always "tissue_general")
      - n_cells_total : int
      - values : list[dict]  # each: {"value": str|None, "count": int}, count > 0 only
    """
    return _census_get_third_dimension_impl(
        organism=organism,
        census_version=census_version,
        uri=uri,
        tissue_general=None,
        cell_type=cell_type,
        disease=disease,
        only_primary_data=only_primary_data,
        candidates_of_interest=candidates_of_interest,
    )


@mcp.tool
def census_list_obs_columns(
    organism: str,
    census_version: str = "2025-11-08",
    uri: Optional[str] = None,
) -> Dict[str, Any]:
    """
    List relevant column names in the Census obs (cell metadata) schema.
    
    Returns a curated list of the most commonly used column names for filtering:
    cell_type, disease, is_primary_data, and tissue_general.
    
    Use this tool to discover the correct column names to use in filter strings
    for census_diffexp, census_get_disease_for_tissue_and_cell_type,
    census_get_cell_type_for_tissue_and_disease, census_get_tissue_for_cell_type_and_disease,
    and other tools.
    
    Parameters
    ----------
    organism : str
        "Homo sapiens" or "Mus musculus" (case-insensitive aliases ok).
        Note: This parameter is kept for API compatibility but is not used.
    census_version : str
        Census version, e.g. "latest", "stable", "2025-01-30", "2024-07-01".
        Note: This parameter is kept for API compatibility but is not used.
    uri : Optional[str]
        Optional SOMA URI for local / S3 Census (overrides census_version).
        Note: This parameter is kept for API compatibility but is not used.
    
    Returns
    -------
    dict with:
      - organism: str (normalized)
      - census_version: str
      - columns: list[str] (relevant column names: cell_type, disease, is_primary_data, tissue_general)
      - n_columns: int (total number of columns)
    """
    org = _normalize_organism(organism)
    columns = _get_relevant_obs_columns()
    
    return {
        "organism": org,
        "census_version": census_version,
        "columns": columns,
        "n_columns": len(columns),
    }


def _test_census_get_filter_options() -> None:
    """
    Quick manual test: run census_get_filter_options.

    Usage:
        python gene_expression_tools.py test-filter-options
    """
    import sys
    
    # Debug: check what's in the summary table
    org = _normalize_organism("Homo sapiens")
    with _open_census(census_version="2025-11-08") as census:
        summary = census["census_info"]["summary_cell_counts"]
        pdf = summary.read().concat().to_pandas()
        
        print("\n=== DEBUG: summary_cell_counts table ===")
        print(f"Total rows: {len(pdf)}")
        print(f"Columns: {list(pdf.columns)}")
        if len(pdf) > 0:
            print(f"\nFirst few rows:")
            print(pdf.head())
            if "organism" in pdf.columns:
                print(f"\nOrganism values: {pdf['organism'].unique()}")
            elif "organism_name" in pdf.columns:
                print(f"\nOrganism values: {pdf['organism_name'].unique()}")
        sys.stdout.flush()
    
    result = _census_get_filter_options_impl(
        organism="Homo sapiens",
        census_version="2025-11-08",
    )

    print("\n=== census_get_filter_options test ===")
    print(f"Organism        : {result['organism']}")
    print(f"Census version  : {result['census_version']}")
    print(f"\nTissue general options ({len(result['filters']['tissue_general'])}):")
    for tissue in result['filters']['tissue_general'][:20]:
        print(f"  - {tissue}")
    if len(result['filters']['tissue_general']) > 20:
        print(f"  ... and {len(result['filters']['tissue_general']) - 20} more")
    
    print(f"\nCell type options ({len(result['filters']['cell_type'])}):")
    for cell_type in result['filters']['cell_type'][:20]:
        print(f"  - {cell_type}")
    if len(result['filters']['cell_type']) > 20:
        print(f"  ... and {len(result['filters']['cell_type']) - 20} more")
    
    print(f"\nDisease options ({len(result['filters']['disease'])}):")
    for disease in result['filters']['disease'][:20]:
        print(f"  - {disease}")
    if len(result['filters']['disease']) > 20:
        print(f"  ... and {len(result['filters']['disease']) - 20} more")


def _test_census_diffexp() -> None:
    """
    Quick manual test: run census_diffexp with compound filters using `in`.

    Usage:
        python gene_expression_tools.py test-diff-expression

    Test cases:
      - Test 1: blood T/NK lineage cells vs B/plasma lineage cells
      - Test 2: brain cerebellar neurons with dementia/Alzheimer vs normal
    """
    import sys
    
    organism = "Homo sapiens"
    census_version = "2025-01-30"
    
    test_cases = [
        {
            "name": "Blood T/NK vs B/plasma cells",
            "group1_filter": (
                "tissue_general == 'blood' and "
                "cell_type in ['T cell', 'NK cell']"
            ),
            "group2_filter": (
                "tissue_general == 'blood' and "
                "cell_type in ['B cell', 'plasma cell']"
            ),
            "genes": ["CD79A"],
        },
        {
            "name": "Brain cerebellar neurons: dementia/Alzheimer vs normal",
            "group1_filter": (
                "tissue_general == 'brain' and "
                "cell_type == 'cerebellar neuron' and "
                "(disease == 'dementia' or disease == 'Alzheimer disease')"
            ),
            "group2_filter": (
                "tissue_general == 'brain' and "
                "cell_type == 'cerebellar neuron' and "
                "disease == 'normal'"
            ),
            "genes": ["CD79A"],
        },
        {
            "name": "Bone marrow AML vs normal (myeloid/progenitors)",
            "group1_filter": (
                "tissue_general == 'bone marrow' and "
                "cell_type in ["
                "'classical monocyte', "
                "'early promyelocyte', "
                "'late promyelocyte', "
                "'myelocyte', "
                "'hematopoietic multipotent progenitor cell'"
                "] and "
                "disease == 'acute myeloid leukemia'"
            ),
            "group2_filter": (
                "tissue_general == 'bone marrow' and "
                "cell_type in ["
                "'classical monocyte', "
                "'early promyelocyte', "
                "'late promyelocyte', "
                "'myelocyte', "
                "'hematopoietic multipotent progenitor cell'"
                "] and "
                "disease == 'normal'"
            ),
            "genes": ["CBL", "MAPK1", "MAPK14", "CDK6", "CDK9"],
        },
    ]
    
    print("\n=== Testing census_diffexp ===\n")
    sys.stdout.flush()
    
    for i, test_case in enumerate(test_cases, 1):
        print(f"\n--- Test {i}: {test_case['name']} ---")
        sys.stdout.flush()
        if i == 1:
            continue
        
        
        try:
            result = _census_diffexp_impl(
                organism=organism,
                group1_filter=test_case["group1_filter"],
                group2_filter=test_case["group2_filter"],
                genes=test_case["genes"],
                census_version=census_version,
                gene_field="feature_name",
                only_primary_data=True,
                lfc_threshold=0.25,
                min_effect_size=0.5,
                top_n=100,
            )
            
            print("Test census_diffexp result:")
            print(result)
            sys.stdout.flush()
            
        except Exception as e:
            print(f"ERROR: {str(e)}")
            import traceback
            traceback.print_exc()
            sys.stdout.flush()
    
    print("\n=== All tests completed ===\n")
    sys.stdout.flush()


def _test_census_get_third_dimension() -> None:
    """
    Test census_get_third_dimension with different combinations of filters.
    
    Usage:
        python gene_expression_tools.py test-third-dimension
    """
    import sys
    
    organism = "Homo sapiens"
    census_version = "2025-11-08"
    
    print("\n=== Testing census_get_third_dimension ===\n")
    sys.stdout.flush()
    
    test_cases = [
        {
            "name": "Given tissue_general and cell_type, find diseases",
            "tissue_general": "blood",
            "cell_type": "T cell",
            "disease": None,
        },
        {
            "name": "Given tissue_general and disease, find cell_types",
            "tissue_general": "brain",
            "cell_type": None,
            "disease": "Alzheimer disease",
        },
        {
            "name": "Given cell_type and disease, find tissue_general",
            "tissue_general": None,
            "cell_type": "B cell",
            "disease": "normal",
        },
    ]
    
    for i, test_case in enumerate(test_cases, 1):
        print(f"\n--- Test {i}: {test_case['name']} ---")
        sys.stdout.flush()
        
        try:
            result = _census_get_third_dimension_impl(
                organism=organism,
                census_version=census_version,
                tissue_general=test_case["tissue_general"],
                cell_type=test_case["cell_type"],
                disease=test_case["disease"],
                only_primary_data=True,
            )
            
            print(f"Organism: {result['organism']}")
            print(f"Census version: {result['census_version']}")
            print(f"Provided filters: {result['provided_filters']}")
            print(f"Missing dimension: {result['missing_dimension']}")
            print(f"Total cells matching filters: {result['n_cells_total']:,}")
            print(f"\nValues for {result['missing_dimension']} ({len(result['values'])} total):")
            
            for j, val_info in enumerate(result['values'][:20], 1):
                val_str = val_info['value'] if val_info['value'] is not None else "null"
                print(f"  {j:3d}. {val_str:40s}  count: {val_info['count']:,}")
            
            if len(result['values']) > 20:
                print(f"  ... and {len(result['values']) - 20} more values")
            
            sys.stdout.flush()
            
        except Exception as e:
            print(f"ERROR: {str(e)}")
            import traceback
            traceback.print_exc()
            sys.stdout.flush()
    
    print("\n=== All tests completed ===\n")
    sys.stdout.flush()



if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "test-filter-options":
        _test_census_get_filter_options()
    elif len(sys.argv) > 1 and sys.argv[1] == "test-diff-expression":
        _test_census_diffexp()
    elif len(sys.argv) > 1 and sys.argv[1] == "test-third-dimension":
        _test_census_get_third_dimension()
    else:
        # Default: run as an MCP server
        mcp.run()