"""Chao2 saturation curves for clustered single-cell or other incidence datasets."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Union, Dict, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import anndata as ad
except ImportError:  # AnnData support is optional
    ad = None


@dataclass
class Chao2SaturationResult:
    """Container returned by Chao2_saturation()."""
    curve: pd.DataFrame
    incidence: pd.DataFrame
    summary: pd.DataFrame
    fig: Optional[plt.Figure]
    ax: Optional[plt.Axes]


def _chao2_stats(X: np.ndarray) -> Dict[str, float]:
    """Calculate observed richness, Chao2, bias-corrected Chao2 and iChao2."""
    X = (np.asarray(X) > 0).astype(int)
    if X.ndim != 2:
        raise ValueError("Incidence matrix must be 2-dimensional: rows=clusters, columns=samples.")

    T = X.shape[1]
    incidence = X.sum(axis=1)
    S_obs = int(np.sum(incidence > 0))
    Q1 = int(np.sum(incidence == 1))
    Q2 = int(np.sum(incidence == 2))
    Q3 = int(np.sum(incidence == 3))
    Q4 = int(np.sum(incidence == 4))

    # Classic Chao2. If Q2=0, use the usual fallback to S_obs.
    S_Chao2 = S_obs + (Q1**2) / (2 * Q2) if Q2 > 0 else float(S_obs)

    # Bias-corrected Chao2; useful for finite sample counts.
    if T > 0 and Q2 > 0:
        S_Chao2_bc = S_obs + ((T - 1) / T) * (Q1**2) / (2 * Q2)
    elif T > 0:
        S_Chao2_bc = S_obs + ((T - 1) / T) * (Q1 * (Q1 - 1)) / 2
    else:
        S_Chao2_bc = float(S_obs)

    # iChao2 extension; defined only sensibly for T >= 4.
    if T >= 4:
        denom = Q4 if Q4 > 0 else 1
        extra = ((T - 3) / (4 * T)) * (Q3 / denom) * max(
            Q1 - ((T - 3) / (2 * (T - 1))) * (Q2 * Q3 / denom),
            0,
        )
    else:
        extra = 0.0

    return {
        "n_samples": int(T),
        "S_obs": float(S_obs),
        "Q1": float(Q1),
        "Q2": float(Q2),
        "Q3": float(Q3),
        "Q4": float(Q4),
        "S_Chao2": float(S_Chao2),
        "S_Chao2_bc": float(S_Chao2_bc),
        "S_iChao2": float(S_Chao2_bc + extra),
    }


def _make_incidence_table(
    data: Union[pd.DataFrame, "ad.AnnData"],
    cluster_col: Optional[str],
    sample_col: Optional[str],
    layer: str = "obs",
) -> pd.DataFrame:
    """Create rows=clusters, columns=samples binary incidence table."""
    if ad is not None and isinstance(data, ad.AnnData):
        if cluster_col is None or sample_col is None:
            raise ValueError("For AnnData input, provide cluster_col and sample_col from adata.obs.")
        obs = data.obs
        if cluster_col not in obs or sample_col not in obs:
            missing = [c for c in [cluster_col, sample_col] if c not in obs]
            raise KeyError(f"Missing column(s) in adata.obs: {missing}")
        df = obs[[cluster_col, sample_col]].dropna().copy()
        ctab = pd.crosstab(df[cluster_col].astype(str), df[sample_col].astype(str))
        return (ctab > 0).astype(int)

    if isinstance(data, pd.DataFrame):
        if cluster_col is None and sample_col is None:
            # Treat as an already-prepared cluster x sample table.
            return (data.fillna(0).astype(float) > 0).astype(int)
        if cluster_col is None or sample_col is None:
            raise ValueError("Provide both cluster_col and sample_col, or neither for a precomputed incidence table.")
        if cluster_col not in data.columns or sample_col not in data.columns:
            missing = [c for c in [cluster_col, sample_col] if c not in data.columns]
            raise KeyError(f"Missing column(s) in DataFrame: {missing}")
        df = data[[cluster_col, sample_col]].dropna().copy()
        ctab = pd.crosstab(df[cluster_col].astype(str), df[sample_col].astype(str))
        return (ctab > 0).astype(int)

    raise TypeError("data must be a pandas DataFrame, an AnnData object, or a precomputed incidence DataFrame.")


def Chao2_saturation(
    data: Union[pd.DataFrame, "ad.AnnData"],
    cluster_col: Optional[str] = None,
    sample_col: Optional[str] = None,
    group_col: Optional[str] = None,
    groups: Optional[Sequence[str]] = None,
    n_permutations: int = 200,
    random_state: Optional[int] = 0,
    estimator: str = "S_Chao2_bc",
    plot: bool = True,
    ax: Optional[plt.Axes] = None,
    title: Optional[str] = None,
) -> Chao2SaturationResult:
    """
    Plot Chao2 saturation curves from clustered data.

    Parameters
    ----------
    data
        pandas DataFrame, AnnData object, or precomputed incidence table.
        For raw metadata, rows are cells/observations. For incidence tables,
        rows are clusters/species and columns are samples/donors.
    cluster_col
        Column containing cluster/species labels, e.g. "leiden" or "cell_type".
        Required unless `data` is already an incidence table.
    sample_col
        Column containing sample/donor labels, e.g. "donor_id".
        Required unless `data` is already an incidence table.
    group_col
        Optional column used to stratify the analysis, e.g. "donor_group".
        Only supported for raw DataFrame or AnnData input.
    groups
        Optional subset/order of groups to include.
    n_permutations
        Number of random sample-order permutations used for the curve.
    random_state
        Seed for reproducible curves.
    estimator
        Which estimate to plot: "S_obs", "S_Chao2", "S_Chao2_bc", or "S_iChao2".
    plot
        Whether to draw a matplotlib curve.
    ax
        Optional matplotlib axes.
    title
        Optional plot title.

    Returns
    -------
    Chao2SaturationResult with:
        curve: long DataFrame of saturation values
        incidence: binary cluster x sample table for all included data
        summary: final estimate per group
        fig, ax: matplotlib objects, or None if plot=False
    """
    valid_estimators = {"S_obs", "S_Chao2", "S_Chao2_bc", "S_iChao2"}
    if estimator not in valid_estimators:
        raise ValueError(f"estimator must be one of {sorted(valid_estimators)}")
    if n_permutations < 1:
        raise ValueError("n_permutations must be >= 1")

    rng = np.random.default_rng(random_state)

    # Convert AnnData to obs DataFrame for grouped runs.
    if ad is not None and isinstance(data, ad.AnnData):
        obs_data = data.obs.copy()
    elif isinstance(data, pd.DataFrame):
        obs_data = data.copy()
    else:
        obs_data = None

    if group_col is not None:
        if obs_data is None:
            raise ValueError("group_col requires raw DataFrame or AnnData input, not a precomputed incidence table.")
        if group_col not in obs_data.columns:
            raise KeyError(f"Missing group_col in data: {group_col}")
        group_values = list(groups) if groups is not None else list(pd.Series(obs_data[group_col]).dropna().astype(str).unique())
        incidence_all = _make_incidence_table(obs_data, cluster_col, sample_col)
        incidence_by_group = {
            str(g): _make_incidence_table(obs_data[obs_data[group_col].astype(str) == str(g)], cluster_col, sample_col)
            for g in group_values
        }
    else:
        incidence_all = _make_incidence_table(data, cluster_col, sample_col)
        incidence_by_group = {"all": incidence_all}

    curve_rows = []
    summary_rows = []
    for group_name, incidence in incidence_by_group.items():
        samples = np.array(incidence.columns)
        if len(samples) == 0:
            continue

        for rep in range(n_permutations):
            ordered = rng.permutation(samples)
            for k in range(1, len(ordered) + 1):
                Xk = incidence.loc[:, ordered[:k]].to_numpy()
                stats = _chao2_stats(Xk)
                curve_rows.append({
                    "group": group_name,
                    "permutation": rep,
                    "n_samples": k,
                    **stats,
                })

        final_stats = _chao2_stats(incidence.to_numpy())
        summary_rows.append({"group": group_name, **final_stats})

    curve = pd.DataFrame(curve_rows)
    summary = pd.DataFrame(summary_rows)

    fig = None
    if plot:
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))
        else:
            fig = ax.figure

        mean_curve = (
            curve.groupby(["group", "n_samples"], as_index=False)[estimator]
            .agg(["mean", "std"])
            .reset_index()
        )

        for group_name, sub in mean_curve.groupby("group"):
            x = sub["n_samples"].to_numpy(dtype=float)
            y = sub["mean"].to_numpy(dtype=float)
            sd = sub["std"].fillna(0).to_numpy(dtype=float)
            ax.plot(x, y, marker="o", label=str(group_name))
            ax.fill_between(x, y - sd, y + sd, alpha=0.2)

        ax.set_xlabel("Number of samples")
        ax.set_ylabel(estimator)
        ax.set_title(title or "Chao2 saturation curve")
        ax.legend(frameon=False)
        fig.tight_layout()

    return Chao2SaturationResult(curve=curve, incidence=incidence_all, summary=summary, fig=fig, ax=ax)
