# Chao2 saturation curves

`Chao2_saturation()` estimates and plots how cluster/species richness saturates as more samples are added. It works with single-cell metadata where clusters have already been defined, for example Leiden clusters or annotated cell types.

## Installation

Copy `chao2_saturation.py` into your project, then import it:

```python
from chao2_saturation import Chao2_saturation
```

Dependencies:

```bash
pip install numpy pandas matplotlib anndata
```

`anndata` is only needed if you want to pass an `AnnData` object directly.

## Basic use with AnnData.

```python
#We suggest users to select either classic (S_Chao2) or bias-corrected Chao2 (S_Chao2_bc).
#adata is an AnnData object obtainable through any Scanpy workflow (https://scanpy.scverse.org/en/stable/)
result = Chao2_saturation(
    adata,
    cluster_col="leiden",
    sample_col="donor_id",
    group_col="assay",
    estimator="S_Chao2_bc",
    min_samples=3,
    n_permutations=200,
    random_state=0,
    title="Comparison of assays"
)

result.summary
result.fig.savefig("chao2_saturation.svg", bbox_inches="tight")
```

## Compare groups, for example fetus vs mother

```python
result = Chao2_saturation(
    adata,
    cluster_col="leiden",
    sample_col="donor_id",
    group_col="donor_group",
    groups=["fetus", "mother"],
    estimator="S_Chao2_bc",
)
```

## Use a pandas DataFrame

```python
result = Chao2_saturation(
    metadata,
    cluster_col="cell_type",
    sample_col="donor_id",
    group_col="donor_group",
)
```

## Use a precomputed incidence table

Rows should be clusters/species and columns should be samples/donors. Values can be counts or 0/1 presence/absence.

```python
result = Chao2_saturation(ctab_binary)
```

## Output

The function returns a `Chao2SaturationResult` object:

- `result.curve`: long-format table used for the saturation curve
- `result.incidence`: binary cluster-by-sample incidence table
- `result.summary`: final observed and estimated richness per group
- `result.fig`, `result.ax`: matplotlib figure and axes

## Notes

Chao2 is incidence-based, so the function first converts cluster-by-sample counts to presence/absence. In single-cell data this means a cluster is counted as present in a donor/sample if at least one cell from that cluster is observed in that donor/sample.
