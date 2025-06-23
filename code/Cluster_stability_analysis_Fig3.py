"""
Cluster-stability benchmark
"""
# Imports 
import os, random, argparse, re
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import linear_sum_assignment

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help=".h5ad with integrated data")
parser.add_argument("--outdir", default="stability_out")
parser.add_argument("--iters",  type=int, default=30)
parser.add_argument("--kperc",  type=float, default=1.0,
                    help="neighbours per-100-cells (k)")
parser.add_argument("--res_mega", type=float, default=1.2)
parser.add_argument("--res_sub",  type=float, default=1.2)
parser.add_argument("--pilot_frac", type=float, default=None,
                    help="0<frac≤1 subsample remaining dataset")
parser.add_argument("--skip", default="",
                    help="Comma-separated substrings; samples whose name "
                         "contains any substring are removed.")
ARGS = parser.parse_args()

random.seed(0); np.random.seed(0)
os.makedirs(ARGS.outdir, exist_ok=True)

# Helper functions

def subsample_grid(n):
    if n <= 150:
        return np.array([n], int)
    pts = min(10, max(2, n // 150))
    return np.unique(np.linspace(150, n, pts, dtype=int))

def summarise_by_2d_groups(sim_mat: sp.csr_matrix, rows, cols):
   
    r_u, r_code = np.unique(rows, return_inverse=True)
    c_u, c_code = np.unique(cols, return_inverse=True)
    sums   = np.zeros((len(r_u), len(c_u)))
    counts = np.zeros_like(sums, int)
    coo = sim_mat.tocoo()
    np.add.at(sums,   (r_code[coo.row], c_code[coo.col]), coo.data)
    np.add.at(counts, (r_code[coo.row], c_code[coo.col]), 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        means = sums / counts
    means[counts == 0] = 0.0
    return means

def similarity_knn(adata, n_neighbors):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    return adata.obsp["connectivities"].copy()

def save_heat(mat, path, title):
    sns.heatmap(mat, cmap="coolwarm", vmin=0, vmax=1,
                linewidths=.5, linecolor="gray")
    plt.title(title); plt.tight_layout(); plt.savefig(path, dpi=300); plt.close()

# Load data
print("Loading data")
adata = sc.read(ARGS.input)

if ARGS.skip:
    patterns = [re.escape(s) for s in re.split(r"[,\s]+", ARGS.skip) if s]
    mask_drop = adata.obs["sample"].str.contains("|".join(patterns))
    n_drop = int(mask_drop.sum())
    if n_drop:
        print(f"Removing {n_drop} cells from skipped samples")
        adata = adata[~mask_drop].copy()
    else:
        print(" No cells matched --skip patterns")

# pilot subsample  

if ARGS.pilot_frac:
    f = ARGS.pilot_frac
    assert 0 < f <= 1, "--pilot_frac must be in (0,1]"
    idx_keep = np.random.choice(adata.n_obs, int(f * adata.n_obs),
                                replace=False)
    adata = adata[idx_keep].copy()
    print(f"Pilot subsample {f:.2%} → {adata.n_obs} cells total")

# PCA  
print("Computing PCA")
sc.pp.pca(adata, n_comps=50,  svd_solver="arpack")

print(" Building global neighbours graph")
k_full = max(100, round(ARGS.kperc * adata.n_obs / 100))
sc.pp.neighbors(adata, n_neighbors=k_full)


sc.tl.leiden(adata, resolution=ARGS.res_mega, key_added="megacluster")
mega_codes = adata.obs["megacluster"].astype("category").cat.codes.values
mega_names = adata.obs["megacluster"].cat.categories.tolist()

# Per-sample stability permutations
for sample in adata.obs["sample"].unique():
    print(f"\n=== SAMPLE {sample} ===")
    sample_dir = os.path.join(ARGS.outdir, sample)
    os.makedirs(sample_dir, exist_ok=True)

    idx_all = np.where(adata.obs["sample"] == sample)[0]
    N_total = idx_all.size
    print(f"{N_total} cells")

    GRID = subsample_grid(N_total)
    print("Grid:", GRID.tolist())

    # accumulator for final ratio heat-map
    mat = np.zeros((len(mega_names), len(GRID)))

    for it in range(ARGS.iters):
        np.random.shuffle(idx_all)
        iter_dir = os.path.join(sample_dir, f"iter_{it+1:02d}")
        os.makedirs(iter_dir, exist_ok=True)

        for j, N in enumerate(GRID):
            idxN = idx_all[:N]
            sub  = adata[idxN].copy()

            sc.pp.pca(sub, n_comps=50, zero_center=True)
            k_i = max(5, round(ARGS.kperc * N / 100))
            sim = similarity_knn(sub, n_neighbors=k_i)
            sc.tl.leiden(sub, resolution=ARGS.res_sub,
                         key_added="subcluster")
            sub_lab = sub.obs["subcluster"].astype("category") \
                                   .cat.codes.values

            full_summ  = summarise_by_2d_groups(sim,
                                                mega_codes[idxN],
                                                mega_codes[idxN])
            inter_summ = summarise_by_2d_groups(sim,
                                                mega_codes[idxN],
                                                sub_lab)
            row_ind, col_ind = linear_sum_assignment(
                -np.nan_to_num(inter_summ, nan=0.0))

            #  save diagnostics
            np.save(os.path.join(iter_dir, f"full_summ_N{N}.npy"),
                    full_summ)
            np.save(os.path.join(iter_dir, f"inter_summ_N{N}.npy"),
                    inter_summ)
            np.savez(os.path.join(iter_dir, f"assign_N{N}.npz"),
                     row_ind=row_ind, col_ind=col_ind)

            if it == ARGS.iters - 1:  # final permutation → heat-maps
                save_heat(full_summ,
                          os.path.join(iter_dir,
                                       f"full_summ_heat_N{N}.png"),
                          f"{sample} full_summ N={N}")
                save_heat(inter_summ,
                          os.path.join(iter_dir,
                                       f"inter_summ_heat_N{N}.png"),
                          f"{sample} inter_summ N={N}")

            #stability ratio
            present_rows = np.unique(mega_codes[idxN])
            ratio_present = np.zeros(len(present_rows))

            for r_local, c_local in zip(row_ind, col_ind):
                within = np.nan_to_num(full_summ[r_local, r_local])
                if within > 0:
                    ratio_present[r_local] = \
                        inter_summ[r_local, c_local] / within

            ratio_full = np.zeros(len(mega_names))
            ratio_full[present_rows] = ratio_present
            mat[:, j] += ratio_full

        print(f"perm {it+1}/{ARGS.iters} done", end="\r")

    #  final heat-map for this sample 
    mat /= ARGS.iters
    df = pd.DataFrame(mat,
                      index=[f"M{m}" for m in mega_names],
                      columns=[str(n) for n in GRID])
    plt.figure(figsize=(12, max(4, len(mega_names) * 0.25)))
    sns.heatmap(df, cmap="coolwarm", vmin=0, vmax=1,
                linewidths=.5, linecolor="gray")
    plt.title(sample)
    plt.xlabel("Subsample size"); plt.ylabel("Megacluster")
    plt.tight_layout()
    final_png = os.path.join(ARGS.outdir, f"heat_{sample}.png")
    plt.savefig(final_png, dpi=300); plt.close()
    print(f"\n Final heat-map {final_png}")
