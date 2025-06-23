
"""
Collapse stability-benchmark diagnostics into ONE Parquet file
but **ignore files older than --since**.
Usage

# timestamp with time component 
python prepare_for_R.py --indir stability_out \
                        --since 2025-05-30T14:00:00 \
                        --outfile stability_long.parquet
"""

import os, argparse, re, glob, numpy as np, pandas as pd, datetime as dt
import pyarrow as pa, pyarrow.parquet as pq
from dateutil import parser as dtparse   

p = argparse.ArgumentParser()
p.add_argument("--indir",   required=True,
               help="Top-level output directory of the benchmark")
p.add_argument("--outfile", default="stability_long.parquet",
               help="Parquet file to write")
p.add_argument("--precision", type=int, default=6,
               help="round matrix values to this many decimals")
p.add_argument("--since",
               help="ignore files last-modified *before* this date/time "
                    "(e.g. 2025-05-30 or 2025-05-30T14:00:00)")
args = p.parse_args()

since_ts = None
if args.since:
    since_dt = dtparse.parse(args.since)
    since_ts = since_dt.timestamp()
    print("Only files modified on/after", since_dt.isoformat())

fmt = f"%.{args.precision}f"
rows = []
pat_iter = re.compile(r"iter_(\d+)")
pat_size = re.compile(r"_N(\d+)\.")

def keep(path):

    return since_ts is None or os.path.getmtime(path) >= since_ts

for sample in sorted(next(os.walk(args.indir))[1]):          
    sdir = os.path.join(args.indir, sample)
    for iter_dir in sorted(glob.glob(os.path.join(sdir, "iter_*"))):
        if not keep(iter_dir):
            continue
        iter_no = int(pat_iter.search(iter_dir).group(1))

        
        for kind in ("full", "inter"):
            for npy in glob.glob(os.path.join(iter_dir,
                                              f"{kind}_summ_N*.npy")):
                if not keep(npy):
                    continue
                N = int(pat_size.search(npy).group(1))
                mat = np.round(np.load(npy), args.precision)
                r, c = mat.shape
                rr, cc = np.meshgrid(np.arange(r), np.arange(c),
                                     indexing="ij")
                rows.append(pd.DataFrame({
                    "sample": sample,
                    "iter"  : iter_no,
                    "N"     : N,
                    "kind"  : kind,
                    "row"   : rr.ravel(),
                    "col"   : cc.ravel(),
                    "value" : mat.ravel()
                }))

    
        for npz in glob.glob(os.path.join(iter_dir, "assign_N*.npz")):
            if not keep(npz):
                continue
            N = int(pat_size.search(npz).group(1))
            data = np.load(npz)
            rows.append(pd.DataFrame({
                "sample": sample,
                "iter"  : iter_no,
                "N"     : N,
                "kind"  : "assign",
                "row"   : data["row_ind"],
                "col"   : data["col_ind"],
                "value" : np.nan
            }))


# ───────────────
if not rows:
    raise SystemExit("No files matched the --since filter; nothing written.")

df = pd.concat(rows, ignore_index=True)
pq.write_table(pa.Table.from_pandas(df), args.outfile,
               compression="zstd")
print(f"Wrote {len(df):,} rows to {args.outfile}")