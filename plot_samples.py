#!/usr/bin/env python3
"""
Generate per-sample spatial QC plots from segmented AnnData.

Usage:
  python plot_samples.py --h5ad visium_hd_sample_split.h5ad --out sample_plots.png
  python plot_samples.py --h5ad visium_hd_sample_split.h5ad --feature nCount_Spatial
"""

import argparse
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_per_sample(adata, feature="total_counts", log_transform=True):
    """One panel per sample, colored by expression feature."""
    sids = sorted([int(x) for x in set(adata.obs["sample_id"]) if int(x) > 0])
    n = len(sids)
    if n == 0:
        raise ValueError("No samples found (all sample_id == 0)")

    fig, axes = plt.subplots(1, n, figsize=(6 * n, 6))
    if n == 1:
        axes = [axes]

    # Compute feature values
    if feature in adata.obs.columns:
        vals_all = adata.obs[feature].values
    elif feature in adata.var_names:
        vals_all = adata[:, feature].X.toarray().flatten() if hasattr(adata[:, feature].X, "toarray") else np.asarray(adata[:, feature].X).flatten()
    else:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        if feature in adata.obs.columns:
            vals_all = adata.obs[feature].values
        else:
            raise ValueError(f"Feature '{feature}' not found in obs or var")

    for i, sid in enumerate(sids):
        mask = adata.obs["sample_id"].astype(int) == sid
        coords = adata.obsm["spatial"][mask]
        vals = vals_all[mask]
        if log_transform:
            vals = np.log1p(vals)

        ax = axes[i]
        scat = ax.scatter(coords[:, 0], coords[:, 1],
                          c=vals, s=1, cmap="viridis", rasterized=True)
        ax.set_title(f"sample_{sid:03d} (n={mask.sum():,})", fontsize=13)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.axis("off")
        label = f"log1p({feature})" if log_transform else feature
        plt.colorbar(scat, ax=ax, label=label, shrink=0.7)

    plt.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser(description="Per-sample spatial QC plots")
    parser.add_argument("--h5ad", required=True, help="Path to visium_hd_sample_split.h5ad")
    parser.add_argument("--feature", default="total_counts",
                        help="Feature to color by (obs column or gene name)")
    parser.add_argument("--out", default="sample_plots.png", help="Output PNG path")
    parser.add_argument("--no-log", action="store_true", help="Skip log1p transform")
    args = parser.parse_args()

    print(f"Loading: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)

    # Ensure QC metrics
    if "total_counts" not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(adata, inplace=True)

    print(f"Samples: {sorted(set(adata.obs['sample_id'].astype(int)))}")
    fig = plot_per_sample(adata, feature=args.feature, log_transform=(not args.no_log))
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
