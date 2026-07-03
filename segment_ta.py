#!/usr/bin/env python3
"""
Visium HD Tissue Array — sample/core 拆分与表达矩阵映射 v3
=====================================================

v3 的核心变化：
  1. 默认不再对 foreground 做 binary_closing / fill holes 再连通域标记。
     closing 很容易把相邻 sample 之间的窄缝补上，导致两个样本被合成一个。
  2. 先在原始 bin foreground 上做 conservative connected components。
  3. 可选用 opening 打断 1-2 个 bin 宽的细桥。
  4. 用“样本面积大致相近 / 预期样本数”的先验，把小孤岛合并到最近的大 sample，
     而不是把小孤岛当成独立 sample。

推荐先跑：
  python segment_ta_modified_v3.py \
      --base /share/project/visium_hd_human_TA \
      --bin 016um \
      --open-radius 0 \
      --expected-samples 4

如果两个样本仍然连在一起：
  python segment_ta_modified_v3.py \
      --base /share/project/visium_hd_human_TA \
      --bin 016um \
      --open-radius 1 \
      --expected-samples 4

如果不知道 expected-samples：
  python segment_ta_modified_v3.py \
      --base /share/project/visium_hd_human_TA \
      --bin 016um \
      --open-radius 1 \
      --size-prior
"""

import argparse
import json
import warnings
from pathlib import Path
from typing import Iterable, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image
from scipy import sparse
from scipy.spatial import cKDTree
from skimage.measure import label, regionprops
from skimage.morphology import binary_opening, remove_small_objects, disk

warnings.filterwarnings("ignore")


# -----------------------------
# 基础工具
# -----------------------------


def _first_existing(paths: Iterable[Path], required: bool = True) -> Optional[Path]:
    for p in paths:
        if p.exists():
            return p
    if required:
        raise FileNotFoundError("未找到文件，检查过这些路径：\n" + "\n".join(str(p) for p in paths))
    return None


def _matrix_row_sum(x) -> np.ndarray:
    s = x.sum(axis=1)
    if sparse.issparse(s):
        return np.asarray(s).ravel()
    if hasattr(s, "A1"):
        return s.A1
    return np.asarray(s).ravel()


def _auto_params(bin_size: str) -> Tuple[int, int, int]:
    """返回 min_bins, rescue_max_dist, small_merge_max_dist，单位均是 bin grid。"""
    if bin_size == "016um":
        return 80, 4, 80
    if bin_size == "008um":
        return 300, 8, 160
    if bin_size == "002um":
        return 3000, 25, 600
    return 80, 4, 80


def ensure_position_columns(pos: pd.DataFrame) -> pd.DataFrame:
    required = [
        "barcode", "in_tissue", "array_row", "array_col",
        "pxl_row_in_fullres", "pxl_col_in_fullres",
    ]
    missing = [c for c in required if c not in pos.columns]
    if missing:
        raise ValueError(
            "tissue_positions.parquet 缺少必需列：" + ", ".join(missing) +
            "\n当前列为：" + ", ".join(pos.columns)
        )
    return pos


def _renumber_left_to_right(sample_id: np.ndarray, array_row: np.ndarray, array_col: np.ndarray) -> np.ndarray:
    old_ids = [int(x) for x in np.unique(sample_id) if int(x) > 0]

    def key_func(old_id: int):
        m = sample_id == old_id
        return (float(np.median(array_col[m])), float(np.median(array_row[m])))

    old_ids = sorted(old_ids, key=key_func)
    remap = {old: new for new, old in enumerate(old_ids, start=1)}
    return np.array([remap.get(int(x), 0) for x in sample_id], dtype=np.int32)


def _simple_kmeans_farthest(coords: np.ndarray, k: int, n_iter: int = 80) -> np.ndarray:
    coords = coords.astype(np.float64, copy=False)
    n = coords.shape[0]
    if k <= 1 or n <= k:
        return np.zeros(n, dtype=np.int32)

    centers = np.empty((k, 2), dtype=np.float64)
    centers[0] = coords.mean(axis=0)
    dist2 = ((coords - centers[0]) ** 2).sum(axis=1)
    for i in range(1, k):
        centers[i] = coords[int(np.argmax(dist2))]
        dist2 = np.minimum(dist2, ((coords - centers[i]) ** 2).sum(axis=1))

    labels = np.zeros(n, dtype=np.int32)
    for _ in range(n_iter):
        d2 = ((coords[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
        new_labels = np.argmin(d2).astype(np.int32) if d2.ndim == 1 else np.argmin(d2, axis=1).astype(np.int32)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels
        for j in range(k):
            m = labels == j
            if m.any():
                centers[j] = coords[m].mean(axis=0)
    return labels


# -----------------------------
# 数据读取
# -----------------------------


def load_visium_hd(base: Path, bin_size: str):
    bin_dir = base / "binned_outputs" / f"square_{bin_size}"
    if not bin_dir.exists():
        raise FileNotFoundError(
            f"找不到 {bin_dir}\n"
            "请确认 --base 指向解压后的 Visium HD 数据根目录。"
        )

    matrix_h5 = _first_existing([
        bin_dir / "filtered_feature_bc_matrix.h5",
        bin_dir / "raw_feature_bc_matrix.h5",
    ])
    pos_path = _first_existing([
        bin_dir / "spatial" / "tissue_positions.parquet",
        base / "spatial" / "tissue_positions.parquet",
    ])
    scalefactors_path = _first_existing([
        bin_dir / "spatial" / "scalefactors_json.json",
        base / "spatial" / "scalefactors_json.json",
    ])
    hires_image_path = _first_existing([
        base / "spatial" / "tissue_hires_image.png",
        bin_dir / "spatial" / "tissue_hires_image.png",
    ], required=False)
    lowres_image_path = _first_existing([
        base / "spatial" / "tissue_lowres_image.png",
        bin_dir / "spatial" / "tissue_lowres_image.png",
    ], required=False)

    print("加载数据：")
    print(f"  base:              {base}")
    print(f"  bin_dir:           {bin_dir}")
    print(f"  matrix_h5:         {matrix_h5}")
    print(f"  tissue_positions:  {pos_path}")
    print(f"  scalefactors:      {scalefactors_path}")
    print(f"  hires_image:       {hires_image_path}")

    adata = sc.read_10x_h5(matrix_h5)
    if not adata.var_names.is_unique:
        adata.var_names_make_unique()

    pos = ensure_position_columns(pd.read_parquet(pos_path))
    with open(scalefactors_path, "r", encoding="utf-8") as f:
        sf = json.load(f)

    img = None
    if hires_image_path is not None:
        img = np.array(Image.open(hires_image_path).convert("RGB"))
    elif lowres_image_path is not None:
        img = np.array(Image.open(lowres_image_path).convert("RGB"))
        print("  注意：未找到 tissue_hires_image.png，使用 tissue_lowres_image.png 做可视化。")

    print(f"  AnnData:           {adata.n_obs:,} bins × {adata.n_vars:,} genes")
    if img is not None:
        print(f"  image:             {img.shape[1]} × {img.shape[0]}")
    return adata, pos, sf, img


# -----------------------------
# 组件后处理
# -----------------------------


def estimate_target_bins(counts: np.ndarray, min_bins: int, target_bins: Optional[int]) -> float:
    if target_bins is not None and target_bins > 0:
        return float(target_bins)
    counts = np.asarray(counts, dtype=float)
    counts = counts[counts >= min_bins]
    if len(counts) == 0:
        return float(min_bins)
    # 去掉明显小碎片后，用偏大的分位数估计正常 sample 面积。
    max_count = float(counts.max())
    candidates = counts[counts >= max(min_bins, 0.15 * max_count)]
    if len(candidates) == 0:
        candidates = counts
    return float(np.median(candidates))


def merge_fragments_to_major_samples(
    sample_id: np.ndarray,
    array_row: np.ndarray,
    array_col: np.ndarray,
    expected_samples: int = 0,
    min_bins: int = 80,
    target_bins: Optional[int] = None,
    small_island_frac: float = 0.35,
    major_frac: float = 0.45,
    large_component_frac: float = 1.8,
    max_split_k: int = 6,
    merge_max_dist: float = 80.0,
    discard_far_islands: bool = True,
) -> np.ndarray:
    """
    小碎片合并策略：
      - 如果给了 expected_samples：最大 expected_samples 个 component 作为 major samples；其他 component 合并到最近 major。
      - 如果没给 expected_samples：根据 target area 识别 major / small islands。
      - 明显过大的 component 可按面积比例用 k-means 拆开，作为兜底。
    """
    out = sample_id.copy().astype(np.int32)
    ids, cnts = np.unique(out[out > 0], return_counts=True)
    if len(ids) == 0:
        return out

    print("\n面积/样本数先验后处理：")

    # 0. 对明显过大的 component 兜底拆分。
    target = estimate_target_bins(cnts, min_bins=min_bins, target_bins=target_bins)
    print(f"  target bins:       {target:.1f}")

    next_id = int(out.max()) + 1
    split_events = []
    ids, cnts = np.unique(out[out > 0], return_counts=True)
    for cid, n in zip(ids.tolist(), cnts.tolist()):
        if n <= large_component_frac * target:
            continue
        k = int(np.round(n / target))
        k = max(2, min(k, max_split_k, int(n)))
        idx = np.where(out == cid)[0]
        coords = np.column_stack([array_row[idx], array_col[idx]])
        labs = _simple_kmeans_farthest(coords, k=k)
        for t, lab_id in enumerate(sorted(np.unique(labs).tolist())):
            new_id = int(cid) if t == 0 else next_id
            if t > 0:
                next_id += 1
            out[idx[labs == lab_id]] = new_id
        split_events.append((int(cid), int(n), int(k)))

    if split_events:
        for cid, n, k in split_events:
            print(f"  split large:       component {cid} ({n:,} bins) -> {k} parts")
    else:
        print("  split large:       none")

    ids, cnts = np.unique(out[out > 0], return_counts=True)
    count_map = {int(i): int(c) for i, c in zip(ids, cnts)}

    if expected_samples and expected_samples > 0:
        sorted_ids = sorted(count_map, key=lambda x: count_map[x], reverse=True)
        major_ids = sorted_ids[:expected_samples]
        fragment_ids = sorted_ids[expected_samples:]
        print(f"  expected samples:  {expected_samples}")
        print("  major components:  " + ", ".join([f"{i}({count_map[i]})" for i in major_ids]))
    else:
        major_ids = [i for i, c in count_map.items() if c >= major_frac * target]
        fragment_ids = [i for i, c in count_map.items() if i not in major_ids]
        print(f"  major threshold:   >= {major_frac * target:.1f} bins")
        print(f"  small threshold:   <  {small_island_frac * target:.1f} bins")
        print(f"  major components:  {len(major_ids)}")

    if len(major_ids) == 0:
        print("  warning:           no major component found; skip fragment merge")
        return out

    major_mask = np.isin(out, major_ids)
    major_xy = np.column_stack([array_row[major_mask], array_col[major_mask]])
    major_sid = out[major_mask]
    tree = cKDTree(major_xy)

    merged = 0
    discarded = 0
    kept_far = 0
    for sid in fragment_ids:
        idx = np.where(out == sid)[0]
        if len(idx) == 0:
            continue
        xy = np.column_stack([array_row[idx], array_col[idx]])
        dists, nn = tree.query(xy, k=1)
        j = int(np.argmin(dists))
        min_dist = float(dists[j])
        target_sid = int(major_sid[nn[j]])
        # expected_samples 模式下，碎片通常就是边缘组织或小岛，默认合并到最近大样本。
        allow_merge = (expected_samples and expected_samples > 0) or (min_dist <= merge_max_dist)
        if allow_merge:
            out[idx] = target_sid
            merged += 1
        elif discard_far_islands:
            out[idx] = 0
            discarded += 1
        else:
            kept_far += 1

    print(f"  merge fragments:   {merged}/{len(fragment_ids)}")
    if discarded:
        print(f"  discard far:       {discarded}")
    if kept_far:
        print(f"  keep far:          {kept_far}")
    return out


# -----------------------------
# 主分割算法
# -----------------------------


def split_samples_by_conservative_components(
    adata: sc.AnnData,
    pos: pd.DataFrame,
    min_bins: int,
    min_counts: float,
    open_radius: int = 0,
    connectivity: int = 1,
    rescue_unassigned: bool = True,
    rescue_max_dist: float = 4.0,
    expected_samples: int = 0,
    use_size_prior: bool = False,
    target_bins: Optional[int] = None,
    small_island_frac: float = 0.35,
    major_frac: float = 0.45,
    large_component_frac: float = 1.8,
    max_split_k: int = 6,
    small_merge_max_dist: float = 80.0,
    discard_far_islands: bool = True,
):
    pos0 = pos.copy().set_index("barcode", drop=False)
    common = [bc for bc in adata.obs_names if bc in pos0.index]
    if len(common) == 0:
        raise ValueError("AnnData.obs_names 与 tissue_positions.parquet 的 barcode 没有交集。")

    adata2 = adata[common].copy()
    pos2 = pos0.loc[common].copy()

    counts = _matrix_row_sum(adata2.X)
    in_tissue = pos2["in_tissue"].to_numpy().astype(int) == 1
    foreground = in_tissue & (counts >= min_counts)

    array_row = pos2["array_row"].astype(int).to_numpy()
    array_col = pos2["array_col"].astype(int).to_numpy()
    r0, c0 = int(array_row.min()), int(array_col.min())
    rr, cc = array_row - r0, array_col - c0
    grid_shape = (int(rr.max()) + 1, int(cc.max()) + 1)

    fg_grid = np.zeros(grid_shape, dtype=bool)
    fg_grid[rr[foreground], cc[foreground]] = True

    print("\n构建 bin foreground：")
    print(f"  grid shape:        {grid_shape[1]} × {grid_shape[0]}  [col × row]")
    print(f"  in_tissue bins:    {int(in_tissue.sum()):,}")
    print(f"  foreground bins:   {int(foreground.sum()):,}  (in_tissue 且 counts >= {min_counts})")
    print(f"  opening radius:    {open_radius}")
    print(f"  connectivity:      {connectivity}  (1=4邻域, 2=8邻域)")

    work_grid = fg_grid.copy()
    if open_radius > 0:
        # opening = erosion + dilation；它会打断细桥，但不会像 closing 那样把样本间空隙填上。
        work_grid = binary_opening(work_grid, footprint=disk(open_radius))

    work_grid = remove_small_objects(work_grid, min_size=min_bins)
    raw_lab = label(work_grid, connectivity=connectivity)

    sample_id = np.zeros(len(pos2), dtype=np.int32)
    sample_id[foreground] = raw_lab[rr[foreground], cc[foreground]]

    # 如果 opening 后某些 foreground bin 没有 label，保守地补回最近 label。
    if rescue_unassigned and np.any(sample_id > 0):
        assigned_mask = sample_id > 0
        rescue_mask = foreground & (sample_id == 0)
        if rescue_mask.any():
            assigned_xy = np.column_stack([array_row[assigned_mask], array_col[assigned_mask]])
            rescue_xy = np.column_stack([array_row[rescue_mask], array_col[rescue_mask]])
            tree = cKDTree(assigned_xy)
            dists, idxs = tree.query(rescue_xy, k=1)
            valid = dists <= rescue_max_dist
            assigned_ids = sample_id[assigned_mask]
            rescue_indices = np.where(rescue_mask)[0]
            sample_id[rescue_indices[valid]] = assigned_ids[idxs[valid]]
            print(f"  rescue opening:    {int(valid.sum()):,}/{int(rescue_mask.sum()):,}")

    # 过滤一次极小 component。
    ids, cnts = np.unique(sample_id[sample_id > 0], return_counts=True)
    tiny = set(ids[cnts < min_bins].tolist())
    if tiny:
        sample_id[np.isin(sample_id, list(tiny))] = 0

    if expected_samples > 0 or use_size_prior:
        sample_id = merge_fragments_to_major_samples(
            sample_id=sample_id,
            array_row=array_row,
            array_col=array_col,
            expected_samples=expected_samples,
            min_bins=min_bins,
            target_bins=target_bins,
            small_island_frac=small_island_frac,
            major_frac=major_frac,
            large_component_frac=large_component_frac,
            max_split_k=max_split_k,
            merge_max_dist=small_merge_max_dist,
            discard_far_islands=discard_far_islands,
        )

    sample_id[~foreground] = 0
    sample_id = _renumber_left_to_right(sample_id, array_row=array_row, array_col=array_col)

    label_grid = np.zeros(grid_shape, dtype=np.int32)
    label_grid[rr[foreground], cc[foreground]] = sample_id[foreground]

    pos2["counts"] = counts
    pos2["is_foreground"] = foreground
    pos2["sample_id"] = sample_id
    pos2["sample_name"] = [f"sample_{x:03d}" if x > 0 else "background" for x in sample_id]

    print("\nSample/core 检测结果：")
    sample_ids = sorted([int(x) for x in np.unique(sample_id) if int(x) > 0])
    print(f"  samples:           {len(sample_ids)}")
    for sid in sample_ids:
        print(f"  sample_{sid:03d}:      {(sample_id == sid).sum():,} bins")

    grid_info = {
        "row_offset": r0,
        "col_offset": c0,
        "array_row": array_row,
        "array_col": array_col,
        "rr": rr,
        "cc": cc,
        "foreground_grid": fg_grid,
        "opened_grid": work_grid,
        "label_grid": label_grid,
    }
    return adata2, pos2, grid_info


# -----------------------------
# 写出与可视化
# -----------------------------


def add_spatial_to_adata(adata: sc.AnnData, pos: pd.DataFrame, sf: dict):
    pos_indexed = pos.set_index("barcode", drop=False)
    common = [bc for bc in adata.obs_names if bc in pos_indexed.index]
    adata = adata[common].copy()
    p = pos_indexed.loc[common]

    for col in ["sample_id", "sample_name", "counts", "is_foreground"]:
        if col in p.columns:
            if col == "sample_name":
                adata.obs["sample"] = p[col].values
            else:
                adata.obs[col] = p[col].values

    adata.obsm["spatial_fullres"] = p[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
    scale = float(sf.get("tissue_hires_scalef", 1.0))
    adata.obsm["spatial"] = adata.obsm["spatial_fullres"] * scale
    return adata


def sample_summary_table(pos: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sid in sorted([int(x) for x in pos["sample_id"].unique() if int(x) > 0]):
        p = pos[pos["sample_id"] == sid]
        rows.append({
            "sample_id": sid,
            "sample_name": f"sample_{sid:03d}",
            "n_bins": int(len(p)),
            "sum_counts": float(p["counts"].sum()) if "counts" in p else np.nan,
            "median_array_row": float(p["array_row"].median()),
            "median_array_col": float(p["array_col"].median()),
            "min_array_row": int(p["array_row"].min()),
            "max_array_row": int(p["array_row"].max()),
            "min_array_col": int(p["array_col"].min()),
            "max_array_col": int(p["array_col"].max()),
        })
    return pd.DataFrame(rows)


def visualize_sample_split(img, pos: pd.DataFrame, grid_info: dict, sf: dict,
                           out_png: Path, overlay_he: bool = True):
    label_grid = grid_info["label_grid"]
    fg_grid = grid_info["foreground_grid"]
    opened_grid = grid_info["opened_grid"]
    sample_ids = sorted([int(x) for x in np.unique(label_grid) if int(x) > 0])
    n_samples = len(sample_ids)

    n_cols = 4 if (img is not None and overlay_he) else 3
    fig, axes = plt.subplots(1, n_cols, figsize=(7 * n_cols, 7))
    if n_cols == 3:
        axes = list(axes)

    axes[0].imshow(fg_grid, cmap="gray")
    axes[0].set_title("Raw foreground bins")
    axes[0].axis("off")

    axes[1].imshow(opened_grid, cmap="gray")
    axes[1].set_title("After optional opening")
    axes[1].axis("off")

    axes[2].imshow(label_grid, interpolation="nearest")
    axes[2].set_title(f"Final labels: {n_samples} samples")
    axes[2].axis("off")

    if img is not None and overlay_he:
        axes[3].imshow(img)
        scale = float(sf.get("tissue_hires_scalef", 1.0))
        foreground = pos["sample_id"].to_numpy() > 0
        p = pos[foreground]
        n_colors = max(n_samples + 1, 20)
        cmap = plt.cm.get_cmap("tab20", n_colors)
        for sid in sorted(p["sample_id"].unique()):
            q = p[p["sample_id"] == sid]
            x = q["pxl_col_in_fullres"].to_numpy(dtype=float) * scale
            y = q["pxl_row_in_fullres"].to_numpy(dtype=float) * scale
            axes[3].scatter(x, y, s=0.8, alpha=0.55, color=cmap((int(sid) - 1) % n_colors), rasterized=True)
        axes[3].set_title("Sample assignment on H&E")
        axes[3].axis("off")

    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure:            {out_png}")


def write_outputs(base: Path, bin_size: str, adata: sc.AnnData, pos: pd.DataFrame, grid_info: dict, img, sf: dict, tag: str):
    out_dir = base / "_sample_split_v3" / f"square_{bin_size}" / tag
    out_dir.mkdir(parents=True, exist_ok=True)

    h5ad_path = out_dir / "visium_hd_sample_split.h5ad"
    pos_parquet = out_dir / "spot_sample_assignment.parquet"
    pos_csv = out_dir / "spot_sample_assignment.csv"
    summary_csv = out_dir / "sample_summary.csv"
    label_npy = out_dir / "sample_label_grid.npy"
    fg_npy = out_dir / "foreground_grid.npy"
    opened_npy = out_dir / "opened_grid.npy"
    fig_path = out_dir / "sample_split_overview.png"

    adata.write_h5ad(h5ad_path)
    pos.reset_index(drop=True).to_parquet(pos_parquet)
    pos.reset_index(drop=True).to_csv(pos_csv, index=False)
    sample_summary_table(pos).to_csv(summary_csv, index=False)
    np.save(label_npy, grid_info["label_grid"])
    np.save(fg_npy, grid_info["foreground_grid"])
    np.save(opened_npy, grid_info["opened_grid"])
    visualize_sample_split(img, pos, grid_info, sf, fig_path,
                           overlay_he=(not args.no_he_image))

    print("\n输出文件：")
    print(f"  AnnData:           {h5ad_path}")
    print(f"  assignment:        {pos_parquet}")
    print(f"  summary:           {summary_csv}")
    print(f"  label_grid:        {label_npy}")
    print(f"  overview:          {fig_path}")
    return out_dir


def main():
    parser = argparse.ArgumentParser(description="Visium HD tissue array sample/core splitting v3")
    parser.add_argument("--base", type=str, required=True, help="解压后的 Visium HD 数据根目录")
    parser.add_argument("--bin", default="016um", choices=["002um", "008um", "016um"])
    parser.add_argument("--min-bins", type=int, default=None, help="小于该 bin 数的初始 component 删除")
    parser.add_argument("--min-counts", type=float, default=1.0, help="bin 总 counts 至少达到该值才作为 foreground")
    parser.add_argument("--open-radius", type=int, default=0, help="用 opening 打断细桥；0=不做。16um 可试 1-2")
    parser.add_argument("--connectivity", type=int, default=1, choices=[1, 2], help="1=4邻域，更保守；2=8邻域")
    parser.add_argument("--no-rescue", action="store_true", help="不补回 opening 后丢失的 foreground bin")
    parser.add_argument("--no-he-image", action="store_true", help="可视化时不叠加 H&E 背景")
    parser.add_argument("--rescue-max-dist", type=float, default=None, help="opening 后补回距离阈值，单位 bin grid")

    parser.add_argument("--expected-samples", type=int, default=0, help="已知样本数时强烈推荐设置，例如 4")
    parser.add_argument("--size-prior", action="store_true", help="不知道样本数时启用面积先验")
    parser.add_argument("--target-bins", type=int, default=None, help="正常 sample 的预期 bin 数；不设则自动估计")
    parser.add_argument("--small-island-frac", type=float, default=0.35)
    parser.add_argument("--major-frac", type=float, default=0.45)
    parser.add_argument("--large-component-frac", type=float, default=1.8)
    parser.add_argument("--max-split-k", type=int, default=6)
    parser.add_argument("--small-merge-max-dist", type=float, default=None, help="小碎片合并到最近 major 的最大距离，单位 bin grid")
    parser.add_argument("--keep-far-islands", action="store_true", help="远离 major 的小岛保留为独立 sample；默认丢弃")
    parser.add_argument("--tag", type=str, default="run")
    args = parser.parse_args()

    auto_min_bins, auto_rescue, auto_merge = _auto_params(args.bin)
    min_bins = args.min_bins if args.min_bins is not None else auto_min_bins
    rescue_max_dist = args.rescue_max_dist if args.rescue_max_dist is not None else auto_rescue
    small_merge_max_dist = args.small_merge_max_dist if args.small_merge_max_dist is not None else auto_merge

    base = Path(args.base)
    adata, pos, sf, img = load_visium_hd(base, args.bin)

    print("\n参数：")
    print(f"  min_bins:          {min_bins}")
    print(f"  min_counts:        {args.min_counts}")
    print(f"  open_radius:       {args.open_radius}")
    print(f"  rescue_max_dist:   {rescue_max_dist}")
    print(f"  expected_samples:  {args.expected_samples}")
    print(f"  size_prior:        {args.size_prior}")
    print(f"  merge_max_dist:    {small_merge_max_dist}")

    adata2, pos2, grid_info = split_samples_by_conservative_components(
        adata=adata,
        pos=pos,
        min_bins=min_bins,
        min_counts=args.min_counts,
        open_radius=args.open_radius,
        connectivity=args.connectivity,
        rescue_unassigned=not args.no_rescue,
        rescue_max_dist=rescue_max_dist,
        expected_samples=args.expected_samples,
        use_size_prior=args.size_prior,
        target_bins=args.target_bins,
        small_island_frac=args.small_island_frac,
        major_frac=args.major_frac,
        large_component_frac=args.large_component_frac,
        max_split_k=args.max_split_k,
        small_merge_max_dist=small_merge_max_dist,
        discard_far_islands=not args.keep_far_islands,
    )

    adata2 = add_spatial_to_adata(adata2, pos2, sf)
    tag = args.tag
    if args.expected_samples > 0:
        tag = f"{tag}_N{args.expected_samples}"
    if args.open_radius > 0:
        tag = f"{tag}_open{args.open_radius}"
    write_outputs(base, args.bin, adata2, pos2, grid_info, img, sf, tag=tag)


if __name__ == "__main__":
    main()
