# Visium HD Sample Split

Automatically split multiple tissue samples/cores in a Visium HD capture area by connected components on the bin grid — no image processing needed.

## What it does

Given Visium HD binned output (from Space Ranger), this tool:

1. Builds a 2D binary grid from `array_row` × `array_col` of foreground bins
2. Uses morphological opening to break thin bridges between adjacent samples
3. Runs connected components to label each tissue piece
4. Merges small fragments into the nearest major sample
5. Optionally splits over-sized components via K-means
6. Outputs an AnnData with a `sample_id` column

## Input requirements

Your `--base` directory must contain:

```
base/
├── spatial/
│   ├── tissue_hires_image.png        # optional, for visualization only
│   └── tissue_lowres_image.png       # optional fallback
└── binned_outputs/
    └── square_<bin>/
        ├── filtered_feature_bc_matrix.h5
        └── spatial/
            ├── tissue_positions.parquet
            └── scalefactors_json.json
```

### Required fields in `tissue_positions.parquet`

| Column | Type | Description |
|--------|------|-------------|
| `barcode` | str | Must match AnnData.obs_names |
| `in_tissue` | int | 1 = foreground, 0 = background |
| `array_row` | int | Bin row index on the 2D grid |
| `array_col` | int | Bin column index on the 2D grid |

These fields are standard Space Ranger output — no modification needed.

### Implicit assumptions

- `array_row`/`array_col` form a regular grid (neighboring bins differ by 1)
- Tissue samples are separated by gaps in the bin grid
- Barcodes are consistent between `.h5` and `.parquet`

## Installation

```bash
pip install -r requirements.txt
```

## Quick start

```bash
# If you know the number of samples (recommended):
python segment_ta.py \
    --base /path/to/visium_hd_data \
    --bin 016um \
    --expected-samples 3 \
    --open-radius 1

# If you don't know the number of samples:
python segment_ta.py \
    --base /path/to/visium_hd_data \
    --bin 016um \
    --open-radius 1 \
    --size-prior

# Try different opening radii if samples remain connected:
python segment_ta.py \
    --base /path/to/visium_hd_data \
    --bin 016um \
    --expected-samples 5 \
    --open-radius 2
```

## Output

```
base/_sample_split_v3/square_<bin>/<tag>/
├── sample_split_overview.png          # 3-panel diagnostic figure
├── visium_hd_sample_split.h5ad        # AnnData with sample_id column
├── spot_sample_assignment.parquet     # per-barcode assignment
├── spot_sample_assignment.csv         # same, CSV format
├── sample_summary.csv                 # per-sample statistics
└── sample_label_grid.npy              # integer label grid (npy)
```

### Per-sample AnnData after loading

```python
import scanpy as sc
adata = sc.read_h5ad("visium_hd_sample_split.h5ad")

# Each spot now has:
adata.obs["sample_id"]    # integer, 0 = background
adata.obs["sample"]       # "sample_001", "sample_002", ...

# Extract a single sample:
sample1 = adata[adata.obs["sample_id"] == 1].copy()
```

## Algorithm

```
① Build foreground mask
   foreground = in_tissue AND total_counts >= min_counts
   → 2D boolean grid (array_row × array_col)

② Optional: morphological opening
   if open_radius > 0: disk(open_radius) erosion then dilation
   → breaks 1-2 bin wide bridges between samples

③ Conservative connected components
   connectivity=1 (4-neighbor) by default
   → each tissue piece gets a unique label

④ Component filtering
   Remove components < min_bins bins
   → eliminates noise / tiny fragments

⑤ Post-processing (fragment merge + large split)
   If --expected-samples N:
     top-N largest components → major samples
     all others → merge to nearest major (by grid distance)
   
   If --size-prior:
     auto-estimate target bin count
     small fragments (< 35% target) → merge to nearest major
     over-sized components (> 1.8x target) → K-means split

⑥ Renumber samples left-to-right by median array_col
```

## Key parameters

| Parameter | Default | When to adjust |
|-----------|---------|----------------|
| `--expected-samples` | 0 (off) | **Set this if you know N.** Small fragments merge to nearest major. |
| `--open-radius` | 0 | Increase (1–2) if samples are connected by thin bridges. |
| `--min-bins` | auto (30 for 16µm) | Increase to filter more noise; decrease to keep small cores. |
| `--min-counts` | 1.0 | Raise to exclude low-quality bins from foreground. |
| `--connectivity` | 1 (4-neighbor) | Use 2 (8-neighbor) for tighter connections. |
| `--size-prior` | off | Enable when you don't know the sample count. |

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Two samples merged into one | Increase `--open-radius` (try 1, then 2) |
| One sample split into pieces | Decrease `--open-radius` or set it to 0 |
| Too many tiny samples | Increase `--min-bins` or set `--expected-samples` |
| Large sample not split | Lower `--large-component-frac` (default 1.8) |
| Small islands dropped | Use `--keep-far-islands` |

## Tested dataset

The algorithm was developed and tested on:

- **10x Visium HD 11mm Human Tissue Array (FFPE)**
- ~328,000 bins at 16µm resolution
- 3–5 distinct tissue samples per capture area
- Download: [10x Genomics Datasets](https://www.10xgenomics.com/datasets)

## License

MIT
