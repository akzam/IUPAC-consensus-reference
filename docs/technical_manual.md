# Technical Manual — calculate_fp_burden.py

**Version:** 5.3.7  
**Language:** Python 3.8+  
**Dependencies:** pandas >= 1.3.0, numpy >= 1.20.0  
**License:** MIT

---

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Data Flow](#data-flow)
3. [Core Algorithms](#core-algorithms)
4. [Percentage-Point Derivation](#percentage-point-derivation)
5. [TP:FP Ratio & Net Benefit](#tpfp-ratio--net-benefit)
6. [Error Handling & Validation](#error-handling--validation)
7. [Aggregation Engine](#aggregation-engine)
8. [Output Format Specifications](#output-format-specifications)
9. [Testing & Validation](#testing--validation)
10. [Performance Considerations](#performance-considerations)
11. [API Reference](#api-reference)

---

## Architecture Overview

```
Input (Manifest or CLI) → Validation → hap.py Loader → compute_burden()
                                                          |
              ┌─────────────────┼─────────────────┐
              ▼                 ▼                 ▼
        Raw CSV output    Aggregate (median)    Format (MD/TEX/CSV)
              │                 │                 │
              └─────────────────┴─────────────────┘
                                │
                    ┌───────────┴───────────┐
                    ▼                       ▼
            Summary CSV          Reproducibility Audit
```

### Design principles

1. **Reproducibility-first:** Every derived metric is traceable to raw hap.py values
2. **Fail-fast validation:** Input checks before computation to avoid silent errors
3. **Defensive aggregation:** Handles missing samples, empty groups, and edge cases
4. **Self-documenting outputs:** Generated tables include interpretation notes
5. **Multiple output formats:** Raw data, summary statistics, human-readable tables, and machine-readable CSV
6. **Edge-case robustness:** Zero counts, negative values, and `inf`/`nan` handled gracefully

---

## Data Flow

### Stage 1: Input ingestion

**hap.py loader** (`load_happy()`)
- Auto-detects delimiter (comma or tab)
- Strips whitespace from column names
- Validates required columns against `REQUIRED_HAPPY_COLUMNS`
- Checks numeric types and value ranges

**Manifest loader** (`load_manifest()`)
- Normalizes column names (lowercase, spaces → underscores)
- Validates against `REQUIRED_MANIFEST_COLUMNS`
- Checks for missing values and duplicate comparisons
- Verifies all referenced files exist and are readable

### Stage 2: Pairwise computation

For each baseline-consensus pair:

```python
df_merged = df_cons.merge(df_base, on=MERGE_KEYS, suffixes=("_cons", "_base"))
```

Merge keys: `Type`, `Subtype`, `Subset`, `Filter`, `Genotype`

**Unmatched row handling:**
- Rows present in only one file are flagged (`in_both = False`)
- Comparative metrics (pp gains, additional counts) are NaN for unmatched rows
- Warning emitted to stderr with count of unmatched rows

### Stage 3: Metric derivation

All metrics computed row-wise. No sample-level aggregation at this stage.

### Stage 4: Aggregation (batch mode only)

Group by: `Caller`, `Condition`, `Type`, `Subset`, `Filter`, `Genotype`

Statistics per group:
- `median`: Central tendency across samples
- `min`, `max`: Inter-sample range
- `count`: Number of samples contributing

### Stage 5: Formatting

Six output files produced:

| # | File | Format | Contents |
|---|------|--------|----------|
| 1 | `_raw_per_sample.csv` | CSV | All columns, all rows, no formatting |
| 2 | `_summary_median_range.csv` | CSV | Median/min/max/count + formatted strings |
| 3 | `_table.txt` | Markdown | Human-readable table with renamed columns |
| 4 | `_table.tex` | LaTeX | Journal-ready table |
| 5 | `_table.csv` | CSV | Machine-readable table with clean column names |
| 6 | `_reproducibility.txt` | Plain text | Formula documentation and data provenance |

Plus optional overall and median-of-medians outputs.

---

## Core Algorithms

### 1. Percentage-point gain

```python
# Raw hap.py outputs (0-1 scale)
recall_base    = METRIC.Recall_base      # e.g., 0.985
recall_cons    = METRIC.Recall_consensus  # e.g., 0.988

# Convert to percentage scale
recall_base_pct = recall_base * 100      # 98.5
recall_cons_pct = recall_cons * 100      # 98.8

# Percentage-point difference
recall_gain_pp  = recall_cons_pct - recall_base_pct   # 0.3 pp

# Verification (must match to machine precision)
recall_gain_pp_check = (recall_cons - recall_base) * 100
assert abs(recall_gain_pp - recall_gain_pp_check) < 1e-10
```

### 2. Absolute burden

```python
additional_tp = QUERY.TP_consensus - QUERY.TP_baseline
additional_fp = QUERY.FP_consensus - QUERY.FP_baseline
additional_fn = (TRUTH.TOTAL_cons - QUERY.TP_cons) - (TRUTH.TOTAL_base - QUERY.TP_base)
```

### 3. TP:FP ratio and Net Benefit

**Raw ratio** (may produce `inf`/`nan`):

```python
with np.errstate(divide='ignore', invalid='ignore'):
    ratio_raw = np.where(
        merged["Additional_FP"] != 0,
        merged["Additional_TP"] / merged["Additional_FP"],
        np.where(merged["Additional_TP"] != 0, np.inf, np.nan)
    )
```

**Adjusted ratio** (conditional Laplace smoothing):

```python
merged["TP_FP_ratio_adj"] = np.where(
    (merged["Additional_TP"] == 0) & (merged["Additional_FP"] == 0),
    np.nan,
    np.where(
        merged["Additional_FP"] >= 0,
        (merged["Additional_TP"] + 1) / (merged["Additional_FP"] + 1),
        merged["Additional_TP"] / merged["Additional_FP"]
    )
)
```

**Net benefit** (signed, handles all cases):

```python
net_benefit = additional_tp - additional_fp
```

### 4. Per-megabase normalization

```python
conf_size_mb = Subset.IS_CONF.Size / 1_000_000
additional_tp_per_mbp = additional_tp / conf_size_mb   # Avoid division by zero
```

### 5. Per-100-truth normalization

```python
truth_total = TRUTH.TOTAL_base
additional_tp_per_100 = (additional_tp / truth_total) * 100
```

---

## Percentage-Point Derivation

### Why explicit derivation matters

hap.py reports `METRIC.Recall` and `METRIC.Precision` as **proportions (0–1)**. Manuscripts typically report changes as **percentage points (pp)**. The conversion is:

```
pp_gain = (proportion_consensus - proportion_baseline) * 100
```

This is mathematically trivial but **operationally critical** for reproducibility. The script preserves the full derivation chain:

| Stage | Column | Value | Scale |
|-------|--------|-------|-------|
| Raw hap.py | `METRIC.Recall` | 0.733 | 0–1 |
| Preserved | `Recall_cons_raw` | 0.733 | 0–1 |
| Converted | `Recall_cons_pct` | 73.3 | 0–100 |
| Derived | `Recall_gain_pp` | 3.3 | pp |
| Verified | `Recall_gain_pp_check` | 3.3 | pp |

### Drift detection

The script computes each pp gain **two ways**:

1. **Primary:** `(cons_pct - base_pct)`
2. **Check:** `(cons_raw - base_raw) * 100`

If `abs(primary - check) > 1e-10`, a warning is emitted. This catches:
- Floating-point anomalies
- Data type corruption (e.g., integer division)
- Manual editing errors in input files

### Audit trail

The `_reproducibility.txt` output documents:
- Exact formulas used
- Column names where raw values can be found
- Verification procedure for independent checking

---

## TP:FP Ratio & Net Benefit

### The problem with raw ratios

When `Additional_FP = 0`, the raw ratio `Additional_TP / Additional_FP` is mathematically undefined (`inf`). This is common when:
- A sample has no false positives in a small region
- The consensus reference strictly improves precision
- Baseline FP count is zero (rare but possible)

### Solution: Three complementary metrics

| Metric | Formula | Handles zeros? | Use case |
|--------|---------|----------------|----------|
| **Raw ratio** | `TP / FP` | No (inf/nan) | Per-sample transparency |
| **Adjusted ratio** | Conditional Laplace | Yes | Summary statistics, plotting |
| **Net benefit** | `TP - FP` | Yes | Clinical decision-making |

### Edge-case handling

| Scenario | Raw ratio | Adjusted ratio | Net benefit | Display string |
|----------|-----------|----------------|-------------|----------------|
| TP > 0, FP = 0 | `inf` | (TP+1)/1 | +TP | **inf (pure gain)** |
| TP = 0, FP > 0 | 0 | 1/(FP+1) | −FP | **0 (cost only)** |
| TP = 0, FP = 0 | `nan` | NaN | 0 | **- (no change)** |
| TP < 0, FP = 0 | `-inf` | (TP+1)/1 | TP | **-inf (pure loss)** |
| TP > 0, FP > 0 | TP/FP | (TP+1)/(FP+1) | TP−FP | Numeric value |

### Clinical interpretation

**Net benefit is preferred for clinical decisions** because:
1. It is always defined (no inf/nan)
2. It is signed (positive = good, negative = bad)
3. It is additive across regions
4. It directly answers: "How many more true variants do I get per genome?"

**Adjusted ratio is preferred for statistical reporting** because:
1. It prevents aggregation failures
2. It preserves ranking (higher = better)
3. It is familiar to genomics researchers

**Raw ratio is preserved for transparency** because:
1. It shows exact computational steps
2. It allows independent verification
3. It reveals when edge cases occur

---

## Error Handling & Validation

### Exception hierarchy

```
BurdenCalculatorError (base)
├── FileAccessError        # Files missing, unreadable, empty
├── DataValidationError    # Wrong columns, bad types, bad values
├── ManifestError          # Wrong format, missing values, duplicates
└── AggregationError       # Empty groups, no successful comparisons
```

### Validation stages

| Stage | Checks | Failure action |
|-------|--------|---------------|
| File access | Exists, is file, readable, non-empty | `FileAccessError` |
| Manifest structure | Required columns, no nulls | `ManifestError` |
| hap.py structure | Required columns, no null merge keys | `DataValidationError` |
| Data types | Numeric columns are numeric | `DataValidationError` |
| Value sanity | No negative FP/TP, recall/precision in [0,1] | Warning or `DataValidationError` |
| Merge compatibility | Matching stratification keys | Warning for unmatched rows |
| PP drift | `gain_pp == (cons_raw - base_raw) * 100` | Warning if > 1e-10 |
| Region-size drift | `Subset.IS_CONF.Size` within 5% | Warning if exceeded |
| Aggregation | Non-empty groups after filtering | `AggregationError` |

### Exit codes

| Code | Meaning | Typical cause |
|------|---------|---------------|
| 0 | Success | — |
| 1 | No successful comparisons | All manifest rows failed |
| 2 | File access error | Missing files, permission denied |
| 3 | Manifest format error | Wrong columns, missing values |
| 4 | Data validation error | Wrong hap.py format, bad values |

---

## Aggregation Engine

### Grouping logic

Default groups: `Caller`, `Condition`, `Type`, `Subset`, `Filter`, `Genotype`

Customizable via `--group-by`.

### Statistics computed

For each metric in each group:

```python
grouped = df.groupby(group_cols).agg({
    "Recall_gain_pp":    ["median", "min", "max", "count"],
    "Additional_TP":     ["median", "min", "max", "count"],
    "Additional_FP":     ["median", "min", "max", "count"],
    "TP_FP_ratio_adj":   ["median", "min", "max", "count"],
    "Net_benefit":       ["median", "min", "max", "count"],
    # ... etc
})
```

### Formatted output

```python
def format_median_range(median, min_val, max_val):
    if abs(median) >= 100:   fmt = ".0f"
    elif abs(median) >= 10:  fmt = ".1f"
    elif abs(median) >= 1:   fmt = ".2f"
    else:                    fmt = ".3f"
    return f"{median:{fmt}} [{min_val:{fmt}}-{max_val:{fmt}}]"
```

Example: `4.75 [2.27-7.11]`

### Edge cases handled

| Scenario | Handling |
|----------|----------|
| Single sample in group | min = max = median (range collapses to point) |
| All NaN in group | Formatted string is empty |
| Empty group after filtering | `AggregationError` with descriptive message |
| `inf` values in raw ratio | Preserved in raw; adjusted ratio used for stats |
| Negative net benefit | Formatted with sign (e.g., `-49 [-172-+74]`) |

### Overall aggregation

For overall summaries (summing across stratification rows), counts are summed within each sample first, then recall/precision/F1 are **recomputed from aggregated counts** rather than averaged:

```python
Recall_overall  = sum(QUERY.TP) / sum(TRUTH.TOTAL)
Precision_overall = sum(QUERY.TP) / sum(QUERY.TP + QUERY.FP)
```

This is mathematically required because truth-set sizes differ across stratifications.

### Median-of-medians

As a sensitivity analysis, the framework computes the median across all stratification rows *within* each sample, then reports the median [min–max] of those per-sample medians. This gives equal weight to each sample regardless of stratification row count.

---

## Output Format Specifications

### 1. Raw per-sample CSV

**Filename:** `{prefix}_raw_per_sample.csv`  
**Rows:** One per sample × stratification × filter  
**Columns:** 40+ including all derivation stages

Key column categories:
- **Identifiers:** `Sample`, `Caller`, `Condition`, `Type`, `Subset`, `Filter`
- **Raw metrics (0-1):** `Recall_base_raw`, `Recall_cons_raw`, etc.
- **Percentage metrics (0-100):** `Recall_base_pct`, `Recall_cons_pct`, etc.
- **PP gains:** `Recall_gain_pp`, `Precision_gain_pp`
- **Verification:** `Recall_gain_pp_check`, `Precision_gain_pp_check`
- **Absolute counts:** `Additional_TP`, `Additional_FP`, `Additional_FN`
- **Ratios:** `TP_FP_ratio` (adjusted), `TP_FP_ratio_raw`, `TP_FP_ratio_adj`
- **Net benefit:** `Net_benefit`
- **Normalized:** `*_per_Mbp`, `*_per_100_truth`
- **Metadata:** `Subset.IS_CONF.Size_base`, `in_both`, `Duplicate_Count`

### 2. Summary CSV

**Filename:** `{prefix}_summary_median_range.csv`  
**Rows:** One per aggregation group  
**Columns:** Group keys + `{metric}_{stat}` + `{metric}_formatted`

Statistics per metric:
- `{metric}_median`
- `{metric}_min`
- `{metric}_max`
- `{metric}_count`
- `{metric}_formatted` (e.g., `4.75 [2.27-7.11]`)

### 3. Markdown table

**Filename:** `{prefix}_table.txt`  
**Format:** Plain text with aligned columns  
**Columns:** Renamed for readability

| Internal column | Display name |
|-----------------|--------------|
| `Recall_gain_pp_formatted` | Recall Gain (pp) |
| `Precision_gain_pp_formatted` | Precision Gain (pp) |
| `Additional_TP_formatted` | Additional TP |
| `Additional_FP_formatted` | Additional FP |
| `TP_FP_ratio_formatted` | TP:FP Ratio |
| `Net_benefit_formatted` | Net Benefit |
| `Additional_TP_per_Mbp_formatted` | Additional TP/Mbp |
| `Additional_FP_per_Mbp_formatted` | Additional FP/Mbp |

Includes interpretation notes at the bottom.

### 4. LaTeX table

**Filename:** `{prefix}_table.tex`  
**Format:** Standard `tabular` environment  
**Features:**
- Caption included
- Compatible with `longtable` if needed

### 5. CSV table

**Filename:** `{prefix}_table.csv`  
**Format:** Plain CSV with clean column names  
**Purpose:** Machine-readable, direct import to Excel/R/Python

| Internal column | CSV column name |
|-----------------|-----------------|
| `Recall_gain_pp_formatted` | `Recall_Gain_pp` |
| `Precision_gain_pp_formatted` | `Precision_Gain_pp` |
| `Additional_TP_formatted` | `Additional_TP` |
| `Additional_FP_formatted` | `Additional_FP` |
| `TP_FP_ratio_formatted` | `TP_FP_Ratio` |
| `Net_benefit_formatted` | `Net_Benefit` |
| `Additional_TP_per_Mbp_formatted` | `Additional_TP_per_Mbp` |
| `Additional_FP_per_Mbp_formatted` | `Additional_FP_per_Mbp` |
| `Additional_TP_per_100_truth_formatted` | `Additional_TP_per_100_truth` |
| `Additional_FP_per_100_truth_formatted` | `Additional_FP_per_100_truth` |

**Key design decisions:**
- No spaces in column names (Excel/R/Python friendly)
- No special characters (underscores only)
- Median-range strings preserved as single cell values
- Identical row ordering to Markdown and LaTeX tables
- Edge cases displayed as descriptive strings (e.g., `inf (pure gain)`)

### 6. Reproducibility audit

**Filename:** `{prefix}_reproducibility.txt`  
**Sections:**
1. Percentage-point derivation formula
2. Verification procedure
3. Absolute burden derivation
4. Edge-case handling
5. Normalization formulas
6. Aggregation method
7. Data provenance (samples, callers, conditions, row counts)
8. Validation checks performed

---

## Testing & Validation

### Unit test checklist

```python
# Test 1: PP derivation correctness
def test_pp_derivation():
    base = pd.DataFrame({"METRIC.Recall": [0.985]})
    cons = pd.DataFrame({"METRIC.Recall": [0.988]})
    result = compute_burden(base, cons, "S", "C", "T")
    assert result["Recall_gain_pp"].iloc[0] == 0.3
    assert result["Recall_gain_pp_check"].iloc[0] == 0.3

# Test 2: TP:FP ratio edge cases
def test_tp_fp_ratio():
    assert calc_tp_fp_ratio_raw(10, 5) == 2.0      # Normal
    assert calc_tp_fp_ratio_raw(10, 0) == np.inf    # Pure gain
    assert calc_tp_fp_ratio_raw(0, 0) == np.nan     # No change
    assert calc_tp_fp_ratio_raw(-10, 0) == -np.inf  # Pure loss

    # Adjusted ratio handles all cases
    assert (10 + 1) / (0 + 1) == 11.0             # Adjusted pure gain
    assert (0 + 1) / (5 + 1) == 1/6               # Adjusted cost only

    # Net benefit
    assert 10 - 5 == 5                             # Net positive
    assert 0 - 5 == -5                             # Net negative
    assert 0 - 0 == 0                              # No change

# Test 3: Empty group handling
def test_empty_group():
    df = pd.DataFrame({"A": [], "Metric": []})
    with pytest.raises(AggregationError):
        aggregate_summary(df, group_cols=["A"])

# Test 4: File validation
def test_missing_file():
    with pytest.raises(FileAccessError):
        validate_file_exists("/nonexistent/file.csv", "test")

# Test 5: CSV table output
def test_csv_table():
    summary = pd.DataFrame({
        "Caller": ["bcftools"],
        "Recall_gain_pp_formatted": ["0.30 [0.25-0.35]"],
        "Net_benefit_formatted": ["-49 [-172-74]"]
    })
    csv_df = to_csv_table(summary)
    assert "Recall_Gain_pp" in csv_df.columns
    assert "Net_Benefit" in csv_df.columns
    assert "Recall_gain_pp_formatted" not in csv_df.columns

# Test 6: Edge-case formatting
def test_ratio_formatting():
    assert format_ratio(100, 0, 101.0) == "inf (pure gain)"
    assert format_ratio(0, 50, 1/51) == "0 (cost only)"
    assert format_ratio(0, 0, np.nan) == "- (no change)"
    assert format_ratio(50, 100, 51/101) == "0.50"
```

### Integration test

```bash
# Generate synthetic data
python -c "
import pandas as pd
import numpy as np
np.random.seed(42)

for sample in ['HG001', 'HG002']:
    for caller in ['bcftools']:
        base = pd.DataFrame({
            'Type': ['SNV', 'INDEL'],
            'Subtype': ['*', '*'],
            'Subset': ['AllAutosomes', 'AllAutosomes'],
            'Filter': ['ALL', 'ALL'],
            'Genotype': ['*', '*'],
            'METRIC.Recall': [0.985, 0.680],
            'METRIC.Precision': [0.995, 0.750],
            'QUERY.FP': [2000, 8000],
            'QUERY.TP': [395000, 54400],
            'TRUTH.TOTAL': [400000, 80000],
            'Subset.IS_CONF.Size': [2500000000, 2500000000]
        })
        base.to_csv(f'{sample}_baseline_{caller}.csv', index=False)

        cons = base.copy()
        cons['METRIC.Recall'] += [0.003, 0.050]
        cons['METRIC.Precision'] -= [0.002, 0.030]
        cons['QUERY.TP'] += [1200, 4000]
        cons['QUERY.FP'] += [1500, 3500]
        cons.to_csv(f'{sample}_consensus_{caller}.csv', index=False)
"

# Create manifest
cat > test_manifest.csv << EOF
sample,caller,condition,baseline_path,consensus_path
HG001,bcftools,AF10,HG001_baseline_bcftools.csv,HG001_consensus_bcftools.csv
HG002,bcftools,AF10,HG002_baseline_bcftools.csv,HG002_consensus_bcftools.csv
EOF

# Run
python scripts/calculate_fp_burden.py \
    --manifest test_manifest.csv \
    --output-prefix test_output/ \
    --verbose

# Verify outputs exist
ls test_output/
```

---

## Performance Considerations

### Memory usage

- Raw DataFrames held in memory: ~2x number of comparisons x stratification rows
- Typical: 18 comparisons x 50 regions = ~900 rows -> < 1 MB
- Even for 100+ samples: well under 100 MB

### Speed

| Operation | Typical time | Bottleneck |
|-----------|-------------|------------|
| File loading | < 100 ms/file | Disk I/O |
| Merge + derivation | < 10 ms/pair | pandas merge |
| Aggregation | < 50 ms | GroupBy median |
| Markdown/LaTeX formatting | < 100 ms | tabulate (if installed) |
| CSV table formatting | < 50 ms | DataFrame rename + to_csv |
| **Total (18 comparisons)** | **< 5 seconds** | — |

### Scalability

- Tested up to 50 samples x 3 callers x 5 conditions = 750 comparisons
- Memory: ~50 MB
- Time: ~30 seconds
- For larger cohorts, consider chunking the manifest or using Dask

---

## API Reference

### Public functions

#### `load_happy(path: str) -> pd.DataFrame`

Load and validate a hap.py summary file.

**Args:**
- `path`: File path to hap.py CSV/TSV

**Returns:**
- Validated DataFrame with cleaned column names

**Raises:**
- `FileAccessError`: File missing or unreadable
- `DataValidationError`: Missing columns, bad types, bad values

---

#### `load_manifest(path: str) -> pd.DataFrame`

Load and validate a batch manifest.

**Args:**
- `path`: File path to manifest CSV/TSV

**Returns:**
- Validated DataFrame with normalized column names

**Raises:**
- `FileAccessError`: File missing or unreadable
- `ManifestError`: Wrong format, missing values, duplicates

---

#### `compute_burden(df_base, df_cons, sample, caller, condition) -> Tuple[pd.DataFrame, List[pd.DataFrame]]`

Compute all burden metrics for a single baseline-consensus pair.

**Args:**
- `df_base`: Baseline hap.py DataFrame
- `df_cons`: Consensus hap.py DataFrame
- `sample`: Sample identifier string
- `caller`: Caller name string
- `condition`: Condition label string

**Returns:**
- `(result_df, mismatch_records)`: Derived metrics and any mismatch reports

**Raises:**
- `DataValidationError`: Merge keys missing or data incompatible

**Key derived columns:**
- `Recall_gain_pp`, `Precision_gain_pp`
- `Additional_TP`, `Additional_FP`, `Net_benefit`
- `TP_FP_ratio` (adjusted), `TP_FP_ratio_raw`, `TP_FP_ratio_adj`
- `Additional_TP_per_Mbp`, `Additional_FP_per_Mbp`

---

#### `aggregate_summary(df_raw, group_cols=None) -> pd.DataFrame`

Aggregate per-sample results across samples.

**Args:**
- `df_raw`: Output from `compute_burden()` for multiple samples
- `group_cols`: List of columns to group by (default: `DEFAULT_GROUP_COLS`)

**Returns:**
- DataFrame with `median`, `min`, `max`, `count`, and `_formatted` columns

**Raises:**
- `AggregationError`: Empty groups or no valid metrics

---

#### `compute_overall_by_caller_region(df_raw, levels=None) -> pd.DataFrame`

Compute overall metrics at configurable aggregation levels.

**Args:**
- `df_raw`: Output from `compute_burden()` for multiple samples
- `levels`: List of level names (e.g., `["caller-region", "caller"]`)

**Returns:**
- DataFrame with overall metrics recomputed from summed counts

---

#### `compute_median_of_medians(df_raw, levels=None) -> pd.DataFrame`

Compute median-of-medians across samples at configurable levels.

**Args:**
- `df_raw`: Output from `compute_burden()` for multiple samples
- `levels`: List of level names (e.g., `["caller-condition", "caller"]`)

**Returns:**
- DataFrame with median-of-medians sensitivity analysis

---

#### `to_csv_table(df_summary, region_order=None, caller_order=None) -> pd.DataFrame`

Convert summary to a clean CSV-ready DataFrame with machine-friendly column names.

**Args:**
- `df_summary`: Output from `aggregate_summary()`
- `region_order`: Optional list to enforce region ordering
- `caller_order`: Optional list to enforce caller ordering

**Returns:**
- DataFrame with clean column names (no spaces, underscores only)

**Column mapping:**

| Input column | Output column |
|--------------|---------------|
| `Recall_gain_pp_formatted` | `Recall_Gain_pp` |
| `Precision_gain_pp_formatted` | `Precision_Gain_pp` |
| `Additional_TP_formatted` | `Additional_TP` |
| `Additional_FP_formatted` | `Additional_FP` |
| `TP_FP_ratio_formatted` | `TP_FP_Ratio` |
| `Net_benefit_formatted` | `Net_Benefit` |
| `Additional_TP_per_Mbp_formatted` | `Additional_TP_per_Mbp` |
| `Additional_FP_per_Mbp_formatted` | `Additional_FP_per_Mbp` |
| `Additional_TP_per_100_truth_formatted` | `Additional_TP_per_100_truth` |
| `Additional_FP_per_100_truth_formatted` | `Additional_FP_per_100_truth` |

---

#### `process_manifest(manifest_path, output_prefix, regions=None, filters=None, group_cols=None, dry_run=False, overall_levels=None, mom_levels=None, no_overall=False, no_mom=False, report_mismatches=False) -> Tuple[pd.DataFrame, pd.DataFrame]`

End-to-end batch processing.

**Args:**
- `manifest_path`: Path to manifest file
- `output_prefix`: Prefix for all output files
- `regions`: Optional list of Subset names to filter
- `filters`: Optional list of Filter values to filter
- `group_cols`: Optional custom grouping columns
- `dry_run`: If True, validate only; do not write outputs
- `overall_levels`: List of overall aggregation level names
- `mom_levels`: List of median-of-medians level names
- `no_overall`: Skip overall aggregation
- `no_mom`: Skip median-of-medians calculation
- `report_mismatches`: Write detailed mismatch report

**Returns:**
- `(df_raw, df_summary)`: Raw and aggregated DataFrames

**Side effects:**
Writes 6+ files to disk (see Output Files section).

**Raises:**
- Any `BurdenCalculatorError` subclass

---

### Internal functions

#### `validate_file_exists(path, label) -> Path`

Check file existence, readability, and non-emptiness.

#### `validate_happy_dataframe(df, source)`

Comprehensive validation of hap.py DataFrame structure and content.

#### `validate_manifest(manifest_df, path)`

Comprehensive validation of manifest DataFrame.

#### `format_median_range(median, min_val, max_val) -> str`

Smart formatting with adaptive precision.

#### `to_markdown_table(df_summary, region_order=None, caller_order=None) -> str`

Convert summary to plain-text table.

#### `to_latex_table(df_summary, region_order=None, caption="") -> str`

Convert summary to LaTeX table.

#### `write_reproducibility_doc(output_prefix, df_raw, df_summary, manifest_path)`

Generate reproducibility audit document.

---

## References

1. Krusche P, et al. **Haplotype comparison tools / hap.py.** 2018.
2. Dwarshuis N, et al. **The GIAB genomic stratifications resource for human reference genomes.** *Nat Commun* 2024.
3. Saidin A, Ricos MG, Dibbens LM. **IUPAC consensus references improve variant detection in clinically challenging genomic regions.** *Cell Genomics*. 2026.
