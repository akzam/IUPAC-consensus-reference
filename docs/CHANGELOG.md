# Changelog

All notable changes to `calculate_fp_burden.py` are documented in this file.

## [5.3.7] – 2026-05-23

### Added
- Explicit percentage-point derivation with raw 0-1 value preservation.
- `_check` columns for PP calculation verification (`Recall_gain_pp_check`, `Precision_gain_pp_check`).
- Drift detection warning (tolerance > 1e-10).
- `Net_benefit` metric (signed TP - FP) for clinical interpretation.
- `TP_FP_ratio_raw` and `TP_FP_ratio_adj` columns with conditional Laplace smoothing.
- Edge-case formatting for ratios ("inf (pure gain)", "0 (cost only)", etc.).
- Overall aggregation by configurable levels (`--overall-levels`).
- Median-of-medians sensitivity analysis (`--mom-levels`).
- Duplicate handling: median-averaging of technical replicates in manifest.
- Comprehensive error hierarchy (`BurdenCalculatorError` subclasses).
- Reproducibility audit document auto-generation (`_reproducibility.txt`).
- Dry-run mode (`--dry-run`).
- Adaptive precision formatting for median [min-max] strings.

### Changed
- `TP_FP_ratio` now uses Laplace-adjusted values by default.
- All output tables include `Net_Benefit` column.
- Per-Mbp normalisation uses baseline `Subset.IS_CONF.Size` consistently.
- Aggregation now recomputes recall/precision from summed counts for overall metrics.

### Fixed
- `inf`/`nan` values in aggregated statistics.
- Handling of zero FP counts in ratio computation.
- Missing-column validation now checks all merge keys.

## [5.1] – 2026-05-23

### Added
- `Net_benefit` metric (signed TP - FP).
- `TP_FP_ratio_raw` and `TP_FP_ratio_adj` columns.
- Edge-case formatting for ratios.

### Changed
- `TP_FP_ratio` defaults to adjusted ratio.

### Fixed
- `inf`/`nan` in aggregated statistics.

## [5.0] – 2026-05-23

### Added
- Explicit PP derivation with raw 0-1 preservation.
- Drift detection and `_check` columns.
- Reproducibility audit auto-generation.
- Comprehensive error hierarchy.
- Dry-run mode.

## [4.0]

### Added
- Comprehensive error handling and logging.
- Manifest validation.
- Progress tracking and exit codes.

## [3.0]

### Added
- TP burden and TP:FP ratio.
- Per-Mbp and per-100-truth normalization.

## [2.0]

### Added
- Batch processing and aggregation.
- Markdown/LaTeX output.

## [1.0]

### Added
- Initial single-file implementation.
