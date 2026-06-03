# Documentation Index

This directory contains all human-readable documentation for the IUPAC consensus reference project.

| Document | Audience | Description |
|----------|----------|-------------|
| [`consensus_generation.md`](consensus_generation.md) | Bioinformaticians | Step-by-step guide to building the IUPAC consensus reference with novoUtil. |
| [`usage_manual.md`](usage_manual.md) | End users | How to run `calculate_fp_burden.py`, interpret outputs, and troubleshoot. |
| [`technical_manual.md`](technical_manual.md) | Developers / Reviewers | Architecture, data flow, API reference, error codes, and testing. |
| [`supplemental_methods.md`](supplemental_methods.md) | Readers / Reviewers | Complete mathematical derivation of all burden framework formulas. |
| [`CHANGELOG.md`](CHANGELOG.md) | All users | Version history for `calculate_fp_burden.py`. |

---

## Quick Reference

### Burden framework formulas

See [`supplemental_methods.md`](supplemental_methods.md) for the full derivation. Key equations:

- **Recall gain (pp):** `(METRIC.Recall_consensus - METRIC.Recall_baseline) * 100`
- **Net benefit:** `Additional_TP - Additional_FP`
- **Per-Mbp normalisation:** `Additional_TP / (Subset.IS_CONF.Size_baseline / 1,000,000)`

### Exit codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | No successful comparisons |
| 2 | File access error |
| 3 | Manifest format error |
| 4 | Data validation error |

### Required hap.py columns

`Type, Subtype, Subset, Filter, Genotype, METRIC.Recall, METRIC.Precision, QUERY.FP, QUERY.TP, TRUTH.TOTAL, Subset.IS_CONF.Size`
