# IUPAC Consensus Reference for Variant Calling

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Reference-bias mitigation using IUPAC degenerate codes and ambiguity-aware alignment (novoAlign).
This repository contains the complete workflow for constructing consensus genomes, benchmarking variant calls, and quantifying absolute TP/FP burden trade-offs.

---

## Quick Links

| I want to... | Go to... |
|--------------|----------|
| Build an IUPAC consensus reference | [`docs/consensus_generation.md`](docs/consensus_generation.md) |
| Run the burden calculation framework | [`docs/usage_manual.md`](docs/usage_manual.md) |
| Understand the math / formulas | [`docs/supplemental_methods.md`](docs/supplemental_methods.md) |
| Read the full technical spec | [`docs/technical_manual.md`](docs/technical_manual.md) |
| Reproduce the paper figures & tables | [`paper/`](paper/) |
| Browse the API / run tests | [`docs/technical_manual.md`](docs/technical_manual.md) |

---

## What is this?

Clinical short-read WGS pipelines struggle with reference bias in low-mappability regions, segmental duplications, the MHC, and challenging medically relevant genes (CMRG). This project evaluates a pragmatic mitigation: substituting common alleles into GRCh38 using IUPAC degenerate codes, then aligning with an ambiguity-aware aligner (novoAlign). The burden framework (`calculate_fp_burden.py`) quantifies absolute true-positive and false-positive trade-offs relative to a BWA-MEM/GRCh38 baseline.

**Key finding:** SNV recall increased by 3.1–3.9 pp in low-mappability regions and 1.8–3.1 pp in segmental duplications; INDEL recall rose by 4.5–5.8 pp and 2.4–3.8 pp, respectively.

> **Note:** These are descriptive, hypothesis-generating findings from three GIAB samples. They warrant validation in larger, ancestrally diverse cohorts before clinical deployment.

---

## Repository Structure

```
IUPAC-consensus-reference/
├── docs/                    # Documentation (methods, usage, technical)
├── scripts/                 # Production code
│   └── calculate_fp_burden.py
├── examples/                # Manifests & toy data for quick start
├── paper/                   # Manuscript figure & table scripts
├── tests/                   # pytest suite
└── .github/                 # Issue templates & CI
```

---

## Installation

```bash
git clone https://github.com/akzam/IUPAC-consensus-reference.git
cd IUPAC-consensus-reference
pip install pandas numpy matplotlib seaborn
```

---

## Quick Start

### 1. Single comparison

```bash
python scripts/calculate_fp_burden.py \
    --baseline examples/happy_outputs/HG001_baseline_bcftools.csv \
    --consensus examples/happy_outputs/HG001_consensus_bcftools.csv \
    --sample HG001 --caller bcftools --condition 1KG_AF10pc \
    --output results.csv
```

### 2. Batch processing (recommended)

```bash
python scripts/calculate_fp_burden.py \
    --manifest examples/manifests/samples_manifest.tsv \
    --output-prefix my_run/ \
    --regions AllAutosomes lowmappabilityall segdups CMRG MHC \
    --filters ALL
```

See [`docs/usage_manual.md`](docs/usage_manual.md) for full options.

---

## Citation

If you use this framework, please cite:

> Saidin A, Ricos MG, Dibbens LM. *IUPAC consensus references improve variant detection in clinically challenging genomic regions.* Cell Genomics. 2026.

```bibtex
@article{saidin2026iupac,
  title={IUPAC consensus references improve variant detection in clinically challenging genomic regions},
  author={Saidin, Akzam and Ricos, Michael G. and Dibbens, Leanne M.},
  journal={Cell Genomics},
  year={2026},
  publisher={Elsevier}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) file.

## Contact

For issues or feature requests, please open a [GitHub Issue](https://github.com/akzam/IUPAC-consensus-reference/issues).

For questions about data access or reproduction, contact the lead author: **Leanne M. Dibbens** (leanne.dibbens@adelaide.edu.au).
