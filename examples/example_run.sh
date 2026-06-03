#!/bin/bash
# Example run script for calculate_fp_burden.py
# This demonstrates batch processing with the example data.

set -e

echo "=== IUPAC Consensus Burden Framework - Example Run ==="
echo ""

# Check Python version
python3 --version

# Install dependencies (if not already installed)
# pip install pandas numpy

# Run batch processing with example manifest
echo "Running batch processing..."
python3 scripts/calculate_fp_burden.py \
    --manifest examples/manifests/samples_manifest.tsv \
    --output-prefix example_output/ \
    --regions AllAutosomes lowmappabilityall segdups CMRG MHC \
    --filters ALL \
    --verbose

echo ""
echo "=== Output files ==="
ls -lh example_output/

echo ""
echo "=== Preview of CSV table ==="
head -5 example_output/example_output_table.csv

echo ""
echo "=== Done ==="
