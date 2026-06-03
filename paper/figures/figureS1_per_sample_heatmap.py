#!/usr/bin/env python3
"""
Figure S1. Heatmap of recall gain by sample, region, caller, and variant type.
Extended version of Figure 3 showing per-sample recall gains as faceted grid.
Three panels (one per sample), each showing region x caller heatmap.

Usage:
    python paper/figures/figureS1_per_sample_heatmap.py --input results_table.csv --output figureS1.pdf
"""

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to burden output CSV")
    parser.add_argument("--output", default="figureS1_per_sample_heatmap.pdf")
    args = parser.parse_args()

    data = pd.read_csv(args.input)

    samples = ['HG001', 'HG002', 'HG005']
    regions = ['LM', 'SD', 'HP', 'TR', 'MHC', 'CMRG', 'AS', 'D']
    callers = ['BCFtools', 'FreeBayes', 'GATK']

    fig, axes = plt.subplots(3, 2, figsize=(12, 18))

    for s_idx, sample in enumerate(samples):
        for v_idx, vtype in enumerate(['SNP', 'INDEL']):
            ax = axes[s_idx, v_idx]
            matrix = np.zeros((len(regions), len(callers)))

            for r_idx, region in enumerate(regions):
                for c_idx, caller in enumerate(callers):
                    subset = data[(data['Type'] == vtype) &
                                 (data['Subset'] == region) &
                                 (data['Caller'] == caller) &
                                 (data['Condition'] == 'IUPAC median')]
                    if len(subset) > 0:
                        col_name = f'{sample}_Recall_Gain_pp'
                        matrix[r_idx, c_idx] = subset[col_name].values[0]

            sns.heatmap(matrix, annot=True, fmt='.2f', cmap='RdBu_r',
                     center=0, vmin=-2, vmax=7, ax=ax,
                     xticklabels=callers, yticklabels=regions,
                     linewidths=0.5)
            ax.set_title(f'{sample} - {vtype}')

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.savefig(args.output.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
