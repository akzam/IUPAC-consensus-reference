#!/usr/bin/env python3
"""
Figure 3. Heatmap of recall gain by region, caller, and variant type.
Heatmap with rows = regions, columns = callers, panels = SNV/INDEL.
Colour scale from blue (loss) to red (gain), white at 0.
Asterisk denotes directional consistency.

Usage:
    python paper/figures/figure3_heatmap.py --input results_table.csv --output figure3.pdf
"""

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to burden output CSV")
    parser.add_argument("--output", default="figure3_heatmap.pdf")
    args = parser.parse_args()

    data = pd.read_csv(args.input)
    consensus_data = data[data['Condition'] == 'IUPAC median']

    regions = ['LM', 'SD', 'MHC', 'CMRG', 'AS']
    callers = ['BCFtools', 'FreeBayes', 'GATK']

    fig, axes = plt.subplots(1, 2, figsize=(10, 8))

    for v_idx, vtype in enumerate(['SNP', 'INDEL']):
        ax = axes[v_idx]
        matrix = np.zeros((len(regions), len(callers)))
        consistency = np.zeros((len(regions), len(callers)), dtype=bool)

        for r_idx, region in enumerate(regions):
            for c_idx, caller in enumerate(callers):
                subset = consensus_data[(consensus_data['Type'] == vtype) &
                                       (consensus_data['Subset'] == region) &
                                       (consensus_data['Caller'] == caller)]
                if len(subset) > 0:
                    val = subset['Recall_Gain_pp'].values[0]
                    matrix[r_idx, c_idx] = val
                    # Directional consistency: all three samples same sign
                    per_sample = subset[['HG001_Recall_Gain', 'HG002_Recall_Gain',
                                        'HG005_Recall_Gain']].values[0]
                    consistency[r_idx, c_idx] = (np.all(per_sample > 0) or
                                                 np.all(per_sample < 0))

        sns.heatmap(matrix, annot=True, fmt='.2f', cmap='RdBu_r',
                   center=0, vmin=-1, vmax=6, ax=ax,
                   xticklabels=callers, yticklabels=regions,
                   linewidths=0.5, linecolor='white')

        for r_idx in range(len(regions)):
            for c_idx in range(len(callers)):
                if consistency[r_idx, c_idx]:
                    ax.text(c_idx + 0.85, r_idx + 0.15, '*',
                           fontsize=16, color='black', weight='bold')

        ax.set_title(f'{vtype} Recall Gain (pp)')

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.savefig(args.output.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
