#!/usr/bin/env python3
"""
Figure 2. Decomposition of recall gains by aligner versus IUPAC contribution.
Stacked bar chart showing total recall gain decomposed into aligner contribution
and IUPAC contribution, by region and variant type.

Usage:
    python paper/figures/figure2_decomposition.py --input results_table.csv --output figure2.pdf
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to burden output CSV")
    parser.add_argument("--output", default="figure2_decomposition.pdf")
    args = parser.parse_args()

    data = pd.read_csv(args.input)

    regions = ['LM', 'SD', 'MHC', 'CMRG']
    callers = ['BCFtools', 'FreeBayes', 'GATK']
    variant_types = ['SNP', 'INDEL']

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

    for v_idx, vtype in enumerate(variant_types):
        ax = axes[v_idx]
        x = np.arange(len(regions))
        width = 0.25

        for c_idx, caller in enumerate(callers):
            aligner_gains = []
            iupac_gains = []

            for region in regions:
                base_val = data[(data['Caller'] == caller) &
                               (data['Condition'] == 'novoAlign GRCh38 (base)') &
                               (data['Type'] == vtype) &
                               (data['Subset'] == region)]['Recall_Gain_pp'].values[0]

                consensus_val = data[(data['Caller'] == caller) &
                                    (data['Condition'] == 'IUPAC median') &
                                    (data['Type'] == vtype) &
                                    (data['Subset'] == region)]['Recall_Gain_pp'].values[0]

                aligner_gains.append(base_val)
                iupac_gains.append(consensus_val - base_val)

            offset = (c_idx - 1) * width
            ax.bar(x + offset, aligner_gains, width, label=f'{caller} (aligner)',
                   alpha=0.7, color=colors[c_idx])
            ax.bar(x + offset, iupac_gains, width, bottom=aligner_gains,
                   label=f'{caller} (IUPAC)', alpha=0.4, color=colors[c_idx])

        ax.set_ylabel('Recall Gain (pp)')
        ax.set_title(f'{vtype} Recall Gain Decomposition')
        ax.set_xticks(x)
        ax.set_xticklabels(regions)
        ax.axhline(y=0, color='black', linewidth=0.5)
        ax.legend(loc='upper left', fontsize=8)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.savefig(args.output.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
