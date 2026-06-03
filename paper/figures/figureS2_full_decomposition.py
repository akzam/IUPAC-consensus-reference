#!/usr/bin/env python3
"""
Figure S2. Stacked bar chart showing total recall gain decomposed into
aligner and IUPAC contribution for all stratifications.
Same as Figure 2 but including all stratifications (AS, D, LM, SD, HP, TR, MHC, CMRG)
and all database/AF conditions.

Usage:
    python paper/figures/figureS2_full_decomposition.py --input results_table.csv --output figureS2.pdf
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to burden output CSV")
    parser.add_argument("--output", default="figureS2_full_decomposition.pdf")
    args = parser.parse_args()

    data = pd.read_csv(args.input)

    # All stratifications
    regions = ['AS', 'D', 'LM', 'SD', 'HP', 'TR', 'MHC', 'CMRG']
    callers = ['BCFtools', 'FreeBayes', 'GATK']
    variant_types = ['SNP', 'INDEL']

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

    for v_idx, vtype in enumerate(variant_types):
        ax = axes[v_idx]
        x = np.arange(len(regions))
        width = 0.25

        for c_idx, caller in enumerate(callers):
            aligner_gains = []
            iupac_gains = []

            for region in regions:
                base_rows = data[(data['Caller'] == caller) &
                                (data['Condition'] == 'novoAlign GRCh38 (base)') &
                                (data['Type'] == vtype) &
                                (data['Subset'] == region)]
                consensus_rows = data[(data['Caller'] == caller) &
                                     (data['Condition'] == 'IUPAC median') &
                                     (data['Type'] == vtype) &
                                     (data['Subset'] == region)]

                if len(base_rows) > 0 and len(consensus_rows) > 0:
                    base_val = base_rows['Recall_Gain_pp'].values[0]
                    consensus_val = consensus_rows['Recall_Gain_pp'].values[0]
                    aligner_gains.append(base_val)
                    iupac_gains.append(consensus_val - base_val)
                else:
                    aligner_gains.append(0)
                    iupac_gains.append(0)

            offset = (c_idx - 1) * width
            ax.bar(x + offset, aligner_gains, width, label=f'{caller} (aligner)',
                   alpha=0.7, color=colors[c_idx])
            ax.bar(x + offset, iupac_gains, width, bottom=aligner_gains,
                   label=f'{caller} (IUPAC)', alpha=0.4, color=colors[c_idx])

        ax.set_ylabel('Recall Gain (pp)')
        ax.set_title(f'{vtype} Recall Gain Decomposition (All Stratifications)')
        ax.set_xticks(x)
        ax.set_xticklabels(regions, rotation=45, ha='right')
        ax.axhline(y=0, color='black', linewidth=0.5)
        ax.legend(loc='upper left', fontsize=7)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.savefig(args.output.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
