#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np


def main():
    #INPUTS- TOUCH ONLY THAT
    #Your working directory
    os.getcwd()

    #Logistic regression result on discrete phenotype
    pheno_pop_logistic = r"hapmap1_POP_GWAS.assoc.logistic"
    #Linear regression result on quantitative phenotype
    pheno_qt_linear = r"hapmap1_QT_GWAS.assoc.linear"
    
    #=====================DON'T TOUCH ANYTHING ELSE=================
    #Unless you really want to :)
    outfile_pop = pheno_pop_logistic.split(".")[0] + ".manhattan_plot.png"
    outfile_qt = pheno_qt_linear.split(".")[0] + ".manhattan_plot.png"
    # Manhattan Plot
    plot_manhattan(pheno_pop_logistic, outfile_pop)
    plot_manhattan(pheno_qt_linear, outfile_qt)
    
def plot_manhattan(infile, outfile):
    # Load assoc/linear file (expecting CHR, BP, P, BETA/OR columns)
    df = pd.read_csv(infile, delim_whitespace=True)
    
    # Determine effect size column
    effect_col = None
    if 'BETA' in df.columns:
        effect_col = 'BETA'
    elif 'OR' in df.columns:
        effect_col = 'OR'
    
    # Ensure needed columns exist
    if not {'CHR', 'BP', 'P'}.issubset(df.columns):
        raise ValueError("Input file must have CHR, BP, and P columns")

    # Drop missing P-values
    df = df.dropna(subset=['P'])

    # Sort by chromosome and base pair
    df = df.sort_values(['CHR', 'BP'])

    # Create cumulative base position column across chromosomes
    chr_sizes = df.groupby('CHR').size()
    chr_offsets = chr_sizes.cumsum() - chr_sizes
    df['SNP_index'] = df.groupby('CHR').cumcount() + 1
    df['pos_cum'] = df.apply(lambda row: row['SNP_index'] + chr_offsets[row['CHR']], axis=1)

    # Plot setup
    sns.set_style("white")
    palette = sns.color_palette("Set2", n_colors=df['CHR'].nunique())

    plt.figure(figsize=(14, 6))
    for i, (chrom, group) in enumerate(df.groupby('CHR')):
        plt.scatter(
            group['pos_cum'], 
            -np.log10(group['P']), 
            c=[palette[i % len(palette)]], 
            s=20, label=f'Chr {chrom}'
        )

    # Significance threshold (Bonferroni)
    threshold = -np.log10(0.05 / 1_000_000)
    plt.axhline(y=threshold, color='red', linestyle='--')

    # Annotate significant SNPs with Beta/OR if available
    if effect_col:
        sig = df[df['P'] < 0.05 / 1_000_000]  # Bonferroni significant
        for _, row in sig.iterrows():
            plt.text(
                row['pos_cum'], 
                -np.log10(row['P']) + 0.2,  # slightly above the point
                f"{row[effect_col]:.2f}", 
                fontsize=10, 
                ha='center',
                rotation=0
            )

    # X-axis ticks at chromosome centers
    chrom_ticks = []
    chrom_labels = []
    for chrom, group in df.groupby('CHR'):
        center = (group['pos_cum'].min() + group['pos_cum'].max()) / 2
        chrom_ticks.append(center)
        chrom_labels.append(str(chrom))

    plt.xticks(chrom_ticks, chrom_labels)
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(P)")
    name = os.path.basename(infile).split(".")[0]
    plt.title("Manhattan Plot - " + name)
    plt.legend(markerscale=2, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()
    plt.savefig(outfile, dpi=300)

if __name__ == "__main__":
    main()
    
    
    
    
    
# commands used in the classwork (tutorial 11):
# 3: 
#   plink --file hapmap1 --make-bed --out hapmap1_binary
#       verify with: ls hapmap1_binary.*
# 4:
#       qt.phe → --linear
#       pop.phe → --logistic
#   plink --bfile hapmap1_binary --pheno qt.phe --linear --out hapmap1_QT_GWAS
#   plink --bfile hapmap1_binary --pheno pop.phe --logistic --out hapmap1_POP_GWAS

