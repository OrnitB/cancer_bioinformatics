# -*- coding: utf-8 -*-
"""
Created on Sat Nov 29 18:35:54 2025

@author: Ornit Bhonkar
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import spearmanr
from statannotations.Annotator import Annotator
import itertools
from matplotlib.backends.backend_pdf import PdfPages

# 1
mrna_seq = pd.read_csv('data_mrna_seq_v2_rsem.txt', delimiter='\t')
log2_cna = pd.read_csv('data_log2_cna.txt', delimiter='\t')
clinical_patient = pd.read_csv('data_clinical_patient.txt', delimiter='\t')

""" sanity checks:
mrna_seq.shape
mrna_seq.head()
mrna_seq.iloc[0]

log2_cna.shape
log2_cna.head()
log2_cna.iloc[0]

clinical_patient.shape
clinical_patient.head()
clinical_patient.iloc[0]
"""

# 2
mrna_seq = mrna_seq.dropna()
log2_cna = log2_cna.dropna()

# 3
sample_cols = mrna_seq.columns[2:]
mrna_seq_filtered = mrna_seq[~(mrna_seq[sample_cols] == 0).all(axis=1)].copy()
mrna_seq_filtered.shape

# 4
log2_cna["Entrez_Gene_Id"] = log2_cna["Entrez_Gene_Id"].astype(int)
mrna_seq_filtered["Hugo_Entrez"] = mrna_seq_filtered["Hugo_Symbol"] + ';' + mrna_seq_filtered["Entrez_Gene_Id"].astype(str)
log2_cna["Hugo_Entrez"] = log2_cna["Hugo_Symbol"] + ';' + log2_cna["Entrez_Gene_Id"].astype(str)

""" sanity checks:
mrna_seq_filtered.head()
log2_cna.head()
"""

# 5
mrna_seq_filtered = mrna_seq_filtered.drop_duplicates(subset="Hugo_Entrez", keep='first')
log2_cna = log2_cna.drop_duplicates(subset="Hugo_Entrez", keep='first')

# 6
mrna_set = set(mrna_seq_filtered["Hugo_Entrez"])
cna_set = set(log2_cna["Hugo_Entrez"])

common_genes = mrna_set & cna_set

mrna_seq_filtered = mrna_seq_filtered[mrna_seq_filtered["Hugo_Entrez"].isin(common_genes)].copy()
cna_set_filtered = log2_cna[log2_cna["Hugo_Entrez"].isin(common_genes)].copy()

mrna_seq_filtered.drop(['Hugo_Symbol', 'Entrez_Gene_Id'], axis=1, inplace=True)
cna_set_filtered.drop(['Hugo_Symbol', 'Entrez_Gene_Id'], axis=1, inplace=True)

# 7
mrna_seq_filtered = mrna_seq_filtered.set_index("Hugo_Entrez")
cna_set_filtered = cna_set_filtered.set_index("Hugo_Entrez")

""" sanity checks:
mrna_seq_filtered.head()
cna_set_filtered.head()
"""

# 8
common_samples = list(set(mrna_seq_filtered.columns) & set(cna_set_filtered.columns))
mrna_seq_filtered = mrna_seq_filtered[common_samples]
cna_set_filtered = cna_set_filtered[common_samples]

""" sanity checks:
mrna_seq_filtered.shape
cna_set_filtered.shape
"""

# 9
cna_transposed = cna_set_filtered.T
cna_transposed.index = cna_transposed.index.str[:12]
cna_transposed.head()

clinical_patient_indexed = clinical_patient.set_index("PATIENT_ID")
common_patients_set = set(cna_transposed.index) & set(clinical_patient_indexed.index)
common_patients_list = list(common_patients_set)

cna_transposed_filtered = cna_transposed.loc[common_patients_list]
clinical_final = clinical_patient_indexed.loc[common_patients_list]

""" sanity checks:
cna_transposed_filtered.head()
clinical_final.head()

cna_transposed_filtered.shape
clinical_final.shape
"""

# 10
cna_transposed_filtered["SUBTYPE"] = clinical_final["SUBTYPE"]

# 11
cna_group_by_subtype = cna_transposed_filtered.groupby("SUBTYPE")

# 12
cna_group_by_subtype = cna_group_by_subtype.median()

# 13
cna_final = cna_group_by_subtype.T

""" sanity checks:
cna_final.head()
cna_final.shape
"""

# 14
# copy number alterations usually affect large segments of a chromosome 
# (amplifying or deleting millions of base pairs at once). therefore neighbors
# on the same chromosome often share the exact same copy number score.

# 15
top_genes_by_subtype = {}
for subtype in cna_final.columns:
    top25 = cna_final[subtype].sort_values(ascending=False).head(25)
    top_genes_by_subtype[subtype] = top25
    
for subtype, genes in top_genes_by_subtype.items():
    print(f'\nTop 25 genes for {subtype}')
    print(genes)
    
# 16
cna_transposed_filtered.head()

def gene_boxplot(gene_name, pdf_pages=None):
    subtypes = sorted(cna_transposed_filtered["SUBTYPE"].dropna().unique())
    plt.figure(figsize=(14,10))
    
    ax = sns.boxplot(data=cna_transposed_filtered, x="SUBTYPE", y=gene_name, order=subtypes)
    pairs = list(itertools.combinations(subtypes, 2))
    
    annotator = Annotator(ax, pairs, data=cna_transposed_filtered, x="SUBTYPE", y=gene_name, order=subtypes)
    annotator.configure(test="Kruskal", text_format="simple", loc="inside", verbose=0)
    annotator.apply_and_annotate()
    
    plt.title(f'CNA Score for {gene_name}')
    
    if pdf_pages:
        pdf_pages.savefig()
        plt.close()
    else:
        plt.show()
    
gene_boxplot("ERBB2;2064")
gene_boxplot("IKZF3;22806")
gene_boxplot("ITPKB;3707")

# 17
basal_high_genes = cna_final[cna_final["BRCA_Basal"] > 0.4].index

mrna_transposed = mrna_seq_filtered.T
mrna_transposed.head()

mrna_transposed.index = mrna_transposed.index.str[:12]
mrna_matched = mrna_transposed.loc[common_patients_list]

basal_patients = clinical_final[clinical_final["SUBTYPE"] == "BRCA_Basal"].index

cna_basal = cna_transposed_filtered.loc[basal_patients]
mrna_basal = mrna_matched.loc[basal_patients]

correlation_results = {}

for gene in basal_high_genes:
    x = cna_basal[gene]
    y = mrna_basal[gene]
    
    corr, p_val = spearmanr(x, y)
    
    correlation_results[gene] = corr

basal_cna_mrna_corr = pd.Series(correlation_results)
basal_cna_mrna_corr = basal_cna_mrna_corr.dropna()
basal_cna_mrna_corr = basal_cna_mrna_corr.sort_values(ascending=False)
print(basal_cna_mrna_corr.sort_index().head())

# 18
top_gene1 = basal_cna_mrna_corr.index[0]
top_gene2 = basal_cna_mrna_corr.index[1]

with PdfPages('HW1_Plots.pdf') as pdf:
    gene_boxplot(top_gene1, pdf_pages=pdf)
    gene_boxplot(top_gene2, pdf_pages=pdf)
    