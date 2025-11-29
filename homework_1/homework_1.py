# -*- coding: utf-8 -*-
"""
Created on Sat Nov 29 18:35:54 2025

@author: Ornit Bhonkar
"""


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

#1
mrna_seq = pd.read_csv('data_mrna_seq_v2_rsem.txt', delimiter='\t')
log2_cna = pd.read_csv('data_log2_cna.txt', delimiter='\t')
clinical_patient = pd.read_csv('data_clinical_patient.txt', delimiter='\t')

mrna_seq.shape
mrna_seq.head()
mrna_seq.iloc[0]
mrna_seq.columns

log2_cna.shape
log2_cna.head()
log2_cna.iloc[0]
log2_cna.columns

clinical_patient.shape
clinical_patient.head()
clinical_patient.iloc[0]
clinical_patient.columns

#2
mrna_seq = mrna_seq.dropna()
log2_cna = log2_cna.dropna()

mrna_seq.shape
log2_cna.shape

#3
sample_cols = mrna_seq.columns[2:]
mrna_seq_filtered = mrna_seq[~(mrna_seq[sample_cols] == 0).all(axis=1)]
mrna_seq_filtered.shape

#4
log2_cna["Entrez_Gene_Id"] = log2_cna["Entrez_Gene_Id"].astype(int)
mrna_seq_filtered["Hugo_Entrez"] = mrna_seq_filtered["Hugo_Symbol"] + ';' + mrna_seq_filtered["Entrez_Gene_Id"].astype(str)
log2_cna["Hugo_Entrez"] = log2_cna["Hugo_Symbol"] + ';' + log2_cna["Entrez_Gene_Id"].astype(str)

mrna_seq_filtered.head()
log2_cna.head()

#5
mrna_seq_filtered = mrna_seq_filtered.drop_duplicates(subset="Hugo_Entrez", keep='first')
log2_cna = log2_cna.drop_duplicates(subset="Hugo_Entrez", keep='first')

#6
mrna_set = set(mrna_seq_filtered["Hugo_Entrez"])
cna_set = set(log2_cna["Hugo_Entrez"])

common_genes = mrna_set & cna_set

mrna_seq_filtered = mrna_seq_filtered[mrna_seq_filtered["Hugo_Entrez"].isin(common_genes)]
cna_set_filtered = log2_cna[log2_cna["Hugo_Entrez"].isin(common_genes)]

mrna_seq_filtered.drop(['Hugo_Symbol', 'Entrez_Gene_Id'], axis=1, inplace=True)
cna_set_filtered.drop(['Hugo_Symbol', 'Entrez_Gene_Id'], axis=1, inplace=True)

#7
mrna_seq_filtered = mrna_seq_filtered.set_index("Hugo_Entrez")
cna_set_filtered = cna_set_filtered.set_index("Hugo_Entrez")

mrna_seq_filtered.head()
cna_set_filtered.head()

#8
common_samples = list(set(mrna_seq_filtered.columns) & set(cna_set_filtered.columns))
mrna_seq_filtered = mrna_seq_filtered[common_samples]
cna_set_filtered = cna_set_filtered[common_samples]

mrna_seq_filtered.shape
cna_set_filtered.shape

#9
cna_transposed = cna_set_filtered.T
cna_transposed.index = cna_transposed.index.str[:12]
cna_transposed.head()

clinical_patient_indexed = clinical_patient.set_index("PATIENT_ID")
common_patients_set = set(cna_transposed.index) & set(clinical_patient_indexed.index)
common_patients_list = list(common_patients_set)

cna_transposed_filtered = cna_transposed.loc[common_patients_list]
clinical_final = clinical_patient_indexed.loc[common_patients_list]

cna_transposed_filtered.head()
clinical_final.head()

cna_transposed_filtered.shape
clinical_final.shape
