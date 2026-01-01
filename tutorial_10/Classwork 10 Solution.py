# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 20:47:22 2025

@author: stavnaky
"""

#Classwork 10 Solution

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import os

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)   # or a large int
pd.set_option('display.expand_frame_repr', True)

#Q1 - Read the files
os.chdir(r"C:\Users\stavnaky\Desktop\Cancer Bioinformatics 2025_6\Tutorial_10_2025_2026 Mutational Signatures")
abb = pd.read_csv("TCGA_Study_Abbreviations.csv")
sig_activity = pd.read_csv("sig_activity.csv",index_col="Sample")



##Q2 - For each ttype, find the top COSMIC signature by average normalized activity (COSMICX.norm).
#Q2.1 - Relevant columns
#Using list comprehension:
cols=["COSMIC"+str(x)+".norm" for x in range(1,31)]
#Using a loop:
cols = []
for x in range(1,31):
    cols.append("COSMIC"+str(x)+".norm")

#Q2.2 - Group by ttype, for each relevant column compute the mean
avg_norm_activity_by_ttype = sig_activity.groupby("ttype")[cols].mean()
#Q2.3 - Make a new column and fill with nan
avg_norm_activity_by_ttype["top_signature"] = np.nan
#Q2.4 - Fill that column with the top signature of each tumor type
for ttype in avg_norm_activity_by_ttype.index:
    avg_norm_activity_by_ttype.loc[avg_norm_activity_by_ttype.index==ttype,"top_signature"] = \
        avg_norm_activity_by_ttype.loc[avg_norm_activity_by_ttype.index==ttype,cols].idxmax(axis=1)

#Q3 - SBS1 and SBS5. They're clock-like signatures. It means that it comes with age.
#It doesn't give us much information.

#Q4 - Remove these two signatures and reanalyze
cols = list(set(cols).difference({"COSMIC1.norm","COSMIC5.norm"}))

avg_norm_activity_by_ttype = sig_activity.groupby("ttype")[cols].mean()
avg_norm_activity_by_ttype["top_signature"] = np.nan
for ttype in avg_norm_activity_by_ttype.index:
    avg_norm_activity_by_ttype.loc[avg_norm_activity_by_ttype.index==ttype,"top_signature"] = \
        avg_norm_activity_by_ttype.loc[avg_norm_activity_by_ttype.index==ttype,cols].idxmax(axis=1)
avg_norm_activity_by_ttype["top_signature"]
