#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  1 10:46:37 2026

@author: ornitb
"""

import pandas as pd
import numpy as np
import os
import scanpy as sc

print(os.getcwd())

sig_activity = pd.read_csv("Documents/Technion/research_and_studies/courses/ביואינפורמטיקה של סרטן/cancer_bioinformatics/tutorial_10/sig_activity.csv")
tcga_tumor_abbr = pd.read_csv("Documents/Technion/research_and_studies/courses/ביואינפורמטיקה של סרטן/cancer_bioinformatics/tutorial_10/TCGA_Study_Abbreviations.csv")

desired_cosmic_columns = sig_activity.columns.tolist()

desired_cosmic_columns = [
    item for item in desired_cosmic_columns
    if item.startswith("COSMIC") and item.endswith(".norm")
]

desired_cosmic_columns

avg_norm_activity_by_ttype = sig_activity.groupby("tterm")[desired_cosmic_columns].mean()

avg_norm_activity_by_ttype["top_signature"] = np.nan

for ttype in avg_norm_activity_by_ttype.index:
    avg_norm_activity_by_ttype.loc[]