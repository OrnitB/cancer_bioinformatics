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

mrna_seq = pd.read_csv('homework_1/data_mrna_seq_v2_rsem.txt', delimiter='\t')
log2_cna = pd.read_csv('homework_1/data_log2_cna.txt', delimiter='\t')
clinimal_patient = pd.read_csv('homework_1/data_mrna_seq_v2_rsem.txt', delimiter='\t')

mrna_seq.shape
mrna_seq.head()
mrna_seq.iloc[0]
mrna_seq.columns