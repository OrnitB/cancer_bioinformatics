# -*- coding: utf-8 -*-
import os
import seaborn as sns
import pandas as pd

os.getcwd()


penguins = sns.load_dataset("penguins")
penguins_extra = pd.read_csv("penguins_extra.csv", index_col="Unnamed: 0")
penguins.head()
penguins_extra.head()

penguins.index = penguins.index.astype(int)

penguins_extra = penguins_extra.T
penguins_extra.index = penguins_extra.index.astype(int)

penguins_merge = pd.merge(penguins, penguins_extra, how="outer", left_index=True, right_index=True)
print(penguins_merge)
