# -*- coding: utf-8 -*-
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

os.getcwd()

# ex 1
penguins = sns.load_dataset("penguins")
penguins_extra = pd.read_csv("penguins_extra.csv", index_col="Unnamed: 0")
penguins.head()
penguins_extra.head()

penguins.index = penguins.index.astype(int)

penguins_extra = penguins_extra.T
penguins_extra.index = penguins_extra.index.astype(int)

penguins_merge = pd.merge(penguins, penguins_extra, how="outer", left_index=True, right_index=True)
print(penguins_merge)

grouped_by_bill_length = penguins_merge.groupby(by=["island"])["bill_length_mm"].median()
grouped_by_bill_length.head()

plt.bar(penguins_merge["island"], height=penguins_merge["bill_length_mm"])

plt.xlabel("Island")
plt.ylabel("Bill Length (mm)")
plt.title("Median Bill Length by Island")

plt.show()

penguins_merge["IQ"] = penguins_merge["IQ"].astype(float)
group_by_iq_sex_preg = penguins_merge.groupby(by=["sex", "is_pregnant"])["IQ"].mean()

print(group_by_iq_sex_preg)

# ex 2
