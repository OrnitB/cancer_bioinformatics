# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 14:50:42 2024

@author: stavn
"""
#Tutorial 7- more data manipulations

import pandas as pd
import seaborn as sns
# import matplotlib.pyplot as plt
# Create DataFrames
df1 = pd.DataFrame({
    'ID': [1, 2, 3],
    'Name': ['Alice', 'Bob', 'Charlie']
})

df2 = pd.DataFrame({
    'ID': [2, 3, 4],
    'Age': [25, 30, 35],
})

# Merge on 'ID' column
result = pd.merge(df1, df2, on='ID')
print(result)

# Merge on 'ID' column for left or right dataframe
result = pd.merge(df1, df2, on='ID', how='left')
print(result)

result = pd.merge(df1, df2, on='ID', how='right')
print(result)
#Outer join
result = pd.merge(df1, df2, on='ID', how='outer')
print(result)
#Inner join
result = pd.merge(df1, df2, on='ID', how='inner')
print(result)

#On different columns
df1 = pd.DataFrame({
    'That_ID': [1, 2, 3],
    'Name': ['Alice', 'Bob', 'Charlie']
})

df2 = pd.DataFrame({
    'Other_ID': [2, 3, 4],
    'Age': [25, 30, 35],
})
#Returns the inner join between That ID and Other ID
result = pd.merge(df1,df2,right_on="Other_ID",left_on="That_ID")
print(result)
#For indices
df1 = df1.set_index("That_ID")
df2 = df2.set_index("Other_ID")
result = pd.merge(df1,df2,right_index = True,left_index = True,how='outer')
print(result)


#Overlapping column names
df1 = pd.DataFrame({
    'ID': [1, 2, 3],
    'Name': ['Alice', 'Bob', 'Charlie'],
    'Score': [85, 92, 88]
})

df2 = pd.DataFrame({
    'ID': [2, 3, 4],
    'Name': ['Bob', 'Charlie', 'David'],
    'Score': [80, 90, 78]
})

result = pd.merge(df1, df2, on='ID', how='inner')
print(result)
#different suffixes
result = pd.merge(df1, df2, on='ID', how='inner', suffixes=('_df1', '_df2'))
print(result)
#indicators
result = pd.merge(df1, df2, on='ID', how='outer', indicator=True)
print(result)

# DataFrames
dept1 = pd.DataFrame({
    'ID': [1, 2, 3],
    'Name': ['Alice', 'Bob', 'Charlie'],
    'Salary': [50000, 60000, 55000]
})

dept2 = pd.DataFrame({
    'ID': [2, 3, 4],
    'Name': ['Bob', 'Charlie', 'David'],
    'Salary': [65000, 58000, 70000]
})

# Merge with suffixes and sort by ID
result = pd.merge(dept1, dept2, on='ID', how='outer', suffixes=('_Dept1', '_Dept2'), sort=False)
print(result)


##Grouping
#Iterating groups
penguins = sns.load_dataset("penguins")
penguins_grouped_island = penguins.groupby("island",dropna = True)
print(penguins.head())

#Applying mean functions
penguins.groupby(by=["island","sex"])["flipper_length_mm"].mean()

penguins.groupby(by=["species","sex"])["bill_length_mm"].median()

penguins.groupby(by=["species"])["bill_length_mm"].median()
penguins.groupby(by=["species"])["bill_length_mm"].mean()


#As index
diamonds = sns.load_dataset("diamonds")
diamonds.groupby(["cut","color"],as_index=True)["price"].median()
diamonds.groupby(["cut","color"],as_index=False)["price"].median()

##Transpose
diamonds
diamonds.T
