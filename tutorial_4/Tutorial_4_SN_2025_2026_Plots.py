# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 14:16:34 2024

@author: stavn
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.stats import pearsonr



#Seaborn
np.random.seed(0)
data = pd.DataFrame({
    'x': np.random.rand(50),
    'y': np.random.rand(50),
    'size': np.random.rand(50) * 100
})


#Seaborn datasets

penguins = sns.load_dataset('penguins')
iris = sns.load_dataset('iris')
diamonds = sns.load_dataset('diamonds')
mpg = sns.load_dataset('mpg')
car_crashes = sns.load_dataset('car_crashes')
#Scatterplot

sns.relplot(data=penguins, x='bill_length_mm', y= 'bill_depth_mm', size= 'flipper_length_mm', hue="island")
sns.relplot(data=penguins, x='bill_length_mm', y= 'flipper_length_mm', hue="island",markers='*',)
sns.scatterplot(
    data=penguins,
    x='flipper_length_mm',
    y='body_mass_g',
    hue='species',
    style='species',
    palette=sns.color_palette("hls",8)
)
sns.set_theme()

font_title = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }

font_lables = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 13,
        }

#Barplot - Comparing numerical single values (mean, for example) among groups
sns.barplot(data = penguins,x='island',y='bill_depth_mm',palette="Set2")
plt.title('Penguins bill depth by island and sex',fontdict = font_title)
plt.xlabel('Island',fontdict=font_lables)
plt.ylabel('Bill depth [MM]',fontdict=font_lables)

#Boxplot - Comparing THE SPREAD of each group with information of median, quarties and outliers
sns.boxplot(data = penguins, x='island', y='flipper_length_mm',palette="Set2")
plt.title('Penguins flipper length by island',fontdict = font_title)
plt.xlabel('Island',fontdict=font_lables)
plt.ylabel('Bill depth [MM]',fontdict=font_lables)

#Violinplot - Like boxplot, but with density as well
sns.violinplot(data=mpg,x="origin",y='mpg',palette="Set2",fill=False,width=.5)
plt.title('Car fuel consumption by country of origin',fontdict = font_title)
plt.xlabel('Country of Origin',fontdict=font_lables)
plt.ylabel('MPG',fontdict=font_lables)

#Scatterplot - Using dots to represent values for two different numeric variables
sns.relplot(data=iris, x='sepal_length', y= 'petal_length', size="petal_width",
            kind='scatter',palette="Set2",hue='species')
plt.title('Petal length as function of sepal length',fontdict = font_title)
plt.xlabel('Sepal Length [mm]',fontdict=font_lables)
plt.ylabel('Petal Length [mm]',fontdict=font_lables)
# size = 'petal_width'

#Histogram - Plots the distribution of a numeric variable's values as a series of bars
sns.histplot(data=penguins,x='flipper_length_mm',hue="island")
plt.title('Flipper Length Distribution Amongst Penguins',fontdict = font_title)
plt.xlabel('Flipper Length [mm]',fontdict=font_lables)
plt.ylabel('Count',fontdict=font_lables)
plt.show()
sns.color_palette()

#Density plot - Like the histogram, but continuous
sns.kdeplot(data=diamonds,x='price',hue='cut')
plt.title('Diamonds Prices Distribution',fontdict = font_title)
plt.xlabel('Diamond price [$]',fontdict=font_lables)
plt.ylabel('Count',fontdict=font_lables)
#Line chart
sns.relplot(data=diamonds,x="carat",y="price",kind="line")

#Linear regression - Correlation between two numerical variables
sns.regplot(
    data=mpg, x="weight", y="horsepower",
    ci=99, marker="x", color=".3", line_kws=dict(color="r")
)
corr, p_value = pearsonr(diamonds["carat"], diamonds["price"])

plt.title(f"Weight as function of Horsepower. r={corr:.2f}, p={p_value:.2e}",fontdict=font_title)
plt.xlabel('Weight [KG]',fontdict=font_lables)
plt.ylabel('Horsepower',fontdict=font_lables)

#Line Chart - To show numerical variables as function of continuous numerical variables (mainly time)
x = [1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14]
y = [10, 30, 70, 90, 100, 105, 110, 113, 117, 120, 121, 122, 122, 122]

plt.plot(x,y,marker='o',linestyle='dashed')
plt.title("OD as function of Time (min)",fontdict=font_title)
plt.xlabel('Time [min]',fontdict=font_lables)
plt.ylabel('OD (600nm)',fontdict=font_lables)


#Pie charts - Distribution of categorical values among data
values= [25,15,45,15]
labels= ["Arabs","Ultra-orthodox","Secular","Religious"]
plt.pie(values,labels=labels)
plt.title("The different sectors in Israel",fontdict=font_title)

# Calculate the frequency of each island
island_counts = penguins['island'].value_counts()

# Create a pie chart
plt.pie(
    island_counts,
    labels=island_counts.index,
    autopct='%1.1f%%',  # Show percentages
    startangle=90,      # Rotate to start at the top
    colors=sns.color_palette('pastel'),  
)
plt.title("Penguin Distribution by Island",fontdict=font_title)
plt.show()

#Subplot

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# Plot on the first subplot
axes[0].plot(x, y, color='blue')
axes[0].set_title('First Plot')

# Plot on the second subplot
axes[1].bar(x, y, color='green')
axes[1].set_title('Second Plot')

