# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 21:52:07 2024

@author: stavn
"""

#Differential Expression using processed and clustered COVID data
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import matplotlib.pyplot as plt
import warnings
import os
import urllib.request
import anndata

adata = sc.read_h5ad("/Users/ornitb/Documents/Technion/research_and_studies/courses/ביואינפורמטיקה של סרטן/cancer_bioinformatics/tutorial_8/scanpy_covid_qc_dr_scanorama_cl.h5ad")
adata

#Get to know Adata
#Expression matrix
expression_matrix = adata.X
#.var- Metadata for variables (genes)
gene_data = adata.var
#.obs- Metadata for observations (cells)
cell_data = adata.obs
#.obsm- Multidimensional data for observation (cells)
cells_multidimensional_data = adata.obsm
#.varm - Multidimensional data for variables (genes)
genes_multidimensional_data = adata.varm
#.uns - Unstructured data for storing results or other arbitrary information
unsturctured = adata.uns
# Visualize clusters and cell types
sc.pl.umap(adata, color=['louvain_0.6'],legend_loc='on data')

#T-test for differential expression
sc.tl.rank_genes_groups(adata, 'louvain_0.6', method='t-test', key_added = "t-test")
#Plot
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=True, key = "t-test")

#Wilcoxon for differential expression
sc.tl.rank_genes_groups(adata, 'louvain_0.6', method='wilcoxon', key_added = "wilcoxon")
#Plot
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key="wilcoxon")

#Differential expression between clusters 4 and 7
sc.tl.rank_genes_groups(adata, 'louvain_0.6', groups=['4'], reference='7', method='wilcoxon')
#Plot
sc.pl.rank_genes_groups(adata, groups=['4'], n_genes=20)

#Differential expression between covid/ctrl cells
sc.tl.rank_genes_groups(adata, 'type', method='wilcoxon',key_added = "wilcoxon_types")
#Plot
sc.pl.rank_genes_groups(adata, n_genes=20,key="wilcoxon_types")
#Visualize types of data
sc.pl.umap(adata, color=['type'],legend_loc='lower right')

#Heatmap
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon", groupby="louvain_0.6", show_gene_labels=True)

#Dotplot
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon", groupby="louvain_0.6")


# Class Work
# 2
adata.obs['batch']
adata.obs['sample']

sc.tl.rank_genes_groups(adata, groupby='batch', method='t-test', key_added='batch_ttest')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=True, key="t-test")
# 3
sc.tl.rank_genes_groups(adata, groupby='sample', method='wilcoxon', key_added='sample_wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=True, key="sample_wilcoxon")
# 4
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, key='wilcoxon', groupby="sample", show_gene_labels=True)


##GSEA

# UMAP
sc.pl.umap(adata,color="louvain_0.6")


#  Find differentially expressed genes 
sc.tl.rank_genes_groups(adata, groupby="louvain_0.6", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

#Get table format and filter 2 relevant columns
cluster = ["4","7"]
ranked_genes = sc.get.rank_genes_groups_df(adata, group=cluster)
gene_list = ranked_genes[["names", "scores"]]

#Run GSEA prerank using KEGG 2019 Human gene set 
pre_res = gp.prerank(
    rnk=gene_list,
    gene_sets="KEGG_2019_Human",
    processes=4,
    permutation_num=1000,  
    outdir=None,          
    seed=42
)

pre_res
# Visualize
summary = pre_res.res2d
gp.dotplot(summary, column= 'NES', title=f"Cluster {cluster} enrichment",
           top_term=20)


# Class Work
# 1
adata = sc.datasets.pbmc3k_processed()

# 2
sc.pl.umap(adata, color=['louvain'])

# 3
sc.tl.rank_genes_groups(adata, groupby="louvain", method="t-test", groups=["CD8 T cells"])
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

# 4
cluster = ["CD8 T cells"]
ranked_genes = sc.get.rank_genes_groups_df(adata, group=cluster)
gene_list = ranked_genes[["names", "scores"]]
gene_list

# 5
pre_res = gp.prerank(
    rnk=gene_list,
    gene_sets="KEGG_2019_Human",
    processes=4,
    permutation_num=1000,  
    outdir=None,          
    seed=42
)

# 6
summary = pre_res.res2d
gp.dotplot(summary, column= 'NES', title=f"Cluster {cluster} enrichment",
           top_term=20)

