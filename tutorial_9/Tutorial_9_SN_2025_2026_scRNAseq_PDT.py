# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 09:50:28 2025

@author: stavnaky
"""

import scanpy as sc
import scvelo as scv
import anndata
import os
# Load a toy dataset
adata = scv.datasets.pancreas()

sc.pl.umap(adata,color=["clusters_coarse","G2M_score"])


sc.tl.diffmap(adata)

# choose a root cell
# pick a progenitor cluster cell index as root
progenitor = 'Ductal'
root_cell = adata.obs[adata.obs['clusters'] == progenitor].index[0]
adata.uns['iroot'] = adata.obs_names.get_loc(root_cell)

# pseudotime - 
sc.tl.dpt(adata)
sc.pl.umap(adata, color=["dpt_pseudotime", "clusters_coarse"])

progenitor = 'Endocrine'
root_cell = adata.obs[adata.obs['clusters_coarse'] == progenitor].index[0]
adata.uns['iroot'] = adata.obs_names.get_loc(root_cell)

sc.tl.dpt(adata)
sc.pl.umap(adata, color=["dpt_pseudotime", "clusters_coarse"])

# So how can we guess the directionality?

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# 3) Velocity (dynamical model recommended)
scv.tl.recover_dynamics(adata,n_top_genes=50,n_jobs=10)           # fits gene-wise kinetics
scv.tl.velocity(adata, mode='dynamical') # direction vectors
scv.tl.velocity_graph(adata)             # transition graph

# 4) Visualize directional flow
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_loc="right")
# or arrows:
scv.pl.velocity_embedding(adata, basis="umap", arrow_length=10, arrow_size=5)

# 5) Latent time (a continuous directed time)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color="latent_time", basis="umap")

# 6) Root / terminal state guesses from velocity
scv.tl.terminal_states(adata)  # optional (velocity-based)
scv.pl.scatter(adata, color=["root_cells", "end_points"], basis="umap")

os.chdir(r"C:\Users\stavnaky\Desktop\Cancer Bioinformatics 2025_6\Tutorial_9_2025_2026 scRNAseq GSEA")
anndata.AnnData.write(adata,"pancreas_adata_with_pseudotime_and_rna_velocity.h5ad")



#Let's try with a smaller gene set
adata = scv.datasets.pancreas()
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata,n_top_genes=50,n_jobs=10)           # fits gene-wise kinetics
scv.tl.velocity(adata, mode='dynamical') # direction vectors
scv.tl.velocity_graph(adata)             # transition graph

# 4) Visualize directional flow
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_loc="right")
# or arrows:
# scv.pl.velocity_embedding(adata, basis="umap", arrow_length=10, arrow_size=5)

# 5) Latent time (a continuous directed time)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color="latent_time", basis="umap")

# 6) Root / terminal state guesses from velocity
scv.tl.terminal_states(adata)  # optional (velocity-based)
scv.pl.scatter(adata, color=["root_cells", "end_points"], basis="umap")

