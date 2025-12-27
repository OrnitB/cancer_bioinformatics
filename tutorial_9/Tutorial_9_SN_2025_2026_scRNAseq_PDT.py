# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 09:50:28 2025

@author: stavnaky
"""

import scanpy as sc
import scvelo as scv
from pathlib import Path

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

out_dir = Path("/Users/ornitb/Documents/Technion/research_and_studies/courses/ביואינפורמטיקה של סרטן/cancer_bioinformatics/tutorial_9")
out_path = out_dir / "pancreas_adata_with_pseudotime_and_rna_velocity.h5ad"

adata.write_h5ad(out_path)


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



# Class work
adata = sc.read("/Users/ornitb/Documents/Technion/research_and_studies/courses/ביואינפורמטיקה של סרטן/cancer_bioinformatics/tutorial_9/data/gastrulation_e75.h5ad")

sc.pl.umap(adata,color="celltype")
adata.shape

sc.pp.neighbors(adata,n_neighbors=15)
sc.tl.diffmap(adata)

# choose a root cell
# pick a progenitor cluster cell index as root
progenitor = 'Blood progenitors 1'
root_cell = adata.obs[adata.obs['celltype'] == progenitor].index[0]
adata.uns['iroot'] = adata.obs_names.get_loc(root_cell)

# pseudotime - 
sc.tl.dpt(adata)
sc.pl.umap(adata, color=["dpt_pseudotime", "celltype"])

progenitor = 'Erythroid1'
root_cell = adata.obs[adata.obs['celltype'] == progenitor].index[0]
adata.uns['iroot'] = adata.obs_names.get_loc(root_cell)

sc.tl.dpt(adata)
sc.pl.umap(adata, color=["dpt_pseudotime", "celltype"])

# So how can we guess the directionality?

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# 3) Velocity (dynamical model recommended)
scv.tl.recover_dynamics(adata,n_top_genes=30,n_jobs=10)           # fits gene-wise kinetics
scv.tl.velocity(adata, mode='dynamical') # direction vectors
scv.tl.velocity_graph(adata)             # transition graph

# 4) Visualize directional flow
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_loc="right")
# or arrows:
scv.pl.velocity_embedding(adata, basis="umap", arrow_length=10, arrow_size=5)
# 5) Latent time (a continuous directed time)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color="latent_time", basis="umap",color_map="viridis")

# 6) Root / terminal state guesses from velocity
scv.tl.terminal_states(adata)  # optional (velocity-based)
scv.pl.scatter(adata, color=["root_cells", "end_points"], basis="umap",color_map="viridis")
