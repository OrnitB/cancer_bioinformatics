#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 18:25:45 2026

@author: ornitb
"""

import scanpy as sc
import scvelo as scv
import gseapy as gp

adata = scv.datasets.dentategyrus()
sc.pp.log1p(adata)


# Q1
print(adata)

adata.obs.head()
adata.var.head()
adata.obs.columns

adata.shape[0]
adata.shape[1]

n_cells = adata.n_obs
n_genes = adata.n_vars

print(f"Number of cells: {n_cells}")
print(f"Number of genes: {n_genes}")


# Q2
print(f"Available columns in adata.obs: {adata.obs.columns.tolist()}")
print("Categorical columns (potential clustering methods):")
for col in adata.obs.columns:
    if adata.obs[col].dtype.name == "category":
        n_clusters = adata.obs[col].nunique()
        print(f"  {col}: {n_clusters} unique values")
        print(f"    Values: {adata.obs[col].cat.categories.tolist()}")
# Number of clustering methods: 2
# for "clusters", there are 14 unique clusters 
# for "clusters_enlarged", there are 22 unique clusters


# Q3
print("Available embeddings:", list(adata.obsm.keys()))

sc.pl.umap(adata, color="clusters",
           title="Clustering Method 1: clusters (14 clusters)",
           frameon=True, legend_loc="right margin")

sc.pl.umap(adata, color="clusters_enlarged",
           title="Clustering Method 2: clusters_enlarged (22 clusters)",
           frameon=True, legend_loc="right margin")

sc.pl.umap(adata, color=["clusters", "clusters_enlarged"],
           title=["clusters (14)", "clusters_enlarged (22)"], frameon=True,
           wspace=0.5)

# age is not a clustering method, but just for fun:
sc.pl.umap(adata, color="age(days)", title="Cells by Age (P12 vs P35)")


# Q4

sc.tl.rank_genes_groups(adata, groupby="clusters", method="wilcoxon",
                        key_added="rank_genes_clusters")
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, fontsize=8,
                        key="rank_genes_clusters")

sc.tl.rank_genes_groups(adata, groupby="clusters_enlarged", method="wilcoxon",
                        key_added="rank_genes_clusters_enlarged")
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, fontsize=6,
                        key="rank_genes_clusters_enlarged")


# Q5
sc.tl.rank_genes_groups(adata, groupby="clusters", groups=["Neuroblast"],
                        reference="rest", method="wilcoxon")

neuroblast_de = sc.get.rank_genes_groups_df(adata, group="Neuroblast")

neuroblast_top10 = neuroblast_de.head(10)[["names", "logfoldchanges", "pvals", "pvals_adj"]]

print("Top 10 genes differentiating Neuroblasts from rest:")
print(neuroblast_top10.to_string(index=False))

neuroblast_de_for_gsea = neuroblast_de.copy()


# Q6
sc.tl.rank_genes_groups(adata, groupby="clusters", groups=["Granule immature"],
                        reference="Granule mature", method="wilcoxon",
                        key_added="rank_genes_granule_immature_vs_mature")

sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False,
                        key="rank_genes_granule_immature_vs_mature")


granule_immature_de = sc.get.rank_genes_groups_df(adata, group="Granule immature",
                                                   key="rank_genes_granule_immature_vs_mature")
granule_immature_top10 = granule_immature_de.head(10)[["names", "logfoldchanges", "pvals", "pvals_adj"]]
granule_immature_de_for_gsea = granule_immature_de.copy()

print("Top 10 genes differentiating Granule Immature from Granule Mature:")
print(granule_immature_top10.to_string(index=False))


# Q7
sc.tl.rank_genes_groups(adata, groupby="age(days)", groups=["12"],
                          reference="35", method="wilcoxon",
                          key_added="rank_genes_12_vs_35")

sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False,
                        key="rank_genes_12_vs_35")


age_de = sc.get.rank_genes_groups_df(adata, group="12", key="rank_genes_12_vs_35")
age_top10 = age_de.head(10)[["names", "logfoldchanges", "pvals", "pvals_adj"]]
age_de_for_gsea = age_de.copy()

print(age_top10.to_string(index=False))


# Q8
neuroblast_ranked = neuroblast_de_for_gsea.set_index("names")["logfoldchanges"]

neuroblast_gsea = gp.prerank(
    rnk=neuroblast_ranked,
    gene_sets="GO_Biological_Process_2021",
    organism="mouse",
    outdir=None,
    seed=42,
    permutation_num=100
)

neuroblast_gsea_results = neuroblast_gsea.res2d.copy()
neuroblast_gsea_results["Adjusted P-value"] = neuroblast_gsea_results["FDR q-val"]

print("Top enriched pathways for Neuroblasts:")
print(neuroblast_gsea_results.sort_values("FDR q-val").head(10)[["Term", "NES", "FDR q-val"]])

gp.dotplot(neuroblast_gsea_results, title="Neuroblast vs Rest - Enriched Pathways",
           figsize=(8,6))


# Q9
granule_ranked = granule_immature_de_for_gsea.set_index("names")["logfoldchanges"]

granule_gsea = gp.prerank(
      rnk=granule_ranked,
      gene_sets="GO_Biological_Process_2021",
      organism="mouse",
      outdir=None,
      seed=42,
      permutation_num=100
      )

granule_gsea_results = granule_gsea.res2d.copy()
granule_gsea_results["Adjusted P-value"] = granule_gsea_results["FDR q-val"]

print("Top enriched pathways for Granule Immature vs Granule Mature:")
print(granule_gsea_results.sort_values("FDR q-val").head(10)[["Term", "NES", "FDR q-val"]])

gp.dotplot(granule_gsea_results,
             title="Granule Immature vs Granule Mature - Enriched Pathways",
             figsize=(8,6),
             top_term=10,
             cutoff=1.0)

# The results show that no pathways are significantly enriched (all FDR > 0.05)


# Q10
age_ranked = age_de_for_gsea.set_index("names")["logfoldchanges"]

age_gsea = gp.prerank(
    rnk=age_ranked,
    gene_sets="GO_Biological_Process_2021",
    organism="mouse",
    outdir=None,
    seed=42,
    permutation_num=100
    )

age_gsea_results = age_gsea.res2d.copy()
age_gsea_results["Adjusted P-value"] = age_gsea_results["FDR q-val"]

print("Top enriched pathways for 12-days vs 35-days:")
print(age_gsea_results.sort_values("FDR q-val").head(10)[["Term", "NES", "FDR q-val"]])

gp.dotplot(age_gsea_results,
           title="12-days vs 35-days - Enriched Pathways",
           figsize=(8,6),
           top_term=10,
           cutoff=1.0)

# Positive NES = pathways enriched in 12-days (P12)
# Negative NES = pathways enriched in 35-days (P35)
# The results show that no pathways are significantly enriched (all FDR > 0.05)


# Q11
# Based on the UMAP plots and the cluster annotations, the Radial Glia-like cells
# seem to be the least differentiated cell type in the lineage. They are positioned
# at the beginning of the trajectory and appear upstream of neuroblasts and granule
# cells. Since radial glia-like cells are known to act as neural progenitors
# and show less specialized gene expression compared to mature neuronal populations,
# they are likely the earliest and least differentiated cells in the dataset.


# Q12
sc.pp.neighbors(adata)

sc.tl.diffmap(adata)

root_cell_candidates = adata.obs[adata.obs["clusters"] == "Radial Glia-like"].index
print(f"Number of Radial Glia-like cells: {len(root_cell_candidates)}")

adata.uns["iroot"] = adata.obs.index.get_loc(root_cell_candidates[0])

sc.tl.dpt(adata)

sc.pl.umap(adata, color=["dpt_pseudotime", "clusters"],
           title=["Diffusion Pseudotime", "Cell Clusters"],
           cmap="viridis")

pseudotime_by_cluster = adata.obs.groupby("clusters")["dpt_pseudotime"].mean().sort_values()
print("Mean pseudotime by cell type (low = least differentiated, high = most differentiated):")
print(pseudotime_by_cluster)

print(f"Least differentiated: {pseudotime_by_cluster.index[0]}")
print(f"Most differentiated: {pseudotime_by_cluster.index[-1]}")


# Q13
print("Layers available:", list(adata.layers.keys()))

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis='umap', color='clusters', title='RNA Velocity')

scv.tl.latent_time(adata)

sc.pl.umap(adata, color=['latent_time', 'clusters'],
           title=['Latent Time', 'Cell Clusters'],
           cmap='viridis')

latent_time_by_cluster = adata.obs.groupby('clusters')['latent_time'].mean().sort_values()
print("\nMean latent time by cell type:")
print(latent_time_by_cluster)

print(f"\nEarliest cells (RNA velocity): {latent_time_by_cluster.index[0]}")
print(f"Latest cells (RNA velocity): {latent_time_by_cluster.index[-1]}")

# Endothelial cells show the lowest latent time (0.08), but these are blood
# vessel cells and not part of the neuronal lineage. Within the neuronal lineage,
# nIPC and Radial Glia-like cells are the earliest (0.36-0.38). Granule immature and
# Granule mature have nearly identical latent times (0.84 vs 0.83), suggesting
# they are at similar stages in the RNA velocity timeline.
# nIPC/Radial Glia-like > Neuroblast > Granule cells


# Q14
# The three methods gave us different results because they're measuring different things.
# Diffusion pseudotime (DPT) is basically measuring how "far" each cell is from
# the root cell I picked based on their overall gene expression. It's like
# measuring distance on a map, granule mature cells are the most different from
# the starting point, so they get the highest pseudotime values. The problem is that this 
# method depends on which cell I choose as the root.
# Rna velocity looks at something completely different, it compares unspliced vs spliced
# mRNA to figure out which genes are being actively transcribed right now. This
# tells me where cells are "going" rather than where they "are". It doesn't need me to
# define a starting point, it figures out the direction from the data itself.
# Latent time comes from RNA velocity and captures the timing of when transcriptional
# changes are happening.
# The interesting finding is that in DPT, Granule immature and Granule mature are clearly
# separated (immature comes before mature), but in latent time they're almost identical
# (~0.83-0.84). This actually makes biological sense - both are terminally differentiated
# neurons that have stopped actively changing. Their gene expression profiles are
# different (which DPT captures), but their transcriptional dynamics are similar
# because neither is actively transitioning to a new state anymore.
# Another thing worth noting is that Endothelial cells showed the lowest latent time,
# but they're blood vessel cells and not part of the neuronal lineage at all.
# RNA velocity just analyzes each cell independently without knowing which lineage
# it belongs to.

# So basically: DPT tells us the differentiation order along the trajectory we defined,
# RNA velocity shows us the direction cells are moving and predicts their future states,
# and latent time shows when cells are actively changing vs when they've stabilized

# Each method gives us a different piece of the puzzle, and together they help us
# understand both the static states and the dynamics of cell differentiation in the
# dentate gyrus.