"""Run pyVIA trajectory analysis for peripheral-blood macrophages."""

import matplotlib.pyplot as plt
import omicverse as ov
import scanpy as sc


ov.utils.ov_plot_set()

adata = sc.read_h5ad("mph_pb_downsample.h5ad")
if adata.raw is not None:
    adata = adata.raw.to_adata()
adata.var_names_make_unique()
adata.obs_names_make_unique()
sc.tl.pca(adata, svd_solver="arpack", n_comps=200)

trajectory = ov.single.pyVIA(
    adata=adata,
    adata_key="X_pca",
    adata_ncomps=80,
    basis="X_umap",
    clusters="cl5",
    knn=30,
    random_seed=4,
    root_user="None",
)
trajectory.run()

figure, axis = plt.subplots(1, 1, figsize=(4, 4))
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=["cl5"],
    frameon=False,
    ncols=1,
    wspace=0.5,
    show=False,
    ax=axis,
)
figure.savefig("via_fig1.png", dpi=300, bbox_inches="tight")

figure, _, _ = trajectory.plot_piechart_graph(
    clusters="cl5", cmap="Reds", dpi=80, show_legend=False,
    ax_text=False, fontsize=4,
)
figure.savefig("via_fig2.png", dpi=300, bbox_inches="tight")

genes = ["CXCR4", "ELL2", "CCL4", "CEACAM8"]
figure, _ = trajectory.plot_clustergraph(gene_list=genes, figsize=(12, 3))
figure.savefig("via_fig2_1.png", dpi=300, bbox_inches="tight")

figure, _ = trajectory.plot_stream(
    basis="X_umap", clusters="cl5", density_grid=0.8,
    scatter_size=30, scatter_alpha=0.3, linewidth=0.5,
)
figure.savefig("via_fig4.png", dpi=300, bbox_inches="tight")

figure, _ = trajectory.plot_stream(
    basis="X_umap", density_grid=0.8, scatter_size=30,
    color_scheme="time", linewidth=0.5, min_mass=1, cutoff_perc=5,
    scatter_alpha=0.3, marker_edgewidth=0.1, density_stream=2,
    smooth_transition=1, smooth_grid=0.5,
)
figure.savefig("via_fig5.png", dpi=300, bbox_inches="tight")
