"""Merge the AnnData objects exported by 1_preprocess.R."""

import scanpy as sc


all_cells = sc.read_h5ad("all.h5ad")
integrated_cells = sc.read_h5ad("C4_inte.h5ad")
merged = sc.concat([all_cells, integrated_cells], merge="same")
merged.write_h5ad("all_merge.h5ad")
