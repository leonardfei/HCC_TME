# External data

Large and controlled-access data are intentionally excluded from this repository. Place the required inputs in this directory, or point `HCC_TME_DATA_DIR` to another directory containing them.

The scripts reference the following external files:

- `allcell_anno.qs`, `orign_tissue.qs`, `C3_seurat.qs`, `allcell_merge.qs`
- `all_subtype_anno2.qs`, `all_downsample.qs`, `myeloid2_anno.qs`, `Mo_Mph.qs`
- `tumor_epi.qs`, `tumor_hcc_infercnv.qs`, `tumor_hcc_noepi.qs`, `monocle2_result.qs`
- `sample.csv`, `patient_cli.csv`, `cluster4.csv`, `Myeloid.csv`, `cancersea_features.csv`
- `TCGA_LIHC_matrix.txt`, `LIHC_cli.csv`
- `GSE14520_data_surv.Rdata`, `GSE40873_data_surv.Rdata`, `GSE116174_data_surv.Rdata`
- `lr_network_human_21122021.rds`, `ligand_target_matrix_nsga2r_final.rds`, `weighted_networks_nsga2r_final.rds`

Do not commit patient-level, controlled-access, or personally identifying data.
