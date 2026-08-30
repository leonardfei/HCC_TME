# HCC_TME

Analysis code accompanying a study of tumor microenvironment subtypes in hepatocellular carcinoma (HCC). The repository contains workflows for single-cell data integration and annotation, comparisons across disease groups, macrophage and tumor-cell analyses, clustering, and cell-cell communication analysis.

## Repository structure

- `1_Data_merge_annotation/`: preprocessing, integration, and major cell-type annotation
- `2_compare/`: cancer-type comparisons and Bhattacharyya-distance analysis
- `3_mph_analysis/`: macrophage subtyping, IL1B-related analyses, and peripheral-blood trajectory analysis
- `4_cluster/`: correlation analysis for TME clusters
- `5_tumor_analysis/`: malignant-cell subtyping and downstream analyses
- `6_cellchat/`: CellChat and NicheNet analyses

## Data configuration

Large analysis objects and patient-level data are not included. By default, scripts look for external inputs in `data/`. To use another location, set an environment variable before starting R:

```bash
export HCC_TME_DATA_DIR=/path/to/HCC_TME_data
```

See `data/README.md` for the expected filenames. Only data that can be redistributed under the applicable consent, repository, and journal policies should be placed in a public repository.

## Software requirements

The analysis was written in R and uses CRAN, Bioconductor, and GitHub packages. Major dependencies include Seurat, harmony, qs, tidyverse, ggplot2, survival, survminer, GSVA, clusterProfiler, CellChat, nichenetr, CytoTRACE2, UCell, irGSEA, IOBR, and Monocle. Additional packages are loaded close to the relevant analysis steps in each script.

Python steps require `scanpy`, `omicverse`, and `matplotlib`:

```bash
python -m pip install -r requirements.txt
```

For a reproducible release, record the R and package versions used for the manuscript with `sessionInfo()` or an `renv.lock` file. Package-installation commands are intentionally not executed inside the analysis scripts.

## Suggested workflow

Run scripts from their own directories because intermediate outputs use relative filenames. The main order is:

1. Run `1_Data_merge_annotation/1_preprocess.R` and the cell-type annotation scripts in numerical order.
2. If AnnData merging is required, run `1_Data_merge_annotation/merge_anndata.py` after the two `.h5ad` files have been exported.
3. Run `1_Data_merge_annotation/6_all_subtype_merge.R` to create the combined annotations.
4. Run the comparison, macrophage, clustering, tumor, and communication workflows as required for the corresponding analyses.
5. For the peripheral-blood trajectory analysis, run `3_mph_analysis/Mph_pb_inter.R` followed by `3_mph_analysis/mph_pb_trajectory.py`.

Several scripts contain independent analysis sections and assume that upstream intermediate objects already exist. Review file paths, sample annotations, factor levels, and resource requirements before running a full workflow. Single-cell objects may require substantial memory and compute time.

## Reproducibility and privacy

- Random seeds are retained where present in the original analyses.
- Local absolute paths have been removed; external inputs are resolved through `HCC_TME_DATA_DIR`.
- Generated data objects, tables, and figures are ignored by Git to reduce the risk of committing large or sensitive files.
- Before publication, add the manuscript citation, data-accession links, exact software versions, and a license approved by all authors and the institution.
