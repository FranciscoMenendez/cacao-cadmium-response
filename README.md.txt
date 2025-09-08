# Molecular and Physiological Mechanisms of the Cadmium Response in Seedlings of Two *Theobroma cacao* L. Genotypes

**Repository purpose.**  
Data, code, and supplementary materials supporting the manuscript:

> Menéndez-Burns, F. M., Delgadillo-Duran, P., Rodríguez-Medina, C., Montenegro, A. C., Istvan, A., Guiltinan, M. J., Yockteng, R., Maximova, S. N. (2025). *Molecular and Physiological Mechanisms of the Cadmium Response in Seedlings of Two Theobroma cacao L. Genotypes*. _Plant Physiology_. (DOI TBA)

This repository reproduces the transcriptomic (RNA-seq) and physiological (gas exchange, ABA) analyses and figures.

---

## Directory layout

data/
rna_seq/
raw/ # links/README pointing to PRJNA943175
processed/ # count matrices, sample metadata
gas_exchange/
licor_raw/ # original LI-6400XT Excel exports (4 sessions)
licor_qc/ # row-removal logs (what was excluded + why)
licor_clean/ # cleaned, analysis-ready CSVs
design/ # experimental design & measurement settings
hormones/
aba_raw/ # LC-MS outputs (core facility)
aba_processed/ # tidy tables for stats/plots
annotations/
kofam/ planttribes2/ go/ kegg/
code/
00_utils/ # helpers (paths, plotting, IO)
01_qc_licor/ # LICOR cleaning pipeline (R)
02_deg_edger/ # DEG contrasts (edgeR/limma)
03_go_enrichment/ # GO enrichment + semantic similarity clustering
04_kegg_pt_mapping/ # KOFAM + PlantTribes2 mapping
05_fig2_upset/ # intersections & regulation proportions
06_fig4_hormones/ # indoles, apocarotenoids, ethylene panels
07_fig5_roots/ # catabolism, ion transport, metal transporters
08_fig6_leaf_metab/ # TCA/CHO/FA/terpenoids panels
09_fig7_phys_aba/ # Gs/A + ABA stats and plots
10_fig8_model/ # inputs/exports used by schematic
90_tables/ # writes Supplementary Tables S1–S9
figures/
main/ # submission PDFs/SVGs for Figs 1–8
supp/ # Figures S1–S2
supplement/
tables_xlsx/ # single Excel workbook with sheets S1–S9
methods_snippets/ # short method blocks referenced by tables
licenses/
LICENSE-code # recommended: MIT (code)
LICENSE-data # recommended: CC BY 4.0 (datasets)
---

## Reproducibility

- **R ≥ 4.2** (recommend using `renv` for a locked environment)
- Suggested packages: `data.table`, `tidyverse`, `readxl`, `ggplot2`,  
  `edgeR`, `limma`, `ComplexHeatmap`, `clusterProfiler`/`fgsea`,  
  `ComplexUpset`/`UpSetR`, `BiocParallel`, `AnnotationDbi`.
- External tools: PlantTribes2, KofamScan/KOfam, DeepLoc2.

### Quick start

```bash
# 1) set up R environment (optional but recommended)
R -q -e 'install.packages("renv"); renv::init(); renv::snapshot()'

# 2) run LICOR cleaning to produce analysis-ready tables
Rscript code/01_qc_licor/clean_LICOR.R

# 3) run RNA-seq DE analysis and write tables/plots
Rscript code/02_deg_edger/run_deg_pipeline.R

# 4) GO & clustering, KEGG mapping, and figure generation
Rscript code/03_go_enrichment/go_semantic_clustering.R
Rscript code/04_kegg_pt_mapping/kofam_pt_mapping.R

# 5) Figures
Rscript code/05_fig2_upset/make_upset.R
Rscript code/06_fig4_hormones/make_hormone_panels.R
Rscript code/07_fig5_roots/make_root_panels.R
Rscript code/08_fig6_leaf_metab/make_leaf_metab.R
Rscript code/09_fig7_phys_aba/make_gs_aba.R
Each script writes CSVs to data/.../processed and figures to figures/main or figures/supp.
Data availability
RNA-seq reads: GenBank BioProject PRJNA943175.
LICOR raw outputs: deposit unaltered session Excel files (4 sessions) in a public repository (Zenodo/Figshare/Dryad) and link the DOI here.
Supplementary tables: provided as a single Excel workbook (Sheets S1–S9) under supplement/tables_xlsx/.
LICOR transparency policy (Table S9)
S9.1 Cleaned dataset used in analyses (tidy, one row per measurement).
S9.2 Removal log listing every excluded row and the reason (misread / error / no data).
S9.3 Session metadata (date, settings, filenames), and a link to the raw files DOI.
Citation
If you use this repository, please cite the manuscript above once DOI is available.
Contact
Open an issue in this repository for questions or reproducibility problems.
