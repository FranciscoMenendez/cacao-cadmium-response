# Molecular and Physiological Mechanisms of the Cadmium Response in Seedlings of Two *Theobroma cacao* L. Genotypes

> **Repository for data, code, and reproducible analyses accompanying the manuscript:**  
> *Molecular and Physiological Mechanisms of the Cadmium Response in Seedlings of Two Theobroma cacao L. Genotypes*

---

## Authors (manuscript order)

- **P. Delgadillo‑Duran**\*¹  
- **F. M. Menéndez‑Burns**\*² ³  
- **C. Rodríguez‑Medina**⁴  
- **A. C. Montenegro**¹  
- **A. Istvan**⁵  
- **M. J. Guiltinan**² ³  
- **R. Yockteng**¹ ⁷ ✉️  
- **S. N. Maximova**² ³ ✉️  

\*Equal contribution

### Affiliations
1. Corporación Colombiana de Investigación Agropecuaria – **AGROSAVIA**, Centro de Investigación Tibaitatá, Mosquera, Colombia  
2. Department of Plant Science, The Pennsylvania State University, University Park, PA, USA  
3. Huck Institutes of the Life Sciences, The Pennsylvania State University, University Park, PA, USA  
4. Corporación Colombiana de Investigación Agropecuaria – **AGROSAVIA**, Centro de Investigación Palmira, Valle del Cauca, Colombia  
5. Department of Biochemistry and Molecular Biology, The Pennsylvania State University, University Park, PA, USA  
6. Department of Statistics, The Pennsylvania State University, University Park, PA, USA  
7. Institut de Systématique, Evolution, Biodiversité (UMR CNRS 7205), Muséum National d’Histoire Naturelle, France

**Corresponding authors:**  
Siela N. Maximova — <snm104@psu.edu>  
Roxana Yockteng — <ryockteng@agrosavia.co>

---

## Overview

This repository contains the cleaned data, analysis code, and figure-generation scripts used to quantify genotypic differences in cadmium (Cd) uptake and to characterize the physiological and transcriptomic responses of cacao seedlings to Cd exposure. We analyze roots and leaves across time points and genotypes, integrate DEG calls with GO/KEGG and pathway-centric gene family heatmaps, and test physiological hypotheses with gas exchange and ABA measurements.

**Highlights**
- Cd accumulation and biomass in contrasting cacao genotypes (Fig. 1)
- DEG counts, intersections, and regulation patterns by tissue, time, and genotype (Fig. 2)
- GO-term clustering and labels via semantic similarity and keyword enrichment (Fig. 3)
- Hormone-related pathways (indoles, ABA/apocarotenoids, ethylene) and signaling components (Fig. 4)
- Root-enriched functions: catabolism/detoxification, ion/macronutrient transport, metal transporters (Fig. 5)
- Leaf metabolic reprogramming: TCA, carbohydrate, fatty acids, terpenoids (Fig. 6)
- Physiology & hormones: stomatal conductance and ABA (Fig. 7)
- Working model integrating molecular and physiological responses (Fig. 8)

---

## Repository structure

```
cacao-cadmium-response/
├── data/
│   ├── raw/                    # Raw inputs (RNA-seq counts, LI-COR session files, etc.)
│   ├── interim/                # Intermediate files (tidy LI-COR exports, DEG tables)
│   └── processed/              # Final analysis-ready tables (DEG, GO/KEGG, metadata)
├── metadata/
│   ├── sample_metadata.csv     # Sample and experimental metadata
│   ├── contrasts.csv           # Linear model contrasts used for DE analysis
│   └── licor_sessions/         # Session-level notes and row-removal logs
├── scripts/
│   ├── R/
│   │   ├── 00_setup.R          # Package install/load, global options
│   │   ├── 01_qc_counts.R      # RNA-seq QC and normalization
│   │   ├── 02_deg_models.R     # Linear modeling and DEG extraction
│   │   ├── 03_go_kegg.R        # GO enrichment and KEGG mapping
│   │   ├── 04_families_maps.R  # Gene-family/orthogroup mapping & heatmaps
│   │   ├── 05_figures.R        # End-to-end figure generation
│   │   └── licor_pipeline.R    # LI-COR read/tidy/merge + stats for Gs/A & ABA
│   └── bash/
│       └── make_project_tree.sh
├── results/
│   ├── figures/                # Final figure panels (PDF/PNG/SVG)
│   ├── tables/                 # Main & supplementary tables (CSV/TSV/XLSX)
│   └── reports/                # Rendered notebooks/reports (HTML/PDF)
├── notebooks/
│   ├── figure_checks.Rmd       # Lightweight figure sanity checks
│   └── exploratory.Rmd         # Optional exploratory analyses
├── LICENSE
└── README.md
```

> **Note:** For reproducibility, all supplemental tables are provided as an Excel workbook with one sheet per table. LI-COR session files are deposited unmodified in `data/raw/` with per-session row-removal logs in `metadata/licor_sessions/` and a scripted tidy pipeline in `scripts/R/licor_pipeline.R`.

---

## Data availability

- **RNA‑seq reads:** NCBI BioProject **PRJNA943175**.  
- **Processed tables:** Provided under `results/tables/` (main and supplementary).  
- **Raw LI‑COR outputs:** Provided under `data/raw/licor/` with session-level readme and row-removal justifications.  
- **Cleaned LI‑COR tables:** Under `data/processed/` and fully reproducible via `scripts/R/licor_pipeline.R`.

---

## Reproducibility

### Requirements
- R (≥ 4.2), data.table, edgeR/DESeq2 (per your analysis), topGO/clusterProfiler, KEGGREST, ggplot2, ComplexHeatmap
- (Optional) Python ≥ 3.10 for small utilities
- LI‑COR parsing depends on readxl/openxlsx

Install core R packages:
```r
source("scripts/R/00_setup.R")
```

### End‑to‑end figures
```r
# Run individual steps or render all figures
source("scripts/R/05_figures.R")
```

### LI‑COR pipeline (Gs, A, ABA)
```r
source("scripts/R/licor_pipeline.R")
# Produces tidy tables in data/processed/ and summary stats in results/tables/
```

---

## Citation

If you use this repository, please cite the manuscript (authors in order as listed above) and this repository. A formal citation will be added upon acceptance.

---

## Funding & acknowledgments

This work was supported by AGROSAVIA (Colombia MinCiencias Agreements TV‑18/TV‑19, Project ID 1000429) and USDA‑NIFA Hatch Project #PEN04707 (Accession #1019863). We thank the Huck Institutes’ Metabolomics Core Facility (RRID: SCR_023864), the USDA‑ARS Tropical Agriculture Research Station (Mayagüez, PR) for cacao pods, Dr. Naomi Altman for statistical guidance, and BioRender® for Figure 8 template support.

---

## License

- **Code:** MIT License (see `LICENSE`)  
- **Data:** CC BY 4.0 (attribution required). Some third‑party data may have separate terms.

---

## Contact

- **Siela N. Maximova** — snm104@psu.edu  
- **Roxana Yockteng** — ryockteng@agrosavia.co
