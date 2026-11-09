# Township-level life expectancy in Taiwan, 2000–2021

Analysis code for:

> Wang J-Z, Liu P-C, Bennett JE, Lo W-C, Lu T-H, Ezzati M, Chiang T-L, Lin H-H.
> *Evolving drivers of geographic inequalities in life expectancy under universal health coverage: a nationwide spatiotemporal analysis of Taiwan, 2000–2021.*

The pipeline fits a Bayesian spatio-temporal Poisson model to death counts by
township (n = 358), sex, 5-year age group and year, builds abridged life tables
for every posterior draw, and decomposes the life-expectancy gap between the
highest- and lowest-LE township deciles by age and cause of death (Arriaga).

## Data availability

The death counts come from the Cause of Death database of the Ministry of
Health and Welfare, accessed through the Health and Welfare Data Science Center
(HWDC) under a data-use agreement. **They are not included in this repository
and cannot be redistributed.** Researchers can apply for access through HWDC.
The pipeline starts from the aggregated counts released by HWDC,
which are git-ignored.

What *is* included:

| Path | Content | Source |
|---|---|---|
| `data/NB_mat.csv` | 358 × 358 township adjacency matrix (queen contiguity) used by the ICAR terms | `code/00_build_adjacency.R` (TWmap boundaries) |
| `data/MOI_ZH.csv` | Township code ↔ Chinese name (Big5 encoded) | Ministry of the Interior |

See [`data/README.md`](data/README.md) for the format of the restricted input
files, so the code can be run on data obtained from HWDC.

## Repository layout

```
R/                  shared functions
  lifeTable.R       abridged life table (Kannisto–Thatcher old-age ax, Coale–Demeny infant ax, Chiang SE)
  arriaga.R         Arriaga decomposition by age and by cause
code/               analysis pipeline, run in numeric order
data/               public inputs; data/restricted/ for HWDC data (not tracked)
output/             generated results (not tracked)
environment/        Python (conda) environment specification
```

## Pipeline

All R scripts assume the **repository root** as working directory. The
notebook assumes `code/` as working directory (Jupyter's default).
Code for drawing the figures is not included.

| Step | Script | Main output | Notes |
|---|---|---|---|
| 0 | `00_build_adjacency.R` | `data/NB_mat.csv` | Already included |
| 1 | `01_icar_poisson.ipynb` | `output/posterior/PPd{sex}_{chain}.csv.gz`, `PPmean_{sex}.csv`, `output/TraceSummary_{m,f}.csv` | PyMC + nutpie (JAX); run once with `SEX = 1` and once with `SEX = 2` |
| 2 | `02_le_ppmean_subgroups.R` | `output/MOIsubgroup.csv`, `LEsmoothed_topbottom.csv`, `LEsmoothed_topbottom_ageDecom.csv` | Defines top/bottom 10% township subgroups |
| 3 | `03_le_draws_server.R` | `output/posterior/LEdraws/` | Life tables for 4 chains × 1000 draws; many-core server |
| 4 | `04_e0_subset.R` | `output/posterior/LEsmoothed_{sex}.csv.gz` | Keeps e0 from every draw |
| 5 | `05_posterior_summary.R` | `output/LEsmoothed_mean.csv`, gap / change summaries | Contrasts computed within draws |
| 6 | `06_cause_model.R` | `output/DCcauseSmoothed.csv` | INLA cause-specific smoothing by subgroup |
| 7 | `07_cause_decompose.R` | `output/LEsmoothed_topbottom_*causeDecom.csv` | Point-estimate cause decomposition |
| 8 | `08_decompose_draws_server.R` | `output/decom_*_summary.csv` | Decomposition for every draw → 95% CrI |
| 9 | `09_le_raw.R` | `output/LEraw.csv` | Life tables from observed counts (validation) |
| 10 | `10_spatial_analysis.R` | `output/MortStand_bySex.csv`, Moran's I | Age-standardised mortality (WHO standard) |

Steps 3 and 8 were suggested to run on a Linux server with ~120 cores (`future::plan(multicore)`);
the per-draw posterior files total several GB.


## Software

- R 4.x: data.table, dplyr, stringr, future, future.apply, ggplot2, sf, spdep,
  INLA, and TWmap (Taiwan administrative boundaries). Versions last used: `environment/R-packages.txt`.
- Python 3.12: PyMC 5.20, nutpie 0.13, JAX 0.5, ArviZ 0.20, pandas 2.2, NumPy 1.26.
  Full conda specification: `environment/conda-spec-list.txt`
  (`conda create --name pymc_env --file environment/conda-spec-list.txt`, linux-64).

## Contact

Hsien-Ho Lin (corresponding author), Institute of Epidemiology and Preventive
Medicine, National Taiwan University.
