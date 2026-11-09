# Data

## Public files (included)

- `NB_mat.csv` — 358 × 358 binary adjacency matrix of the townships analysed
  (header `V1 … V358`). Rows and columns follow ascending township code
  (`MOI`) after excluding the offshore counties (Kinmen, Lienchiang); this must
  match the township order in `DeathCount_v5.csv`.
- `MOI_ZH.csv` — township code (`MOI`, 8-digit Ministry of the Interior code),
  Chinese county + township name. **Big5 encoded** (read with
  `fileEncoding = "BIG-5"`).

## Restricted files (`data/restricted/`, not included)

Obtained from the Health and Welfare Data Science Center (HWDC), Ministry of
Health and Welfare, Taiwan. They contain death counts in small cells and cannot
be redistributed. Their structure is documented here so that the pipeline can
be rerun on data obtained through HWDC.

### `DeathCount_v5.csv` — all-cause deaths

One row per year × sex × township × age group (2000–2021, 2 sexes, 368
townships, 19 age groups = 307,648 rows), **sorted by `YEAR`, `SEX`, `MOI`,
`AGP`**. The model and the per-draw life tables match posterior draws to rows by
position, so the order matters.

| Column | Description |
|---|---|
| `YEAR` | 2000–2021 |
| `SEX` | 1 = male, 2 = female |
| `MOI` | 8-digit township code |
| `AGP` | Lower bound of age group: 0, 1, 5, 10, …, 80, 85 (85+) |
| `POP` | Mid-year population |
| `COUNT` | Number of deaths |
| `County_ZH`, `Township_ZH` | Chinese names |
| `Urbanicity_AST2023` | Township urbanisation cluster 1 (most urban) – 8; 9 = offshore islands (excluded from all analyses) |
| `MOH` | Ministry of Health and Welfare area code |

### `DeathCount_cause.csv` — cause-specific deaths

| Column | Description |
|---|---|
| `YEAR`, `SEX`, `AGP`, `MOI` | As above |
| `CAUSE` | One of 27 underlying-cause groups (see below) |
| `COUNT` | Number of deaths |
| `POP` | Mid-year population |
| `Proportion` | `COUNT / POP` |

Cause groups: Breast Cancer, Chronic Respiratory Diseases, Cirrhosis, CKD,
Colon Cancer, Diabetes, Esophageal Cancer, Injuries, Ischemic Heart Disease,
Leukemia, Lip Cancer, Liver Cancer, Lower Respiratory Infection, Lung Cancer,
Neuro Dementia, Other Cancer, Other Causes, Other Communicable Diseases, Other
Digestive Diseases, Other Non-Communicable Diseases, Pancreatic Cancer,
Prostate Cancer, Self Harm, Stomach Cancer, Stroke & Other CVD, TB, Upper
Respiratory Infection.
