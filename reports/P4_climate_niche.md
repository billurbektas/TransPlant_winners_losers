# P4: Climate and niche

Experiment climate for the 43 included experiments from the cached CHELSA V2.1 monthly files (2008 to 2021; read, not downloaded): T = destination June to August mean temperature (day-weighted), W = destination minus origin JJA temperature, P_summer and P_annual = log destination JJA and annual precipitation, each averaged over the experiment's years. Niche optima from P6's GBIF niche (`prep/gbif_niche.R`; GBIF.org (28 October 2026) GBIF Occurrence Download https://doi.org/10.15468/dl.7p28y9): T_opt = mean CHELSA 1981-2010 bio10 and P_opt = log mean bio18 over one value per species x CHELSA cell. Lynn's niche file is no longer read (decided 28 Oct).

## Climate

| variable | min | median | max |
| --- | --- | --- | --- |
| T (destination JJA, degC) | 8.791 | 13.281 | 21.018 |
| W (destination - origin JJA, degC) | 0.357 | 2.314 | 11.089 |
| P_summer (log mm) | 5.059 | 5.931 | 6.978 |
| P_annual (log mm) | 6.171 | 6.979 | 8.377 |

Correlation of P_annual and P_summer across the 43 included experiments: r = 0.776.

Partial coverage (5 experiments; CHELSA has no 2022 to 2023), averaged over the years available:

| experiment_id | destination_years_available | destination_years_expected | missing_years |
| --- | --- | --- | --- |
| FR_AlpeHuez::HIGH::LOW | 7 | 9 | 2022;2023 |
| FR_Lautaret::G::L | 5 | 6 | 2022 |
| US_Colorado::mo::pf | 4 | 6 | 2022;2023 |
| US_Colorado::mo::um | 4 | 6 | 2022;2023 |
| US_Colorado::pf::um | 4 | 6 | 2022;2023 |

The comparison with the frozen legacy climate and niche of step 00 is in `reports/P5_compare_to_frozen.md`.

### The niche's climatology at the destination sites

CHELSA 1981-2010 bio10 (mean temperature of the warmest quarter) and bio18 (precipitation of the warmest quarter) at the 43 destination sites, against the ladder's destination JJA T and JJA P averaged over the experiment years: T against bio10 r = 0.979, mean difference (ours minus climatology) 0.95 degC, range -0.27 to 1.97; log JJA P against log bio18 r = 0.944, median ratio JJA P / bio18 1.03. Table: `data/derived/prep/site_bioclim_check.csv`.

## Niche

Study species with a niche: 901 of 903; with n_obs_sub >= 30: 880. Counts at every cleaning step, the key ledger and the comparison with Lynn's niche are in `reports/P6_gbif_niche.md`.

## Warnings captured

None.

## Deviations from the plan

Niche from P6's GBIF download instead of Lynn's file, P_opt from bio18 instead of bio12 (decided 28 Oct); Alopecurus gerardii -> Alopecurus gerardi in `prep/config/accepted_renames.csv` (decided 26 Oct).
