# Step 00: Freeze inputs

The script read 8 objects from the legacy targets store and Lynn's niche file and wrote them unchanged in content to `data/derived/00_inputs`, one `.rds` per object, with `manifest.csv`. The only changes are: the four site coordinates joined onto `experiments`; T, W and log precipitation derived in `experiment_climate`; the legacy `relative_cover` dropped from `community`, because PLAN.md binds that name to the cover step 01 recomputes after the deep-lineage removal; traits pivoted to one row per species; and niche columns renamed to `T_opt` and `P_opt`. No filter changes the analysis set in this step.

Run with R version 4.6.1 (2026-06-24), targets 1.12.0.

## Counts

| object | step | rows |
| --- | --- | --- |
| sites | read and written | 59 |
| experiments | read | 44 |
| experiments | written (all rows) | 44 |
| experiments | of which included | 42 |
| experiment_climate | read and written | 42 |
| community | read | 50260 |
| community | written | 50260 |
| community | of which cover NA (US_Arizona, supplied relative cover only) | 244 |
| community | of which year and relative_year NA (FR_AlpeHuez) | 1024 |
| species_taxonomy | read and written | 797 |
| traits | read (species x trait rows) | 1411 |
| traits | after dropping blank-name rows | 1410 |
| traits | written (species) | 558 |
| functional_groups | read and written | 782 |
| phylogeny_axes | read and written (species) | 526 |
| niche_optima | read (turf x species rows) | 17112 |
| niche_optima | written (species) | 785 |

Experiments excluded by the legacy pipeline (kept in `experiments.rds` with `included == FALSE`): FR_AlpeHuez::HIGH::LOW (warmed observations have missing Year); IN_Kashmir::HIGH::LOW (insufficient duration).

## Key results

### Manifest

| file | rows | columns | sha256 | source | target_data_hash | target_time_built | source_sha256 | source_modified | units |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| experiments.rds | 44 | 17 | 119f89386e398b5833718edb7b40f1a75dd4eea578e3948828f8656c3dd21cfe | experiments; sites | 55b8b0b64789ea5d; db900525ccc1b622 | 2026-10-01 10:35:45; 2026-10-01 10:35:45 |  |  |  |
| sites.rds | 59 | 11 | 3a0c42f2919855e998ec4730658fa67f29a7e3256a57740dfdabdf5ea2a732e5 | sites | db900525ccc1b622 | 2026-10-01 10:35:45 |  |  |  |
| experiment_climate.rds | 42 | 18 | a4d0af8f7d01915619bbd345526ab602a44c9f1a484e342a59cd3fe667ed9db3 | experiment_climate | ab00a7383fe8e62d | 2026-10-01 10:36:12 |  |  |  |
| community.rds | 50260 | 18 | b74f3946cb8320ea03b765bccf7489739bcd82a07f252c01243c89a702888fe9 | analysis_community | b7e64c95af08fc74 | 2026-10-01 10:36:48 |  |  |  |
| species_taxonomy.rds | 797 | 21 | 371ae3f0761342d08f4a555f30ff2702c9800fe5dc981637c5c465fefe312091 | community_taxonomy_species | 2b93fbc67207ddac | 2026-10-01 10:36:12 |  |  |  |
| traits.rds | 558 | 7 | 3ee1dbd0c3e6e66d8a808b1b7f143e2d59c231dd1f90ad6378c816c13a1806d3 | try_trait_outputs | da4cc7a4839d8a1e | 2026-10-01 10:38:47 |  |  | LDMC: g/g; LMA: mg mm-2; height: m |
| functional_groups.rds | 782 | 5 | 15fa57167fb3520346eaa1bfb2291355a256dc10649ea114f879808dad1590bc | functional_groups | 7ee4a31bf910981a | 2026-10-01 10:36:47 |  |  |  |
| phylogeny_axes.rds | 526 | 521 | d8183ceccd664bce6e8596b06ebae5b548ffafedf743efd0804223efd94838ce | resident_pcoa_candidate | fb7025fab8bf9d17 | 2026-10-01 10:48:00 |  |  |  |
| niche_optima.rds | 785 | 5 | f3d6e905962c3b353ecc462f8feca51d2be2e9782f298f4669db11cdb2ac30f2 | data/MacroContextTransplantDat_2.csv |  |  | 51167f1172b213fefc49ec813f45c6a75159cace26ecb612285ab3a7b099614c | 2026-10-23 12:01:10 |  |

### Community (`analysis_community`, legacy `relative_cover` dropped)

| regions | plots | plot_years | species | relative_year | cover_NA_rows | cover_NA_plot_years | year_NA_rows |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 22 | 763 | 2961 | 797 | 1 to 9 | 244 | 84 | 1024 |

### Climate of the 42 included experiments

T = destination JJA temperature (degC), W = destination minus origin JJA temperature, P_summer = log destination JJA precipitation, P_annual = log destination annual precipitation (legacy sum of monthly CHELSA precipitation per year, averaged over the experiment's years; `legacy/R/data/climate.R:346`, not recomputed here). Not standardised.

| variable | min | median | max |
| --- | --- | --- | --- |
| T | 8.791 | 13.224 | 21.018 |
| W | 0.357 | 2.287 | 11.089 |
| P_summer | 5.059 | 5.922 | 6.978 |
| P_annual | 6.171 | 6.864 | 8.377 |

Correlation of P_annual and P_summer across the 42 included experiments: r = 0.778 (PLAN.md expects about 0.78).

Experiments with partial climate coverage (4), averaged over the years available in the CHELSA cache:

| experiment_id | origin_years_available | origin_years_expected | destination_years_available | destination_years_expected | missing_years |
| --- | --- | --- | --- | --- | --- |
| FR_Lautaret::G::L | 5 | 6 | 5 | 6 | 2022 |
| US_Colorado::mo::pf | 4 | 6 | 4 | 6 | 2022;2023 |
| US_Colorado::mo::um | 4 | 6 | 4 | 6 | 2022;2023 |
| US_Colorado::pf::um | 4 | 6 | 4 | 6 | 2022;2023 |

### Species properties

- Traits: 558 species with at least one of LDMC, LMA, height; 434 with LDMC, 499 with LMA, 477 with height. Legacy per-species medians; 1 blank-name LMA row dropped.
- Functional groups: 782 species, 142 graminoid, 35 woody.
- Phylogeny: 526 resident species, 520 PCoA axes with eigenvalues, all exported. 20 axes exceed the broken-stick expectation; PEV1 to PEV3 carry 43.3% of the positive variation.
- Niche optima: 785 species in Lynn's file, keyed by Lynn's name (no join to community names here); 2 without a finite T_opt or P_opt; 761 with `n_obs_sub >= 30`.

Niche provenance: `MacroContextTransplantDat_2.csv` was built by J. Lynn (11.2024) with `DataProcessing_MacroContext.R` in TransPlant-Updates_oct2021, following the method of Lynn et al. 2021 (Global Change Biology). Species niche summaries are CHELSA 1981 to 2010 climatologies over GBIF occurrences; the `_sub` columns use occurrences thinned to one per CHELSA (30 arc-second) cell, and `n_obs_sub` is the number of thinned occurrences behind them. None of Lynn's processing is rerun; only columns are taken. Source of this description: `legacy/R/additional_checks/thermal_niche/random_forest_thermal_niche.R`, lines 31 to 37.

`T_opt` = `mean_sum_sub`, taken as summer temperature (the bio10 analogue); this reading of the `sum` columns still needs confirming. `P_opt` = log(`mean_ap_sub`), annual precipitation (bio12). `mean_mat_sub` (bio1) is carried only for the bio10 against JJA check of Fig S5.

All 100 columns of the CSV: `Region`, `year_start`, `originSiteID`, `originBlockID`, `destSiteID`, `destBlockID`, `destPlotID`, `Treatment`, `turfID`, `SpeciesName`, `cover_start`, `rel_cov_start`, `Elevation`, `year_end`, `cover_end`, `rel_cov_end`, `species`, `mean_mat_full`, `mean_ap_full`, `mean_sum_full`, `mean_pet_full`, `mean_gsl_full`, `median_mat_full`, `median_ap_full`, `median_sum_full`, `median_pet_full`, `median_gsl_full`, `min_mat_full`, `min5_mat_full`, `max_mat_full`, `max95_mat_full`, `min_ap_full`, `min5_ap_full`, `max_ap_full`, `max95_ap_full`, `min_sum_full`, `min5_sum_full`, `max_sum_full`, `max95_sum_full`, `min_pet_full`, `min5_pet_full`, `max_pet_full`, `max95_pet_full`, `min_gsl_full`, `min5_gsl_full`, `max_gsl_full`, `max95_gsl_full`, `n_obs_full`, `mean_mat_sub`, `mean_ap_sub`, `mean_sum_sub`, `mean_pet_sub`, `mean_gsl_sub`, `median_mat_sub`, `median_ap_sub`, `median_sum_sub`, `median_pet_sub`, `median_gsl_sub`, `min_mat_sub`, `min5_mat_sub`, `max_mat_sub`, `max95_mat_sub`, `min_ap_sub`, `min5_ap_sub`, `max_ap_sub`, `max95_ap_sub`, `min_sum_sub`, `min5_sum_sub`, `max_sum_sub`, `max95_sum_sub`, `min_pet_sub`, `min5_pet_sub`, `max_pet_sub`, `max95_pet_sub`, `min_gsl_sub`, `min5_gsl_sub`, `max_gsl_sub`, `max95_gsl_sub`, `n_obs_sub`, `ahull`, `aoo_2`, `aoo`, `dup`, `dif_sub_rel`, `simp_dif_rel`, `exp_dur`, `extinct`, `colonization`, `Longitude`, `Latitude`, `mat_site`, `sum_temp_site`, `ap_site`, `pet_site`, `gsl_site`, `mat_meandiff`, `sum_meandiff`, `ap_meandiff`, `gsl_meandiff`, `pet_meandiff`.

## Diagnostic figures

None in this step.

## Warnings captured

None.

## Deviations from PLAN.md

None in the analysis. Open question for step 07: the PCoA covers resident species only, so whether colonizers get a phylogeny property, and from which PCoA, is not settled.
