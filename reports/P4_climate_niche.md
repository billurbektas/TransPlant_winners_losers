# P4: Climate and niche

Experiment climate for the 43 included experiments from the cached CHELSA V2.1 monthly files (2008 to 2021; read, not downloaded): T = destination June to August mean temperature (day-weighted), W = destination minus origin JJA temperature, P_summer and P_annual = log destination JJA and annual precipitation, each averaged over the experiment's years. Lynn's niche optima joined to the P1 accepted species by Lynn's name: T_opt = mean_sum_sub, P_opt = log(mean_ap_sub).

## Climate

| variable | min | median | max |
| --- | --- | --- | --- |
| T (destination JJA, degC) | 8.791 | 13.281 | 21.018 |
| W (destination - origin JJA, degC) | 0.357 | 2.314 | 11.089 |
| P_summer (log mm) | 5.059 | 5.931 | 6.978 |
| P_annual (log mm) | 6.171 | 6.979 | 8.377 |

Correlation of P_annual and P_summer across the 43 included experiments: r = 0.776 (PLAN.md: about 0.78 on the earlier 41).

Partial coverage (5 experiments; CHELSA has no 2022 to 2023), averaged over the years available:

| experiment_id | destination_years_available | destination_years_expected | missing_years |
| --- | --- | --- | --- |
| FR_AlpeHuez::HIGH::LOW | 7 | 9 | 2022;2023 |
| FR_Lautaret::G::L | 5 | 6 | 2022 |
| US_Colorado::mo::pf | 4 | 6 | 2022;2023 |
| US_Colorado::mo::um | 4 | 6 | 2022;2023 |
| US_Colorado::pf::um | 4 | 6 | 2022;2023 |

Against step 00: 41 experiments with the same years agree exactly (largest difference 0). Experiments whose years changed: 1. New experiments: FR_AlpeHuez::HIGH::LOW.

| experiment_id | n_years_00 | n_years | T_00 | T | W_00 | W | P_annual_00 | P_annual |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CH_Calanda2::Cal::Nes | 3 | 4 | 14.02 | 13.699 | 3.234 | 3.234 | 7.248 | 7.194 |

### Summer temperature: CHELSA 1981-2010 climatology against the experiment JJA T

Lynn's site summer temperature (`sum_temp_site`, the bio10 analogue behind T_opt) against our destination JJA T averaged over the experiment years: r = 0.979; mean difference (ours minus climatology) 0.95 degC, range -0.27 to 1.97. Table: `data/derived/prep/site_bio10_check.csv`.

## Niche

Study species with a niche: 767 of 904; with n_obs_sub >= 30: 745. Lynn names: 785, of which 778 used (11 species from more than one Lynn name). Step 00 had 785 Lynn names, keyed by Lynn's own names.

Lynn names not used:

| reason | lynn_names |
| --- | --- |
| Lynn's niche values are NA | 2 |
| accepted species is not a study species | 4 |
| name not resolved in P1 | 1 |

### Species with more than one Lynn name (11): the name with most occurrences is used (decided 26 Oct)

The second name is a synonym or an infraspecific name recorded in one region; Lynn computed a niche per recorded name, so the same plant has two samples of its range. The values differ by little except for Solidago virgaurea, whose alpine subsp. minuta (recorded in CH_Lavey) has a niche about 5 degC colder than the species, which is used.

| species | lynn_name | regions | n_obs_sub | T_opt | AP_mm | used |
| --- | --- | --- | --- | --- | --- | --- |
| Hieracium lachenalii | Hieracium lachenalii | CH_Lavey | 21752 | 16.49 | 963 | TRUE |
| Hieracium lachenalii | Hieracium vulgatum | NO_Gudmedalen; NO_Skjellingahaugen; NO_Ulvhaugen | 1960 | 15.37 | 1034 | FALSE |
| Noccaea caerulescens | Noccaea caerulescens | NO_Ulvhaugen | 17576 | 15.21 | 906 | TRUE |
| Noccaea caerulescens | Noccaea brachypetala | FR_AlpeHuez | 2464 | 14.54 | 1099 | FALSE |
| Omalotheca norvegica | Omalotheca norvegica | NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen | 9887 | 10.58 | 1437 | TRUE |
| Omalotheca norvegica | Gnaphalium norvegicum | SE_Abisko | 357 | 10.49 | 1514 | FALSE |
| Omalotheca sylvatica | Omalotheca sylvatica | NO_Skjellingahaugen | 39892 | 15.07 | 966 | TRUE |
| Omalotheca sylvatica | Gnaphalium sylvaticum | CH_Lavey | 2061 | 15.4 | 972 | FALSE |
| Pilosella piloselloides | Pilosella piloselloides | DE_Susalps; DE_TransAlps | 8590 | 16.75 | 1118 | TRUE |
| Pilosella piloselloides | Hieracium piloselloides | CH_Lavey | 366 | 17.82 | 929 | FALSE |
| Prunella vulgaris | Prunella vulgaris | CH_Calanda2; CH_Lavey; DE_Grainau; DE_Susalps; DE_TransAlps; NO_Gudmedalen; NO_Skjellingahaugen; NO_Ulvhaugen; US_Montana | 302067 | 17.09 | 1031 | TRUE |
| Prunella vulgaris | Prunella vulgaris subsp. hispida | CN_Gongga | 47 | 18.44 | 1139 | FALSE |
| Ranunculus acris | Ranunculus acris | CH_Calanda; CH_Calanda2; DE_Grainau; DE_Susalps; DE_TransAlps; FR_AlpeHuez; IT_MatschMazia1; IT_MatschMazia2; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen; SE_Abisko | 266688 | 15.94 | 978 | TRUE |
| Ranunculus acris | Ranunculus acris subsp. friesianus | CH_Lavey | 11864 | 17.69 | 1115 | FALSE |
| Ranunculus brotherusii | Ranunculus brotherusii | CN_Heibei | 320 | 13.44 | 820 | TRUE |
| Ranunculus brotherusii | Ranunculus brotherusii var. tanguticus | CN_Gongga | 228 | 13.6 | 752 | FALSE |
| Scorzoneroides autumnalis | Scorzoneroides autumnalis | IT_MatschMazia2 | 157246 | 15.56 | 959 | TRUE |
| Scorzoneroides autumnalis | Leontodon autumnalis | CH_Lavey; DE_Grainau; DE_Susalps; DE_TransAlps; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen | 2440 | 14.5 | 1183 | FALSE |
| Solidago virgaurea | Solidago virgaurea | CH_Lavey; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; SE_Abisko | 168869 | 15.24 | 1037 | TRUE |
| Solidago virgaurea | Solidago virgaurea subsp. minuta | CH_Lavey | 2771 | 10.42 | 1719 | FALSE |
| Viola biflora | Viola biflora | CN_Heibei; NO_Gudmedalen; NO_Lavisdalen; NO_Ulvhaugen; SE_Abisko | 20341 | 11.36 | 1441 | TRUE |
| Viola biflora | Viola rockiana | CN_Gongga | 83 | 12.48 | 890 | FALSE |

### Rows where Lynn's own code-to-name mapping disagrees with ours (16)

Lynn's niche for these rows belongs to the name Lynn assigned, not to our species; joining by name keeps it away from our species.

| region | local_species_code | lynn_name | lynn_accepted | our_species |
| --- | --- | --- | --- | --- |
| CH_Lavey | Ranunculus tuberosus | Ranunculus montanus | Ranunculus montanus | Ranunculus tuberosus |
| CH_Lavey | Festuca pratensis | Festuca ovina | Festuca ovina | Lolium pratense |
| CH_Lavey | Festuca pratense | Festuca ovina | Festuca ovina | Lolium pratense |
| CH_Calanda | Festuca arundinacea | Festuca rubra | Festuca rubra | Lolium arundinaceum |
| CN_Damxung | Anemone obtusiloba | Anthoxanthum nipponicum | Anthoxanthum nipponicum | Anemonastrum obtusilobum |
| DE_Grainau | Festuca pratensis | Festuca ovina | Festuca ovina | Lolium pratense |
| DE_Susalps | Festuca pratensis | Festuca ovina | Festuca ovina | Lolium pratense |
| DE_TransAlps | Festuca pratensis | Festuca ovina | Festuca ovina | Lolium pratense |
| IT_MatschMazia1 | Festuca pratensis | Festuca ovina | Festuca ovina | Lolium pratense |
| IT_MatschMazia2 | Festuca pratensis | Festuca ovina | Festuca ovina | Lolium pratense |
| CH_Lavey | Cerastium fontanum subsp. Vulgare | Cerastium fontanum | Cerastium fontanum | Cerastium holosteoides |
| US_Montana | Hieracium_gracile | Hieracium graciliusculum | Hieracium graciliusculum | Hieracium triste |
| IN_Kashmir | Helianthus anus | Helianthus annuus | Helianthus annuus | Helianthus petiolaris |
| CN_Heibei | Stellaria arenaria | Stellaria arenarioides | Stellaria arenarioides | Arenaria hispanica |
| DE_Susalps | Vicia angustifolia | Vicia lathyroides | Vicia lathyroides | Vicia sativa |
| CN_Heibei | Koeleria cristata | Koeleria cristata |  | Koeleria macrantha |

## Warnings captured

None.

## Deviations from the plan

Alopecurus gerardii -> Alopecurus gerardi added to `prep/config/accepted_renames.csv` (decided 26 Oct), so that Lynn's and TRY's records of the species join it.
