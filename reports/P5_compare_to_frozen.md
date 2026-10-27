# P5: Comparison with the frozen legacy inputs

Every file the analysis reads from `data/derived/prep/` was compared with its counterpart in the frozen legacy inputs of step 00, row by row, after mapping legacy names to the P1 species. Each of the 13797 differences has a named cause, classed as a correction (a rule, ledger or decision of P1 to P4) or new data (the network SQLite against the RData); the script stops on any difference without one. The differences are listed in full in `reports/P5_differences/`. Step 00 is retired with this step: its script and report remain in the git history, its data are archived outside the repository, and steps 01 to 10 read `data/derived/prep/` only. This report carries every comparison with step 00 that the P1 to P4 reports made before.

## Summary

| file | rows_step00 | rows_prep | columns_added | columns_dropped | same_name | renamed_or_merged | lost | gained | differences |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| community | 50260 | 58803 (54980 analysis rows) | plot_id, name_original, n_records, drop_reason, is_vegetation | turf_id, supplied_relative_cover, elevation, stable_plot_id, year_established, p1_species | 721 | 42 | 1 | 130 | 10869 |
| sites | 59 | 59 | gradient |  |  |  |  |  | 6 |
| experiments | 44 | 44 |  | origin_longitude, origin_latitude, destination_longitude, destination_latitude |  |  |  |  | 7 |
| species_taxonomy | 797 | 1126 | accepted_species, accepted_family, taxonomic_status, decided_by, needs_review_reason, n_records, regions | submitted_name, matched_name, score, name_normalized, accepted_normalized, decision, accepted_name_override, family, note, accepted_name_reviewed, needs_review, original_word_count, accepted_word_count, explicit_unspecified_species, family_name, include_species_level, resolution_decision |  |  |  |  | 376 |
| traits | 558 | 674 | n_datasets_LMA, n_datasets_LDMC, n_datasets_height |  | 522 | 23 | 13 | 133 | 1355 |
| functional_groups | 782 | 904 | genus, woody_basis | functional_group | 712 | 69 | 1 | 147 | 224 |
| phylogeny_axes | 526 | 904 | PEV521 to PEV903 |  | 476 | 49 | 1 | 389 | 915 |
| experiment_climate | 42 | 43 |  |  |  |  |  |  | 2 |
| niche_optima | 785 | 767 | lynn_name_used, lynn_names_merged, n_lynn_names |  |  |  |  |  | 43 |

Differences by kind:

| file | correction | new data | correction + new data |
| --- | --- | --- | --- |
| community | 4939 | 5930 | 0 |
| experiment_climate | 1 | 1 | 0 |
| experiments | 2 | 5 | 0 |
| functional_groups | 224 | 0 | 0 |
| niche_optima | 43 | 0 | 0 |
| phylogeny_axes | 915 | 0 | 0 |
| sites | 0 | 6 | 0 |
| species_taxonomy | 365 | 10 | 1 |
| traits | 1355 | 0 | 0 |

## Community

Frozen rows (recorded names) were re-keyed to the P1 species through region x local code and summed per plot-year; prep rows are one per plot-year x species, plus every other vegetation record with its `drop_reason`. Rows of each cause (they add up to both tables):

| kind | cause | rows_step00 | rows_prep |
| --- | --- | --- | --- |
|  | same row, same cover | 48978 | 48963 |
| correction | vegetation record kept with its drop_reason (relative-cover denominator, step 01) | 0 | 3174 |
| correction | species in one plot-year of its experiment x treatment: legacy filter not applied (legacy/R/data/community.R:242-245) | 0 | 964 |
| correction | CH_Calanda2 2020 kept (legacy drop C6 removed, decided 26 Oct) | 0 | 798 |
| correction | cover summed over names merged into one species | 1 | 1 |
| correction | excluded by decision (prep/config/taxonomy_decisions.csv) | 2 | 0 |
| new data | FR_AlpeHuez: survey years in the SQLite only (the RData had no year) | 1024 | 4654 |
| new data | US_Arizona: absolute cover from the SQLite (step 00 had supplied relative cover only) | 244 | 244 |
| new data | SQLite value differs from the RData (data cleaning repository) | 5 | 5 |
| new data | SQLite recodes the record (FR_Lautaret Pilosella officinarum -> Pilosella) | 6 | 0 |
|  | total | 50260 | 58803 |

Record by record, the SQLite equals the RData outside FR_AlpeHuez except 355 US_Arizona covers (absent from the RData), 5 CN_Heibei covers and 6 FR_Lautaret records recoded.

| region | year | destination_plot | local_species_code | cover_rdata | cover_sqlite |
| --- | --- | --- | --- | --- | --- |
| CN_Heibei | 2013 | 3200_3200_2 | Lancea tibetica | 1 | 2 |
| CN_Heibei | 2014 | 3200_3200_1 | Potentilla bifurca | 1 | 2 |
| CN_Heibei | 2014 | 3200_3200_1 | Potentilla nivea | 15 | 30 |
| CN_Heibei | 2014 | 3200_3200_2 | Potentilla nivea | 9 | 18 |
| CN_Heibei | 2014 | 3200_3200_3 | Potentilla anserina | 5 | 10 |

Species of the analysis rows: 764 legacy species, 887 P1 species; same name 721, renamed or merged 42, lost 1, gained 130. Lost and gained species by cause:

| status | cause | species |
| --- | --- | --- |
| gained | FR_AlpeHuez: survey years in the SQLite only (the RData had no year) | 9 |
| gained | FR_AlpeHuez: survey years in the SQLite only (the RData had no year) + species in one plot-year of its experiment x treatment: legacy filter not applied (legacy/R/data/community.R:242-245) | 1 |
| gained | species in one plot-year of its experiment x treatment: legacy filter not applied (legacy/R/data/community.R:242-245) | 120 |
| lost | excluded by decision (prep/config/taxonomy_decisions.csv) | 1 |

## Sites and experiments

| key | difference | cause |
| --- | --- | --- |
| FR_AlpeHuez::HIGH | survey_year_min: NA -> 2015 | FR_AlpeHuez: survey years in the SQLite only (the RData had no year) |
| FR_AlpeHuez::LOW | survey_year_min: NA -> 2015 | FR_AlpeHuez: survey years in the SQLite only (the RData had no year) |
| FR_AlpeHuez::HIGH | survey_year_max: NA -> 2023 | FR_AlpeHuez: survey years in the SQLite only (the RData had no year) |
| FR_AlpeHuez::LOW | survey_year_max: NA -> 2023 | FR_AlpeHuez: survey years in the SQLite only (the RData had no year) |
| FR_AlpeHuez::HIGH | supplied_year_range: NA -> 9 | FR_AlpeHuez: survey years in the SQLite only (the RData had no year) |
| FR_AlpeHuez::LOW | supplied_year_range: NA -> 9 | FR_AlpeHuez: survey years in the SQLite only (the RData had no year) |

| key | difference | cause |
| --- | --- | --- |
| FR_AlpeHuez::HIGH::LOW | first_year: NA -> 2015 | experiment included: FR_AlpeHuez years in the SQLite |
| CH_Calanda2::Cal::Nes | last_year: 2019 -> 2020 | experiment years changed (CH_Calanda2 2020 kept) |
| FR_AlpeHuez::HIGH::LOW | last_year: NA -> 2023 | experiment included: FR_AlpeHuez years in the SQLite |
| CH_Calanda2::Cal::Nes | year_range: 3 -> 4 | experiment years changed (CH_Calanda2 2020 kept) |
| FR_AlpeHuez::HIGH::LOW | year_range: NA -> 9 | experiment included: FR_AlpeHuez years in the SQLite |
| FR_AlpeHuez::HIGH::LOW | exclusion_reason: warmed observations have missing Year -> NA | experiment included: FR_AlpeHuez years in the SQLite |
| FR_AlpeHuez::HIGH::LOW | included: FALSE -> TRUE | experiment included: FR_AlpeHuez years in the SQLite |

## Species taxonomy

Recorded names: 797 in step 00 (species level only), 1126 in prep (every recorded name). By cause:

| kind | cause | differences |
| --- | --- | --- |
| correction | species in one plot-year of its experiment x treatment: legacy filter not applied (legacy/R/data/community.R:242-245) | 134 |
| correction | name below species level, now kept with its identification level | 131 |
| correction | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 39 |
| correction | unidentified morphotaxon, now listed (prep/config/unidentified_taxa.csv) | 35 |
| correction | name recorded only in establishment years (dropped by both pipelines; P1 resolves every recorded name) | 17 |
| correction | taxonomy decision (prep/config/taxonomy_decisions.csv) | 5 |
| correction | vegetation record kept with its drop_reason (relative-cover denominator, step 01) | 2 |
| correction | cover summed over names merged into one species | 1 |
| correction | excluded by decision (prep/config/taxonomy_decisions.csv) | 1 |
| correction + new data | FR_AlpeHuez: survey years in the SQLite only (the RData had no year) + species in one plot-year of its experiment x treatment: legacy filter not applied (legacy/R/data/community.R:242-245) | 1 |
| new data | name new in the SQLite (FR_AlpeHuez) | 9 |
| new data | SQLite recodes the record (FR_Lautaret Pilosella officinarum -> Pilosella) | 1 |

### Accepted species that differ from the frozen legacy taxonomy (45 names)

| name_original | legacy_accepted | accepted_species | cause | n_records |
| --- | --- | --- | --- | --- |
| Festuca pratensis | Festuca ovina | Lolium pratense | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 365 |
| Potentilla leuconota | Potentilla leuconota | Argentina leuconota | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 199 |
| Agrostis alpina | Agrostis rupestris subsp. rupestris | Alpagrostis alpina | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 194 |
| Helictotrichon pubescens | Helictotrichon pubescens | Avenula pubescens | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 173 |
| Gentiana campestris | Gentianella germanica | Gentianella campestris | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 154 |
| Polygonum macrophyllum | Polygonum macrophyllum | Bistorta macrophylla | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 154 |
| Minuartia verna | Minuartia verna | Sabulina verna | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 118 |
| Potentilla anserina | Potentilla anserina | Argentina anserina | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 111 |
| Helictotrichon versicolor | Helictotrichon versicolor | Helictochloa versicolor | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 108 |
| Minuartia sedoides | Minuartia sedoides | Cherleria sedoides | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 98 |
| Persicaria vivipara | Persicaria vivipara | Bistorta vivipara | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 95 |
| Potentilla stenophylla | Potentilla stenophylla | Argentina stenophylla | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 88 |
| Alopecurus gerardii | Alopecurus gerardii | Alopecurus gerardi | taxonomy decision (prep/config/taxonomy_decisions.csv) | 65 |
| Polygonum runcinatum | Polygonum runcinatum | Persicaria runcinata | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 51 |
| Artemisia flaccida | Artemisia flaccida | Artemisia fulgens | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 48 |
| Festuca arundinacea | Scolochloa festucacea | Lolium arundinaceum | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 37 |
| Polygonum sibiricum | Polygonum sibiricum | Knorringia sibirica | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 37 |
| Sedum integrifolium | Sedum integrifolium | Rhodiola integrifolia | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 35 |
| Aquilegia caerulea | Aquilegia caerulea | Aquilegia coerulea | taxonomy decision (prep/config/taxonomy_decisions.csv) | 31 |
| Cerastium cerastoides | Cerastium cerastoides | Dichodon cerastoides | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 28 |
| Anemone nemorosa | Anemonoides altaica | Anemonoides nemorosa | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 26 |
| Salix serpyllifolia | Salix serpyllifolia | Salix serpillifolia | taxonomy decision (prep/config/taxonomy_decisions.csv) | 25 |
| Oxytropis glacialis | Oxytropis cinerascens | Oxytropis humifusa | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 24 |
| Sibbaldia pentaphylla | Sibbaldia pentaphylla | Potentilla clandestina | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 24 |
| Polygonum cyanandrum | Polygonum cyanandrum | Koenigia cyanandra | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 21 |
| Perideridia montana | Perideridia montana | Perideridia gairdneri | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 19 |
| Koeleria cristata | Rostraria cristata | Koeleria macrantha | taxonomy decision (prep/config/taxonomy_decisions.csv) | 17 |
| Rumex alpestris | Rumex alpestris | Rumex scutatus | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 17 |
| Anemone obtusiloba | Anemonastrum polyanthes | Anemonastrum obtusilobum | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 14 |
| Vaccinium globulare | Vaccinium globulare | Vaccinium membranaceum | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 13 |
| Acinos alpinus | Acinos alpinus | Clinopodium alpinum | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 10 |
| Arenaria kansuensis | Arenaria kansuensis | Eremogone kansuensis | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 10 |
| Astragalus membranaceus | Astragalus trimestris | Astragalus mongholicus | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 8 |
| Minuartia biflora | Minuartia biflora | Cherleria biflora | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 7 |
| Arctostaphylos alpinus | Arctostaphylos uva-ursi | Arctous alpina | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 5 |
| Hieracium hoppeanum | Hieracium hoppeanum | Pilosella hoppeana | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 5 |
| Arenaria roborowskii | Arenaria roborowskii | Eremogone roborowskii | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 4 |
| Hieracium gracile | Hieracium nigrescens subsp. graciliusculum | Hieracium triste | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 4 |
| Vicia angustifolia | Vicia setifolia var. setifolia | Vicia sativa | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 4 |
| Anemone narcissiflora | Anemone narcissiflora | Anemonastrum narcissiflorum | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 3 |
| Osmorhiza chilensis | Osmorhiza chilensis | Osmorhiza berteroi | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 3 |
| Gnaphalium supinum | Gnaphalium supinum | Omalotheca supina | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 2 |
| Primula serratifolia | Primula serratifolia | Primula auricula | legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match | 2 |
| Senecio corymbosus | Kleinia grandiflora | excluded by decision | excluded by decision (prep/config/taxonomy_decisions.csv) | 2 |
| Thlaspi montanum | Lepidium campestre | Noccaea fendleri | taxonomy decision (prep/config/taxonomy_decisions.csv) | 2 |

## Traits

Species with traits now: 674; in step 00: 558. By name, as P2 counted: gained 152, lost 36 (names that changed count as both); species x trait medians changed by more than 5%: 572. After mapping legacy names: same name 522, renamed or merged 23, lost 13, gained 133.

| kind | cause | differences |
| --- | --- | --- |
| correction | P2 cleaning rules (experiments removed per observation, ErrorRisk < 3, non-positive values dropped) | 1135 |
| correction | study species new to the community (P1 taxonomy or FR_AlpeHuez) | 90 |
| correction | records attached by the contributor's own name, not TRY's AccSpeciesName | 60 |
| correction | species renamed or merged by the P1 taxonomy | 57 |
| correction | not a P1 study species | 11 |
| correction | no TRY record left after the P2 rules | 2 |

The 10 largest relative differences of a median:

| species_00 | species_prep | trait | step00 | n_00 | prep | n_prep | relative_difference | cause |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Euphrasia officinalis | Euphrasia officinalis | height | 0.0395 | 22 | 0.232 | 2 | 4.886 | records attached by the contributor's own name, not TRY's AccSpeciesName |
| Phleum alpinum | Phleum alpinum | height | 0.0694 | 9 | 0.292 | 30 | 3.215 | P2 cleaning rules (experiments removed per observation, ErrorRisk < 3, non-positive values dropped) |
| Helianthemum oelandicum subsp. alpestre | Helianthemum oelandicum | height | 0.0315 | 2 | 0.118 | 6 | 2.762 | species renamed or merged by the P1 taxonomy |
| Cerastium fontanum | Cerastium fontanum | height | 0.0683 | 22 | 0.25 | 26 | 2.66 | P2 cleaning rules (experiments removed per observation, ErrorRisk < 3, non-positive values dropped) |
| Galium anisophyllon | Galium anisophyllon | height | 0.0642 | 1 | 0.232 | 2 | 2.617 | P2 cleaning rules (experiments removed per observation, ErrorRisk < 3, non-positive values dropped) |
| Trifolium hybridum | Trifolium hybridum | height | 0.152 | 1 | 0.494 | 11 | 2.25 | P2 cleaning rules (experiments removed per observation, ErrorRisk < 3, non-positive values dropped) |
| Geranium sylvaticum | Geranium sylvaticum | height | 0.115 | 14 | 0.37 | 24 | 2.217 | P2 cleaning rules (experiments removed per observation, ErrorRisk < 3, non-positive values dropped) |
| Pyrola minor | Pyrola minor | height | 0.0316 | 6 | 0.0995 | 22 | 2.149 | P2 cleaning rules (experiments removed per observation, ErrorRisk < 3, non-positive values dropped) |
| Trollius europaeus | Trollius europaeus | height | 0.198 | 11 | 0.528 | 26 | 1.664 | P2 cleaning rules (experiments removed per observation, ErrorRisk < 3, non-positive values dropped) |
| Euphrasia officinalis | Euphrasia officinalis | LMA | 0.065 | 11 | 0.173 | 4 | 1.661 | records attached by the contributor's own name, not TRY's AccSpeciesName |

Legacy species without traits now:

| species_00 | difference | cause |
| --- | --- | --- |
| Argyrotegium poliochlorum | no traits now; P1 species: none | not a P1 study species |
| Astragalus trimestris | no traits now; P1 species: Astragalus mongholicus | no TRY record left after the P2 rules |
| Cheirolophus sempervirens | no traits now; P1 species: none | not a P1 study species |
| Erigeron andicola | no traits now; P1 species: Erigeron andicola | not a P1 study species |
| Euphorbia serpyllifolia | no traits now; P1 species: Euphorbia serpyllifolia | not a P1 study species |
| Hieracium graciliusculum | no traits now; P1 species: Hieracium graciliusculum | not a P1 study species |
| Humulus lupulus | no traits now; P1 species: Humulus lupulus | not a P1 study species |
| Koeleria pyramidata | no traits now; P1 species: Koeleria pyramidata | not a P1 study species |
| Koenigia alpina | no traits now; P1 species: Koenigia alpina | not a P1 study species |
| Salix serpyllifolia | no traits now; P1 species: Salix serpillifolia | no TRY record left after the P2 rules |
| Saussurea hieracioides | no traits now; P1 species: Saussurea hieracioides | not a P1 study species |
| Taraxacum campylodes | no traits now; P1 species: Taraxacum campylodes | not a P1 study species |
| Trichophorum pumilum | no traits now; P1 species: Trichophorum pumilum | not a P1 study species |

## Functional groups

Species: 782 in step 00, 904 in prep. By cause:

| kind | cause | differences |
| --- | --- | --- |
| correction | study species new to the community (P1 taxonomy or FR_AlpeHuez) | 147 |
| correction | species renamed or merged by the P1 taxonomy | 69 |
| correction | woodiness ledger (prep/config/woodiness.csv) | 7 |
| correction | excluded by decision (prep/config/taxonomy_decisions.csv) | 1 |

Species whose graminoid or woody status differs from step 00 (7):

| species_00 | species_prep | graminoid_00 | graminoid | woody_00 | woody | cause |
| --- | --- | --- | --- | --- | --- | --- |
| Acer pseudoplatanus | Acer pseudoplatanus | FALSE | FALSE | FALSE | TRUE | woodiness ledger (prep/config/woodiness.csv) |
| Artemisia tridentata | Artemisia tridentata | FALSE | FALSE | FALSE | TRUE | woodiness ledger (prep/config/woodiness.csv) |
| Fargesia spathacea | Fargesia spathacea | TRUE | TRUE | FALSE | TRUE | woodiness ledger (prep/config/woodiness.csv) |
| Genista sagittalis | Genista sagittalis | FALSE | FALSE | FALSE | TRUE | woodiness ledger (prep/config/woodiness.csv) |
| Larix decidua | Larix decidua | FALSE | FALSE | FALSE | TRUE | woodiness ledger (prep/config/woodiness.csv) |
| Polygaloides chamaebuxus | Polygaloides chamaebuxus | FALSE | FALSE | FALSE | TRUE | woodiness ledger (prep/config/woodiness.csv) |
| Rubus idaeus | Rubus idaeus | FALSE | FALSE | FALSE | TRUE | woodiness ledger (prep/config/woodiness.csv) |

## Phylogeny axes

Species: 526 in step 00 (residents), 904 in prep (both pools). Against step 00 (476 species under the same name): correlation of scores PEV1 0.998, PEV2 0.794, PEV3 0.028; the destination-pool species rotate the axes after the first, so PEV2 and PEV3 are read again in step 10 (Table S6). Eigenvalues: 520 axes in step 00, 903 in prep.

| kind | cause | differences |
| --- | --- | --- |
| correction | axes recomputed on the tree of all 904 study species | 476 |
| correction | destination-pool species, now included (P3 covers both pools) | 242 |
| correction | study species new to the community (P1 taxonomy or FR_AlpeHuez) | 147 |
| correction | species renamed or merged by the P1 taxonomy | 49 |
| correction | excluded by decision (prep/config/taxonomy_decisions.csv) | 1 |

## Experiment climate

41 experiments with the same years agree exactly (largest difference 0). Differences:

| key | difference | cause |
| --- | --- | --- |
| CH_Calanda2::Cal::Nes | years 3 -> 4; T 14.02 -> 13.699; W 3.234 -> 3.234; P_annual 7.248 -> 7.194 | experiment years changed (CH_Calanda2 2020 kept) |
| FR_AlpeHuez::HIGH::LOW | years NA -> 9; T NA -> 14.487; W NA -> 2.381; P_annual NA -> 7.411 | experiment included: FR_AlpeHuez years in the SQLite |

## Niche optima

Step 00 had 785 Lynn names, keyed by Lynn's own names; prep has 767 study species, keyed by the P1 species through Lynn's name. Unchanged (same name, same values): 742. Differences:

| kind | cause | differences |
| --- | --- | --- |
| correction | Lynn name resolves to another accepted species | 25 |
| correction | one of several Lynn names of one species; the one with most occurrences is used | 11 |
| correction | Lynn name is not a P1 study species | 4 |
| correction | Lynn's niche values are NA | 2 |
| correction | Lynn name not resolved by the P1 rule (needs review; community decisions are not applied to Lynn's names) | 1 |

| key | difference | cause |
| --- | --- | --- |
| Trommsdorffia maculata | Trommsdorffia maculata -> Hypochaeris maculata (values unchanged) | Lynn name resolves to another accepted species |
| Mutellina adonidifolia | Mutellina adonidifolia -> Mutellina purpurea (values unchanged) | Lynn name resolves to another accepted species |
| Bellidiastrum michelii | Bellidiastrum michelii -> Aster bellidiastrum (values unchanged) | Lynn name resolves to another accepted species |
| Pachypleurum mutellinoides | Pachypleurum mutellinoides -> Neogaya simplex (values unchanged) | Lynn name resolves to another accepted species |
| Thymus praecox subsp. polytrichus | Thymus praecox subsp. polytrichus -> Thymus praecox (values unchanged) | Lynn name resolves to another accepted species |
| Viola rockiana | not used (P1 species Viola biflora uses Viola biflora) | one of several Lynn names of one species; the one with most occurrences is used |
| Rhinanthus serotinus | Rhinanthus serotinus -> Rhinanthus major (values unchanged) | Lynn name resolves to another accepted species |
| Chamaenerion angustifolium | Chamaenerion angustifolium -> Epilobium angustifolium (values unchanged) | Lynn name resolves to another accepted species |
| Hieracium graciliusculum | not used (P1 species Hieracium graciliusculum uses none) | Lynn's niche values are NA |
| Phlox gracilis | Phlox gracilis -> Microsteris gracilis (values unchanged) | Lynn name resolves to another accepted species |
| Helianthus annuus | not used (P1 species Helianthus annuus uses none) | Lynn name is not a P1 study species |
| Prunella vulgaris subsp. hispida | not used (P1 species Prunella vulgaris uses Prunella vulgaris) | one of several Lynn names of one species; the one with most occurrences is used |
| Gentiana prostrata var. crenulatotruncata | not used (P1 species Gentiana prostrata uses none) | Lynn's niche values are NA |
| Cryptothladia kokonorica | Cryptothladia kokonorica -> Morina kokonorica (values unchanged) | Lynn name resolves to another accepted species |
| Geum longifolium | Geum longifolium -> Coluria longifolia (values unchanged) | Lynn name resolves to another accepted species |
| Stellaria arenarioides | not used (P1 species Stellaria arenarioides uses none) | Lynn name is not a P1 study species |
| Vicia lathyroides | not used (P1 species Vicia lathyroides uses none) | Lynn name is not a P1 study species |
| Leontodon autumnalis | not used (P1 species Scorzoneroides autumnalis uses Scorzoneroides autumnalis) | one of several Lynn names of one species; the one with most occurrences is used |
| Hieracium vulgatum | not used (P1 species Hieracium lachenalii uses Hieracium lachenalii) | one of several Lynn names of one species; the one with most occurrences is used |
| Empetrum hermaphroditum | Empetrum hermaphroditum -> Empetrum nigrum (values unchanged) | Lynn name resolves to another accepted species |
| Solidago virgaurea subsp. minuta | not used (P1 species Solidago virgaurea uses Solidago virgaurea) | one of several Lynn names of one species; the one with most occurrences is used |
| Hieracium piloselloides | not used (P1 species Pilosella piloselloides uses Pilosella piloselloides) | one of several Lynn names of one species; the one with most occurrences is used |
| Agrostis alpina | Agrostis alpina -> Alpagrostis alpina (values unchanged) | Lynn name resolves to another accepted species |
| Arctostaphylos alpinus | Arctostaphylos alpinus -> Arctous alpina (values unchanged) | Lynn name resolves to another accepted species |
| Helianthemum alpestre | Helianthemum alpestre -> Helianthemum oelandicum (values unchanged) | Lynn name resolves to another accepted species |
| Gnaphalium sylvaticum | not used (P1 species Omalotheca sylvatica uses Omalotheca sylvatica) | one of several Lynn names of one species; the one with most occurrences is used |
| Hieracium hoppeanum | Hieracium hoppeanum -> Pilosella hoppeana (values unchanged) | Lynn name resolves to another accepted species |
| Mahonia repens | Mahonia repens -> Berberis repens (values unchanged) | Lynn name resolves to another accepted species |
| Blepharoneuron tricholepis | Blepharoneuron tricholepis -> Muhlenbergia tricholepis (values unchanged) | Lynn name resolves to another accepted species |
| Dysphania graveolens | Dysphania graveolens -> Dysphania incisa (values unchanged) | Lynn name resolves to another accepted species |
| Pleurospermum hookeri | Pleurospermum hookeri -> Hymenidium hookeri (values unchanged) | Lynn name resolves to another accepted species |
| Oxytropis glacialis | Oxytropis glacialis -> Oxytropis humifusa (values unchanged) | Lynn name resolves to another accepted species |
| Ranunculus brotherusii var. tanguticus | not used (P1 species Ranunculus brotherusii uses Ranunculus brotherusii) | one of several Lynn names of one species; the one with most occurrences is used |
| Polygonum runcinatum | Polygonum runcinatum -> Persicaria runcinata (values unchanged) | Lynn name resolves to another accepted species |
| Equisetum pratense | not used (P1 species Equisetum pratense uses none) | Lynn name is not a P1 study species |
| Gnaphalium norvegicum | not used (P1 species Omalotheca norvegica uses Omalotheca norvegica) | one of several Lynn names of one species; the one with most occurrences is used |
| Alopecurus gerardii | Alopecurus gerardii -> Alopecurus gerardi (values unchanged) | Lynn name resolves to another accepted species |
| Centaurea pseudophrygia | Centaurea pseudophrygia -> Centaurea phrygia (values unchanged) | Lynn name resolves to another accepted species |
| Anemone nemorosa | Anemone nemorosa -> Anemonoides nemorosa (values unchanged) | Lynn name resolves to another accepted species |
| Ranunculus acris subsp. friesianus | not used (P1 species Ranunculus acris uses Ranunculus acris) | one of several Lynn names of one species; the one with most occurrences is used |
| Conyza sumatrensis | Conyza sumatrensis -> Erigeron sumatrensis (values unchanged) | Lynn name resolves to another accepted species |
| Koeleria cristata | not resolved | Lynn name not resolved by the P1 rule (needs review; community decisions are not applied to Lynn's names) |
| Noccaea brachypetala | not used (P1 species Noccaea caerulescens uses Noccaea caerulescens) | one of several Lynn names of one species; the one with most occurrences is used |

## Found while comparing

Lynn names that are decided community names (`prep/config/taxonomy_decisions.csv`) but resolve otherwise in P4, because community decisions are not applied to Lynn's names (1):

| lynn_name | regions | n_obs_sub | lynn_resolves_to | decision_species | decision_species_has_niche |
| --- | --- | --- | --- | --- | --- |
| Koeleria cristata | CN_Heibei | 443 | unresolved | Koeleria macrantha | FALSE |

Lynn's names come from the same community records, so the decision applies to them as well; changing this is a P4 change (open question for Bills).

## Columns

Columns dropped from step 00 that a later step (MILESTONES steps 01 to 10, PLAN.md) names:

None.

Other dropped columns (named by no later step): community `turf_id`; community `supplied_relative_cover`; community `elevation`; community `stable_plot_id`; community `year_established`; experiments `origin_longitude`; experiments `origin_latitude`; experiments `destination_longitude`; experiments `destination_latitude`; species_taxonomy `submitted_name`; species_taxonomy `matched_name`; species_taxonomy `score`; species_taxonomy `name_normalized`; species_taxonomy `accepted_normalized`; species_taxonomy `decision`; species_taxonomy `accepted_name_override`; species_taxonomy `family`; species_taxonomy `note`; species_taxonomy `accepted_name_reviewed`; species_taxonomy `needs_review`; species_taxonomy `original_word_count`; species_taxonomy `accepted_word_count`; species_taxonomy `explicit_unspecified_species`; species_taxonomy `family_name`; species_taxonomy `include_species_level`; species_taxonomy `resolution_decision`; functional_groups `functional_group`.

## Counts to confirm for steps 01 to 10

| count | expected | prep |
| --- | --- | --- |
| included experiments |  | 43 |
| included experiments without US_Arizona | 42 | 42 |
| gradients of the included experiments |  | 21 |
| gradients without US_Arizona | 20 | 20 |

The gradient count uses the SQLite's `gradient` per destination site; step 01 defines `gradient_id`. The number of effects (about 1,596) is confirmed in step 02.

## Switch-over

`grep -rnE "_targets|legacy/" analysis/` outside comments: no hits (`analysis/` holds no script yet). Steps 01 to 10 run the same check on their own script.
Code in `prep/` (outside this script) that reads step 00 files, the legacy store or `legacy/`: none.
The P3 placement draft reuses legacy labels through `prep/config/legacy_taxonomy.csv` (step 00's submitted name and legacy accepted name, copied by this script). `legacy/` is read by nothing; prep comments cite its lines, and it stays until step 15 is closed (also on the tag `targets-pipeline-2026-10`).

## Warnings captured

None.

## Deviations from the plan

(1) The comparison is a one-off script (`prep/compare_to_frozen.R`), not a target: it needs the step 00 files, which are archived after this run. (2) The approved P1 to P4 reports no longer carry their step 00 sections (decided 27 Oct); their content is in this report.
