# P1: Community and taxonomy

The prep pipeline read the network SQLite (64876 records, 22 regions) and the local-code dictionary of the RData, resolved every community and TRY name once with tidyTRY 0.2.1 (TNRS 0.3.6, 2026-10-28; legitimate names only, WCVP before WFO; the answers for 960 study and 33908 TRY names are read from `prep/tnrs/`, not from a new TNRS call), and wrote the community table with one row per plot-year x accepted species for the analysis and every other vegetation record kept with its `drop_reason`. No relative cover is computed here: step 01 computes it over all vegetation records of a plot-year. 44 warmed comparisons, 43 included.

## Counts

| step | rows |
| --- | --- |
| SQLite community rows | 64876 |
| dropped: establishment year | 6053 |
| dropped: year edit | 0 |
| rows kept (all vegetation records) | 58823 |
|   of which set aside: ambiguous code excluded | 14 |
|   of which set aside: community edit C1 | 56 |
|   of which set aside: community edit C2 | 21 |
|   of which set aside: community edit C3 | 7 |
|   of which set aside: community edit C4 | 4 |
|   of which set aside: community edit C5 | 2 |
|   of which set aside: excluded by decision | 2 |
|   of which set aside: family_only | 24 |
|   of which set aside: genus_only | 3339 |
|   of which set aside: unidentified | 354 |
| analysis records before collapsing | 55000 |
| analysis rows (plot-year x species) | 54980 |

## Reconciliation with the RData

Rows per region: SQLite 64876, RData `dat` 61430. Region-years where the two differ:

| region | year | rows_sqlite | rows_rdata | difference |
| --- | --- | --- | --- | --- |
| FR_AlpeHuez | 2015 | 465 | 0 | 465 |
| FR_AlpeHuez | 2016 | 564 | 0 | 564 |
| FR_AlpeHuez | 2017 | 553 | 0 | 553 |
| FR_AlpeHuez | 2018 | 528 | 0 | 528 |
| FR_AlpeHuez | 2019 | 529 | 0 | 529 |
| FR_AlpeHuez | 2020 | 563 | 0 | 563 |
| FR_AlpeHuez | 2021 | 560 | 0 | 560 |
| FR_AlpeHuez | 2022 | 560 | 223 | 337 |
| FR_AlpeHuez | 2023 | 549 | 0 | 549 |
| FR_AlpeHuez |  | 0 | 1202 | -1202 |

## How names were decided

Every community name gets one of three origins (`decided_by`) and one identification level.

- **rule**: chosen automatically by `tidyTRY::resolve_species()` (v 0.2.1 ) from all TNRS candidates: only legitimate names (taxonomic status Accepted or Synonym, never Illegitimate, Invalid or No opinion); the highest match score; WCVP before WFO; Accepted before Synonym. A match into another genus is never chosen, and a name is flagged for review (not chosen automatically) when it comes only from WFO, matches with a score below 0.9, ties between different accepted names, or has no legitimate candidate.
- **decision**: set by Bills in `prep/config/taxonomy_decisions.csv`, for every community name the rule flagged and for names the rule resolved to a plant from the wrong region (with the reason in the ledger).
- **unidentified**: an unidentified field morphotaxon (`prep/config/unidentified_taxa.csv`, from the legacy overrides); its records stay as vegetation but never enter the analysis.

Identification level (legacy rule, `legacy/R/data/community.R:252-285`): **species_or_lower** when the recorded and the accepted name both have at least two words; **genus_only** when either has one word or the recorded name ends in sp./spp./spec./indet.; **family_only** when the name ends in -aceae and has one word. Only species_or_lower names enter the analysis; the others keep their records with a `drop_reason`.

After the rule and the decisions, these accepted names are replaced in every table (`prep/config/accepted_renames.csv`): Taraxacum sect. Taraxacum -> Taraxacum officinale; Alopecurus gerardii -> Alopecurus gerardi.

**needs_review** in the table below counts the names the rule could not settle on its own (the reasons above). For community names every one of them is decided in the ledger before the pipeline continues; for TRY names they are counted here and handled where those names are used (P2).

## Name resolution

TRY names are resolved here, in the same run as the community names, because traits (P2) are attached to a study species through the TRY record's own name, not through TRY's `AccSpeciesName`, which merges species (TRY files Festuca pratensis records under Festuca ovina). One run, one rule and one TNRS date serve every join.

TRY names: 149997 distinct names in the export; 33908 in a genus of a study name or of its accepted name were resolved (the rest cannot attach to a study species unless filed under a third genus).

| names | submitted | resolved | needs_review |
| --- | --- | --- | --- |
| community | 960 | 958 | 9 |
| TRY | 33908 | 33604 | 4082 |

Community names by identification level and how they were decided:

| identification_level | decided_by | names |
| --- | --- | --- |
| excluded by decision | decision | 1 |
| family_only | rule | 9 |
| genus_only | rule | 122 |
| species_or_lower | decision | 13 |
| species_or_lower | rule | 946 |
| unidentified | unidentified | 35 |

TRY names flagged for review: 4082. Those whose accepted name, or one of the tied names, is a study species could carry trait records for the analysis (`could_be_study`); their handling is decided in P2.

| reason | could_be_study | names | try_rows |
| --- | --- | --- | --- |
| fallback source wfo | TRUE | 164 | 8588 |
| fallback source wfo; fuzzy match | TRUE | 30 | 484 |
| fuzzy match | TRUE | 377 | 7268 |
| tie | TRUE | 13 | 1097 |
| fallback source wfo | FALSE | 1068 | 54623 |
| fallback source wfo; fuzzy match | FALSE | 212 | 1471 |
| fuzzy match | FALSE | 1927 | 54751 |
| no legitimate candidate | FALSE | 255 | 31507 |
| tie | FALSE | 36 | 2231 |

The comparison with the frozen legacy taxonomy of step 00 (P5, 27 Oct, before the GBIF niche and the trait screening) is kept in the git history: `reports/P5_compare_to_frozen.md` at commit 04e4b70.

### Names merged into one species (106 names in 50 species; full table in `data/derived/prep/species_merges.csv`)

| merge_type | names |
| --- | --- |
| accepted name | 44 |
| infraspecific to species | 20 |
| spelling | 3 |
| synonym | 39 |

| accepted_species | name_original | accepted_name | merge_type | regions | n_records |
| --- | --- | --- | --- | --- | --- |
| Anthoxanthum nipponicum | Anthoxanthum alpinum | Anthoxanthum nipponicum | synonym | CH_Calanda | 1 |
| Anthoxanthum nipponicum | Anthoxanthum odoratum subsp. nipponicum | Anthoxanthum nipponicum | synonym | CH_Calanda; CH_Lavey; DE_TransAlps | 330 |
| Anthyllis vulneraria | Anthyllis vulneraria | Anthyllis vulneraria | accepted name | CH_Calanda; CH_Calanda2; CH_Lavey; DE_Grainau; FR_Lautaret | 258 |
| Anthyllis vulneraria | Anthyllis vulneraria var. alpestris | Anthyllis vulneraria subsp. alpestris | infraspecific to species | CH_Lavey | 37 |
| Arctous alpina | Arctostaphylos alpinus | Arctous alpina | synonym | CH_Lavey | 5 |
| Arctous alpina | Arctous alpinus | Arctous alpina | spelling | NO_Lavisdalen | 2 |
| Arenaria ciliata | Arenaria ciliata | Arenaria ciliata | accepted name | CH_Calanda | 37 |
| Arenaria ciliata | Arenaria multicaulis | Arenaria ciliata | synonym | CH_Calanda | 5 |
| Avenella flexuosa | Avenella flexuosa | Avenella flexuosa | accepted name | CH_Lavey; FR_AlpeHuez; FR_Lautaret; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen | 408 |
| Avenella flexuosa | Deschampsia flexuosa | Avenella flexuosa | synonym | DE_Susalps; DE_TransAlps; SE_Abisko | 27 |
| Avenula pubescens | Avenula pubescens | Avenula pubescens | accepted name | FR_Lautaret; NO_Ulvhaugen | 19 |
| Avenula pubescens | Helictotrichon pubescens | Avenula pubescens | synonym | CH_Calanda; DE_Grainau; DE_Susalps; DE_TransAlps | 173 |
| Avenula pubescens | Homalotrichon pubescens | Avenula pubescens | synonym | IT_MatschMazia1; IT_MatschMazia2 | 7 |
| Bistorta vivipara | Bistorta vivipara | Bistorta vivipara | accepted name | FR_Lautaret; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen; SE_Abisko | 517 |
| Bistorta vivipara | Persicaria vivipara | Bistorta vivipara | synonym | CN_Gongga; IT_MatschMazia2 | 95 |
| Bistorta vivipara | Polygonum viviparum | Bistorta vivipara | synonym | CH_Calanda; CH_Calanda2; CH_Lavey; CN_Damxung; CN_Heibei | 507 |
| Bromus erectus | Bromopsis erecta | Bromus erectus | synonym | FR_Lautaret | 7 |
| Bromus erectus | Bromus erectus | Bromus erectus | accepted name | CH_Calanda; CH_Lavey | 41 |
| Carex myosuroides | Carex myosuroides | Carex myosuroides | accepted name | FR_Lautaret | 3 |
| Carex myosuroides | Elyna myosuroides | Carex myosuroides | synonym | CH_Calanda; CH_Lavey | 103 |
| Cerastium holosteoides | Cerastium fontanum subsp. vulgare | Cerastium holosteoides | synonym | CH_Lavey | 3 |
| Cerastium holosteoides | Cerastium holosteoides | Cerastium holosteoides | accepted name | DE_Grainau; DE_Susalps; DE_TransAlps; IT_MatschMazia1; IT_MatschMazia2 | 197 |
| Clinopodium alpinum | Acinos alpinus | Clinopodium alpinum | synonym | CH_Calanda; CH_Calanda2 | 10 |
| Clinopodium alpinum | Ziziphora granatensis | Clinopodium alpinum subsp. meridionale | infraspecific to species | FR_AlpeHuez | 23 |
| Comastoma tenellum | Comastoma tenellum | Comastoma tenellum | accepted name | NO_Gudmedalen; NO_Ulvhaugen | 3 |
| Comastoma tenellum | Gentianella tenella | Comastoma tenellum | synonym | CH_Calanda | 86 |
| Crepis pontana | Crepis bocconei | Crepis pontana | synonym | FR_Lautaret | 2 |
| Crepis pontana | Crepis pontana | Crepis pontana | accepted name | FR_Lautaret | 6 |
| Crocus vernus | Crocus albiflorus | Crocus vernus | synonym | CH_Calanda; CH_Calanda2; CH_Lavey; DE_Grainau; IT_MatschMazia1; IT_MatschMazia2 | 146 |
| Crocus vernus | Crocus vernus | Crocus vernus | accepted name | DE_Susalps; DE_TransAlps; FR_AlpeHuez; FR_Lautaret | 105 |
| Deschampsia cespitosa | Deschampsia alpina | Deschampsia cespitosa subsp. cespitosa | infraspecific to species | NO_Lavisdalen | 1 |
| Deschampsia cespitosa | Deschampsia caespitosa | Deschampsia cespitosa | spelling | CN_Heibei | 2 |
| Deschampsia cespitosa | Deschampsia cespitosa | Deschampsia cespitosa | accepted name | CH_Lavey; DE_Grainau; DE_TransAlps; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen | 377 |
| Euphrasia officinalis | Euphrasia officinalis | Euphrasia officinalis | accepted name | DE_Grainau | 3 |
| Euphrasia officinalis | Euphrasia rostkoviana | Euphrasia officinalis subsp. pratensis | infraspecific to species | CH_Lavey | 21 |
| Euphrasia officinalis | Euphrasia rostkoviana subsp. montana | Euphrasia officinalis subsp. monticola | infraspecific to species | CH_Calanda; CH_Calanda2 | 72 |
| Festuca violacea | Festuca melanopsis | Festuca violacea subsp. violacea | infraspecific to species | CH_Lavey | 4 |
| Festuca violacea | Festuca violacea | Festuca violacea | accepted name | CH_Lavey; FR_Lautaret | 101 |
| Gentiana prostrata | Gentiana crenulatotruncata | Gentiana prostrata var. crenulatotruncata | infraspecific to species | CN_Heibei | 19 |
| Gentiana prostrata | Gentiana pudica | Gentiana prostrata var. pudica | infraspecific to species | CN_Heibei | 9 |
| Gentianella campestris | Gentiana campestris | Gentianella campestris | synonym | CH_Calanda; CH_Calanda2; CH_Lavey | 154 |
| Gentianella campestris | Gentianella campestris | Gentianella campestris | accepted name | FR_Lautaret; NO_Ulvhaugen | 4 |
| Helictochloa pratensis | Avenula pratensis | Helictochloa pratensis | synonym | FR_AlpeHuez | 55 |
| Helictochloa pratensis | Helictochloa pratensis | Helictochloa pratensis | accepted name | FR_AlpeHuez | 5 |
| Helictochloa pratensis | Helictotrichon alpinum | Helictochloa pratensis subsp. pratensis | infraspecific to species | CH_Calanda | 8 |
| Helictochloa versicolor | Helictochloa versicolor | Helictochloa versicolor | accepted name | FR_AlpeHuez; FR_Lautaret | 54 |
| Helictochloa versicolor | Helictotrichon versicolor | Helictochloa versicolor | synonym | CH_Lavey; DE_TransAlps | 108 |
| Hieracium lachenalii | Hieracium lachenalii | Hieracium lachenalii | accepted name | CH_Lavey | 4 |
| Hieracium lachenalii | Hieracium vulgatum | Hieracium lachenalii subsp. lachenalii | infraspecific to species | NO_Gudmedalen; NO_Skjellingahaugen; NO_Ulvhaugen | 19 |
| Hieracium piliferum | Hieracium glanduliferum | Hieracium piliferum subsp. glanduliferum | infraspecific to species | FR_Lautaret | 6 |
| Hieracium piliferum | Hieracium piliferum | Hieracium piliferum | accepted name | FR_Lautaret | 1 |
| Hypericum maculatum | Hypericum maculatum | Hypericum maculatum | accepted name | CH_Calanda; CH_Lavey; DE_Grainau; FR_Lautaret; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen | 94 |
| Hypericum maculatum | Hypericum maculatum maculatum | Hypericum maculatum subsp. maculatum | infraspecific to species | FR_Lautaret | 4 |
| Leontodon hispidus | Leontodon hispidus | Leontodon hispidus | accepted name | CH_Calanda; CH_Lavey; DE_Grainau; DE_Susalps; DE_TransAlps; FR_AlpeHuez; FR_Lautaret; IT_MatschMazia1; IT_MatschMazia2 | 365 |
| Leontodon hispidus | Leontodon hispidus subsp. hispidus | Leontodon hispidus subsp. hispidus | infraspecific to species | CH_Lavey | 32 |
| Mutellina purpurea | Ligusticum mutellina | Mutellina purpurea | synonym | CH_Calanda2; CH_Lavey | 138 |
| Mutellina purpurea | Mutellina adonidifolia | Mutellina purpurea | synonym | DE_TransAlps | 49 |
| Neogaya simplex | Ligusticum mutellinoides | Neogaya simplex | synonym | CH_Calanda; CH_Lavey | 125 |
| Neogaya simplex | Pachypleurum mutellinoides | Neogaya simplex | synonym | FR_Lautaret | 1 |
| Noccaea caerulescens | Noccaea brachypetala | Noccaea caerulescens subsp. brachypetala | infraspecific to species | FR_AlpeHuez | 44 |
| Noccaea caerulescens | Noccaea caerulescens | Noccaea caerulescens | accepted name | NO_Lavisdalen; NO_Ulvhaugen | 51 |
| Noccaea fendleri | Noccaea fendleri | Noccaea fendleri | accepted name | US_Colorado | 41 |
| Noccaea fendleri | Thlaspi montanum | Noccaea fendleri | synonym | US_Montana | 2 |
| Omalotheca norvegica | Gnaphalium norvegicum | Omalotheca norvegica | synonym | SE_Abisko | 1 |
| Omalotheca norvegica | Omalotheca norvegica | Omalotheca norvegica | accepted name | NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen | 7 |
| Omalotheca supina | Gnaphalium supinum | Omalotheca supina | synonym | DE_TransAlps | 2 |
| Omalotheca supina | Omalotheca supina | Omalotheca supina | accepted name | NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen | 76 |
| Omalotheca sylvatica | Gnaphalium sylvaticum | Omalotheca sylvatica | synonym | CH_Lavey | 4 |
| Omalotheca sylvatica | Omalotheca sylvatica | Omalotheca sylvatica | accepted name | NO_Skjellingahaugen | 3 |
| Osmorhiza berteroi | Osmorhiza berteroi | Osmorhiza berteroi | accepted name | US_Montana | 4 |
| Osmorhiza berteroi | Osmorhiza chilensis | Osmorhiza berteroi | synonym | US_Montana | 3 |
| Oxalis acetosella | Oxalis acetocella | Oxalis acetosella | spelling | IN_Kashmir | 6 |
| Oxalis acetosella | Oxalis acetosella | Oxalis acetosella | accepted name | NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen | 25 |
| Phleum alpinum | Phleum alpinum | Phleum alpinum | accepted name | CH_Calanda; FR_Lautaret; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen; SE_Abisko | 157 |
| Phleum alpinum | Phleum alpinum subsp. rhaeticum | Phleum alpinum | synonym | CH_Lavey | 4 |
| Phleum alpinum | Phleum rhaeticum | Phleum alpinum | synonym | CH_Lavey; DE_Grainau; IT_MatschMazia2 | 63 |
| Pilosella lactucella | Hieracium lactucella | Pilosella lactucella | synonym | CH_Calanda; CH_Calanda2; CH_Lavey | 121 |
| Pilosella lactucella | Pilosella lactucella | Pilosella lactucella | accepted name | CH_Calanda | 1 |
| Pilosella piloselloides | Hieracium piloselloides | Pilosella piloselloides | synonym | CH_Lavey | 52 |
| Pilosella piloselloides | Pilosella piloselloides | Pilosella piloselloides | accepted name | DE_Susalps; DE_TransAlps | 100 |
| Plantago maritima | Plantago maritima | Plantago maritima | accepted name | FR_AlpeHuez | 109 |
| Plantago maritima | Plantago maritima subsp. serpentina | Plantago maritima subsp. maritima | infraspecific to species | FR_Lautaret | 39 |
| Prunella vulgaris | Prunella hispida | Prunella vulgaris subsp. hispida | infraspecific to species | CN_Gongga | 15 |
| Prunella vulgaris | Prunella vulgaris | Prunella vulgaris | accepted name | CH_Calanda; CH_Calanda2; CH_Lavey; DE_Grainau; DE_Susalps; DE_TransAlps; FR_AlpeHuez; NO_Gudmedalen; NO_Skjellingahaugen; NO_Ulvhaugen; US_Montana | 338 |
| Pulsatilla alpina | Anemone alpina | Pulsatilla alpina | synonym | FR_Lautaret | 1 |
| Pulsatilla alpina | Pulsatilla alpina | Pulsatilla alpina | accepted name | CH_Lavey; FR_AlpeHuez | 64 |
| Ranunculus acris | Ranunculus acris | Ranunculus acris | accepted name | CH_Calanda; CH_Calanda2; DE_Grainau; DE_Susalps; DE_TransAlps; FR_AlpeHuez; IT_MatschMazia1; IT_MatschMazia2; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen; SE_Abisko | 911 |
| Ranunculus acris | Ranunculus acris subsp. friesianus | Ranunculus acris subsp. friesianus | infraspecific to species | CH_Lavey | 8 |
| Ranunculus brotherusii | Ranunculus brotherusii | Ranunculus brotherusii | accepted name | CN_Heibei | 36 |
| Ranunculus brotherusii | Ranunculus brotherusii var. tanguticus | Ranunculus brotherusii var. tanguticus | infraspecific to species | CN_Gongga | 55 |
| Scorzoneroides autumnalis | Leontodon autumnalis | Scorzoneroides autumnalis | synonym | CH_Lavey; DE_Grainau; DE_Susalps; DE_TransAlps; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen | 301 |
| Scorzoneroides autumnalis | Scorzoneroides autumnalis | Scorzoneroides autumnalis | accepted name | IT_MatschMazia2 | 2 |
| Scorzoneroides helvetica | Leontodon helveticus | Scorzoneroides helvetica | synonym | CH_Lavey | 9 |
| Scorzoneroides helvetica | Scorzoneroides helvetica | Scorzoneroides helvetica | accepted name | DE_Susalps; DE_TransAlps | 138 |
| Scorzoneroides pyrenaica | Leontodon pyrenaicus | Scorzoneroides pyrenaica | synonym | FR_Lautaret | 1 |
| Scorzoneroides pyrenaica | Scorzoneroides pyrenaica | Scorzoneroides pyrenaica | accepted name | FR_AlpeHuez; FR_Lautaret | 149 |
| Silene flos-cuculi | Lychnis flos-cuculi | Silene flos-cuculi | synonym | DE_Susalps | 6 |
| Silene flos-cuculi | Silene flos-cuculi | Silene flos-cuculi | accepted name | DE_Susalps | 1 |
| Solidago virgaurea | Solidago virgaurea | Solidago virgaurea | accepted name | CH_Lavey; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; SE_Abisko | 72 |
| Solidago virgaurea | Solidago virgaurea subsp. minuta | Solidago virgaurea subsp. minuta | infraspecific to species | CH_Lavey | 19 |
| Symphyotrichum ascendens | Aster adscendens | Symphyotrichum ascendens | synonym | US_Montana | 1 |
| Symphyotrichum ascendens | Symphyotrichum ascendens | Symphyotrichum ascendens | accepted name | US_Colorado; US_Montana | 11 |
| Vaccinium membranaceum | Vaccinium globulare | Vaccinium membranaceum | synonym | US_Montana | 13 |
| Vaccinium membranaceum | Vaccinium membranaceum | Vaccinium membranaceum | accepted name | US_Montana | 3 |
| Viola biflora | Viola biflora | Viola biflora | accepted name | CN_Heibei; DE_TransAlps; NO_Gudmedalen; NO_Lavisdalen; NO_Skjellingahaugen; NO_Ulvhaugen; SE_Abisko | 318 |
| Viola biflora | Viola biflora var. rockiana | Viola biflora var. rockiana | infraspecific to species | CN_Gongga | 130 |

## Community edits (legacy hard codes, kept and documented)

| edit_id | region | field | value | rows | legacy_source |
| --- | --- | --- | --- | --- | --- |
| C1 | CN_Damxung | name_original | Stipa capillacea | 56 | legacy/R/data/community.R:210-223 |
| C2 | CN_Damxung | name_original | Stipa capillata | 21 | legacy/R/data/community.R:210-223 |
| C3 | CN_Damxung | local_species_code | Polygonaceae | 7 | legacy/R/data/community.R:210-223 |
| C4 | DE_Susalps | name_original | Festuca pratensis | 4 | legacy/R/data/community.R:210-223 |
| C5 | DE_Susalps | name_original | Elymus repens | 2 | legacy/R/data/community.R:210-223 |

## Codes not in the dictionary (recorded name kept)

| region | destination_site | local_species_code | rows |
| --- | --- | --- | --- |
| FR_Lautaret | L | Pilosella | 6 |

## Experiments

| experiment_id | first_year | last_year | year_range | included | exclusion_reason |
| --- | --- | --- | --- | --- | --- |
| CH_Calanda::Cal::Nes | 2013 | 2018 | 6 | TRUE |  |
| CH_Calanda::Pea::Cal | 2013 | 2018 | 6 | TRUE |  |
| CH_Calanda2::Cal::Nes | 2017 | 2020 | 4 | TRUE |  |
| CH_Lavey::RIO::CRE | 2017 | 2019 | 3 | TRUE |  |
| CH_Lavey::RIO::MAR | 2017 | 2019 | 3 | TRUE |  |
| CH_Lavey::RIO::PRA | 2017 | 2019 | 3 | TRUE |  |
| CN_Damxung::HIGH::LOW | 2014 | 2018 | 5 | TRUE |  |
| CN_Gongga::A::M | 2013 | 2016 | 4 | TRUE |  |
| CN_Gongga::H::A | 2013 | 2016 | 4 | TRUE |  |
| CN_Gongga::H::L | 2013 | 2016 | 4 | TRUE |  |
| CN_Gongga::M::L | 2013 | 2016 | 4 | TRUE |  |
| CN_Heibei::3400::3200 | 2008 | 2016 | 9 | TRUE |  |
| CN_Heibei::3800::3200 | 2008 | 2016 | 9 | TRUE |  |
| CN_Heibei::3800::3400 | 2008 | 2016 | 9 | TRUE |  |
| DE_Grainau::HIGH::LOW | 2014 | 2017 | 4 | TRUE |  |
| DE_Susalps::EB::BT | 2017 | 2020 | 4 | TRUE |  |
| DE_Susalps::EB::FE | 2017 | 2020 | 4 | TRUE |  |
| DE_Susalps::EB::GW | 2017 | 2020 | 4 | TRUE |  |
| DE_Susalps::FE::BT | 2017 | 2020 | 4 | TRUE |  |
| DE_Susalps::GW::BT | 2017 | 2020 | 4 | TRUE |  |
| DE_Susalps::GW::FE | 2017 | 2020 | 4 | TRUE |  |
| DE_TransAlps::FP::BT | 2017 | 2020 | 4 | TRUE |  |
| DE_TransAlps::SP::BT | 2017 | 2020 | 4 | TRUE |  |
| FR_AlpeHuez::HIGH::LOW | 2015 | 2023 | 9 | TRUE |  |
| FR_Lautaret::G::L | 2017 | 2022 | 6 | TRUE |  |
| IN_Kashmir::HIGH::LOW | 2014 | 2015 | 2 | FALSE | insufficient duration |
| IT_MatschMazia1::High::Low | 2011 | 2013 | 3 | TRUE |  |
| IT_MatschMazia2::High::Low | 2011 | 2013 | 3 | TRUE |  |
| NO_Gudmedalen::Gudmedalen::Rambera | 2011 | 2017 | 7 | TRUE |  |
| NO_Gudmedalen::Rambera::Arhelleren | 2011 | 2017 | 7 | TRUE |  |
| NO_Lavisdalen::Hogsete::Vikesland | 2011 | 2017 | 7 | TRUE |  |
| NO_Lavisdalen::Lavisdalen::Hogsete | 2011 | 2017 | 7 | TRUE |  |
| NO_Skjellingahaugen::Skjellingahaugen::Veskre | 2011 | 2017 | 7 | TRUE |  |
| NO_Skjellingahaugen::Veskre::Ovstedal | 2011 | 2017 | 7 | TRUE |  |
| NO_Ulvhaugen::Alrust::Fauske | 2011 | 2017 | 7 | TRUE |  |
| NO_Ulvhaugen::Ulvhaugen::Alrust | 2011 | 2017 | 7 | TRUE |  |
| SE_Abisko::High::Low | 2013 | 2015 | 3 | TRUE |  |
| SE_Abisko::High::Mid | 2013 | 2015 | 3 | TRUE |  |
| SE_Abisko::Mid::Low | 2013 | 2015 | 3 | TRUE |  |
| US_Arizona::MC::PP | 2015 | 2018 | 4 | TRUE |  |
| US_Colorado::mo::pf | 2018 | 2023 | 6 | TRUE |  |
| US_Colorado::mo::um | 2018 | 2023 | 6 | TRUE |  |
| US_Colorado::pf::um | 2018 | 2023 | 6 | TRUE |  |
| US_Montana::High::Low | 2014 | 2018 | 5 | TRUE |  |

## Warnings captured

- TRY names: TNRS returned 40 candidate rows of an unrelated name (another genus with a perfect name score); dropped. ℹ Affected: "Pimpinella peregrina", "Euphorbia ruficeps", "Myosotis scorpioides   agg.", "Vaccinium vitis-idaea L. Orthodox", "Rosa chinensis Jacq.", "KALMIA   ANGUSTIFOLIA L.", "ROSA MICRANTHA Borrer ex Sm.", "Arabis canadensis L.",   "Senecio piptocoma O.Hoffm.", "Vicia lutea L. subsp. vestita (Boiss.) Rouy",   "Atriplex leucoclada Boiss.", "Symphyotrichum ciliatum", "Trifolium pratense   subsp. nivale", "Senecio latifolius", "Castilleja arctica", "Centaurea   antennata", "Rheum palaestinum", "Urtica bianorii", and "Carex planostachys".

## Deviations from the plan

(1) The TNRS answers are saved in `prep/tnrs/` (one row per name with its date and the tidyTRY and TNRS versions; decided 30 Oct) and the pipeline never calls TNRS, because TNRS can answer differently between runs. A name missing from the saved answers stops the pipeline; `prep/update_tnrs_cache.R` resolves only the new names and never resends a saved one. (2) TRY names are resolved in a second call limited to the genera of the study names and their accepted names, because the TRY export lists about 67,000 species; TRY records filed under a third genus are not attached. (3) This run used tidyTRY 0.2.1 installed from the branch `fix/tnrs-id-check` (TNRS rows checked against their submitted name, rows of unrelated names dropped, no match into another genus, retries); `renv.lock` pins the release and moves to it once it is tagged. (4) FR_AlpeHuez::HIGH::LOW is included: the SQLite has its years 2015 to 2023, so the legacy reason "warmed observations have missing Year" no longer applies (decision at review).
