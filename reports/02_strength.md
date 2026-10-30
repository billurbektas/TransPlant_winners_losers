# Step 02: Strength (Step A, residents)

One `lmer(relative_cover ~ treatment * year * species + (1 | plot_id:species))` was fitted per experiment on the resident species' rows in origin controls and warmed plots (108,884 rows, 43 experiments), year as a factor. `emmeans` gave the warmed minus origin-control contrast per species and survey year (`by_year`, for step 03) and averaged over each experiment's years (`strength`), with asymptotic degrees of freedom in every experiment. Each effect carries the experiment's W, T and P (summer, log) standardised across the 43 experiments, TP = T x P, A0 at origin and destination standardised across the 1,655 effects, and the species properties standardised across the residents that have them. Strength is in relative-cover units (proportion of plot cover).

## Counts

| rule | n |
| --- | --- |
| resident rows, origin controls and warmed plots (step 01) | 108,884 |
| experiments fitted | 43 |
| gradients | 21 |
| species | 528 |
| species x experiment (effects) | 1,655 |
| by_year contrasts (species x survey year) | 7,571 |
| strength_df rows | 1,655 |

## Strength

Quantiles (min / 25% / median / 75% / max). strength: -0.2965 / -0.0090 / -0.0016 / 0.0020 / 0.1885; strength_se: 0.0070 / 0.0091 / 0.0138 / 0.0191 / 0.0667. Effects with the contrast's 95% interval excluding zero: 232 of 1,655.

![strength and strength_se](../figures/diagnostics/02_strength.png)

Within an experiment the SE does not differ between species: the model has one residual variance and the plots are balanced, so `distinct_se` (distinct SEs to 6 significant digits) is 1 in 43 of 43 experiments. V in the later `rma.mv` calls therefore weights experiments rather than species within an experiment.

## Fits per experiment

| experiment_id | rows | species | years | singular | distinct_se | warnings |
| --- | --- | --- | --- | --- | --- | --- |
| CH_Calanda::Cal::Nes | 7739 | 71 | 6 | FALSE | 1 | 1 |
| CH_Calanda::Pea::Cal | 4080 | 34 | 6 | FALSE | 1 | 1 |
| CH_Calanda2::Cal::Nes | 4240 | 53 | 4 | FALSE | 1 | 1 |
| CH_Lavey::RIO::CRE | 3780 | 63 | 3 | FALSE | 1 | 1 |
| CH_Lavey::RIO::MAR | 3780 | 63 | 3 | FALSE | 1 | 1 |
| CH_Lavey::RIO::PRA | 3717 | 63 | 3 | FALSE | 1 | 1 |
| CN_Damxung::HIGH::LOW | 2560 | 32 | 4 | TRUE | 1 | 2 |
| CN_Gongga::A::M | 1848 | 33 | 4 | FALSE | 1 | 1 |
| CN_Gongga::H::A | 1800 | 36 | 4 | FALSE | 1 | 1 |
| CN_Gongga::H::L | 1800 | 36 | 4 | FALSE | 1 | 1 |
| CN_Gongga::M::L | 1568 | 28 | 4 | FALSE | 1 | 1 |
| CN_Heibei::3400::3200 | 2310 | 55 | 7 | FALSE | 1 | 1 |
| CN_Heibei::3800::3200 | 1360 | 34 | 7 | FALSE | 1 | 1 |
| CN_Heibei::3800::3400 | 1360 | 34 | 7 | FALSE | 1 | 1 |
| DE_Grainau::HIGH::LOW | 4160 | 52 | 4 | FALSE | 1 | 1 |
| DE_Susalps::EB::BT | 3240 | 45 | 4 | FALSE | 1 | 1 |
| DE_Susalps::EB::FE | 3240 | 45 | 4 | FALSE | 1 | 1 |
| DE_Susalps::EB::GW | 3240 | 45 | 4 | FALSE | 1 | 1 |
| DE_Susalps::FE::BT | 864 | 12 | 4 | FALSE | 1 | 1 |
| DE_Susalps::GW::BT | 3168 | 44 | 4 | FALSE | 1 | 1 |
| DE_Susalps::GW::FE | 3168 | 44 | 4 | FALSE | 1 | 1 |
| DE_TransAlps::FP::BT | 1512 | 21 | 4 | FALSE | 1 | 1 |
| DE_TransAlps::SP::BT | 2520 | 35 | 4 | FALSE | 1 | 1 |
| FR_AlpeHuez::HIGH::LOW | 9720 | 54 | 9 | FALSE | 1 | 1 |
| FR_Lautaret::G::L | 4960 | 62 | 4 | FALSE | 1 | 1 |
| IT_MatschMazia1::High::Low | 1122 | 34 | 3 | FALSE | 1 | 1 |
| IT_MatschMazia2::High::Low | 720 | 24 | 3 | FALSE | 1 | 1 |
| NO_Gudmedalen::Gudmedalen::Rambera | 2450 | 49 | 5 | FALSE | 1 | 1 |
| NO_Gudmedalen::Rambera::Arhelleren | 2000 | 40 | 5 | FALSE | 1 | 1 |
| NO_Lavisdalen::Hogsete::Vikesland | 1300 | 26 | 5 | FALSE | 1 | 1 |
| NO_Lavisdalen::Lavisdalen::Hogsete | 2700 | 54 | 5 | FALSE | 1 | 1 |
| NO_Skjellingahaugen::Skjellingahaugen::Veskre | 2050 | 41 | 5 | FALSE | 1 | 1 |
| NO_Skjellingahaugen::Veskre::Ovstedal | 2400 | 48 | 5 | FALSE | 1 | 1 |
| NO_Ulvhaugen::Alrust::Fauske | 1911 | 39 | 5 | FALSE | 1 | 1 |
| NO_Ulvhaugen::Ulvhaugen::Alrust | 2250 | 45 | 5 | FALSE | 1 | 1 |
| SE_Abisko::High::Low | 480 | 16 | 3 | FALSE | 1 | 1 |
| SE_Abisko::High::Mid | 480 | 16 | 3 | FALSE | 1 | 1 |
| SE_Abisko::Mid::Low | 390 | 13 | 3 | FALSE | 1 | 1 |
| US_Arizona::MC::PP | 448 | 8 | 4 | FALSE | 1 | 1 |
| US_Colorado::mo::pf | 1593 | 27 | 5 | FALSE | 1 | 1 |
| US_Colorado::mo::um | 1728 | 27 | 5 | FALSE | 1 | 1 |
| US_Colorado::pf::um | 1560 | 26 | 5 | FALSE | 1 | 1 |
| US_Montana::High::Low | 1568 | 28 | 3 | FALSE | 1 | 1 |

## Property coverage on the residents (added 26 Oct)

Species and effects with each property, after steps 01 and 02:

| property | species | effects |
| --- | --- | --- |
| LDMC | 352 | 1350 |
| LMA | 380 | 1386 |
| height | 334 | 1308 |
| any of the three traits | 397 | 1420 |
| growth form (graminoid, woody) | 528 | 1655 |
| PEV1 to PEV3 | 528 | 1655 |
| niche optima (any GBIF cell) | 527 | 1652 |
| niche optima (n_obs_sub >= 30, used) | 508 | 1618 |

Residents without any of the three traits (131): Agrostis nervosa, Alchemilla conjuncta, Alchemilla flabellata, Aletris pauciflora, Allium beesianum, Allium prattii, Alopecurus aequalis, Anaphalis aureopunctata, Anaphalis flavescens, Androsace minor, Androsace tapete, Anemonastrum demissum, Antennaria racemosa, Argentina leuconota, Argentina stenophylla, Arisaema parvum, Artemisia fulgens, Aster asteroides, Aster diplostephioides, Astragalus mongholicus, Berberis jamesiana, Betula pubescens, Boechera stricta, Calamagrostis rubescens, Calamagrostis tripilifera, Cardamine tangutorum, Carex cercostachys, Carex flava, Carex kansuensis, Carex sargentiana, Carex scabrirostris, Carpesium humile, Cerastium szechuense, Chamaesium viridiflorum, Chrysosplenium pilosum, Clinopodium polycephalum, Codonopsis foetens, Comastoma pulmonarium, Comastoma tenellum, Cyananthus hookeri, Cyananthus incanus, Elymus bakeri, Epilobium fangii, Eremogone roborowskii, Erigeron coulteri, Erigeron neglectus, Euphorbia stracheyi, Euphrasia regelii, Eurybia conspicua, Fargesia spathacea, Galearis spathulata, Gentiana crassuloides, Gentiana grumii, Gentiana prostrata, Gentiana squarrosa, Geranium pylzowianum, Halenia corniculata, Halerpestes tricuspis, Hansenia forbesii, Hansenia weberbaueriana, Hedysarum algidum, Helictotrichon parlatorei, Helictotrichon tibeticum, Hemiphragma heterophyllum, Heterotheca pumila, Hypericum wightianum, Iris potaninii, Juncus allioides, Juncus himalensis, Juncus leucanthus, Koenigia cyanandra, Ligularia sagitta, Ligusticopsis scapiformis, Lomatogonium macranthum, Lomatogonium rotatum, Lycoris radiata, Maianthemum japonicum, Melica spectabilis, Micranthes atrata, Microula sikkimensis, Morina chinensis, Myosotis decumbens, Osmorhiza berteroi, Oxytropis deflexa, Oxytropis densa, Oxytropis humifusa, Oxytropis yunnanensis, Packera cymbalaria, Parasenecio latipes, Pedicularis davidii, Pedicularis mussotii, Pedicularis rhynchodonta, Pedicularis sima, Persicaria runcinata, Poa interior, Poa reflexa, Polygonatum cirrhifolium, Polypogon viridis, Potentilla clandestina, Potentilla pulcherrima, Primula malacoides, Primula serratifolia, Ranunculus inamoenus, Rheum pumilum, Rhodiola fastigiata, Rhodiola yunnanensis, Rhododendron websterianum, Sagina glabra, Sagina japonica, Salix serpillifolia, Saussurea ceterach, Saussurea globosa, Saussurea graminea, Saussurea katochaete, Saussurea pachyneura, Saussurea pulchra, Solidago multiradiata, Sorbus aucuparia, Soroseris hookeriana, Swertia macrosperma, Symphyotrichum campestre, Tanacetum tatsienense, Thalictrum occidentale, Thesium alpinum, Vaccinium membranaceum, Veronica allionii, Veronica szechuanica, Viola nuttallii, Viola orbiculata, Viola szetschwanensis, Youngia racemifera.

## Covariates

Raw values before standardisation:

| variable | range |
| --- | --- |
| W (destination - origin JJA T, degC) | 0.3571 / 1.5285 / 2.3140 / 3.3592 / 11.0886 |
| T (destination JJA T, degC) | 8.7911 / 11.3936 / 13.2811 / 16.4864 / 21.0185 |
| P (log destination JJA precipitation, mm) | 5.0594 / 5.5447 / 5.9308 / 6.1053 / 6.9784 |
| A0_origin (log(% + 1)) | 0.0000 / 0.0998 / 0.5494 / 1.2698 / 3.8012 |
| A0_dest (log(% + 1)) | 0.0000 / 0.0000 / 0.0000 / 0.4121 / 3.8040 |

## Warnings captured

- CH_Calanda::Cal::Nes: NOTE: Results may be misleading due to involvement in interactions
- CH_Calanda::Pea::Cal: NOTE: Results may be misleading due to involvement in interactions
- CH_Calanda2::Cal::Nes: NOTE: Results may be misleading due to involvement in interactions
- CH_Lavey::RIO::CRE: NOTE: Results may be misleading due to involvement in interactions
- CH_Lavey::RIO::MAR: NOTE: Results may be misleading due to involvement in interactions
- CH_Lavey::RIO::PRA: NOTE: Results may be misleading due to involvement in interactions
- CN_Damxung::HIGH::LOW: boundary (singular) fit: see help('isSingular')
- CN_Damxung::HIGH::LOW: NOTE: Results may be misleading due to involvement in interactions
- CN_Gongga::A::M: NOTE: Results may be misleading due to involvement in interactions
- CN_Gongga::H::A: NOTE: Results may be misleading due to involvement in interactions
- CN_Gongga::H::L: NOTE: Results may be misleading due to involvement in interactions
- CN_Gongga::M::L: NOTE: Results may be misleading due to involvement in interactions
- CN_Heibei::3400::3200: NOTE: Results may be misleading due to involvement in interactions
- CN_Heibei::3800::3200: NOTE: Results may be misleading due to involvement in interactions
- CN_Heibei::3800::3400: NOTE: Results may be misleading due to involvement in interactions
- DE_Grainau::HIGH::LOW: NOTE: Results may be misleading due to involvement in interactions
- DE_Susalps::EB::BT: NOTE: Results may be misleading due to involvement in interactions
- DE_Susalps::EB::FE: NOTE: Results may be misleading due to involvement in interactions
- DE_Susalps::EB::GW: NOTE: Results may be misleading due to involvement in interactions
- DE_Susalps::FE::BT: NOTE: Results may be misleading due to involvement in interactions
- DE_Susalps::GW::BT: NOTE: Results may be misleading due to involvement in interactions
- DE_Susalps::GW::FE: NOTE: Results may be misleading due to involvement in interactions
- DE_TransAlps::FP::BT: NOTE: Results may be misleading due to involvement in interactions
- DE_TransAlps::SP::BT: NOTE: Results may be misleading due to involvement in interactions
- FR_AlpeHuez::HIGH::LOW: NOTE: Results may be misleading due to involvement in interactions
- FR_Lautaret::G::L: NOTE: Results may be misleading due to involvement in interactions
- IT_MatschMazia1::High::Low: NOTE: Results may be misleading due to involvement in interactions
- IT_MatschMazia2::High::Low: NOTE: Results may be misleading due to involvement in interactions
- NO_Gudmedalen::Gudmedalen::Rambera: NOTE: Results may be misleading due to involvement in interactions
- NO_Gudmedalen::Rambera::Arhelleren: NOTE: Results may be misleading due to involvement in interactions
- NO_Lavisdalen::Hogsete::Vikesland: NOTE: Results may be misleading due to involvement in interactions
- NO_Lavisdalen::Lavisdalen::Hogsete: NOTE: Results may be misleading due to involvement in interactions
- NO_Skjellingahaugen::Skjellingahaugen::Veskre: NOTE: Results may be misleading due to involvement in interactions
- NO_Skjellingahaugen::Veskre::Ovstedal: NOTE: Results may be misleading due to involvement in interactions
- NO_Ulvhaugen::Alrust::Fauske: NOTE: Results may be misleading due to involvement in interactions
- NO_Ulvhaugen::Ulvhaugen::Alrust: NOTE: Results may be misleading due to involvement in interactions
- SE_Abisko::High::Low: NOTE: Results may be misleading due to involvement in interactions
- SE_Abisko::High::Mid: NOTE: Results may be misleading due to involvement in interactions
- SE_Abisko::Mid::Low: NOTE: Results may be misleading due to involvement in interactions
- US_Arizona::MC::PP: NOTE: Results may be misleading due to involvement in interactions
- US_Colorado::mo::pf: NOTE: Results may be misleading due to involvement in interactions
- US_Colorado::mo::um: NOTE: Results may be misleading due to involvement in interactions
- US_Colorado::pf::um: NOTE: Results may be misleading due to involvement in interactions
- US_Montana::High::Low: NOTE: Results may be misleading due to involvement in interactions

## Deviations from PLAN.md

None in the model. Decisions taken where PLAN.md is silent (30 Oct, Bills): asymptotic degrees of freedom for every experiment; PEV1 to PEV3 carried; properties standardised across the residents that have them, traits logged first; graminoid and woody as 0/1. MILESTONES' expected count (about 1,596 effects in 41 experiments) is 1,655 in 43, after FR_AlpeHuez (26 Oct) and US_Arizona (28 Oct) were included.
