# P2: TRY traits

TRY records were read by TRY species ID (1340 IDs whose contributor or TRY name resolves to a study species in P1) and each record was attached to a species through its own resolved name (the contributor's `SpeciesName` first, TRY's `AccSpeciesName` only as fallback). Observations under experimental or ambiguous treatment values were removed, then duplicates, ErrorRisk >= 4, observations without coordinates and those in arid and tropical climate zones, as in the legacy pipeline. 674 of 904 study species have at least one trait. LDMC in g/g, LMA in mg mm-2 (1 / SLA per observation), height in m.

## Counts

| stage | rows | observations |
| --- | --- | --- |
| read by TRY species ID (1340 IDs) | 2877931 | 179998 |
| rows of study species (by our name or TRY's) | 2832174 | 177075 |
| experimental observations removed | 2516912 | 156231 |
| TRY duplicates removed | 2457036 | 154333 |
| quantitative, ErrorRisk < 3 | 145916 | 96782 |
| LDMC, SLA and vegetative height | 89236 | 61254 |
| with coordinates and a climate zone | 77869 | 50106 |
| outside the excluded zones Af, Am, Aw, BWh, BWk, BSh, BSk | 77385 | 49694 |

## Treatment values

Counted as natural (`prep/config/try_natural_values.csv`, compared after lower-casing and dropping trailing punctuation): "none", "no", "no treatment", "untreated", "natural", "natural environment", "natural env", "nat env", "natural conditions", "natural vegetation", "natural habitat", "natural forest environment", "natural grassland", "field", "field conditions", "in situ", "wild", "natural envireonment", "Natural environment (sun leaves)", "natural environment, sun exposed".

Every other value marks the observation as experimental: 437 distinct values on 44891 rows. The 40 most frequent:

| value | n_rows |
| --- | --- |
| natural environment (biodiversity experiment | 7927 |
| greenhouse | 1245 |
| low | 1019 |
| intermediate | 880 |
| high | 678 |
| meadows (m) and pastures (p) on south east to south west exposed slopes | 669 |
| 16 | 668 |
| common garden in growth containers with soil corresponding to seed origin | 597 |
| growth chamber | 585 |
| 1 | 561 |
| glasshouse | 536 |
| control | 526 |
| pm | 511 |
| outdoor cultivation in 30x30cm pots, standardized substrate | 486 |
| open-top chamber | 457 |
| 350 | 452 |
| 300 | 446 |
| 09/28/2016 | 400 |
| 2 | 372 |
| 1.12.2008 | 363 |
| 180 | 363 |
| 360 | 363 |
| 941 | 363 |
| well watered | 363 |
| glasshouse with ambient condition (retractable roof | 344 |
| natural_environment | 343 |
| 33209 | 331 |
| n | 331 |
| 25 | 321 |
| 0 | 314 |
| 65 | 310 |
| pot | 303 |
| conspecific | 302 |
| non-limiting | 300 |
| semi-controlled | 300 |
| growth chambers | 297 |
| heterospecific | 295 |
| sand | 269 |
| mesocosm | 263 |
| rehydr | 246 |

## Values not used

Non-positive values dropped before the medians (decided 26 Oct): 17 observations.

| trait | DatasetID | Dataset | observations |
| --- | --- | --- | --- |
| LDMC | 291 | MARGINS - leaf traits database | 1 |
| SLA | 291 | MARGINS - leaf traits database | 3 |
| vegetative_height | 186 | Growth and Herbivory of Juvenil Trees | 5 |
| vegetative_height | 247 | Traits of Halohytic Species in North-West-Germany | 6 |
| vegetative_height | 468 | Tundra Trait Team | 2 |

Text fields repaired from windows-1252 to UTF-8 after reading: LastName (50133 values), FirstName (116262 values), Dataset (22224 values), SpeciesName (2381 values), DataName (370 values), OriglName (8131 values), OrigValueStr (37743 values), OrigUnitStr (13190 values), Reference (432047 values), Comment (7855 values), StdValueStr (70 values).

## Species with traits

| trait | species |
| --- | --- |
| LDMC | 567 |
| LMA | 629 |
| height | 568 |
| any of the three | 674 |

Study species without any of the three traits (230): Agastache urticifolia, Agrostis nervosa, Ajania przewalskii, Alchemilla alpigena, Alchemilla conjuncta, Alchemilla flabellata, Aletris pauciflora, Allium beesianum, Allium brevistylum, Allium prattii, Alopecurus aequalis, Anaphalis aureopunctata, Anaphalis flavescens, Anaphalis xylorhiza, Androsace graminifolia, Androsace mariae, Androsace minor, Androsace occidentalis, Androsace tapete, Anemonastrum demissum, Antennaria microphylla, Antennaria racemosa, Arenaria hispanica, Argentina leuconota, Argentina stenophylla, Arisaema parvum, Artemisia fulgens, Artemisia hedinii, Aster asteroides, Aster diplostephioides, Aster tataricus, Astragalus australis, Astragalus mongholicus, Athyrium spinulosum, Atriplex argentea, Berberis jamesiana, Boechera stricta, Bromus porteri, Calamagrostis rubescens, Calamagrostis tripilifera, Cardamine tangutorum, Carex cercostachys, Carex flava, Carex foetida, Carex kansuensis, Carex lehmannii, Carex nubigena, Carex phaeocephala, Carex sargentiana, Carex scabrirostris, Carex thibetica, Carpesium humile, Cerastium szechuense, Chamaesium viridiflorum, Chamorchis alpina, Chenopodium atrovirens, Chenopodium leptophyllum, Chrysosplenium pilosum, Claytonia lanceolata, Clinopodium polycephalum, Codonopsis foetens, Collomia linearis, Coluria longifolia, Comastoma pulmonarium, Comastoma tenellum, Corydalis pallida, Crepis pontana, Cryptogramma crispa, Cyananthus hookeri, Cyananthus incanus, Delphinium occidentale, Delphinium trichophorum, Draba spectabilis, Drymocallis glandulosa, Elsholtzia ciliata, Elsholtzia densa, Elymus bakeri, Epilobium fangii, Eremogone roborowskii, Erigeron coulteri, Erigeron elatior, Erigeron neglectus, Eritrichium hemisphaericum, Euphorbia pekinensis, Euphorbia prolifera, Euphorbia serrula, Euphorbia stracheyi, Euphrasia pectinata, Euphrasia regelii, Eurybia conspicua, Eurybia integrifolia, Fargesia spathacea, Fragaria nubicola, Galearis spathulata, Galium bifolium, Gentiana aristata, Gentiana crassuloides, Gentiana grumii, Gentiana parryi, Gentiana prostrata, Gentiana scabra, Gentiana squarrosa, Gentiana trichotoma, Gentianella azurea, Geranium pylzowianum, Goodyera oblongifolia, Hackelia micrantha, Halenia corniculata, Halerpestes tricuspis, Hansenia forbesii, Hansenia weberbaueriana, Hedysarum algidum, Helianthella uniflora, Helianthus petiolaris, Helictotrichon parlatorei, Helictotrichon tibeticum, Hemiphragma heterophyllum, Heterotheca pumila, Hypericum wightianum, Iliamna rivularis, Iris potaninii, Juncus allioides, Juncus himalensis, Juncus leucanthus, Juncus prismatocarpus, Juncus setchuensis, Koenigia cyanandra, Lappula occidentalis, Leontodon hyoseroides, Leontopodium japonicum, Ligularia pleurocaulis, Ligularia sagitta, Ligusticopsis scapiformis, Lithophragma glabrum, Lithophragma parviflorum, Lomatogonium macranthum, Lomatogonium rotatum, Luzula alpina, Lycoris radiata, Maianthemum henryi, Maianthemum japonicum, Malaxis monophyllos, Meconopsis horridula, Melica spectabilis, Micranthes atrata, Micranthes occidentalis, Microsteris gracilis, Microula sikkimensis, Morina chinensis, Morina kokonorica, Morina nepalensis, Myosotis decumbens, Myriactis nepalensis, Nemophila breviflora, Notholirion bulbuliferum, Noticastrum adscendens, Orobanche fasciculata, Osmorhiza berteroi, Oxytropis deflexa, Oxytropis densa, Oxytropis humifusa, Oxytropis yunnanensis, Packera cymbalaria, Parasenecio latipes, Pedicularis davidii, Pedicularis mussotii, Pedicularis rhynchodonta, Pedicularis sima, Perideridia gairdneri, Persicaria runcinata, Poa interior, Poa reflexa, Polygonatum cirrhifolium, Polypogon viridis, Pomatosace filicula, Potentilla brauneana, Potentilla clandestina, Potentilla pulcherrima, Primula deflexa, Primula malacoides, Ranunculus acriformis, Ranunculus inamoenus, Ranunculus japonicus, Rheum pumilum, Rhodiola fastigiata, Rhodiola quadrifida, Rhodiola yunnanensis, Rhododendron websterianum, Rudbeckia occidentalis, Rumex nepalensis, Sagina glabra, Sagina japonica, Salix serpillifolia, Saussurea ceterach, Saussurea globosa, Saussurea graminea, Saussurea katochaete, Saussurea pachyneura, Saussurea pulchra, Saussurea stella, Saxifraga cespitosa, Sedum triactina, Sibbaldia cuneata, Solidago multiradiata, Soroseris hookeriana, Sphaeralcea ambigua, Stellaria irrigua, Stellaria longifolia, Stipa capillacea, Swertia macrosperma, Symphyotrichum campestre, Tanacetum tatsienense, Thalictrum javanicum, Thalictrum occidentale, Thermopsis lanceolata, Thesium alpinum, Townsendia smithii, Vaccinium membranaceum, Vaccinium scoparium, Valeriana excelsa, Veronica allionii, Veronica ciliata, Veronica praecox, Veronica szechuanica, Viola nuttallii, Viola orbiculata, Viola philippica, Viola szetschwanensis, Youngia racemifera, Youngia simulatrix.

## Trusting TRY's accepted name against attaching by the record's own name

Records of species x trait whose count or median differs between the two (median differing by more than 5%): 151 of 1778 species x trait combinations. Festuca ovina now holds 0 records submitted as Festuca pratensis.

| species | trait | n_ours | median_ours | n_trusting_try | median_trusting_try | relative_difference |
| --- | --- | --- | --- | --- | --- | --- |
| Festuca ovina | LMA | 115 | 0.086095566078347 | 146 | 0.0788982259570495 | -0.084 |
| Festuca ovina | height | 91 | 0.07 | 116 | 0.0915 | 0.307 |
| Lolium pratense | LDMC | 245 | 0.27 | 228 | 0.269615384615384 | -0.001 |
| Lolium pratense | LMA | 188 | 0.0521770260927943 | 157 | 0.0535045478865704 | 0.025 |
| Lolium pratense | height | 258 | 0.31 | 233 | 0.28 | -0.097 |
| Festuca ovina | LDMC | 88 | 0.450356911 | 105 | 0.438461538 | -0.026 |

The 30 largest differences:

| species | trait | n_ours | median_ours | n_trusting_try | median_trusting_try | relative_difference |
| --- | --- | --- | --- | --- | --- | --- |
| Taraxacum officinale | LMA | 295 | 0.0369763729103778 |  |  |  |
| Taraxacum officinale | height | 178 | 0.205 |  |  |  |
| Taraxacum officinale | LDMC | 211 | 0.17674517293 |  |  |  |
| Anthoxanthum nipponicum | height | 17 | 0.093 |  |  |  |
| Anthoxanthum nipponicum | LDMC | 29 | 0.31958 |  |  |  |
| Aster bellidiastrum | LDMC | 15 | 0.134859011033919 |  |  |  |
| Trichophorum cespitosum | LDMC | 5 | 0.28929 |  |  |  |
| Anthoxanthum nipponicum | LMA | 32 | 0.051138333107998 |  |  |  |
| Aster bellidiastrum | LMA | 19 | 0.0355346634296281 |  |  |  |
| Centaurea montana | height | 2 | 0.41 |  |  |  |
| Centaurea montana | LMA | 3 | 0.0435016884812633 |  |  |  |
| Ranunculus tuberosus | height | 1 | 0.4325 |  |  |  |
| Ranunculus tuberosus | LMA | 4 | 0.0529513801683058 |  |  |  |
| Aster bellidiastrum | height | 13 | 0.08 |  |  |  |
| Aquilegia coerulea | LMA | 16 | 0.0511165599954287 |  |  |  |
| Euphorbia serpillifolia | LMA | 36 | 0.0320482168089139 |  |  |  |
| Dysphania incisa | height | 15 | 0.17 |  |  |  |
| Dysphania incisa | LMA | 29 | 0.0301617768028517 |  |  |  |
| Euphrasia minima | LDMC | 21 | 0.1809066221 |  |  |  |
| Knorringia sibirica | height | 1 | 0.25 |  |  |  |
| Aquilegia coerulea | height | 7 | 0.3072 |  |  |  |
| Hieracium triste | height | 7 | 0.1535 |  |  |  |
| Hieracium triste | LMA | 13 | 0.0462247341533242 |  |  |  |
| Cynoglossum officinale | height | 1 | 0.4 |  |  |  |
| Myosotis arvensis | height | 1 | 0.2 |  |  |  |
| Lappula squarrosa | height | 1 | 0.4 |  |  |  |
| Saussurea superba | LMA | 2 | 0.105225030902826 |  |  |  |
| Poa alberti | LMA | 9 | 0.0633457830369215 |  |  |  |
| Poa alberti | LDMC | 10 | 0.4825 |  |  |  |
| Trichophorum distigmaticum | LDMC | 12 | 0.419797096 |  |  |  |

## Against the frozen legacy traits (step 00)

Species with traits now: 674; in step 00: 558; species gained: 152; lost: 36 (names that changed count as both); species x trait medians changed by more than 5%: 572. The 30 largest:

| species | trait | new | legacy | relative_difference |
| --- | --- | --- | --- | --- |
| Euphrasia officinalis | height | 0.2325 | 0.0395 | 4.886 |
| Phleum alpinum | height | 0.2925 | 0.0694 | 3.215 |
| Cerastium fontanum | height | 0.25 | 0.068305 | 2.66 |
| Galium anisophyllon | height | 0.23208 | 0.06416 | 2.617 |
| Trifolium hybridum | height | 0.494 | 0.152 | 2.25 |
| Geranium sylvaticum | height | 0.37 | 0.115 | 2.217 |
| Pyrola minor | height | 0.0995 | 0.0316 | 2.149 |
| Trollius europaeus | height | 0.5275 | 0.198 | 1.664 |
| Euphrasia officinalis | LMA | 0.1730821395 | 0.0650364203954214 | 1.661 |
| Trifolium dubium | height | 0.225350459362423 | 0.09 | 1.504 |
| Poa annua | height | 0.278 | 0.1165 | 1.386 |
| Agrostis mertensii | height | 0.14 | 0.06 | 1.333 |
| Hieracium alpinum | height | 0.138 | 0.06 | 1.3 |
| Lactuca serriola | height | 0.8 | 0.36 | 1.222 |
| Ranunculus repens | height | 0.265909090909091 | 0.12 | 1.216 |
| Luzula spicata | height | 0.1218 | 0.0554 | 1.199 |
| Astragalus danicus | height | 0.10975 | 0.05 | 1.195 |
| Gentianella campestris | height | 0.109 | 0.05 | 1.18 |
| Euphrasia salisburgensis | LMA | 0.114625942468789 | 0.055005500550055 | 1.084 |
| Carex pallescens | height | 0.6535 | 0.315 | 1.075 |
| Salix herbacea | height | 0.049 | 0.0237 | 1.068 |
| Betula pubescens | height | 0.27 | 5.27 | -0.949 |
| Hedysarum hedysaroides | height | 0.38 | 0.2 | 0.9 |
| Echium vulgare | height | 0.0454 | 0.4227 | -0.893 |
| Stellaria graminea | height | 0.5 | 0.269 | 0.859 |
| Dactylorhiza viridis | LMA | 0.0461039768439346 | 0.0249149221203853 | 0.85 |
| Androsace septentrionalis | height | 0.02 | 0.011 | 0.818 |
| Noccaea montana | height | 0.0335 | 0.0185 | 0.811 |
| Festuca vivipara | height | 0.11 | 0.061 | 0.803 |
| Festuca ovina | height | 0.07 | 0.3504 | -0.8 |

## Flagged TRY names of study species

A flagged TRY name is used when its accepted species is a study species and the matched name has the submitted genus and epithet (decided 26 Oct).

| used | reason | names |
| --- | --- | --- |
| FALSE | fallback source wfo | 3 |
| FALSE | fuzzy match | 11 |
| TRUE | fallback source wfo | 167 |
| TRUE | fallback source wfo; fuzzy match | 30 |
| TRUE | fuzzy match | 369 |

## Warnings captured

None.

## Deviations from the plan

(1) Community name decisions (e.g. Thlaspi montanum -> Noccaea fendleri for US_Montana) are not carried over to TRY names, because they concern the study records; TRY names are resolved by the rule alone, with the accepted-name renames of P1 (Taraxacum officinale) applied. (2) ErrorRisk < 3 instead of the legacy < 4 (decided 26 Oct): with experiments removed per observation, datasets that legacy dropped wholesale now contribute, and some of their records 3 to 4 SD from TRY's species mean gave implausible medians (Viola canina height 1.16 m, Euphrasia officinalis LMA 0.43 mg mm-2). (3) Non-positive values dropped (decided 26 Oct). (4) TRY text fields repaired to UTF-8 after reading (legacy repaired names only).
