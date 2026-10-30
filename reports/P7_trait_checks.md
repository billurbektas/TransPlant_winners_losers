# P7: Trait checks

Checks of the 73,809 TRY records behind P2's medians after P2's screening (666 species; LDMC, LMA, height), within TRY only and on the log scale. Nothing is changed here. The rules Bills accepted on 29-30 Oct (trees, two reviewed values, record outliers and unit slips) are P2 cleaning targets, so this report shows what remains after them; new record flags can appear because the species medians moved. Thresholds: robust z 3.5; datasets off by 1.5x on at least 3 shared species; unit signatures at 10, 100 and 1000x within 15%. Full lists in `reports/P7_trait_checks/`.

## Counts

| check | n |
| --- | --- |
| records: outlier (|z| > 3.5) | 175 |
| records: unit signature | 1 |
| species x trait medians moving > 10% without flagged records | 0 |
| datasets off by >= 1.5x | 19 |
| species flagged against relatives or plausibility | 46 |
| fragile medians (< 3 records or one dataset) | 532 |

Species x trait groups with at least 5 records and zero spread (z undefined): 10.
Species check: spread = the larger of the relatives' own spread and the typical within-genus spread (log scale; LMA 0.18, LDMC 0.13, height 0.38).

Known case (P2 notes): Euphrasia officinalis LMA from dataset 524 flagged at record level: FALSE; dataset level: TRUE; species level: FALSE; still in P2's records: FALSE (set to NA by the review ledger, 30 Oct).

## Datasets off by 1.5x or more against the other datasets

| trait | dataset_id | offset | shared_species | records | species_fed | species_mostly_this_dataset |
| --- | --- | --- | --- | --- | --- | --- |
| height | 158 | 4.12 | 7 | 18 | 7 | 0 |
| height | 474 | 0.334 | 8 | 8 | 8 | 1 |
| height | 428 | 0.346 | 20 | 2394 | 20 | 16 |
| LMA | 524 | 2.58 | 18 | 73 | 18 | 0 |
| height | 115 | 2.5 | 37 | 37 | 38 | 5 |
| height | 272 | 0.411 | 5 | 5 | 5 | 1 |
| height | 65 | 2.39 | 14 | 27 | 14 | 1 |
| height | 192 | 2.08 | 18 | 18 | 19 | 2 |
| height | 273 | 1.94 | 7 | 7 | 7 | 0 |
| height | 90 | 1.92 | 9 | 21 | 9 | 2 |
| height | 76 | 1.71 | 98 | 169 | 104 | 21 |
| height | 247 | 1.67 | 3 | 38 | 3 | 2 |
| height | 50 | 1.65 | 42 | 241 | 45 | 8 |
| height | 1 | 0.647 | 6 | 6 | 6 | 0 |
| height | 617 | 1.54 | 10 | 20 | 11 | 2 |
| height | 86 | 1.54 | 92 | 106 | 96 | 13 |
| height | 193 | 0.656 | 22 | 304 | 43 | 32 |
| height | 547 | 0.661 | 18 | 19 | 21 | 5 |
| height | 320 | 0.663 | 42 | 601 | 45 | 16 |

## Records: the 20 largest departures from their species median

| species | trait | dataset_id | observation_id | SpeciesName | value | species_median | ratio | z | reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Lotus corniculatus | height | 320 | 2474870 | Lotus corniculatus | 0.03 | 0.27 | 0.111 | -2.67 | unit signature |
| Rumex acetosa | height | 76 | 940965 | Rumex acetosa | 0.546 | 0.07 | 7.8 | 3.68 | outlier |
| Rumex acetosa | height | 192 | 1396716 | RUMEX ACETOSA L. | 0.51 | 0.07 | 7.29 | 3.56 | outlier |
| Geranium sylvaticum | height | 25 | 175293 | Geranium sylvaticum | 0.06 | 0.378 | 0.159 | -3.89 | outlier |
| Cerastium fontanum | height | 519 | 5995524 | Cerastium fontanum | 0.043 | 0.265 | 0.162 | -4.65 | outlier |
| Cerastium fontanum | height | 519 | 5995533 | Cerastium fontanum | 0.044 | 0.265 | 0.166 | -4.59 | outlier |
| Cerastium fontanum | height | 519 | 5995528 | Cerastium fontanum | 0.045 | 0.265 | 0.17 | -4.53 | outlier |
| Geranium sylvaticum | height | 25 | 1117784 | Geranium sylvaticum | 0.07 | 0.378 | 0.185 | -3.57 | outlier |
| Poa pratensis | height | 428 | 2951250 | Poa pratensis L. | 0.054 | 0.29 | 0.186 | -3.53 | outlier |
| Carum carvi | height | 519 | 5995511 | Carum carvi | 0.056 | 0.2759 | 0.203 | -3.79 | outlier |
| Carum carvi | height | 519 | 5995507 | Carum carvi | 0.057 | 0.2759 | 0.207 | -3.75 | outlier |
| Poa alpina | height | 468 | 5815804 | Poa alpina | 0.057 | 0.2716 | 0.21 | -4.02 | outlier |
| Poa alpina | height | 168 | 1339254 | Poa alpina | 0.06 | 0.2716 | 0.221 | -3.89 | outlier |
| Poa alpina | height | 168 | 1339279 | Poa alpina | 0.06 | 0.2716 | 0.221 | -3.89 | outlier |
| Poa alpina | height | 468 | 5815803 | Poa alpina | 0.061 | 0.2716 | 0.225 | -3.85 | outlier |
| Cerastium fontanum | height | 519 | 5995529 | Cerastium fontanum | 0.06 | 0.265 | 0.226 | -3.8 | outlier |
| Poa alpina | height | 468 | 5815819 | Poa alpina | 0.064 | 0.2716 | 0.236 | -3.73 | outlier |
| Luzula spicata | height | 468 | 5814041 | Luzula spicata | 0.032 | 0.127 | 0.252 | -4.14 | outlier |
| Poa alpina | height | 468 | 5815817 | Poa alpina | 0.069 | 0.2716 | 0.254 | -3.53 | outlier |
| Silene vulgaris | height | 519 | 5996683 | Silene vulgaris | 0.099 | 0.3853 | 0.257 | -5.76 | outlier |

## Species medians that move most without their flagged records

| species | trait | median_all | median_without | n_flagged | shift_pct |
| --- | --- | --- | --- | --- | --- |
| Cerastium fontanum | height | 0.265 | 0.284 | 4 | 7 |
| Geranium sylvaticum | height | 0.378 | 0.4 | 2 | 5.8 |
| Carum carvi | height | 0.276 | 0.29 | 2 | 5.1 |
| Myosotis sylvatica | LMA | 0.019 | 0.0181 | 1 | -4.3 |
| Galium album | LMA | 0.0456 | 0.0441 | 1 | -3.5 |
| Silene vulgaris | height | 0.386 | 0.397 | 1 | 3 |
| Lappula squarrosa | LDMC | 0.232 | 0.226 | 2 | -2.5 |
| Trifolium hybridum | height | 0.504 | 0.514 | 1 | 2 |
| Polygala vulgaris | height | 0.25 | 0.255 | 1 | 2 |
| Lotus corniculatus | height | 0.27 | 0.275 | 1 | 1.9 |
| Galium pumilum | LMA | 0.0438 | 0.0446 | 1 | 1.9 |
| Hippophae rhamnoides | LMA | 0.0784 | 0.077 | 1 | -1.8 |
| Trollius europaeus | height | 0.597 | 0.607 | 1 | 1.6 |
| Carex curvula | LDMC | 0.335 | 0.329 | 5 | -1.6 |
| Hippocrepis comosa | height | 0.031 | 0.0305 | 1 | -1.6 |
| Erigeron canadensis | height | 0.72 | 0.73 | 1 | 1.4 |
| Alchemilla pentaphyllea | LMA | 0.0494 | 0.0487 | 1 | -1.3 |
| Erigeron uniflorus | LMA | 0.0445 | 0.0439 | 1 | -1.3 |
| Lilium martagon | LDMC | 0.109 | 0.108 | 1 | -1.3 |
| Luzula spicata | height | 0.127 | 0.128 | 1 | 1.2 |

## Species flagged against their relatives or by plausibility rules

| species | trait | median | reason |
| --- | --- | --- | --- |
| Arnica angustifolia | height | 0.05 | far from its genus (robust z -3.6, 3 relatives) |
| Bromus inermis | height | 0.72 | far from its genus (robust z 3.8, 3 relatives) |
| Campanula rapunculoides | LDMC | 0.142 | far from its genus (robust z -4.3, 3 relatives) |
| Campanula rapunculoides | LMA | 0.0223 | far from its genus (robust z -4.8, 4 relatives) |
| Cirsium acaule | height | 0.085 | far from its genus (robust z -6.1, 4 relatives) |
| Cirsium oleraceum | LMA | 0.0303 | far from its genus (robust z -3.7, 4 relatives) |
| Crepis aurea | height | 0.08 | far from its genus (robust z -5.1, 4 relatives) |
| Draba aizoides | LDMC | 0.28 | far from its family (robust z 4.8, 8 relatives) |
| Drymaria molluginea | LDMC | 0.138 | LMA and LDMC at opposite extremes |
| Drymaria molluginea | LMA | 0.121 | LMA and LDMC at opposite extremes |
| Epilobium anagallidifolium | height | 0.05 | far from its genus (robust z -7.3, 3 relatives) |
| Epilobium angustifolium | height | 0.87 | far from its genus (robust z 7.6, 3 relatives) |
| Epilobium lactiflorum | height | 0.05 | far from its genus (robust z -7.3, 3 relatives) |
| Epilobium obscurum | height | 0.8 | far from its genus (robust z 7.3, 3 relatives) |
| Erigeron uniflorus | height | 0.0325 | far from its genus (robust z -3.9, 6 relatives) |
| Festuca thurberi | height | 0.572 | far from its genus (robust z 4, 6 relatives) |
| Galium aparine | LMA | 0.0204 | far from its genus (robust z -4, 9 relatives) |
| Gentiana nivalis | LDMC | 0.146 | far from its genus (robust z -3.7, 9 relatives) |
| Geum aleppicum | LDMC | 0.165 | far from its genus (robust z -3.9, 4 relatives) |
| Geum montanum | height | 0.049 | far from its genus (robust z -5.1, 4 relatives) |
| Hieracium umbellatum | height | 1.03 | far from its genus (robust z 4.7, 4 relatives) |
| Holcus lanatus | LDMC | 0.21 | graminoid LDMC below the herb median |
| Juncus alpinoarticulatus | LMA | 0.113 | far from its family (robust z 5, 12 relatives) |
| Kalmia procumbens | height | 0.0135 | far from its family (robust z -3.9, 9 relatives) |
| Kalmia procumbens | height | 0.0135 | woody, shorter than 95% of herbs |
| Lolium multiflorum | LDMC | 0.186 | graminoid LDMC below the herb median |
| Luzula pilosa | LDMC | 0.196 | graminoid LDMC below the herb median |
| Luzula sudetica | LDMC | 0.195 | graminoid LDMC below the herb median |
| Mutellina purpurea | height | 0.0435 | far from its family (robust z -3.7, 23 relatives) |
| Myosotis alpestris | LMA | 0.0462 | far from its genus (robust z 4.7, 3 relatives) |
| Neogaya simplex | height | 0.035 | far from its family (robust z -4.1, 23 relatives) |
| Neottia ovata | LMA | 0.0259 | far from its family (robust z -3.5, 5 relatives) |
| Noccaea montana | height | 0.0335 | far from its family (robust z -4.3, 15 relatives) |
| Pedicularis resupinata | height | 0.4 | far from its genus (robust z 4.1, 3 relatives) |
| Phyteuma hemisphaericum | LMA | 0.0781 | far from its genus (robust z 3.9, 3 relatives) |
| Phyteuma orbiculare | height | 0.283 | far from its genus (robust z 4.6, 3 relatives) |
| Poa alberti | LDMC | 0.482 | far from its genus (robust z 4, 7 relatives) |
| Poa annua | LDMC | 0.187 | far from its genus (robust z -3.8, 7 relatives) |
| Poa annua | LDMC | 0.187 | graminoid LDMC below the herb median |
| Poa supina | LDMC | 0.187 | far from its genus (robust z -3.8, 7 relatives) |
... 9 more in `flagged_species.csv`.

## What to decide

1. Datasets: exclude a flagged dataset for that trait if its offset reflects a unit or protocol difference (the table above lists how many species medians it feeds and for how many it is the main source).
2. Records: drop outliers and unit-signature records as a rule, or only those that move a species median by more than 10%.
3. Species: species flagged against their relatives are inspected by name; the plausibility flags are soft and not exclusions by themselves.

## Warnings captured

None.

## Deviations from PLAN.md

None; a check of P2's inputs, requested 29 Oct.
