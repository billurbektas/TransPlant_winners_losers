# Step 02: Strength (Step A, residents)

One Gaussian `glmmTMB(relative_cover ~ treatment * year * species + (1 | plot_id:species), dispformula = ~ species, REML = TRUE)` (one residual variance per species) was fitted per experiment on the resident species' rows in origin controls and warmed plots (66,412 rows, 43 experiments), year as a factor. `emmeans` gave the warmed minus origin-control contrast per species and survey year (`by_year`, for step 03) and averaged over each experiment's years (`strength`), with degrees of freedom = origin-control plots + warmed plots - 2 of the experiment. Each effect carries the experiment's W, T and P (summer, log) standardised across the 43 experiments, TP = T x P, A0 at origin and destination standardised across the 995 effects, and the species properties standardised across the residents that have them. Strength is in relative-cover units (proportion of plot cover).

## Model as fitted (for reviewers)

The call, as stored in the fitted model objects, and the contrasts (`R/fit_experiment.R`), run once per experiment:

```r
m <- glmmTMB::glmmTMB(formula = relative_cover ~ treatment * year * species + (1 | plot_id:species), data = d, dispformula = ~species, REML = TRUE, ziformula = ~0)
# d: one experiment's resident rows in origin controls and warmed plots;
# treatment (origin_control, warmed), year (factor), species (factor), plot_id
df_plots <- n_distinct(plot_id[origin_control]) + n_distinct(plot_id[warmed]) - 2
by_year  <- contrast(emmeans(m, ~ treatment | species * year, df = df_plots), method = "revpairwise")
strength <- contrast(emmeans(m, ~ treatment | species, df = df_plots), method = "revpairwise")
# keep "warmed - origin_control"; emmeans averages the years with equal weights
```

- Fixed part: one mean per treatment x year x species cell. Residual variance: one per species (`dispformula = ~ species`, decided 30 Oct). REML; Gaussian on relative cover (bounded, zero-heavy: accepted approximation, PLAN.md Known weaknesses).
- Degrees of freedom: origin-control plots + warmed plots - 2 of the experiment (decided 30 Oct), because treatment is applied to whole plots and the yearly surveys of a plot are repeated measures; 4 to 18 (median 12), so the two-sided 5% cut-off is t = 2.10 to 2.78 instead of z = 1.96. They change p-values only, not `strength` or `strength_se`.
- Random part, a known limitation kept by decision (30 Oct): `(1 | plot_id:species)` has one plot x species variance shared by all species of an experiment, while the residual variance is per species. This shared term is more than half of the SE variance in 427 of 995 effects (median share 0.33), and a species' own plot-to-plot spread grows with its residual SD (Spearman 0.93): the shared term overstates the SE of rare species and understates that of dominant ones.

## Counts

| rule | n |
| --- | --- |
| resident rows, origin controls and warmed plots (step 01) | 66,412 |
| experiments fitted | 43 |
| gradients | 21 |
| species | 348 |
| species x experiment (effects) | 995 |
| by_year contrasts (species x survey year) | 4,610 |
| strength_df rows | 995 |

## Strength

Quantiles (min / 25% / median / 75% / max). strength: -0.2965 / -0.0150 / -0.0034 / 0.0030 / 0.1885; strength_se: 0.0004 / 0.0047 / 0.0083 / 0.0140 / 0.0896. Effects significant at p < 0.05 (|strength| / SE above the experiment's t cut-off): 310 of 995.

![strength and strength_se](../figures/diagnostics/02_strength.png)

With one residual variance per species, the SE follows each species' own scatter: within an experiment the largest SE is a median of 6 times the smallest (range 1 to 91), and the Spearman correlation of SE with the species' mean relative cover in the experiment is 0.82. Rare species therefore get the smallest SEs on the absolute cover scale, and the largest weights in the later `rma.mv` calls.

## Fits per experiment

| experiment_id | rows | species | years | df | convergence | pdHess | sd_plot_species | se_min | se_median | se_max | warnings |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CH_Calanda::Cal::Nes | 5341 | 49 | 6 | 18 | 0 | TRUE | 0.00255 | 0.00117 | 0.00211 | 0.0274 | 1 |
| CH_Calanda::Pea::Cal | 3000 | 25 | 6 | 18 | 0 | TRUE | 0.0093 | 0.00418 | 0.00472 | 0.0243 | 1 |
| CH_Calanda2::Cal::Nes | 2880 | 36 | 4 | 18 | 0 | TRUE | 0.0132 | 0.00589 | 0.00598 | 0.0145 | 1 |
| CH_Lavey::RIO::CRE | 2340 | 39 | 3 | 18 | 0 | TRUE | 0.00144 | 0.000707 | 0.00298 | 0.0178 | 1 |
| CH_Lavey::RIO::MAR | 2340 | 39 | 3 | 18 | 0 | TRUE | 0.00216 | 0.00106 | 0.00332 | 0.019 | 1 |
| CH_Lavey::RIO::PRA | 2301 | 39 | 3 | 18 | 0 | TRUE | 0.000678 | 0.000428 | 0.00355 | 0.0207 | 1 |
| CN_Damxung::HIGH::LOW | 1600 | 20 | 4 | 18 | 0 | TRUE | 3.07e-07 | 0.000372 | 0.00472 | 0.0285 | 1 |
| CN_Gongga::A::M | 1344 | 24 | 4 | 12 | 0 | TRUE | 0.01 | 0.00545 | 0.00601 | 0.0328 | 1 |
| CN_Gongga::H::A | 800 | 16 | 4 | 12 | 0 | TRUE | 0.0141 | 0.00771 | 0.0103 | 0.0305 | 1 |
| CN_Gongga::H::L | 800 | 16 | 4 | 12 | 0 | TRUE | 0.015 | 0.00813 | 0.00996 | 0.0244 | 1 |
| CN_Gongga::M::L | 1008 | 18 | 4 | 12 | 0 | TRUE | 0.00125 | 0.00169 | 0.00892 | 0.0212 | 1 |
| CN_Heibei::3400::3200 | 1764 | 42 | 7 | 4 | 0 | TRUE | 0.000853 | 0.000924 | 0.00396 | 0.0178 | 1 |
| CN_Heibei::3800::3200 | 760 | 19 | 7 | 4 | 0 | TRUE | 0.000519 | 0.00129 | 0.00683 | 0.0214 | 1 |
| CN_Heibei::3800::3400 | 760 | 19 | 7 | 4 | 0 | TRUE | 0.0014 | 0.00158 | 0.00606 | 0.0196 | 1 |
| DE_Grainau::HIGH::LOW | 2480 | 31 | 4 | 18 | 0 | TRUE | 0.0129 | 0.00587 | 0.00684 | 0.0194 | 1 |
| DE_Susalps::EB::BT | 1872 | 26 | 4 | 16 | 0 | TRUE | 0.00193 | 0.00118 | 0.0113 | 0.0353 | 1 |
| DE_Susalps::EB::FE | 1872 | 26 | 4 | 16 | 0 | TRUE | 0.00136 | 0.000905 | 0.00691 | 0.0408 | 1 |
| DE_Susalps::EB::GW | 1872 | 26 | 4 | 16 | 0 | TRUE | 0.00177 | 0.00104 | 0.00681 | 0.0425 | 1 |
| DE_Susalps::FE::BT | 504 | 7 | 4 | 16 | 0 | TRUE | 0.118 | 0.0564 | 0.0639 | 0.0741 | 1 |
| DE_Susalps::GW::BT | 1368 | 19 | 4 | 16 | 0 | TRUE | 0.000433 | 0.000635 | 0.00974 | 0.0576 | 1 |
| DE_Susalps::GW::FE | 1368 | 19 | 4 | 16 | 0 | TRUE | 0.000571 | 0.000538 | 0.0132 | 0.0487 | 1 |
| DE_TransAlps::FP::BT | 1224 | 17 | 4 | 16 | 0 | TRUE | 0.0164 | 0.00781 | 0.00934 | 0.0643 | 1 |
| DE_TransAlps::SP::BT | 1656 | 23 | 4 | 16 | 0 | TRUE | 0.00509 | 0.00249 | 0.00783 | 0.0413 | 1 |
| FR_AlpeHuez::HIGH::LOW | 6840 | 38 | 9 | 18 | 0 | TRUE | 0.0126 | 0.00565 | 0.00621 | 0.0183 | 1 |
| FR_Lautaret::G::L | 2640 | 33 | 4 | 18 | 0 | TRUE | 0.00127 | 0.000762 | 0.00261 | 0.0218 | 1 |
| IT_MatschMazia1::High::Low | 825 | 25 | 3 | 10 | 0 | TRUE | 0.00989 | 0.00679 | 0.0119 | 0.0204 | 1 |
| IT_MatschMazia2::High::Low | 510 | 17 | 3 | 10 | 0 | TRUE | 0.00563 | 0.00388 | 0.00891 | 0.0376 | 1 |
| NO_Gudmedalen::Gudmedalen::Rambera | 1350 | 27 | 5 | 8 | 0 | TRUE | 0.0131 | 0.00836 | 0.00929 | 0.018 | 1 |
| NO_Gudmedalen::Rambera::Arhelleren | 1050 | 21 | 5 | 8 | 0 | TRUE | 0.0199 | 0.0127 | 0.0138 | 0.0306 | 1 |
| NO_Lavisdalen::Hogsete::Vikesland | 900 | 18 | 5 | 8 | 0 | TRUE | 0.0289 | 0.0183 | 0.0204 | 0.0317 | 1 |
| NO_Lavisdalen::Lavisdalen::Hogsete | 1550 | 31 | 5 | 8 | 0 | TRUE | 0.0134 | 0.00847 | 0.00895 | 0.0271 | 1 |
| NO_Skjellingahaugen::Skjellingahaugen::Veskre | 1000 | 20 | 5 | 8 | 0 | TRUE | 0.022 | 0.0139 | 0.0148 | 0.0422 | 1 |
| NO_Skjellingahaugen::Veskre::Ovstedal | 1600 | 32 | 5 | 8 | 0 | TRUE | 0.0121 | 0.00773 | 0.00878 | 0.0324 | 1 |
| NO_Ulvhaugen::Alrust::Fauske | 1029 | 21 | 5 | 8 | 0 | TRUE | 0.0127 | 0.00815 | 0.00937 | 0.0278 | 1 |
| NO_Ulvhaugen::Ulvhaugen::Alrust | 1450 | 29 | 5 | 8 | 0 | TRUE | 0.012 | 0.0077 | 0.00858 | 0.0306 | 1 |
| SE_Abisko::High::Low | 270 | 9 | 3 | 8 | 0 | TRUE | 0.0581 | 0.0376 | 0.0399 | 0.0486 | 1 |
| SE_Abisko::High::Mid | 270 | 9 | 3 | 8 | 0 | TRUE | 0.0448 | 0.0291 | 0.0313 | 0.0416 | 1 |
| SE_Abisko::Mid::Low | 330 | 11 | 3 | 8 | 0 | TRUE | 0.0484 | 0.031 | 0.0356 | 0.0636 | 1 |
| US_Arizona::MC::PP | 224 | 4 | 4 | 12 | 0 | TRUE | 0.137 | 0.0781 | 0.0834 | 0.0896 | 1 |
| US_Colorado::mo::pf | 944 | 16 | 5 | 10 | 0 | TRUE | 0.00394 | 0.00299 | 0.009 | 0.0286 | 1 |
| US_Colorado::mo::um | 1024 | 16 | 5 | 11 | 0 | TRUE | 0.00512 | 0.00318 | 0.00704 | 0.0466 | 1 |
| US_Colorado::pf::um | 960 | 16 | 5 | 10 | 0 | TRUE | 0.0293 | 0.0172 | 0.0201 | 0.0412 | 1 |
| US_Montana::High::Low | 392 | 7 | 3 | 18 | 0 | TRUE | 0.0277 | 0.0135 | 0.0177 | 0.0329 | 1 |

## Property coverage on the residents (added 26 Oct)

Species and effects with each property, after steps 01 and 02:

| property | species | effects |
| --- | --- | --- |
| LDMC | 246 | 835 |
| LMA | 260 | 855 |
| height | 231 | 800 |
| any of the three traits | 271 | 874 |
| growth form (graminoid, woody) | 348 | 995 |
| PEV1 to PEV3 | 348 | 995 |
| niche optima (any GBIF cell) | 348 | 995 |
| niche optima (n_obs_sub >= 30, used) | 339 | 979 |

Residents without any of the three traits (77): Agrostis nervosa, Alchemilla conjuncta, Alchemilla flabellata, Aletris pauciflora, Allium prattii, Androsace minor, Androsace tapete, Anemonastrum demissum, Argentina leuconota, Argentina stenophylla, Artemisia fulgens, Aster asteroides, Aster diplostephioides, Boechera stricta, Calamagrostis rubescens, Calamagrostis tripilifera, Carex cercostachys, Carex sargentiana, Carex scabrirostris, Carpesium humile, Cerastium szechuense, Chrysosplenium pilosum, Clinopodium polycephalum, Comastoma pulmonarium, Comastoma tenellum, Cyananthus incanus, Epilobium fangii, Erigeron neglectus, Eurybia conspicua, Fargesia spathacea, Galearis spathulata, Gentiana crassuloides, Gentiana prostrata, Geranium pylzowianum, Halenia corniculata, Hansenia weberbaueriana, Hedysarum algidum, Helictotrichon parlatorei, Helictotrichon tibeticum, Hemiphragma heterophyllum, Hypericum wightianum, Iris potaninii, Juncus himalensis, Juncus leucanthus, Koenigia cyanandra, Ligularia sagitta, Ligusticopsis scapiformis, Maianthemum japonicum, Melica spectabilis, Microula sikkimensis, Morina chinensis, Myosotis decumbens, Oxytropis densa, Oxytropis humifusa, Parasenecio latipes, Pedicularis davidii, Pedicularis sima, Persicaria runcinata, Poa reflexa, Polygonatum cirrhifolium, Potentilla pulcherrima, Ranunculus inamoenus, Rheum pumilum, Rhododendron websterianum, Sagina glabra, Salix serpillifolia, Saussurea ceterach, Saussurea globosa, Saussurea pulchra, Swertia macrosperma, Tanacetum tatsienense, Thesium alpinum, Vaccinium membranaceum, Veronica allionii, Veronica szechuanica, Viola nuttallii, Viola szetschwanensis.

## Covariates

Raw values before standardisation:

| variable | range |
| --- | --- |
| W (destination - origin JJA T, degC) | 0.3571 / 1.5285 / 2.3140 / 3.3592 / 11.0886 |
| T (destination JJA T, degC) | 8.7911 / 11.3936 / 13.2811 / 16.4864 / 21.0185 |
| P (log destination JJA precipitation, mm) | 5.0594 / 5.5447 / 5.9308 / 6.1053 / 6.9784 |
| A0_origin (log(% + 1)) | 0.0000 / 0.3617 / 0.9307 / 1.6384 / 3.8012 |
| A0_dest (log(% + 1)) | 0.0000 / 0.0000 / 0.0000 / 0.6917 / 3.8040 |

## Warnings captured

- CH_Calanda::Cal::Nes: NOTE: Results may be misleading due to involvement in interactions
- CH_Calanda::Pea::Cal: NOTE: Results may be misleading due to involvement in interactions
- CH_Calanda2::Cal::Nes: NOTE: Results may be misleading due to involvement in interactions
- CH_Lavey::RIO::CRE: NOTE: Results may be misleading due to involvement in interactions
- CH_Lavey::RIO::MAR: NOTE: Results may be misleading due to involvement in interactions
- CH_Lavey::RIO::PRA: NOTE: Results may be misleading due to involvement in interactions
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

None in the model: PLAN.md Step A was revised on 30 Oct (Bills) to one residual variance per species (`glmmTMB`, `dispformula = ~ species`). Decisions taken where PLAN.md is silent (30 Oct, Bills): degrees of freedom = origin-control plots + warmed plots - 2 per experiment; the plot x species variance kept shared (see Model as fitted); PEV1 to PEV3 carried; properties standardised across the residents that have them, traits logged first; graminoid and woody as 0/1. Pools follow the rule of 30 Oct (step 01), so the effects are 995 in 43 experiments rather than MILESTONES' original estimate of about 1,596 in 41.
