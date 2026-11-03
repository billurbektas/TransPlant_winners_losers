# Step 07a: Colonizer strength and presence flows (Steps A and B, colonizers)

Each of the 28 experiments with destination-pool species got one mixed model of relative cover on its warmed plots only (Bills, 2 Nov), with the residents' model terms for one group: glmmTMB, warmed cover by year and species, plot x species random intercept, one residual variance per species, degrees of freedom warmed plots - 1. Colonizer strength is the warmed marginal mean, zeros kept, averaged over all years and, separately, over years 3 and later; A0_dest, the propagule pressure, is the species' destination-control cover averaged over every plot and year (no model), as log(100 x cover + 1) like the residents' A0 (Bills, 3 Nov). Both strength windows are kept until the choice is made (Bills, 2 Nov). Each colonizer was classified at years 2, 3, 4 and 6 as absent (no cover in any warmed plot that year, from the raw data), present (warmed mean above zero at one-sided p < 0.05) or present (ns), and presence agreement and the concordance correlation of warmed cover were computed between adjacent anchors.

## Counts

| rule | n |
| --- | --- |
| colonizer effects in the destination pool (species x experiment) | 402 |
| experiments dropped by the pool-size rule (fewer than 5 colonizer species) | 11 |
| their effects | 23 |
| colonizer effects kept | 379 |
| species | 229 |
| experiments | 28 |
| gradients | 19 |
| effects with a strength from years 3 and later | 379 |
| species never present in a warmed plot (any gradient) | 44 |

Dropped by the pool-size rule, as for residents (Bills, 2 Nov): CN_Damxung::HIGH::LOW (3), CN_Gongga::M::L (1), DE_Susalps::EB::BT (1), DE_Susalps::EB::FE (2), DE_Susalps::EB::GW (3), DE_Susalps::GW::FE (1), DE_TransAlps::SP::BT (2), NO_Lavisdalen::Hogsete::Vikesland (2), NO_Skjellingahaugen::Veskre::Ovstedal (1), SE_Abisko::Mid::Low (3), US_Colorado::mo::pf (4); CN_Damxung leaves with its only experiment.

Where each colonizer's strength SE comes from (all years): measured absence 107, species 272. 'species': the warmed-only model with one residual variance per species; 'measured absence': zero cover in every warmed plot, strength 0 and SE 0 from the data; 'experiment': the species' variance could not be told apart from its plot variance (), so its estimate and SE come from the same model with one residual variance for the experiment.

0 effects have no strength from years 3 and later (their experiment has no survey in year 3 or later).

## Fits

| experiment_id | species | years | df_plots | warnings |
| --- | --- | --- | --- | --- |
| CH_Calanda::Cal::Nes | 23 | 1, 2, 3, 4, 5, 6 | 9 | 0 |
| CH_Calanda::Pea::Cal | 30 | 1, 2, 3, 4, 5, 6 | 9 | 0 |
| CH_Calanda2::Cal::Nes | 10 | 1, 2, 3, 4 | 9 | 0 |
| CH_Lavey::RIO::CRE | 2 | 1, 2, 3 | 9 | 0 |
| CH_Lavey::RIO::MAR | 9 | 1, 2, 3 | 9 | 0 |
| CH_Lavey::RIO::PRA | 13 | 1, 2, 3 | 9 | 0 |
| CN_Gongga::A::M | 7 | 1, 2, 3, 4 | 6 | 0 |
| CN_Gongga::H::A | 8 | 1, 2, 3, 4 | 6 | 0 |
| CN_Gongga::H::L | 7 | 1, 2, 3, 4 | 6 | 0 |
| CN_Heibei::3800::3200 | 8 | 1, 2, 3, 6, 7, 8, 9 | 2 | 0 |
| CN_Heibei::3800::3400 | 10 | 1, 2, 3, 6, 7, 8, 9 | 2 | 0 |
| DE_Grainau::HIGH::LOW | 11 | 1, 2, 3, 4 | 9 | 0 |
| DE_Susalps::FE::BT | 6 | 1, 2, 3, 4 | 8 | 0 |
| DE_TransAlps::FP::BT | 8 | 1, 2, 3, 4 | 8 | 0 |
| FR_AlpeHuez::HIGH::LOW | 13 | 1, 2, 3, 4, 5, 6, 7, 8, 9 | 9 | 0 |
| FR_Lautaret::G::L | 9 | 1, 2, 5, 6 | 9 | 0 |
| IT_MatschMazia1::High::Low | 4 | 1, 2, 3 | 2 | 0 |
| IT_MatschMazia2::High::Low | 15 | 1, 2, 3 | 8 | 0 |
| NO_Gudmedalen::Rambera::Arhelleren | 8 | 2, 3, 4, 6, 8 | 4 | 0 |
| NO_Lavisdalen::Lavisdalen::Hogsete | 7 | 2, 3, 4, 6, 8 | 4 | 0 |
| NO_Skjellingahaugen::Skjellingahaugen::Veskre | 7 | 2, 3, 4, 6, 8 | 4 | 0 |
| NO_Ulvhaugen::Alrust::Fauske | 6 | 2, 3, 4, 6, 8 | 4 | 0 |
| NO_Ulvhaugen::Ulvhaugen::Alrust | 12 | 2, 3, 4, 6, 8 | 4 | 0 |
| SE_Abisko::High::Low | 4 | 1, 2, 3 | 4 | 0 |
| SE_Abisko::High::Mid | 3 | 1, 2, 3 | 4 | 0 |
| US_Colorado::mo::um | 10 | 1, 2, 4, 5, 6 | 7 | 0 |
| US_Colorado::pf::um | 5 | 1, 2, 4, 5, 6 | 6 | 0 |
| US_Montana::High::Low | 17 | 1, 2, 5 | 9 | 0 |

## Presence classes per anchor

| year | combinations | Present | Present (ns) | Absent |
| --- | --- | --- | --- | --- |
| 2 | 379 | 39 | 138 | 202 |
| 3 | 331 | 59 | 119 | 153 |
| 4 | 205 | 57 | 92 | 56 |
| 6 | 166 | 56 | 70 | 40 |

## Stability between anchors (`stability.csv`)

| pair | n | presence_agreement | ccc |
| --- | --- | --- | --- |
| 2 to 3 | 331 | 0.770 | 0.514 |
| 3 to 4 | 190 | 0.784 | 0.772 |
| 4 to 6 | 130 | 0.831 | 0.575 |

Residents (step 03), for the year-3 read: 2 to 3 sign agreement 0.75, CCC 0.54; 3 to 4 sign agreement 0.83, CCC 0.78; 4 to 6 sign agreement 0.79, CCC 0.71.

All years against years 3 and later, on the 379 effects with both: Pearson r 0.917, concordance correlation 0.861.

![Stability](../figures/appendix/figS15.png)

## Fig 1

![Fig 1](../figures/main/fig1.png)

Panel a is step 03's, redrawn from its saved classes and flows with the same code; panel b the colonizers below it (Bills, 2 Nov). Ribbons between absent and present are labelled where they enter the next anchor, with their count and their share of the class they started from; those that start or end in Present are drawn darker.

## Never present

44 colonizer species never appear in a warmed plot of any experiment where they are in the pool; 42 of them occur in one gradient and none in three or more, so none enters Step C (`never_present.csv`).

## Warnings captured

- NOTE: Results may be misleading due to involvement in interactions
- The `label.size` argument of `geom_label()` is deprecated as of ggplot2 3.5.0.
ℹ Please use the `linewidth` argument instead.

## Deviations from PLAN.md

- Strength from a model of the warmed plots only, and A0_dest the raw destination-control cover averaged over years, not the destination-control mean of a joint fit; species never recorded in a warmed plot are measured absences (strength 0, SE 0), and species whose variance the model cannot separate from their plot variance take the experiment's error (Bills, 2 Nov; PLAN.md updated).
- The pool-size rule applied to colonizers too, at least 5 colonizer species per experiment (Bills, 2 Nov; PLAN.md updated).
- Strength kept twice, all years and years 3 and later, until the stability is read (Bills, 2 Nov).
- Presence is tested one-sided (warmed mean above zero), as PLAN.md's colonizer code now writes it.
- Figures written to `figures/main/` and `figures/appendix/` so that the report can link them.
