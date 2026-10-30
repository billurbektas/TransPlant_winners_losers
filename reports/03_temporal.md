# Step 03: Temporal consistency (Step B, residents)

Each resident species x experiment was classified at years 2, 3, 4 and 6 after transplanting from its year-specific warmed minus origin-control contrast (step 02, p-values with the plot-based df): Winner or Loser at p < 0.05 by the sign of the estimate, otherwise Winner (ns) or Loser (ns). A combination with no cover in any origin-control or warmed plot in that year was not classified there (decided 30 Oct). Between adjacent anchors, for the combinations present at both, the share keeping its sign and Lin's concordance correlation of the two contrasts were computed, with the 4 x 4 flow counts that draw Fig 1a. Year 1 is excluded.

## Counts per anchor

| year | combinations | absent_both_groups | classified | Loser | Loser (ns) | Winner (ns) | Winner |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 1655 | 232 | 1423 | 128 | 765 | 484 | 46 |
| 3 | 1453 | 174 | 1279 | 161 | 662 | 404 | 52 |
| 4 | 1150 | 125 | 1025 | 131 | 570 | 279 | 45 |
| 6 | 766 | 132 | 634 | 61 | 349 | 202 | 22 |

## Stability between anchors (`stability.csv`)

| pair | n | sign_agreement | ccc |
| --- | --- | --- | --- |
| 2 to 3 | 1182 | 0.740 | 0.591 |
| 3 to 4 | 896 | 0.826 | 0.778 |
| 4 to 6 | 486 | 0.802 | 0.710 |

## Flows (`flows.rds`)

**Year 2 to year 3** (rows: class at the first anchor; columns: class at the second)

| from | Loser | Loser (ns) | Winner (ns) | Winner |
| --- | --- | --- | --- | --- |
| Loser | 53 | 40 | 9 | 0 |
| Loser (ns) | 84 | 409 | 131 | 8 |
| Winner (ns) | 13 | 140 | 224 | 30 |
| Winner | 0 | 6 | 22 | 13 |

**Year 3 to year 4** (rows: class at the first anchor; columns: class at the second)

| from | Loser | Loser (ns) | Winner (ns) | Winner |
| --- | --- | --- | --- | --- |
| Loser | 58 | 52 | 1 | 0 |
| Loser (ns) | 52 | 356 | 51 | 3 |
| Winner (ns) | 4 | 96 | 161 | 25 |
| Winner | 0 | 1 | 24 | 12 |

**Year 4 to year 6** (rows: class at the first anchor; columns: class at the second)

| from | Loser | Loser (ns) | Winner (ns) | Winner |
| --- | --- | --- | --- | --- |
| Loser | 16 | 22 | 1 | 0 |
| Loser (ns) | 24 | 201 | 32 | 1 |
| Winner (ns) | 1 | 57 | 102 | 8 |
| Winner | 0 | 4 | 13 | 4 |

## Fig 1 panel a (draft)

![Fig 1a draft](../figures/diagnostics/03_fig1_residents_draft.png)

Strata are drawn to one absolute scale, so columns shrink as experiments end and as species are absent from both groups. The final figure, with panel b (colonizers, step 07), is drawn in step 10. Columns are evenly spaced: the gap from year 4 to year 6 spans two years, the others one.

Switches between a significant winner and a significant loser, in either direction: 2 to 3: 0; 3 to 4: 0; 4 to 6: 0. Ribbons that cross zero with a significant class at either end are drawn darker; every ribbon that crosses zero is labelled where it enters the next anchor with its count and its share of the class it started from.

## Warnings captured

- The `label.size` argument of `geom_label()` is deprecated as of ggplot2 3.5.0.
ℹ Please use the `linewidth` argument instead.

## Deviations from PLAN.md

Combinations absent from both groups in an anchor year ( 663  of  5,024  anchor contrasts) are not classified there, where PLAN.md's rule would have put their zero contrast in Loser (ns) (decided 30 Oct). Anchors are years 2, 3, 4 and 6 rather than 2, 4 and 6 (decided 30 Oct: agreement between consecutive years rises from year 3 on; PLAN.md updated), labelled Year 2, 3, 4, 6 rather than the mock's First 2, 4, 6 years. The draft figure is written to figures/diagnostics/ rather than data/derived/03_temporal/, so that it is committed with the report.
