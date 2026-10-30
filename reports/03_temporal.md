# Step 03: Temporal consistency (Step B, residents)

Each resident species x experiment was classified at years 2, 3, 4 and 6 after transplanting from its year-specific warmed minus origin-control contrast (step 02, p-values with the plot-based df): Winner or Loser at p < 0.05 by the sign of the estimate, otherwise Winner (ns) or Loser (ns). A combination with no cover in any origin-control or warmed plot in that year was not classified there (decided 30 Oct). Between adjacent anchors, for the combinations present at both, the share keeping its sign and Lin's concordance correlation of the two contrasts were computed, with the 4 x 4 flow counts that draw Fig 1a. Year 1 is excluded.

## Counts per anchor

| year | combinations | absent_both_groups | classified | Loser | Loser (ns) | Winner (ns) | Winner |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 995 | 48 | 947 | 104 | 489 | 323 | 31 |
| 3 | 887 | 38 | 849 | 137 | 401 | 272 | 39 |
| 4 | 687 | 25 | 662 | 111 | 339 | 180 | 32 |
| 6 | 472 | 36 | 436 | 44 | 246 | 132 | 14 |

## Stability between anchors (`stability.csv`)

| pair | n | sign_agreement | ccc |
| --- | --- | --- | --- |
| 2 to 3 | 822 | 0.729 | 0.573 |
| 3 to 4 | 596 | 0.817 | 0.787 |
| 4 to 6 | 336 | 0.783 | 0.734 |

## Flows (`flows.rds`)

**Year 2 to year 3** (rows: class at the first anchor; columns: class at the second)

| from | Loser | Loser (ns) | Winner (ns) | Winner |
| --- | --- | --- | --- | --- |
| Loser | 48 | 30 | 9 | 0 |
| Loser (ns) | 76 | 250 | 94 | 7 |
| Winner (ns) | 10 | 101 | 147 | 22 |
| Winner | 0 | 2 | 17 | 9 |

**Year 3 to year 4** (rows: class at the first anchor; columns: class at the second)

| from | Loser | Loser (ns) | Winner (ns) | Winner |
| --- | --- | --- | --- | --- |
| Loser | 55 | 42 | 1 | 0 |
| Loser (ns) | 42 | 200 | 32 | 3 |
| Winner (ns) | 4 | 68 | 106 | 16 |
| Winner | 0 | 1 | 16 | 10 |

**Year 4 to year 6** (rows: class at the first anchor; columns: class at the second)

| from | Loser | Loser (ns) | Winner (ns) | Winner |
| --- | --- | --- | --- | --- |
| Loser | 12 | 21 | 1 | 0 |
| Loser (ns) | 17 | 132 | 24 | 0 |
| Winner (ns) | 1 | 45 | 64 | 6 |
| Winner | 0 | 2 | 8 | 3 |

## Fig 1 panel a (draft)

![Fig 1a draft](../figures/diagnostics/03_fig1_residents_draft.png)

Strata are drawn to one absolute scale, so columns shrink as experiments end and as species are absent from both groups. The final figure, with panel b (colonizers, step 07), is drawn in step 10. Columns are evenly spaced: the gap from year 4 to year 6 spans two years, the others one.

Switches between a significant winner and a significant loser, in either direction: 2 to 3: 0; 3 to 4: 0; 4 to 6: 0. Ribbons that cross zero with a significant class at either end are drawn darker; every ribbon that crosses zero is labelled where it enters the next anchor with its count and its share of the class it started from.

## Warnings captured

- The `label.size` argument of `geom_label()` is deprecated as of ggplot2 3.5.0.
ℹ Please use the `linewidth` argument instead.

## Deviations from PLAN.md

Combinations absent from both groups in an anchor year ( 147  of  3,041  anchor contrasts) are not classified there, where PLAN.md's rule would have put their zero contrast in Loser (ns) (decided 30 Oct). Anchors are years 2, 3, 4 and 6 rather than 2, 4 and 6 (decided 30 Oct: agreement between consecutive years rises from year 3 on; PLAN.md updated), labelled Year 2, 3, 4, 6 rather than the mock's First 2, 4, 6 years. The draft figure is written to figures/diagnostics/ rather than data/derived/03_temporal/, so that it is committed with the report.
