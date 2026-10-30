# Step 03: Temporal consistency (Step B, residents)

Each resident species x experiment was classified at years 2, 4 and 6 after transplanting from its year-specific warmed minus origin-control contrast (step 02, p-values with the plot-based df): Winner or Loser at p < 0.05 by the sign of the estimate, otherwise Winner (ns) or Loser (ns). A combination with no cover in any origin-control or warmed plot in that year was not classified there (decided 30 Oct). Between adjacent anchors, for the combinations present at both, the share keeping its sign and Lin's concordance correlation of the two contrasts were computed, with the 4 x 4 flow counts that draw Fig 1a. Year 1 is excluded.

## Counts per anchor

| year | combinations | absent_both_groups | classified | Loser | Loser (ns) | Winner (ns) | Winner |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 1655 | 232 | 1423 | 128 | 765 | 484 | 46 |
| 4 | 1150 | 125 | 1025 | 131 | 570 | 279 | 45 |
| 6 | 766 | 132 | 634 | 61 | 349 | 202 | 22 |

## Stability between anchors (`stability.csv`)

| pair | n | sign_agreement | ccc |
| --- | --- | --- | --- |
| 2 to 4 | 951 | 0.752 | 0.454 |
| 4 to 6 | 486 | 0.802 | 0.710 |

## Flows (`flows.rds`)

**Year 2 to year 4** (rows: class at the first anchor; columns: class at the second)

| from | Loser | Loser (ns) | Winner (ns) | Winner |
| --- | --- | --- | --- | --- |
| Loser | 36 | 42 | 1 | 0 |
| Loser (ns) | 73 | 350 | 87 | 9 |
| Winner (ns) | 10 | 122 | 164 | 29 |
| Winner | 1 | 6 | 16 | 5 |

**Year 4 to year 6** (rows: class at the first anchor; columns: class at the second)

| from | Loser | Loser (ns) | Winner (ns) | Winner |
| --- | --- | --- | --- | --- |
| Loser | 16 | 22 | 1 | 0 |
| Loser (ns) | 24 | 201 | 32 | 1 |
| Winner (ns) | 1 | 57 | 102 | 8 |
| Winner | 0 | 4 | 13 | 4 |

## Fig 1 panel a (draft)

![Fig 1a draft](../figures/diagnostics/03_fig1_residents_draft.png)

Strata are drawn to one absolute scale, so columns shrink as experiments end and as species are absent from both groups. Ribbon labels give the combinations that crossed from the loser side (Loser, Loser (ns)) to the winner side or back, with their share of the starting side. The final figure, with panel b (colonizers, step 07), is drawn in step 10.

## Warnings captured

- The `label.size` argument of `geom_label()` is deprecated as of ggplot2 3.5.0.
ℹ Please use the `linewidth` argument instead.

## Deviations from PLAN.md

Combinations absent from both groups in an anchor year ( 489  of  3,571  anchor contrasts) are not classified there, where PLAN.md's rule would have put their zero contrast in Loser (ns) (decided 30 Oct). Anchors are labelled Year 2, 4, 6 rather than the mock's First 2, 4, 6 years. The draft figure is written to figures/diagnostics/ rather than data/derived/03_temporal/, so that it is committed with the report.
