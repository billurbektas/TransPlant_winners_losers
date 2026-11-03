# Step 07c: The colonizer ladder, fitted (Step D route 1, colonizers)

The residents' ladder (step 05) was fitted to the colonizers (Bills, 2 Nov): warmed cover from years 3 and later (step 07a) for the 110 effects of the 27 colonizer species found in at least 3 gradients that carry every property (24 experiments, 18 gradients). The same 12 models, with A0_origin dropped (a colonizer has none) so that Start is the destination abundance A0_dest and its interaction with W. Properties are on the residents' step 02 scale (checked), A0_dest (log(100 x mean destination-control cover + 1)) standardised across all 379 colonizer effects of 07a. 24 effects are measured absences (SE 0).

## Counts

| rule | n |
| --- | --- |
| colonizer effects (07a, years 3 and later) | 379 |
| species in >= 3 gradients | 29 |
| of which lacking a property | 2 |
| ladder effects | 110 |
| species | 27 |
| experiments | 24 |
| gradients | 18 |
| measured absences among them | 24 |

Species lacking a property: Euphrasia officinalis, Pilosella piloselloides.

Growth form on these rows: 1 woody and 6 graminoid species of 27; with one woody species metafor drops the redundant term TP:woody from the Woody model, which rests on that species alone.

## Heterogeneity per model (×10⁻⁴)

| model | τ² | σ²_experiment | H intercept | H slope | optimizer | minutes |
| --- | --- | --- | --- | --- | --- | --- |
| Null | 0.207 | 0.380 | 0.040 | 0.000 | nlminb | 0.0 |
| Start | 0.205 | 0.307 | 0.028 | 0.000 | nlminb | 0.0 |
| Climate | 0.211 | 0.842 | 0.014 | 0.000 | nlminb | 0.0 |
| Climate + start | 0.213 | 0.594 | 0.004 | 0.000 | nlminb | 0.0 |
| Species-dependent | 0.121 | 0.616 | 0.003 | 0.303 | nlminb | 0.1 |
| LDMC | 0.109 | 0.759 | 0.000 | 0.244 | nlminb | 0.1 |
| LMA | 0.121 | 0.613 | 0.000 | 0.378 | nlminb | 0.1 |
| Height | 0.059 | 0.948 | 0.000 | 0.255 | nlminb | 0.1 |
| Graminoid | 0.114 | 0.709 | 0.002 | 0.339 | nlminb | 0.1 |
| Woody | 0.138 | 0.625 | 0.000 | 0.316 | nlminb | 0.1 |
| Phylogeny | 0.129 | 0.688 | 0.000 | 0.448 | nlminb | 0.1 |
| Niche | 0.116 | 0.738 | 0.017 | 0.345 | nlminb | 0.0 |

## Step 05's checks, reported

τ² at Climate against Null: +2.2% (step 05 stops beyond 5%).

| property | τ² move | σ²_experiment move |
| --- | --- | --- |
| LDMC | -1.5% | +19.5% |
| LMA | +0.0% | -0.3% |
| Height | -8.4% | +45.1% |
| Graminoid | -0.9% | +12.6% |
| Woody | +2.3% | +1.2% |
| Phylogeny | +1.1% | +9.8% |
| Niche | -0.6% | +16.6% |

Moves beyond step 05's 5% limit: 6 of 14. On 110 effects the variance components are not pinned down well enough for that check to hold; they are reported, not enforced (Bills, 2 Nov).

## Coefficients whose 95% interval excludes zero

| model | term | estimate | ci |
| --- | --- | --- | --- |
| Start | A0_dest | 0.0039 | 0.0023 to 0.0055 |
| Start | A0_dest:W | -0.0009 | -0.0018 to -0.0000 |
| Climate | P | 0.0068 | 0.0001 to 0.0136 |
| Climate + start | P | 0.0083 | 0.0023 to 0.0143 |
| Climate + start | A0_dest | 0.0044 | 0.0027 to 0.0061 |
| Climate + start | W:A0_dest | -0.0011 | -0.0020 to -0.0002 |
| Species-dependent | P | 0.0082 | 0.0016 to 0.0148 |
| Species-dependent | A0_dest | 0.0043 | 0.0027 to 0.0059 |
| Species-dependent | W:A0_dest | -0.0013 | -0.0024 to -0.0002 |
| LDMC | T | 0.0115 | 0.0014 to 0.0217 |
| LDMC | P | 0.0098 | 0.0027 to 0.0169 |
| LDMC | A0_dest | 0.0044 | 0.0028 to 0.0059 |
| LDMC | W:A0_dest | -0.0014 | -0.0024 to -0.0004 |
| LDMC | W:ldmc_z | -0.0015 | -0.0029 to -0.0001 |
| LMA | A0_dest | 0.0045 | 0.0029 to 0.0062 |
| LMA | W:A0_dest | -0.0016 | -0.0030 to -0.0002 |
| Height | A0_dest | 0.0029 | 0.0016 to 0.0041 |
| Height | W:height_z | -0.0036 | -0.0058 to -0.0013 |
| Height | T:height_z | 0.0080 | 0.0044 to 0.0116 |
| Height | P:height_z | 0.0051 | 0.0004 to 0.0097 |
| Height | height_z:TP | -0.0062 | -0.0113 to -0.0010 |
| Graminoid | P | 0.0074 | 0.0000 to 0.0148 |
| Graminoid | A0_dest | 0.0040 | 0.0024 to 0.0056 |
| Woody | P | 0.0087 | 0.0020 to 0.0155 |
| Woody | A0_dest | 0.0050 | 0.0032 to 0.0068 |
| Woody | W:A0_dest | -0.0015 | -0.0027 to -0.0003 |
| Phylogeny | A0_dest | 0.0056 | 0.0036 to 0.0076 |
| Phylogeny | W:A0_dest | -0.0026 | -0.0044 to -0.0008 |
| Phylogeny | W:PEV1 | -0.0061 | -0.0105 to -0.0018 |
| Phylogeny | T:PEV1 | 0.0103 | 0.0016 to 0.0190 |
| Niche | A0_dest | 0.0046 | 0.0029 to 0.0064 |

## Figures

![Fig S16, colonizer coefficients](../figures/appendix/figS16.png)

## Warnings captured

- Null: Model does not contain an '~ inner | outer' term, so 'struct' argument is disregaded.
- Null: There are outcomes with non-positive sampling variances.
- Null: 'V' appears to be not positive definite.
- Start: Model does not contain an '~ inner | outer' term, so 'struct' argument is disregaded.
- Start: There are outcomes with non-positive sampling variances.
- Start: 'V' appears to be not positive definite.
- Climate: Model does not contain an '~ inner | outer' term, so 'struct' argument is disregaded.
- Climate: There are outcomes with non-positive sampling variances.
- Climate: 'V' appears to be not positive definite.
- Climate + start: Model does not contain an '~ inner | outer' term, so 'struct' argument is disregaded.
- Climate + start: There are outcomes with non-positive sampling variances.
- Climate + start: 'V' appears to be not positive definite.
- Species-dependent: There are outcomes with non-positive sampling variances.
- Species-dependent: 'V' appears to be not positive definite.
- LDMC: There are outcomes with non-positive sampling variances.
- LDMC: 'V' appears to be not positive definite.
- LMA: There are outcomes with non-positive sampling variances.
- LMA: 'V' appears to be not positive definite.
- Height: There are outcomes with non-positive sampling variances.
- Height: 'V' appears to be not positive definite.
- Graminoid: There are outcomes with non-positive sampling variances.
- Graminoid: 'V' appears to be not positive definite.
- Woody: There are outcomes with non-positive sampling variances.
- Woody: 'V' appears to be not positive definite.
- Woody: Redundant predictors dropped from the model.
- Phylogeny: There are outcomes with non-positive sampling variances.
- Phylogeny: 'V' appears to be not positive definite.
- Niche: There are outcomes with non-positive sampling variances.
- Niche: 'V' appears to be not positive definite.
- `geom_errorbarh()` was deprecated in ggplot2 4.0.0.
ℹ Please use the `orientation` argument of `geom_errorbar()` instead.

## Deviations from PLAN.md

- The colonizer ladder on the residents' reduction and models (Bills, 2 Nov); properties on the residents' scale.
- Step 05's tau2 and property checks are reported, not enforced (Bills, 2 Nov: run without stopping).
- The colonizer coefficients are their own appendix figure (Fig S16), beside the residents' Fig S12.
