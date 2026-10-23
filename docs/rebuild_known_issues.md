# Deferred diagnostic from the October 2026 rebuild

The LMA full-model profile-likelihood interval for the species-specific
precipitation-slope variance has a fitted estimate and an upper 95% bound, but
`metafor::confint()` could not obtain its lower bound because the profile
refits did not converge. The affected row is `LMA` / `species_P` in
`output/tables/model_random_effect_profile_intervals.csv`.

This is **not** evidence that the interval includes zero or that it excludes
zero. Its interval-based support is undetermined. The point estimate, fitted
model, and other 34 profile intervals remain available; the appendix figure
shows this point without a confidence bar and identifies the incomplete
interval. The release check allows only this exact, warning-bearing incomplete
row. Another missing or invalid profile interval still fails the check.

The cause and scientific value of these predictor variance profiles are deferred
for review. This exception does not repair the interval or authorize a claim
about support for the LMA precipitation-slope variance. It should be revisited
before treating that appendix panel as a fully resolved inferential result.
