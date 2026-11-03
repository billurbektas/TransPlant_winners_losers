# Step 12: final review of code, runs, documentation, open issues and statistics.
#
# PLAN.md: the whole document (consistency); CLAUDE.md (rules).
# Reads:  data/derived/ (every step's outputs and cached fits, read only, nothing refitted);
#         data/derived/prep/ (prep outputs); prep/_targets (store, for the TRY name map)
# Writes: reports/12_final_review.md
# Checks: none that stop; every number in the report is computed here from the saved outputs.
#
# Decided 3 Nov (Bills): a critical final pass; statistical choices are not changed, only
# reported. The code findings come from a line-by-line read of every script (three read-only
# reviews, each finding re-checked before it was kept); the quantities below are the
# evidence. Runs that cannot be repeated inside one script (the fresh prep rebuild and the
# two analysis runs in a clean clone without legacy/) are summarised in the report text.

suppressPackageStartupMessages({library(dplyr); library(metafor)})
source("R/md_table.R")
warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
nn  <- function(x) format(x, big.mark = ",")

res   <- readRDS("data/derived/03b_strength_y3/strength_df.rds")
col   <- readRDS("data/derived/07a_colonizers/strength_df_y3.rds")
s04   <- read.csv("data/derived/04_classification/species_summary.csv")
s07b  <- read.csv("data/derived/07b_classification/species_summary.csv")
cv    <- readRDS("data/derived/01_relative_cover/cover_long.rds") |> mutate(year = as.integer(as.character(year)))
a0    <- readRDS("data/derived/01_relative_cover/a0.rds")
exps  <- readRDS("data/derived/prep/experiments.rds") |> filter(experiment_id %in% res$experiment_id)
lad   <- readRDS("data/derived/05_ladder_fit/models/all/null.rds"); lad <- if (inherits(lad, "rma")) lad else lad$m
p06   <- read.csv("data/derived/06_ladder_predict/pred_metrics.csv")
p07   <- read.csv("data/derived/07d_ladder_predict/pred_metrics.csv")
f07   <- read.csv("data/derived/07e_forest/forest_vs_metafor.csv")

# ---- Cached fits against today's inputs (no refit) ----------------------------------
rma_ok <- function(dir, rows) sapply(list.files(dir, recursive = TRUE, full.names = TRUE), function(f) {
  m <- readRDS(f); m <- if (inherits(m, "rma")) m else m$m
  cur <- rows[match(m$data$effect_id, rows$effect_id), ]
  !anyNA(cur$effect_id) && max(abs(as.numeric(m$yi) - cur$strength)) == 0 && max(abs(as.numeric(m$vi) - cur$strength_se^2)) == 0 })
ok05 <- rma_ok("data/derived/05_ladder_fit/models", res)
ok07 <- rma_ok("data/derived/07c_ladder_fit/models", readRDS("data/derived/07c_ladder_fit/ladder_rows.rds"))
fold_ok <- function(dir, rows) sapply(list.files(dir, full.names = TRUE, recursive = TRUE), function(f) {
  p <- readRDS(f)$pred; if (is.null(p)) return(NA); cur <- rows$strength[match(p$effect_id, rows$effect_id)]; !anyNA(cur) && max(abs(cur - p$observed)) == 0 })
fo06 <- fold_ok("data/derived/06_ladder_predict/folds", res)
fo07 <- fold_ok("data/derived/07d_ladder_predict/folds", readRDS("data/derived/07c_ladder_fit/ladder_rows.rds"))

# ---- Colonizer Step C with measured absences (SE 0) ---------------------------------
zero_se <- col |> filter(species %in% s07b$species[!is.na(s07b$mu)]) |> group_by(species) |> summarise(n_zero = sum(strength_se == 0), .groups = "drop") # the modelled colonizers
reml_gap <- bind_rows(lapply(c("Lotus corniculatus", "Phleum pratense"), function(sp) {
  d <- col[col$species == sp, ]
  m  <- rma.mv(strength, strength_se^2, random = ~ 1 | effect_id, data = d, method = "REML")
  m0 <- rma.mv(strength, strength_se^2, random = ~ 1 | effect_id, data = d, method = "REML", sigma2 = 1e-12)
  tibble::tibble(species = sp, `effects with SE 0` = sum(d$strength_se == 0), `reported tau2` = signif(m$sigma2, 3),
                 `REML log-lik at reported tau2` = fmt(as.numeric(logLik(m))), `REML log-lik at tau2 = 1e-12` = fmt(as.numeric(logLik(m0))))
}))
tiny_mu <- s07b |> filter(!is.na(mu), mu > 0, mu < 1e-6)

# ---- Data handling ----------------------------------------------------------------------
cal2 <- cv |> filter(experiment_id == "CH_Calanda2::Cal::Nes") |> summarise(years = paste(sort(unique(year)), collapse = ", "), .by = species_pool)
survey_years <- cv |> distinct(experiment_id, year) |> count(experiment_id, name = "surveyed") |>
  left_join(transmute(exps, experiment_id, span = last_year - first_year + 1L), by = "experiment_id")
gaps <- filter(survey_years, surveyed < span)
nm <- targets::tar_read(try_name_map, store = "prep/_targets")
ts <- targets::tar_read(trait_summary, store = "prep/_targets")
fb <- ts$long[!is.na(ts$long$name_source) & ts$long$name_source == "AccSpeciesName", ]
source("prep/R/taxonomy.R"); source("prep/R/traits.R")
fb_rejected <- sum(try_key(fb$SpeciesName) %in% try_key(nm$name[!nm$used]))

# ---- Statistical quantities -----------------------------------------------------------
shared <- exps |> mutate(n_o = n(), .by = c(region, origin_site))
lad_shared <- lad$data |> left_join(select(exps, experiment_id, region, origin_site), by = "experiment_id") |>
  count(species, region, origin_site) |> filter(n > 1)
oc3 <- cv |> filter(species_pool == "resident", treatment == "origin_control", year >= 3) |> summarise(oc3 = mean(relative_cover), .by = c(experiment_id, species))
cpl <- res |> select(experiment_id, species, strength) |> left_join(oc3, by = c("experiment_id", "species")) |> left_join(a0, by = c("experiment_id", "species"))
clim <- readRDS("data/derived/prep/experiment_climate.rds") |> filter(experiment_id %in% res$experiment_id) |> left_join(distinct(res, experiment_id, gradient_id), by = "experiment_id")
leak <- sapply(c(W = "W", T = "T", P = "P_summer"), function(v) { x <- clim[[v]]
  max(sapply(unique(clim$gradient_id), function(g) abs(mean(x[clim$gradient_id != g]) - mean(x)) / sd(x))) })
q_t <- qt(0.975, s04$n_experiments - 1); hw <- q_t * sqrt(s04$tau^2 + s04$se_mu^2)
cls_t <- ifelse(s04$mu + hw < 0, "Consistent loser", ifelse(s04$mu - hw > 0, "Consistent winner", "Context-dependent"))
w_c <- 1 / (res$strength_se[res$species %in% s04$species])^2; w_c <- sort(w_c, decreasing = TRUE)
by_year02 <- readRDS("data/derived/02_strength/by_year.rds") |> mutate(year = as.integer(as.character(year)))
ny <- by_year02 |> filter(year >= 3) |> distinct(experiment_id, year) |> count(experiment_id, name = "n_years")
win <- res |> left_join(ny, by = "experiment_id") |> mutate(window = ifelse(n_years >= 4, "4 or more", as.character(n_years))) |>
  summarise(effects = n(), `median SE` = signif(median(strength_se), 2), `mean |strength|` = signif(mean(abs(strength)), 2), .by = window) |> arrange(window)
spread <- function(p, lv) { r <- p$r2[p$level == lv]; paste0(fmt(median(r)), " (", fmt(quantile(r, .25)), " to ", fmt(quantile(r, .75)), "; above 0 in ", sum(r > 0), " of ", length(r), ")") }
woody7 <- p07$r2[p07$level == "woody"]
col_null <- median(f07$forest[f07$level == "null" & f07$metric == "r2"])
col_clim <- readRDS("data/derived/07c_ladder_fit/ladder_rows.rds")

# ---- Report ---------------------------------------------------------------------------------
report <- c(
  "# Step 12: Final review", "",
  "A critical final pass over the whole repository (Bills, 3 Nov): every script of the prep pipeline and the analysis read line by line; everything run that does not fit a model or download data; the documentation checked against the code; every open item re-classified; and statistical weak points reported without changing any choice. `legacy/` was retired. What was fixed changes no result; what would change a result is listed for Bills under *Decisions*.",
  "", "## Runs", "",
  md_table(tibble::tribble(
    ~check, ~outcome,
    "prep pipeline rebuilt from an empty store in a clean clone without `legacy/`", "all 75 targets ran; every file in `data/derived/prep/` and every P report byte-identical to the working copy",
    "all 15 analysis steps run from the cached fits in the clean clone", "every step finished; every tracked output byte-identical, apart from figures whose two legends swapped order between runs (Figs 2, 3, 4, S10; fixed here)",
    "two consecutive full runs in the clean clone, then the figure steps (04, 07b, 10) three more times", "identical outputs; Figs 2, 3, 4, S3 and S10 byte-identical to the versions reviewed",
    "cached ladder fits (05, 07c) against today's rows", paste0(sum(ok05), " of ", length(ok05), " and ", sum(ok07), " of ", length(ok07), " identical in strength and variance"),
    "cached prediction folds (06, 07d) against today's rows", paste0(sum(fo06, na.rm = TRUE), " of ", sum(!is.na(fo06)), " and ", sum(fo07, na.rm = TRUE), " of ", sum(!is.na(fo07)), " identical"),
    "package versions against renv.lock", "every package the code loads is in the lock at the locked version (R 4.6.1); renv itself is not activated (no project library)")),
  "", "## Fixed in this step (no result changes)", "",
  "- **Nondeterministic figures.** Figs 2, 3, 4 and S10 swapped the order of their two legends between runs (ggplot2 orders them arbitrarily); the order is now fixed to the one in the reviewed figures. Label placement in Figs 4 and S3 now ends on an iteration count, not on ggrepel's half-second clock.",
  "- **Table S4.** 30 resident rows were the five floor and ceiling fits repeated under each property's row set (all row sets are the ladder rows since 2 Nov); dropped, so Table S4 has one row per model and pool. The slope structure is labelled the same way for both pools.",
  "- **Reports.** The 07a fit table showed 0 species for every experiment (`nlevels` on a character column). Steps 05 and 06 said the P and TP slope variances sit at the boundary; only P does (TP is of the size of T). Steps 06 and 07d dropped metafor's 'struct disregarded' warning from their reports, against CLAUDE.md; it is now listed with its count. Old figure numbers (Fig 4, S4, S7, S9, S11, S2) replaced by the final ones in headers, comments and report text. Step 03's report said step 10 draws Fig 1b; step 07a does.",
  "- **README.** The forest sentence compared medians over 19 levels including the Null (the forest's Null is the training mean) and disagreed with the 06b report; it now uses the 06b definition over the 18 model levels. A sentence now states that the colonizer Null is a weak reference (see Statistics). The Fig S12 caption, the ranger snippet and the repository tree were corrected.",
  "- **Prep text.** P2 said 'ErrorRisk >= 4' (the pipeline uses 3, written now from the counts table); P4 called experiments 'destination sites'; the `clean_try()` default threshold and the `gbif_niche.R` header count now match the pipeline. Rerunning the pipeline changed only the reports.",
  "- **Dead code.** `R/yardstick.R` (the shared R² yardstick, replaced 2 Nov) removed.",
  "- **`legacy/` retired.** Nothing read from it: the clean clone ran without it. Every comment citing `legacy/<file>:<lines>` now points to commit 7fed5b1, the last commit holding it; the old pipeline at its original paths is the tag `targets-pipeline-2026-10`. The four provenance files `data/README.md` and `data/source_manifest.csv` refer to moved to `prep/provenance/`. The `.gitignore` lines and the README tree entry removed. 67 git-ignored legacy figures remain on disk (see Decisions).",
  "", "## Decisions for Bills (would change a result)", "",
  paste0("1. **Colonizer Step C with measured absences.** Effects with strength 0 and SE 0 enter `rma.mv` with V = 0. For a species with two or more of them the REML likelihood keeps rising as τ² goes to 0, so the reported τ² is not the maximum (table below). ",
         sum(zero_se$n_zero >= 2), " of ", nrow(zero_se), " modelled colonizers have two or more; ", nrow(tiny_mu), " species (", paste(tiny_mu$species, collapse = ", "), ") have pooled means below 1e-6, which Fig 5 counts as 'warmed cover above zero' (", sum(s07b$mu > 0, na.rm = TRUE), " colonizers; ", sum(s07b$mu > 1e-6, na.rm = TRUE), " with a tolerance). Options: a floor on SE, fitting the absences apart, or a tolerance in Fig 5."),
  "", md_table(reml_gap), "",
  paste0("2. **CH_Calanda2 colonizer years.** The colonizer rows of CH_Calanda2::Cal::Nes have years ", cal2$years[cal2$species_pool == "destination_pool"], " and the residents ", cal2$years[cal2$species_pool == "resident"],
         ": step 01 drops a year when one of the two compared groups was not surveyed, which made sense for the destination-control contrast but not for the warmed-only colonizer model (2 Nov). The colonizer strength from year 3 there rests on one year."),
  paste0("3. **TRY names.** ", nn(nrow(fb)), " kept TRY records reach a species through TRY's accepted name because their own name was not used; for ", nn(fb_rejected), " of them the own name had been resolved and rejected (`prep/R/traits.R` `attach_try_names()`). Some are right (Veratrum album subsp. lobelianum), some move records between species. This is the open step 15 finding 3."),
  paste0("4. **Climate years.** Experiment climate averages every calendar year from the first to the last survey, also years without a survey: ", nrow(gaps), " of ", nrow(survey_years), " analysed experiments have such gaps (",
         paste(head(gaps$experiment_id, 4), collapse = ", "), ", ...). PLAN says 'averaged over the experiment's years', which reads either way."),
  "5. **A clone cannot run P6.** `data/derived/prep/gbif_niche.rds` is a pipeline input but git-ignored; a fresh clone stops at `gbif_niche_file`. The same table is tracked as `gbif_niche.csv`: reading the csv, or tracking the rds, would close this.",
  "6. **Warnings in parallel workers are lost.** `parallel::mclapply` workers muffle warnings into their own copy of `warns`, so the forests' (06b, 07e) 'Warnings: None.' is unverified. Fixing it needs the forests refitted once.",
  paste0("7. **R² trimming in medians.** Taking out an R² below −100 moves the colonizer Woody held-out median from ", fmt(median(woody7)), " to ", fmt(median(woody7[woody7 > -100])), ", and 07e compares metafor over 7 gradients with the forest over 8."),
  "8. **Caches keyed by file existence.** Every model and fold cache is reused when its file exists; today all match their inputs (Runs), but a change upstream would be ignored silently. A hash of the inputs stored with each cache would guard this.",
  "9. **Smaller points.** Fig 2's panels b and c show plain mean cover, which differs in sign from the model strength for a few effects; 30 morphotaxon codes (`*_NID_*`) are labelled genus-level instead of unidentified (same treatment, different label); the IN_Kashmir exclusion is a region name, not a stated rule; Step C reports name the first optimizer's failure cause without reading the error.",
  "10. **Environment and files.** renv is not activated: the code runs on the user library, whose versions match the lock. Nine locked packages are used only by the old pipeline (clubSandwich, crew, geometry, ggalluvial, here, mclogit, tarchetypes, testthat, TNRS). `data/README.md` and `data/source_manifest.csv` (read only for me) still call the RData the network source and name tidyTRY revision 3352199 and `config/external_sources.csv`; installed is tidyTRY 0.2.1 (0169b2f), and the provenance files are now in `prep/provenance/`. CLAUDE.md still states rules about `legacy/` and `figures/diagnostics/`. Remove the leftover ignored files with `rm -rf legacy`.",
  "", "## Open items", "",
  md_table(tibble::tribble(
    ~item, ~status,
    "Response scale (PLAN, open since 30 Oct)", paste0("still open: in Step C the top 10% of effects carry ", fmt(100 * sum(w_c[seq_len(ceiling(length(w_c) / 10))]) / sum(w_c), 0), "% of the weight 1/SE²; the ladder's effective weights are even (step 05)"),
    "GBIF download citation (deadline 28 Apr 2027)", "open: DOI in README; goes into methods and data statement",
    "[P1] reference in the README motivation", "open: needs the network reference",
    "Step 15, upstream-correction ledger", "open (planned): finding 3 still live (decision 3 above) and finding 10 (Kashmir rule) still live; the others are marked fixed or decided in the ledger, not re-checked one by one here",
    "Steps 10 and 11", "awaiting review",
    "GEN slope check (additional_checks/05_gen_check.R)", "obsolete: DIAG taken ladder-wide on 1 Nov",
    "Fig S8 surfaces, Fig 5 per-gradient, results-dependent Fig 6, family rerun", "closed: dropped 2 Nov",
    "Colonizer A0_dest untransformed; SE 0 for absences; colonizer W-slope variance at the boundary", "open (step 11); SE 0 is now decision 1")),
  "", "## Statistics: weak points not discussed before (reported, not changed)", "",
  paste0("- **Shared control plots across experiments.** ", sum(shared$n_o > 1), " of ", nrow(shared), " experiments share their origin controls with another experiment of the gradient; ",
         nrow(lad_shared), " species × origin-site groups hold ", sum(lad_shared$n), " of the ", nrow(lad$data), " ladder effects, whose strengths subtract the same control cover. Their sampling errors are correlated, V is diagonal and the ladder has no gradient term, so the fixed-effect intervals (Figs 7, S11, S15) are too narrow. Held-out prediction is not affected (whole gradients are held out)."),
  paste0("- **The origin abundance shares data with the response.** A0_origin comes from the origin controls whose later cover is subtracted in strength: its correlation with the origin-control cover from year 3 is ", fmt(cor(cpl$A0_origin, log(100 * cpl$oc3 + 1))),
         ", and with strength ", fmt(cor(cpl$A0_origin, cpl$strength)), ". Part of the strong negative A0_origin effect is regression to the mean, on top of 'cannot lose more than you have'."),
  paste0("- **Standardisation includes the held-out gradient.** W, T and P are z-scored on all 42 experiments; leaving one gradient out would move their means by at most ", fmt(leak[["W"]]), ", ", fmt(leak[["T"]]), " and ", fmt(leak[["P"]]), " SD. Small, but the held-out scores use information from the held-out gradient."),
  paste0("- **Fig 4 classes are signs of noisy estimates.** ", nn(sum(abs(res$strength) < res$strength_se)), " of ", nn(nrow(res)), " resident effects have |strength| below one SE; the loser and winner slices partly sort noise."),
  paste0("- **Consistent classes rest on τ² = 0.** ", sum(s04$tau < 1e-4 & s04$class != "Context-dependent"), " of ", sum(s04$class != "Context-dependent"), " consistent residents have τ²_s at its boundary (", sum(s04$tau < 1e-4), " of ", nrow(s04), " species overall). With a t quantile (df = experiments − 1) in place of 1.96, ",
         sum(cls_t != s04$class), " of the ", sum(s04$class != "Context-dependent"), " consistent species become context-dependent."),
  "", md_table(win), "",
  "- **Unequal averaging windows** (table above): experiments contribute different numbers of years from year 3, so precision and the size of strength differ with the window.",
  paste0("- **Held-out medians are fragile.** Residents, median (interquartile range; gradients above 0): initial abundances ", spread(p06, "start"), "; Climate + start ", spread(p06, "climate_start"), "; Species-dependent ", spread(p06, "species"),
         ". Colonizers on 8 gradients: destination abundance ", spread(p07, "start"), "."),
  paste0("- **The colonizer Null is a biased reference.** The forest's plain training mean reaches a median R² of ", fmt(col_null), " against the metafor Null on the colonizer gradients: the measured absences (strength and SE 0) pull the weighted intercept towards zero, so every colonizer R² against the Null is measured from a low baseline."),
  paste0("- **Colonizer climate scale.** The colonizer ladder carries W, T, P on the residents' 42-experiment scale; on its ", nrow(col_clim), " rows T has mean ", fmt(mean(col_clim$T)), " and SD ", fmt(sd(col_clim$T)), "; the Woody colonizer model rests on one woody species."),
  "", "## Known and not reopened", "",
  "Gaussian model on bounded cover (23 Oct); compositional covariance within an experiment (24 Oct); diagonal T, P, TP slopes (1 Nov); no uncertainty on τ_s (25 Oct); forest tuning by out-of-bag R² on all rows (2 Nov); multiple testing in Fig 7 (raised 2 Nov); colonizer A0_dest scale and the W-slope boundary (step 11).",
  "", "## Warnings captured", "", if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "", "None: this step changes no analysis choice.")
writeLines(report, "reports/12_final_review.md")
message("Step 12 written")
