# Display names for the ladder's models and fixed-effect terms, shared by the coefficient figures
# (steps 05 and 07c) and Tables S5 and S6 (step 10), in the words Fig 7 uses. The fitted keys ("Start", "A0_origin",
# "W:height_z") stay as they are in data/derived/; only what a reader sees is renamed (Bills, 4 Nov).

model_label <- function(model, pool) {
  start <- ifelse(pool == "colonizers", "Destination abundance", "Initial abundances") # as Figs 6 and S17
  ifelse(model == "Start", start, ifelse(model == "Climate + start", paste("Climate +", tolower(start)), model))
}

term_label <- function(term, pool) {
  word <- c(intrcpt = "intercept", W = "warming", T = "temperature", P = "precipitation", TP = "temperature × precipitation",
            A0_origin = "initial abundance at origin", ldmc_z = "LDMC", lma_z = "LMA", height_z = "height",
            graminoid = "graminoid", woody = "woody", PEV1 = "PEV1", PEV2 = "PEV2", PEV3 = "PEV3",
            T_opt = "thermal optimum", P_opt = "precipitation optimum")
  climate <- c("W", "T", "P", "TP")
  lab <- mapply(function(t, p) {
    parts <- strsplit(t, ":", fixed = TRUE)[[1]]
    w <- unname(ifelse(parts == "A0_dest", if (p == "colonizers") "mean abundance at destination" else "initial abundance at destination",
                       word[parts]))
    clim <- w[parts %in% climate][order(match(parts[parts %in% climate], climate))]
    other <- w[!parts %in% climate]
    # a property leads its interactions ("Height × warming"); an abundance follows the climate
    # ("Warming × initial abundance at origin"), as Fig 7 names them
    paste(if (any(grepl("^A0", parts))) c(clim, other) else c(other, clim), collapse = " × ")
  }, term, pool, USE.NAMES = FALSE)
  stopifnot(!anyNA(lab), !any(grepl("NA", lab, fixed = TRUE))) # every code has a word
  paste0(toupper(substring(lab, 1, 1)), substring(lab, 2))
}
