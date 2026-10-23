# =============================================================================
# get_chelsa() : robust CHELSA climate data extractor for point coordinates
# -----------------------------------------------------------------------------
# One function to pull CHELSA data for three time domains:
#
#   source = "present" -> CHELSA V2.1 daily        (1979 - ~present, ~1 km)
#   source = "cruts"   -> CHELSAcruts monthly      (1901 - 2016, the "pre-1980"
#                                                    long timeseries, ~1 km)
#   source = "paleo"   -> CHELSA-TraCE21k monthly  (21 kyr BP - present, in
#                                                    centennial time slices)
#
# It works in two modes:
#   * POINT mode (default)  -> pass `coords`; returns a tidy long data.frame of
#                              values extracted at your points.
#   * AREA  mode            -> pass `extent = c(xmin, xmax, ymin, ymax)` in
#                              lon/lat; returns a terra SpatRaster (one layer per
#                              date/slice/variable) cropped to that window, and
#                              optionally writes a GeoTIFF.
#
# It uses the same GDAL /vsicurl/ streaming mechanism as Rchelsa::getChelsa: the
# global file is opened lazily and only the bytes you need are pulled - the
# pixels at your points, or just the blocks inside your extent. You never
# download the whole global grid. Two problems in that package that were causing
# the errors you hit are fixed here:
#
#   1. getChelsa() stacks the whole date range in a single terra::rast() call,
#      so if ONE file is missing on the server the entire call errors out.
#      Here every file is read in its own tryCatch(): missing/unavailable files
#      are skipped with a warning and the rest of the download still completes.
#
#   2. Rchelsa::getBCdate() (used internally for paleo dates) returns the wrong
#      month/slice for negative TraCE21k time slices in any month except
#      January. Here the TraCE21k filename is built directly, so every
#      month/slice combination is correct.
#
# Requires: terra. (Rchelsa is NOT required.)
#
# Author: added 2026-07-14
# =============================================================================

# --- internal: build the download URL for one file -----------------------
.chelsa_build_url <- function(source, var, month = NULL, year = NULL, slice = NULL) {
  if (source == "present") {
    # filename order is day_month_year; day is fixed by the Date, handled by caller
    stop(".chelsa_build_url present handled inline")  # never called for present
  } else if (source == "cruts") {
    sprintf(
      "https://os.zhdk.cloud.switch.ch/chelsav1/chelsa_cruts/%s/CHELSAcruts_%s_%d_%d_V.1.0.tif",
      var, var, as.integer(month), as.integer(year)
    )
  } else if (source == "paleo") {
    sprintf(
      "https://os.zhdk.cloud.switch.ch/chelsa01/chelsa_trace21k/global/centennial/%s/CHELSA_TraCE21k_%s_%02d_%04d_V.1.0.tif",
      var, var, as.integer(month), as.integer(slice)
    )
  }
}

# --- internal: is this variable a temperature variable (needs conversion)? ---
.chelsa_is_temp <- function(var) {
  tolower(var) %in% c("tas", "tasmax", "tasmin", "tmax", "tmin")
}

# --- internal: normalise the user's coordinates to lon/lat/id ------------
.chelsa_prep_coords <- function(coords, id_col = NULL) {
  coords <- as.data.frame(coords)
  nms <- tolower(names(coords))

  lon_i <- which(nms %in% c("lon", "long", "longitude", "x"))[1]
  lat_i <- which(nms %in% c("lat", "latitude", "y"))[1]

  if (is.na(lon_i) || is.na(lat_i)) {
    # fall back to first two columns, assumed lon, lat
    if (ncol(coords) < 2) stop("`coords` must have at least two columns (lon, lat).")
    lon_i <- 1L; lat_i <- 2L
    warning("Could not find lon/lat columns by name; using the first two columns as lon, lat.")
  }

  if (!is.null(id_col)) {
    if (!id_col %in% names(coords)) stop("`id_col` '", id_col, "' not found in coords.")
    id <- coords[[id_col]]
  } else if (any(nms %in% c("id"))) {
    id <- coords[[which(nms == "id")[1]]]
  } else {
    id <- seq_len(nrow(coords))
  }

  data.frame(id = id,
             lon = as.numeric(coords[[lon_i]]),
             lat = as.numeric(coords[[lat_i]]),
             stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------------
# get_chelsa()
# -------------------------------------------------------------------------
# Arguments
#   coords     POINT mode: data.frame / matrix of points. lon & lat are detected
#              by column name (lon/long/longitude/x, lat/latitude/y) or, failing
#              that, taken as the first two columns. An "id" column (or `id_col`)
#              is carried through; otherwise a row index is used. Leave NULL when
#              using `extent`.
#   extent     AREA mode: numeric c(xmin, xmax, ymin, ymax) in decimal degrees
#              (lon/lon/lat/lat). When supplied, each file is cropped to this
#              window and the function returns a terra SpatRaster instead of a
#              data.frame. `coords` and `extent` are mutually exclusive.
#   vars       character vector of CHELSA variables to fetch. Valid values:
#                present : "tas", "tasmax", "tasmin", "pr"
#                cruts   : "tmax", "tmin", "prec"
#                paleo   : "pr", "tasmax", "tasmin", "tz"
#   source     "present", "cruts", or "paleo" (see file header).
#   start,end  Meaning depends on `source`:
#                present : Date or "YYYY-MM-DD"; a daily sequence is built.
#                cruts   : year (numeric) or Date; a yearly x `months` grid.
#                paleo   : TraCE21k time-slice id (integer, -200 .. 20; each
#                          step is ~100 yr, 20 = present, 0 = year 0,
#                          -200 = 20 kyr BP). A slice sequence x `months` grid.
#   months     which months to fetch. Ignored for daily unless you want to keep
#              only certain months of the daily series; used as the month grid
#              for cruts and paleo. Default 1:12.
#   to_celsius convert temperature variables to degrees Celsius (default TRUE).
#              Precipitation is always returned in its native units (mm).
#   output_dir if not NULL, a tidy CSV is written here (created if needed).
#   id_col     name of the id column in `coords` (optional).
#   verbose    print progress (default TRUE).
#
# Value
#   POINT mode: a long-format data.frame with columns
#     id, lon, lat, source, variable, year, month, date, slice, value
#   (`date` is filled for present/cruts, `slice` for paleo; the other is NA).
#   AREA mode:  a terra SpatRaster, one layer per file, layer names encoding
#     variable + date (present/cruts) or variable + slice + month (paleo). A
#     data.frame mapping layer -> variable/year/month/slice is attached as
#     attr(x, "layers").
#   Files missing on the server are omitted and reported either way.
# -------------------------------------------------------------------------
get_chelsa <- function(coords = NULL,
                       vars,
                       source = c("present", "cruts", "paleo"),
                       start,
                       end,
                       months = 1:12,
                       extent = NULL,
                       to_celsius = TRUE,
                       output_dir = NULL,
                       id_col = NULL,
                       verbose = TRUE) {

  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required. install.packages('terra')")
  }
  source <- match.arg(source)

  # decide mode -------------------------------------------------------------
  if (!is.null(coords) && !is.null(extent))
    stop("Pass either `coords` (point mode) or `extent` (area mode), not both.")
  mode <- if (!is.null(extent)) "area" else "points"

  if (mode == "points") {
    if (is.null(coords)) stop("Provide `coords` (point mode) or `extent` (area mode).")
    pts <- .chelsa_prep_coords(coords, id_col = id_col)
    xy  <- pts[, c("lon", "lat")]
  } else {
    if (length(extent) != 4 || !is.numeric(extent))
      stop("`extent` must be numeric c(xmin, xmax, ymin, ymax) in lon/lat degrees.")
    ext_obj <- terra::ext(extent)
  }

  # sanity check on variable names ------------------------------------------
  valid <- switch(source,
                  present = c("tas", "tasmax", "tasmin", "pr"),
                  cruts   = c("tmax", "tmin", "prec"),
                  paleo   = c("pr", "tasmax", "tasmin", "tz"))
  bad <- setdiff(vars, valid)
  if (length(bad)) {
    warning("For source='", source, "' the usual variables are: ",
            paste(valid, collapse = ", "),
            ". Unrecognised variable(s) passed (will still be attempted): ",
            paste(bad, collapse = ", "))
  }

  # build the job list (one row per file to read) ---------------------------
  jobs <- list()

  if (source == "present") {
    d0 <- as.Date(start); d1 <- as.Date(end)
    if (is.na(d0) || is.na(d1)) stop("present: `start`/`end` must be Dates or 'YYYY-MM-DD'.")
    dates <- seq(d0, d1, by = "day")
    dates <- dates[as.integer(format(dates, "%m")) %in% months]
    for (v in vars) for (dt in as.character(dates)) {
      dt <- as.Date(dt)
      url <- sprintf(
        "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/daily/%s/%s/CHELSA_%s_%s_%s_%s_V.2.1.tif",
        v, format(dt, "%Y"), v, format(dt, "%d"), format(dt, "%m"), format(dt, "%Y"))
      jobs[[length(jobs) + 1]] <- list(var = v, url = url,
                                       year = as.integer(format(dt, "%Y")),
                                       month = as.integer(format(dt, "%m")),
                                       date = dt, slice = NA_integer_)
    }

  } else if (source == "cruts") {
    # accept either a bare year (1975) or a Date/"YYYY-MM-DD"
    as_year <- function(x) {
      y <- suppressWarnings(as.integer(x))
      if (is.na(y)) y <- as.integer(format(as.Date(x), "%Y"))
      y
    }
    y0 <- as_year(start); y1 <- as_year(end)
    if (is.na(y0) || is.na(y1)) stop("cruts: `start`/`end` must be years or Dates.")
    if (y0 < 1901 || y1 > 2016)
      warning("CHELSAcruts covers 1901-2016; requested ", y0, "-", y1,
              " will be clamped to that range.")
    years <- seq(max(1901L, y0), min(2016L, y1))
    for (v in vars) for (yr in years) for (m in months) {
      jobs[[length(jobs) + 1]] <- list(var = v,
                                       url = .chelsa_build_url("cruts", v, month = m, year = yr),
                                       year = yr, month = as.integer(m),
                                       date = as.Date(sprintf("%d-%02d-01", yr, m)),
                                       slice = NA_integer_)
    }

  } else { # paleo
    s0 <- as.integer(start); s1 <- as.integer(end)
    if (is.na(s0) || is.na(s1)) stop("paleo: `start`/`end` must be integer TraCE21k slice ids (-200..20).")
    if (s0 < -200 || s1 > 20)
      warning("TraCE21k slices run -200..20; requested ", s0, "..", s1,
              " will be clamped.")
    slices <- seq(max(-200L, s0), min(20L, s1))
    for (v in vars) for (sl in slices) for (m in months) {
      jobs[[length(jobs) + 1]] <- list(var = v,
                                       url = .chelsa_build_url("paleo", v, month = m, slice = sl),
                                       year = NA_integer_, month = as.integer(m),
                                       date = as.Date(NA), slice = as.integer(sl))
    }
  }

  # a short, unique label per job (used to name raster layers in area mode) --
  lab <- function(j) {
    if (source == "paleo") sprintf("%s_sl%d_m%02d", j$var, j$slice, j$month)
    else                   sprintf("%s_%s", j$var, format(j$date, "%Y-%m-%d"))
  }

  if (verbose) cat(sprintf("[get_chelsa] source=%s | mode=%s | %s | %d file(s) to read\n",
                           source, mode,
                           if (mode == "points") paste(nrow(pts), "coordinate(s)")
                           else paste0("extent ", paste(extent, collapse = ",")),
                           length(jobs)))

  # read each file on its own, skipping any that are unavailable -------------
  out    <- vector("list", length(jobs))   # point mode: data.frames
  layers <- vector("list", length(jobs))   # area  mode: cropped SpatRasters
  n_ok <- 0L; n_missing <- 0L; missing_urls <- character(0)

  for (i in seq_along(jobs)) {
    j <- jobs[[i]]
    obj <- tryCatch({
      r <- terra::rast(paste0("/vsicurl/", j$url))
      if (mode == "points") as.numeric(terra::extract(r, xy)[, 2])
      else                  terra::crop(r, ext_obj)          # pulls only this window
    }, error = function(e) NULL, warning = function(w) NULL)

    if (is.null(obj)) {
      n_missing <- n_missing + 1L
      missing_urls <- c(missing_urls, j$url)
      if (verbose) cat(sprintf("  ! missing/unavailable: %s\n", basename(j$url)))
      next
    }
    n_ok <- n_ok + 1L

    # unit conversion (works on both numeric vectors and SpatRasters)
    if (.chelsa_is_temp(j$var) && to_celsius) {
      obj <- if (source == "cruts") obj / 10 else obj - 273.15
    }

    if (mode == "points") {
      out[[i]] <- data.frame(
        id = pts$id, lon = pts$lon, lat = pts$lat,
        source = source, variable = j$var,
        year = j$year, month = j$month, date = j$date, slice = j$slice,
        value = obj, stringsAsFactors = FALSE)
    } else {
      names(obj) <- lab(j)
      layers[[i]] <- obj
    }

    if (verbose && (i %% 25 == 0 || i == length(jobs)))
      cat(sprintf("  ... %d/%d files read (%d ok, %d missing)\n",
                  i, length(jobs), n_ok, n_missing))
  }

  # assemble ----------------------------------------------------------------
  if (mode == "points") {
    result <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
    if (is.null(result)) result <- data.frame()
    n_out <- nrow(result)
  } else {
    layers <- layers[!vapply(layers, is.null, logical(1))]
    result <- if (length(layers)) terra::rast(layers) else NULL
    if (!is.null(result)) {
      meta <- do.call(rbind, lapply(jobs, function(j) data.frame(
        layer = lab(j), variable = j$var, year = j$year,
        month = j$month, slice = j$slice, date = j$date,
        stringsAsFactors = FALSE)))
      attr(result, "layers") <- meta[meta$layer %in% names(result), ]
    }
    n_out <- if (is.null(result)) 0L else terra::nlyr(result)
  }

  if (verbose) {
    cat(sprintf("[get_chelsa] done: %d file(s) read, %d skipped as unavailable, %s.\n",
                n_ok, n_missing,
                if (mode == "points") paste(n_out, "rows returned")
                else paste(n_out, "raster layer(s) returned")))
  }
  if (n_missing > 0) {
    if (!is.null(result)) attr(result, "missing_urls") <- missing_urls
    warning(n_missing, " file(s) were unavailable on the server and were skipped. ",
            "See attr(result, 'missing_urls') for the list.")
  }

  # optional output ---------------------------------------------------------
  if (!is.null(output_dir) && n_out > 0) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (mode == "points") {
      fn <- file.path(output_dir,
                      sprintf("chelsa_%s_%s.csv", source, paste(vars, collapse = "-")))
      write.csv(result, fn, row.names = FALSE)
    } else {
      fn <- file.path(output_dir,
                      sprintf("chelsa_%s_%s.tif", source, paste(vars, collapse = "-")))
      terra::writeRaster(result, fn, overwrite = TRUE)
    }
    if (verbose) cat("[get_chelsa] written: ", fn, "\n", sep = "")
  }

  result
}
