# Data inputs

This directory contains immutable or externally obtained scientific inputs.
Large source files and caches remain at their established paths and are not
duplicated for presentation purposes.

- `TransPlantNetwork_060824.RData`: current network vegetation source.
- `try_47163.txt`: large TRY export (ignored by Git).
- `cache/chelsa_monthly/`: required resumable CHELSA cache (ignored by Git).
- `source_manifest.csv`: checksums and provenance for the inputs used by the
  current workflow.

TRY observations are assigned Köppen–Geiger codes from their coordinates by
`tidyTRY::extract_climate_zones()`, using the raster bundled with pinned tidyTRY
revision `3352199fc99f72170a8983c522d63aa9029d0e9f`. The registered
raster SHA-256 is checked before filtering. Codes are not taken from the TRY
export or the old local raster. That unused local raster, other unused data
files, and the obsolete daily CHELSA cache are recoverable from the verified
ZIP listed in `docs/pre_rebuild_backup_receipts.csv`. The 336 monthly CHELSA
files remain at their current paths and are also copied into that ZIP.

The raw RData and TRY export paths are explicit file targets. CHELSA branches
use only their prescribed year-month-variable checkpoint paths. The bundled
tidyTRY raster is a package input verified at the TRY climate-filtering step.
