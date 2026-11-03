# Data inputs

Raw and external inputs of the analysis. They are read only: nothing in the code
writes here, except `data/derived/` (the outputs of every step) and the caches
under `data/cache/` that `prep/gbif_niche.R` fills when it is run by hand.
Checksums of the inputs are in `source_manifest.csv`.

- `transplant_network_clean.sqlite`: the TransPlant network vegetation data
  (TransPlantNetwork/transplant_network_data_cleaning), the community source of
  P1. Its taxonomy columns are not used; names are resolved once in prep.
- `TransPlantNetwork_060824.RData`: the network's earlier export, read for its
  dictionary of local species codes (`cleanedspecies`) and for the P1
  reconciliation with the earlier data; not a source of cover values.
- `try_47163.txt`: the TRY export (large, ignored by Git), read by P2.
- `cache/chelsa_monthly/`: monthly CHELSA V2.1 temperature and precipitation
  per site and year (P4; ignored by Git).
- `cache/chelsa_bioclim/`: CHELSA 1981-2010 bio10 and bio18 rasters, for the
  niche optima (P6) and the niche climatology check (P4).
- `cache/gbif/`, `cache/naturalearth/`: the GBIF occurrence download
  (https://doi.org/10.15468/dl.7p28y9) and the land polygons used to clean its
  coordinates; read by `prep/gbif_niche.R`, which writes
  `data/derived/prep/gbif_niche.rds` (tracked) and `gbif_download.csv`.
- `MacroContextTransplantDat_2.csv`: Lynn's niche table, read by
  `prep/gbif_niche.R` for comparison only.

TRY observations are assigned Köppen-Geiger zones from their coordinates by
`tidyTRY::extract_climate_zones()`, with the raster bundled in tidyTRY 0.2.1
(revision 0169b2f, the version in `renv.lock`); its SHA-256 is checked before
filtering. Provenance records of the external code and of the backups made
before the rebuild are in `prep/provenance/`.
