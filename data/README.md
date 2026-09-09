# Data inputs

This directory contains immutable or externally obtained scientific inputs.
Large source files and caches remain at their established paths and are not
duplicated for presentation purposes.

- `TransPlantNetwork_060824.RData`: current network vegetation source.
- `try_47163.txt`: large TRY export (ignored by Git).
- `cache/chelsa_monthly/`: required resumable CHELSA cache (ignored by Git).
- `source_manifest.csv`: checksums and provenance for the inputs used by the
  canonical targets lineage.

No target searches this directory. Every input path is supplied explicitly by
configuration and represented by a file target.
