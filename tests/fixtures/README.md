# Frozen regression fixtures

`step22_reference_results.rds` is a byte-identical copy of the immutable Step
22 numerical authority. Its SHA-256 is
`8d5d7544a6e6bcb757827263763ffae95115fa10b52b187669eaa5627af6931f`.

`canonical_transition_reference.rds` is a compact, version-3 R serialization
of the validated Step 19/20/22 interfaces needed by branch-level regression
tests. It contains tables, covariance matrices, two compact model signatures,
and the phylogeny/PCoA reference interface. It excludes fitted-model
environments, raw source data, caches, figures, and reports. Its SHA-256 is
`6aa52701c604bf23f7366970fbf80e4e3b05f7253ae1961e1f7251ff08c962e6`.

The transition fixture can be recreated before legacy outputs are archived by
running:

```sh
Rscript scripts/freeze_transition_reference.R --write
```

Replacing either fixture requires explicit review and a corresponding update
to `config/reference_artifacts.csv`; ordinary pipeline execution never writes
these paths.
