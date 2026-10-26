# P3: Functional groups and phylogeny

Graminoid and woody status and a phylogeny with its PCoA axes for all 904 study species (residents and destination pool), on the P1 accepted species. The tree is cut from the V.PhyloMaker2 mega-tree GBOTB.extended.LCVP (scenario S3); each species' place is set by the reviewed ledger `prep/config/phylogeny_placements.csv`, which reuses the legacy hand-checked labels for the same species (and for renamed species whose epithet is unchanged), and places new species by rule.

## Functional groups

Graminoid (Poaceae, Cyperaceae, Juncaceae): 160 species. Woody (`prep/config/woodiness.csv`): 46 species. Family taken from the backbone genus (no family in P1): Jacobaea alpina, Koeleria macrantha.

| species | family | woody_basis |
| --- | --- | --- |
| Acer pseudoplatanus | Sapindaceae | genus row: tree (added 26 Oct; legacy genus list lacked it) |
| Arctous alpina | Ericaceae | genus row: dwarf shrub (legacy listed it as Arctostaphylos) |
| Artemisia tridentata | Asteraceae | species row: shrub; Artemisia is mixed |
| Berberis jamesiana | Berberidaceae | genus row: shrub (legacy list) |
| Berberis repens | Berberidaceae | genus row: shrub (legacy list) |
| Betula nana | Betulaceae | genus row: tree or shrub (legacy list) |
| Betula pubescens | Betulaceae | genus row: tree or shrub (legacy list) |
| Calluna vulgaris | Ericaceae | genus row: dwarf shrub (legacy list) |
| Daphne striata | Thymelaeaceae | genus row: shrub (legacy list) |
| Dasiphora fruticosa | Rosaceae | genus row: shrub (legacy list) |
| Dryas octopetala | Rosaceae | genus row: dwarf shrub (legacy list) |
| Empetrum nigrum | Ericaceae | genus row: dwarf shrub (legacy list) |
| Fargesia spathacea | Poaceae | genus row: woody bamboo; also graminoid (added 26 Oct) |
| Genista sagittalis | Fabaceae | genus row: dwarf shrub (added 26 Oct) |
| Helianthemum nummularium | Cistaceae | genus row: dwarf subshrub (legacy list) |
| Helianthemum oelandicum | Cistaceae | genus row: dwarf subshrub (legacy list) |
| Hippophae rhamnoides | Elaeagnaceae | genus row: shrub (legacy list) |
| Juniperus communis | Cupressaceae | genus row: shrub (legacy list) |
| Kalmia procumbens | Ericaceae | genus row: dwarf shrub (legacy list) |
| Larix decidua | Pinaceae | genus row: tree (added 26 Oct) |
| Phyllodoce caerulea | Ericaceae | genus row: dwarf shrub (legacy list) |
| Picea abies | Pinaceae | genus row: tree (added 26 Oct) |
| Pinus sylvestris | Pinaceae | genus row: tree (added 26 Oct) |
| Polygaloides chamaebuxus | Polygalaceae | genus row: dwarf shrub (added 26 Oct) |
| Quercus robur | Fagaceae | genus row: tree (added 26 Oct) |
| Rhododendron websterianum | Ericaceae | genus row: shrub (legacy list) |
| Rosa woodsii | Rosaceae | genus row: shrub (legacy list) |
| Rubus idaeus | Rosaceae | species row: shrub; Rubus is mixed |
| Rubus parviflorus | Rosaceae | species row: shrub; Rubus is mixed |
| Salix herbacea | Salicaceae | genus row: tree or shrub (legacy list) |
| Salix lindleyana | Salicaceae | genus row: tree or shrub (legacy list) |
| Salix polaris | Salicaceae | genus row: tree or shrub (legacy list) |
| Salix reticulata | Salicaceae | genus row: tree or shrub (legacy list) |
| Salix retusa | Salicaceae | genus row: tree or shrub (legacy list) |
| Salix serpillifolia | Salicaceae | genus row: tree or shrub (legacy list) |
| Sorbus aucuparia | Rosaceae | genus row: tree or shrub (legacy list) |
| Spiraea betulifolia | Rosaceae | genus row: shrub (legacy list) |
| Symphoricarpos albus | Caprifoliaceae | genus row: shrub (legacy list) |
| Thymus praecox | Lamiaceae | genus row: dwarf subshrub (legacy list) |
| Thymus pulegioides | Lamiaceae | genus row: dwarf subshrub (legacy list) |
| Thymus serpyllum | Lamiaceae | genus row: dwarf subshrub (legacy list) |
| Vaccinium membranaceum | Ericaceae | genus row: dwarf shrub (legacy list) |
| Vaccinium myrtillus | Ericaceae | genus row: dwarf shrub (legacy list) |
| Vaccinium scoparium | Ericaceae | genus row: dwarf shrub (legacy list) |
| Vaccinium uliginosum | Ericaceae | genus row: dwarf shrub (legacy list) |
| Vaccinium vitis-idaea | Ericaceae | genus row: dwarf shrub (legacy list) |

Species whose graminoid or woody status differs from step 00 (7):

| species | graminoid | woody | graminoid_00 | woody_00 |
| --- | --- | --- | --- | --- |
| Acer pseudoplatanus | FALSE | TRUE | FALSE | FALSE |
| Artemisia tridentata | FALSE | TRUE | FALSE | FALSE |
| Fargesia spathacea | TRUE | TRUE | TRUE | FALSE |
| Genista sagittalis | FALSE | TRUE | FALSE | FALSE |
| Larix decidua | FALSE | TRUE | FALSE | FALSE |
| Polygaloides chamaebuxus | FALSE | TRUE | FALSE | FALSE |
| Rubus idaeus | FALSE | TRUE | FALSE | FALSE |

## Placements

| placement | species |
| --- | --- |
| exact tip | 643 |
| exact tip under an older genus | 4 |
| exact tip under an older name (legacy reviewed) | 13 |
| grafted within family | 12 |
| grafted within genus | 217 |
| grafted within older genus | 15 |

Grafted at family level, their genus being absent from the mega-tree (12): Chamorchis alpina, Hemiphragma heterophyllum, Iliamna rivularis, Meum athamanticum, Rudbeckia occidentalis, Rumex acetosa, Rumex acetosella, Rumex arifolius, Rumex nepalensis, Rumex sanguineus, Rumex scutatus, Thesium alpinum.

Ledger rows edited by hand relative to the rule-based draft: 0.

## Tree and axes

Tips: 904, rooted and ultrametric. PCoA on cophenetic distances: 903 positive axes; no correction needed (negative/positive eigenvalue ratio 0). 36 axes exceed the broken-stick expectation; PEV1 to PEV3 carry 38.2% of the positive variation. Orientation anchors: Poa trivialis +PEV1, Erigeron formosissimus -PEV2, Euphrasia minima +PEV3.

| PEV | axis_number | eigenvalue | percent_positive_variation | cumulative_percent | broken_stick_expectation_percent | exceeds_broken_stick | correction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PEV1 | 1 | 5677393.28 | 17.52 | 17.52 | 0.82 | TRUE | none |
| PEV2 | 2 | 3705887.18 | 11.44 | 28.96 | 0.71 | TRUE | none |
| PEV3 | 3 | 2993180.28 | 9.24 | 38.19 | 0.65 | TRUE | none |
| PEV4 | 4 | 2092437.89 | 6.46 | 44.65 | 0.61 | TRUE | none |
| PEV5 | 5 | 1562096.07 | 4.82 | 49.47 | 0.59 | TRUE | none |
| PEV6 | 6 | 1219541.32 | 3.76 | 53.24 | 0.56 | TRUE | none |

Against step 00 (476 resident species in both): correlation of scores PEV1 0.998, PEV2 0.794, PEV3 0.028.

## Warnings captured

None.

## Deviations from the plan

The genus map gained three rows found while checking the draft (Eriocoma -> Stipa instead of Achnatherum, Ptilagrostis -> Stipa, Cherleria -> Minuartia): V.PhyloMaker2 binds only into genera listed in its node table, and Achnatherum, Ptilagrostis and Cherleria have tips but no node there.
