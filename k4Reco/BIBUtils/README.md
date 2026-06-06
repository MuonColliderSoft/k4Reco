# BIBUtils

Gaudi-native ports of the beam-induced-background (BIB) cleaning processors used
in the Muon Collider reconstruction. They were originally Marlin processors:

| Gaudi algorithm   | Original Marlin processor | Original package |
|-------------------|---------------------------|------------------|
| `FilterConeHits`  | `FilterConeHits`          | [MarlinTrkProcessors](https://github.com/MuonColliderSoft/MarlinTrkProcessors) |
| `CaloConer`       | `CaloConer`               | [MyBIBUtils](https://github.com/madbaron/MyBIBUtils) |
| `CaloHitSelector` | `CaloHitSelector`         | [MyBIBUtils](https://github.com/madbaron/MyBIBUtils) |

All three are functional `k4FWCore::MultiTransformer`s. Selected hits are written
to **subset** collections that reference the original hits, accompanied by a
freshly built reco-to-sim link collection. The simulated hit of a given
reconstructed hit is resolved through the input link collection rather than by
relying on positional alignment between the hit and relation collections.

## FilterConeHits

Keeps the tracker hits that lie inside a cone around the trajectory of a
generator-level MC particle. For each selected MC particle a helix is built from
its production vertex, momentum and charge in the detector field (taken from the
`GeoSvc`). A hit is kept when its angular distance to the helix is below
`DeltaRCut` and/or its 3D distance to the helix is below `Dist3DCut`.

The helix math lives in the self-contained, header-only `TrackHelix` (see
`include/TrackHelix.h`); it reproduces the point-to-helix distance and
cylinder-crossing operations of MarlinUtil's `HelixClass` so the package carries
no dependency on MarlinUtil or any other Marlin-era package.

Each instance handles a single tracker subdetector — configure one instance per
collection (vertex/inner/outer × barrel/endcap), as in the original steering.

| Property | Default | Description |
|---|---|---|
| `MCParticleCollection` | `MCParticle` | input MC particles |
| `TrackerHitInputCollections` | `VBTrackerHits` | input reco tracker hits |
| `TrackerHitInputRelations` | `VBTrackerHitsRelations` | input reco→sim links |
| `TrackerHitOutputCollections` | `VBTrackerHitsConed` | output reco hits (subset) |
| `TrackerSimHitOutputCollections` | `VertexBarrelCollectionConed` | output sim hits (subset) |
| `TrackerHitOutputRelations` | `VBTrackerHitsRelationsConed` | output reco→sim links |
| `DeltaRCut` | `-1` | max angular distance to the helix [rad] (disabled if ≤ 0) |
| `Dist3DCut` | `-1` | max 3D distance to the helix [mm] (disabled if ≤ 0) |
| `ConeAroundStatus` | `[1]` | MC generator statuses to cone around |
| `FillHistograms` | `false` | fill diagnostic histograms |
| `TrackerOuterRadius` | `1500` | tracker barrel outer radius used to clip the helix [mm] |

## CaloConer

Keeps the calorimeter hits within a fixed angular cone (`ConeWidth`, in radians)
around the direction of any generator-level (`generatorStatus == 1`) MC particle.

| Property | Default | Description |
|---|---|---|
| `MCParticleCollectionName` | `MCParticle` | input MC particles |
| `CaloHitCollectionName` | `EcalBarrelCollectionRec` | input reco calo hits |
| `CaloRelationCollectionName` | `EcalBarrelRelationsSimRec` | input reco→sim links |
| `GoodHitCollection` | `EcalBarrelCollectionConed` | output reco hits (subset) |
| `GoodRelationCollection` | `EcalBarrelRelationsSimConed` | output reco→sim links |
| `ConeWidth` | `0.2` | half-opening angle of the cone [rad] |

## CaloHitSelector

Applies a per-cell energy threshold and a time-window selection to calorimeter
hits. The threshold is read as a function of polar angle and layer from two ROOT
histograms (`th_2dmode_sym`, `stddev_sym`) in `ThresholdsFilePath`:
`threshold = mode + Nsigma * stddev`. A constant `FlatThreshold` (GeV) can be
used instead. Surviving hits must fall inside `[TimeWindowMin, TimeWindowMax]`
after a time-of-flight correction. The `GeoSvc` is used to decode the layer from
the cellID.

| Property | Default | Description |
|---|---|---|
| `CaloHitCollectionName` | `EcalBarrelCollectionRec` | input reco calo hits |
| `CaloRelationCollectionName` | `EcalBarrelRelationsSimRec` | input reco→sim links |
| `GoodHitCollection` | `EcalBarrelCollectionSel` | output reco hits (subset) |
| `GoodRelationCollection` | `EcalBarrelRelationsSimSel` | output reco→sim links |
| `ThresholdsFilePath` | `""` | ROOT file with the threshold maps |
| `Nsigma` | `3` | number of BIB-energy sigmas above the modal threshold |
| `FlatThreshold` | `0` | constant threshold [GeV]; overrides the maps if > 0 |
| `TimeWindowMin` / `TimeWindowMax` | `-0.5` / `10` | TOF-corrected time window [ns] |
| `DoBIBsubtraction` | `false` | subtract the modal BIB energy from each cell |

See `options/` for runnable example steering files.
