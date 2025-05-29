# v00.02.00

* 2025-05-28 jmcarcell ([PR#32](https://github.com/key4hep/k4Reco/pull/32))
  - Use ubuntu-latest for running pre-commit

* 2025-05-28 jmcarcell ([PR#31](https://github.com/key4hep/k4Reco/pull/31))
  - Follow https://github.com/iLCSoft/ConformalTracking/pull/63 and fix an out of bounds erase and add an exception if trying to erase out of bounds

* 2025-04-29 jmcarcell ([PR#23](https://github.com/key4hep/k4Reco/pull/23))
  - Add ported processors that are used in tracking: `ConformalTracking`, `ClonesAndSplitTracksFinder`, `RefitFinal` and `TruthTrackFinder`
  - Add ported tracking utilities and tools that are needed for the algorithms
  - Add tests to run the new Gaudi algorithms and compare their output to the output obtained from the Marlin wrapper
  - Add some README files with comments about the stuff that is ported, for example, which parts have not been ported.

* 2025-03-27 jmcarcell ([PR#29](https://github.com/key4hep/k4Reco/pull/29))
  - Rename SimTrackerHitCollectionName to SimTrackHitCollectionName in DDPlanarDigi, because this is the original input name

* 2025-03-25 Mateusz Jakub Fila ([PR#28](https://github.com/key4hep/k4Reco/pull/28))
  - Change PRNGs from `thread_local static` members to execution-local variables
  - Remove `mutable` and don't mutate "PhysicsBX" property during execution

* 2025-03-19 Mateusz Jakub Fila ([PR#26](https://github.com/key4hep/k4Reco/pull/26))
  - Use the properties `Input` and `Output` with `IOSvc` instead of the deprecated `input` and `output`.

* 2025-03-18 jmcarcell ([PR#25](https://github.com/key4hep/k4Reco/pull/25))
  - DDPlanarDigi: Add the option "CellIDBits" to select only some number of bits in the Cell IDs. It can be used for detectors that have segmentation. Note that for this to work, the segmentation bits have to be at the end of the cell ID number, since only the first N bits will be used.

* 2025-03-17 jmcarcell ([PR#24](https://github.com/key4hep/k4Reco/pull/24))
  - Use auto std::uint64_t with cellIDs instead of converting to int in DDPlanarDigi

* 2025-02-12 jmcarcell ([PR#20](https://github.com/key4hep/k4Reco/pull/20))
  - Improve failing when retrieving GeoSvc in DDPlanarDigi

* 2025-02-12 jmcarcell ([PR#19](https://github.com/key4hep/k4Reco/pull/19))
  - Fix setting the time for tracker hits in Overlay

* 2025-02-10 jmcarcell ([PR#18](https://github.com/key4hep/k4Reco/pull/18))
  - Use the usual naming `m_<name>` instead of keeping the name they had in OverlayTiming. Makes it easier to read not to have the two mixed.

* 2024-12-06 jmcarcell ([PR#16](https://github.com/key4hep/k4Reco/pull/16))
  - Add const where possible and improve the names of temporary variables in the Overlay algorithm

* 2024-12-06 jmcarcell ([PR#15](https://github.com/key4hep/k4Reco/pull/15))
  - Set the time for SimTrackerHits and CaloHitContributions from background in the Overlay algorithm

* 2024-11-20 Thomas Madlener ([PR#13](https://github.com/key4hep/k4Reco/pull/13))
  - Make sure to find DD4hep first during CMake configuration to avoid running into issues with different python versions.

* 2024-10-22 jmcarcell ([PR#12](https://github.com/key4hep/k4Reco/pull/12))
  - Add an option to copy CellID metadata when overlaying, so that the new collections can be equivalent to the old ones.

* 2024-10-15 jmcarcell ([PR#11](https://github.com/key4hep/k4Reco/pull/11))
  - Move the check for the end of file before getting the next event, making it possible to run with a bacground file with the same number of events as the signal file (with the option of reusing disabled). Previously the check was happening one event too soon.

* 2024-10-15 jmcarcell ([PR#10](https://github.com/key4hep/k4Reco/pull/10))
  - Fix the steering file for Overlay Timing: Make sure to pass a string for the MCParticle name in the background file and not a list

# v00-01-00

* 2024-09-30 jmcarcell ([PR#9](https://github.com/key4hep/k4Reco/pull/9))
  - Fix the steering file for DDPlanarDigi; remove an unused property

* 2024-09-24 jmcarcell ([PR#8](https://github.com/key4hep/k4Reco/pull/8))
  - Fix histograms after Gaudi v39. RootHistogram has been renamed to StaticHistogram, see https://github.com/key4hep/k4FWCore/pull/239

* 2024-09-19 jmcarcell ([PR#2](https://github.com/key4hep/k4Reco/pull/2))
  - Add an Overlay algorithm, ported from https://github.com/iLCSoft/Overlay/blob/master/include/OverlayTiming.h

* 2024-09-09 jmcarcell ([PR#7](https://github.com/key4hep/k4Reco/pull/7))
  - Use the Key4hepConfig flag to set the standard, compiler flags and rpath magic.

* 2024-08-09 jmcarcell ([PR#6](https://github.com/key4hep/k4Reco/pull/6))
  - Change Associations to Links

* 2024-07-12 jmcarcell ([PR#5](https://github.com/key4hep/k4Reco/pull/5))
  - Don't link to GaudiAlgLib

* 2024-07-07 jmcarcell ([PR#4](https://github.com/key4hep/k4Reco/pull/4))
  - Change the default value for the time resolution to -1, which will mean no time smearing by default.
  - Use the histogram service to save histograms instead of saving our own files (which meant that if the algorithm runs multiple time there will be multiple files, or a single one if the name is the name for all of them with the content of the last one that was saved).
  - Other minor changes like removing an unnecessary header or changing `and` to `&&`.

* 2024-06-27 tmadlener ([PR#3](https://github.com/key4hep/k4Reco/pull/3))
  - Switch from TrackerHit Plane asociation to TrackerHit association, because the former has been removed from EDM4hep in [EDM4hep#331](https://github.com/key4hep/EDM4hep/pull/331)

* 2024-06-24 jmcarcell ([PR#1](https://github.com/key4hep/k4Reco/pull/1))
  - Add the DDPlanarDigi processor as a Gaudi algorithm. It should have exactly the same logic as the original DDPlanarDigiProcessor. There has been a bit of clean up, some parameters have been changed and some variable names have been changed to be in the style of the rest of the key4hep repositories. Ported from https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Digitisers/include/DDPlanarDigiProcessor.h

# v00-01

* This file is also automatically populated by the tagging script