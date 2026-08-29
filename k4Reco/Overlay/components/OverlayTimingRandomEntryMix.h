/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/** Background overlay algorithm with random entry sampling.

    This algorithm overlays background events on top of signal events. Each
    background group is treated as a flat collection of EDM4hep event entries,
    independent of how those entries are distributed over files.

    The MCParticleCollection in signal and background are overlaid into one
    collection. SimTrackerHit collections are cropped and overlaid if they are
    in the time window. SimCalorimeterHit collections are overlaid based on the
    cellID. If a signal hit has the same cellID as a background hit, they are
    combined into a single hit. Only hits that have CaloHitContributions in the
    time range are considered.

**/

#include "podio/Frame.h"
#include "podio/Reader.h"

#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

#include "k4FWCore/Transformer.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "Gaudi/Parsers/Factory.h"
#include "Gaudi/Property.h"

#include <condition_variable>
#include <cstddef>
#include <exception>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace OverlayTimingRandomEntryMixNS {
// Entry-sampling infrastructure unique to this component.
struct EventReader {
  struct Request {
    size_t groupIndex;
    size_t entryIndex;
    std::promise<podio::Frame> promise;
  };

  explicit EventReader(std::vector<std::vector<std::string>> fileNames);
  ~EventReader();

  podio::Frame read(size_t groupIndex, size_t entryIndex);
  size_t numberOfGroups() const { return m_fileNames.size(); }
  size_t numberOfEntries(size_t groupIndex) const { return m_numberOfEntries.at(groupIndex); }

private:
  void run(std::promise<void> ready);

  std::vector<std::vector<std::string>> m_fileNames;
  std::vector<podio::Reader> m_readers;
  std::vector<size_t> m_numberOfEntries;

  std::queue<Request> m_requests;
  std::mutex m_queueMutex;
  std::condition_variable m_queueCondition;
  std::thread m_worker;
  bool m_stop{false};
};
} // namespace OverlayTimingRandomEntryMixNS

using OverlayTimingRandomEntryMixReturn =
    std::tuple<edm4hep::MCParticleCollection, std::vector<edm4hep::SimTrackerHitCollection>,
               std::vector<edm4hep::SimCalorimeterHitCollection>, std::vector<edm4hep::CaloHitContributionCollection>>;

struct OverlayTimingRandomEntryMix
    : public k4FWCore::MultiTransformer<OverlayTimingRandomEntryMixReturn(
          const edm4hep::EventHeaderCollection& headers, const edm4hep::MCParticleCollection&,
          const std::vector<const edm4hep::SimTrackerHitCollection*>&,
          const std::vector<const edm4hep::SimCalorimeterHitCollection*>&)> {
  OverlayTimingRandomEntryMix(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc,
                         {KeyValues("EventHeader", {"EventHeader"}), KeyValues("MCParticles", {"DefaultMCParticles"}),
                          KeyValues("SimTrackerHits", {"DefaultSimTrackerHits"}),
                          KeyValues("SimCalorimeterHits", {"DefaultSimCalorimeterHits"})},
                         {KeyValues("OutputMCParticles", {"NewMCParticles"}),
                          KeyValues("OutputSimTrackerHits", {"NewSimTrackerHits"}),
                          KeyValues("OutputSimCalorimeterHits", {"NewSimCalorimeterHits"}),
                          KeyValues("OutputCaloHitContributions", {"OverlayCaloHitContributions"})}) {}

  StatusCode initialize() final;
  StatusCode finalize() final;

  OverlayTimingRandomEntryMixReturn
  operator()(const edm4hep::EventHeaderCollection& headers, const edm4hep::MCParticleCollection& mcParticles,
             const std::vector<const edm4hep::SimTrackerHitCollection*>& simTrackerHits,
             const std::vector<const edm4hep::SimCalorimeterHitCollection*>& simCalorimeterHits) const final;

  std::pair<float, float> defineTimeWindow(const std::string& collectionName) const;

private:
  constexpr static int SIMTRACKERHIT_INDEX_POSITION = 2;
  constexpr static int SIMCALOHIT_INDEX_POSITION = 3;

  // Overlay configuration shared with OverlayTimingRandomMix.
  Gaudi::Property<bool> m_randomBX{this, "RandomBx", false,
                                   "Place the physics event at a random position in the bunch train"};
  Gaudi::Property<int> m_physBX{this, "PhysicsBX", 1, "Bunch crossing containing the physics event"};
  Gaudi::Property<int> m_numberOfBunches{this, "NBunchtrain", 1, "Number of bunches in the bunch train"};

  // Entry-sampling configuration unique to this component.
  Gaudi::Property<std::vector<std::vector<std::string>>> m_inputFileNames{
      this, "BackgroundFileNames", {}, "EDM4hep background files or directories, grouped for independent sampling"};
  Gaudi::Property<std::vector<double>> m_numberBackground{
      this, "NumberBackground", {}, "Fixed number of background entries, or Poisson mean, for each group"};
  Gaudi::Property<std::vector<bool>> m_poisson{
      this, "Poisson_random_NOverlay", {}, "Draw the number of background entries from a Poisson distribution"};
  Gaudi::Property<std::vector<bool>> m_allowReusingEntries{
      this,
      "AllowReusingBackgroundEntries",
      {},
      "Allow entries in each group to be selected more than once in one bunch crossing"};

  // Collection handling shared with OverlayTimingRandomMix.
  Gaudi::Property<std::string> m_MCParticleCollectionName{this, "BackgroundMCParticleCollectionName", "MCParticle",
                                                          "Name of the MCParticle collection in the background files"};
  Gaudi::Property<float> m_deltaT{this, "Delta_t", 0.5f, "Time difference between bunch crossings"};

  std::unique_ptr<OverlayTimingRandomEntryMixNS::EventReader> m_backgroundEvents;

  Gaudi::Property<std::map<std::string, std::vector<float>>> m_timeWindows{
      this, "TimeWindows", {}, "Time windows for the different collections"};
  Gaudi::Property<bool> m_copyCellIDMetadata{this, "CopyCellIDMetadata", false,
                                             "Copy cellID metadata from signal input to overlay output collections"};
  Gaudi::Property<bool> m_mergeMCParticles{this, "MergeMCParticles", true, "Merge the MCParticle collections"};

  SmartIF<IUniqueIDGenSvc> m_uidSvc;
};
