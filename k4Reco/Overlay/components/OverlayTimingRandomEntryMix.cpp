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

#include "OverlayTimingRandomEntryMix.h"

#include "podio/FrameCategories.h"

#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/Constants.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

#include "k4FWCore/MetadataUtils.h"

#include <GaudiKernel/MsgStream.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {
// Entry-sampling implementation unique to this component.
std::vector<std::string> resolveFiles(const std::vector<std::string>& inputs) {
  std::vector<std::string> files;
  for (const auto& input : inputs) {
    std::error_code error;
    if (!fs::is_directory(input, error)) {
      files.push_back(input);
      continue;
    }

    for (fs::directory_iterator entry(input, error), end; entry != end; entry.increment(error)) {
      if (error) {
        throw std::runtime_error("Could not read background directory '" + input + "': " + error.message());
      }
      if (entry->is_regular_file() && entry->path().extension() == ".root") {
        files.push_back(entry->path().string());
      }
    }
    if (error) {
      throw std::runtime_error("Could not read background directory '" + input + "': " + error.message());
    }
  }

  std::sort(files.begin(), files.end());
  return files;
}

// the same as the time-of-flight calculation in OverlayTimingRandomMix.
template <typename T>
float timeOfFlight(const T& position) {
  return std::sqrt(position[0] * position[0] + position[1] * position[1] + position[2] * position[2]) / TMath::C() *
         1e6;
}

// Entry-sampling implementation unique to this component.
std::vector<size_t> selectEntries(size_t poolSize, size_t count, bool withReplacement, std::mt19937& engine) {
  std::vector<size_t> entries;
  entries.reserve(count);

  if (withReplacement) {
    std::uniform_int_distribution<size_t> distribution(0, poolSize - 1);
    for (size_t i = 0; i < count; ++i) {
      entries.push_back(distribution(engine));
    }
    return entries;
  }

  if (count > poolSize) {
    throw std::out_of_range("Requested more background entries than are available without replacement");
  }

  entries.resize(poolSize);
  std::iota(entries.begin(), entries.end(), 0);
  std::shuffle(entries.begin(), entries.end(), engine);
  entries.resize(count);
  return entries;
}
} // namespace

namespace OverlayTimingRandomEntryMixNS {
// Entry-aware reader; all podio and ROOT I/O remains on this worker thread.
EventReader::EventReader(std::vector<std::vector<std::string>> fileNames, std::vector<bool> oneEntryPerFile)
    : m_fileNames(std::move(fileNames)), m_oneEntryPerFile(std::move(oneEntryPerFile)) {
  std::promise<void> ready;
  auto readyFuture = ready.get_future();
  m_worker = std::thread([this, ready = std::move(ready)]() mutable { run(std::move(ready)); });

  try {
    readyFuture.get();
  } catch (...) {
    m_worker.join();
    throw;
  }
}

EventReader::~EventReader() {
  {
    std::lock_guard<std::mutex> lock(m_queueMutex);
    m_stop = true;
  }
  m_queueCondition.notify_one();
  if (m_worker.joinable()) {
    m_worker.join();
  }
}

void EventReader::run(std::promise<void> ready) {
  try {
    m_readers.reserve(m_fileNames.size());
    m_numberOfEntries.reserve(m_fileNames.size());
    for (size_t group = 0; group < m_fileNames.size(); ++group) {
      const auto& names = m_fileNames[group];
      if (m_oneEntryPerFile[group]) {
        m_readers.emplace_back();
        m_numberOfEntries.push_back(names.size());
      } else {
        auto reader = std::make_unique<podio::Reader>(podio::makeReader(names));
        m_numberOfEntries.push_back(reader->getEntries("events"));
        m_readers.push_back(std::move(reader));
      }
    }
    ready.set_value();
  } catch (...) {
    ready.set_exception(std::current_exception());
    return;
  }

  while (true) {
    Request request;
    {
      std::unique_lock<std::mutex> lock(m_queueMutex);
      m_queueCondition.wait(lock, [this]() { return m_stop || !m_requests.empty(); });
      if (m_stop && m_requests.empty()) {
        return;
      }
      request = std::move(m_requests.front());
      m_requests.pop();
    }

    try {
      if (m_oneEntryPerFile.at(request.groupIndex)) {
        const auto& filename = m_fileNames.at(request.groupIndex).at(request.entryIndex);
        auto reader = podio::makeReader(filename);
        request.promise.set_value(reader.readEvent(0));
      } else {
        request.promise.set_value(m_readers.at(request.groupIndex)->readEvent(request.entryIndex));
      }
    } catch (...) {
      request.promise.set_exception(std::current_exception());
    }
  }
}

podio::Frame EventReader::read(size_t groupIndex, size_t entryIndex) {
  Request request{groupIndex, entryIndex, {}};
  auto result = request.promise.get_future();
  {
    std::lock_guard<std::mutex> lock(m_queueMutex);
    m_requests.push(std::move(request));
  }
  m_queueCondition.notify_one();
  return result.get();
}
} // namespace OverlayTimingRandomEntryMixNS

// the same as OverlayTimingRandomMix::define_time_windows.
std::pair<float, float> OverlayTimingRandomEntryMix::defineTimeWindow(const std::string& collectionName) const {
  const auto& window = m_timeWindows.value().at(collectionName);
  return {window[0], window[1]};
}

// Entry-aware input discovery and validation unique to this component.
StatusCode OverlayTimingRandomEntryMix::initialize() {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_inputFileNames.empty()) {
    error() << "BackgroundFileNames must contain at least one group" << endmsg;
    return StatusCode::FAILURE;
  }

  const size_t numberOfGroups = m_inputFileNames.size();
  if (m_numberBackground.empty()) {
    m_numberBackground = std::vector<double>(numberOfGroups, 1.0);
  }
  if (m_poisson.empty()) {
    m_poisson = std::vector<bool>(numberOfGroups, false);
  }
  if (m_allowReusingEntries.empty()) {
    m_allowReusingEntries = std::vector<bool>(numberOfGroups, false);
  }
  if (m_oneEntryPerFile.empty()) {
    m_oneEntryPerFile = std::vector<bool>(numberOfGroups, false);
  }

  if (m_numberBackground.size() != numberOfGroups || m_poisson.size() != numberOfGroups ||
      m_allowReusingEntries.size() != numberOfGroups || m_oneEntryPerFile.size() != numberOfGroups) {
    error() << "NumberBackground, Poisson_random_NOverlay, AllowReusingBackgroundEntries and "
               "OneEntryPerFile must each have one value per background group"
            << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_numberOfBunches.value() <= 0) {
    error() << "NBunchtrain must be positive" << endmsg;
    return StatusCode::FAILURE;
  }
  if (!m_randomBX.value() && (m_physBX.value() < 1 || m_physBX.value() > m_numberOfBunches.value())) {
    error() << "PhysicsBX must be between 1 and NBunchtrain" << endmsg;
    return StatusCode::FAILURE;
  }

  for (const auto& [collection, window] : m_timeWindows.value()) {
    if (window.size() != 2 || !std::isfinite(window[0]) || !std::isfinite(window[1]) || window[0] >= window[1]) {
      error() << "Time window for " << collection << " must contain two finite, increasing values" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  for (const auto& collection : inputLocations(SIMTRACKERHIT_INDEX_POSITION)) {
    if (m_timeWindows.value().find(collection) == m_timeWindows.value().end()) {
      error() << "No time window defined for collection " << collection << endmsg;
      return StatusCode::FAILURE;
    }
  }
  for (const auto& collection : inputLocations(SIMCALOHIT_INDEX_POSITION)) {
    if (m_timeWindows.value().find(collection) == m_timeWindows.value().end()) {
      error() << "No time window defined for collection " << collection << endmsg;
      return StatusCode::FAILURE;
    }
  }

  std::vector<std::vector<std::string>> inputFiles;
  inputFiles.reserve(numberOfGroups);
  try {
    for (size_t group = 0; group < numberOfGroups; ++group) {
      auto files = resolveFiles(m_inputFileNames.value()[group]);
      if (files.empty()) {
        error() << "Background group " << group << " contains no ROOT files" << endmsg;
        return StatusCode::FAILURE;
      }
      inputFiles.push_back(std::move(files));
    }
    m_backgroundEvents = std::make_unique<OverlayTimingRandomEntryMixNS::EventReader>(
        std::move(inputFiles), m_oneEntryPerFile.value());
  } catch (const std::exception& exception) {
    error() << "Could not open background input: " << exception.what() << endmsg;
    return StatusCode::FAILURE;
  }

  for (size_t group = 0; group < numberOfGroups; ++group) {
    const double requested = m_numberBackground.value()[group];
    if (!std::isfinite(requested) || requested < 0.0) {
      error() << "NumberBackground for group " << group << " must be finite and non-negative" << endmsg;
      return StatusCode::FAILURE;
    }
    if (requested > static_cast<double>(std::numeric_limits<size_t>::max())) {
      error() << "NumberBackground for group " << group << " is too large" << endmsg;
      return StatusCode::FAILURE;
    }
    if (!m_poisson.value()[group] && std::floor(requested) != requested) {
      error() << "Fixed NumberBackground for group " << group << " must be an integer" << endmsg;
      return StatusCode::FAILURE;
    }

    const size_t available = m_backgroundEvents->numberOfEntries(group);
    if (available == 0) {
      error() << "Background group " << group << " contains no events entries" << endmsg;
      return StatusCode::FAILURE;
    }
    if (!m_poisson.value()[group] && !m_allowReusingEntries.value()[group] && requested > available) {
      error() << "Background group " << group << " requests " << requested << " entries, but only " << available
              << " are available without replacement" << endmsg;
      return StatusCode::FAILURE;
    }
    info() << "Background group " << group << " contains " << available << " entries" << endmsg;
  }

  return StatusCode::SUCCESS;
}

OverlayTimingRandomEntryMixReturn OverlayTimingRandomEntryMix::operator()(
    const edm4hep::EventHeaderCollection& headers, const edm4hep::MCParticleCollection& particles,
    const std::vector<const edm4hep::SimTrackerHitCollection*>& simTrackerHits,
    const std::vector<const edm4hep::SimCalorimeterHitCollection*>& simCaloHits) const {
  const auto seed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), name());
  auto engine = std::mt19937(seed);

  // Signal collection handling is the same as OverlayTimingRandomMix.
  auto outputParticles = edm4hep::MCParticleCollection();
  auto outputTrackerHits = std::vector<edm4hep::SimTrackerHitCollection>();
  auto outputCaloHits = std::vector<edm4hep::SimCalorimeterHitCollection>();
  auto outputCaloContributions = std::vector<edm4hep::CaloHitContributionCollection>(simCaloHits.size());

  for (const auto particle : particles) {
    outputParticles.push_back(particle.clone(false));
  }
  for (size_t i = 0; i < particles.size(); ++i) {
    for (const auto& parent : particles[i].getParents()) {
      outputParticles[i].addToParents(outputParticles[parent.getObjectID().index]);
    }
    for (const auto& daughter : particles[i].getDaughters()) {
      outputParticles[i].addToDaughters(outputParticles[daughter.getObjectID().index]);
    }
  }

  for (size_t i = 0; i < simTrackerHits.size(); ++i) {
    const auto& collection = simTrackerHits[i];
    const auto collectionName = inputLocations(SIMTRACKERHIT_INDEX_POSITION)[i];
    const auto [windowStart, windowStop] = defineTimeWindow(collectionName);
    auto outputCollection = edm4hep::SimTrackerHitCollection();
    for (const auto hit : *collection) {
      const float flightTime = timeOfFlight(hit.getPosition());
      if (hit.getTime() > windowStart + flightTime && hit.getTime() < windowStop + flightTime) {
        auto outputHit = hit.clone(false);
        if (hit.getParticle().getObjectID().index != -1) {
          outputHit.setParticle(outputParticles[hit.getParticle().getObjectID().index]);
        }
        outputCollection.push_back(outputHit);
      }
    }
    outputTrackerHits.emplace_back(std::move(outputCollection));
  }

  std::map<int, std::map<uint64_t, edm4hep::MutableSimCalorimeterHit>> cellIDs;
  for (size_t i = 0; i < simCaloHits.size(); ++i) {
    const auto& collection = simCaloHits[i];
    const auto collectionName = inputLocations(SIMCALOHIT_INDEX_POSITION)[i];
    const auto [windowStart, windowStop] = defineTimeWindow(collectionName);
    auto& calorimeterHits = cellIDs[i];
    auto& contributions = outputCaloContributions[i];
    for (const auto hit : *collection) {
      const float flightTime = timeOfFlight(hit.getPosition());
      bool withinTimeWindow = false;
      std::vector<int> contributionIndices;
      for (const auto& contribution : hit.getContributions()) {
        if (contribution.getTime() <= windowStart + flightTime || contribution.getTime() >= windowStop + flightTime) {
          continue;
        }
        withinTimeWindow = true;
        auto outputContribution = contribution.clone(false);
        outputContribution.setParticle(outputParticles[contribution.getParticle().getObjectID().index]);
        contributionIndices.push_back(contributions.size());
        contributions.push_back(std::move(outputContribution));
      }
      if (withinTimeWindow) {
        auto outputHit = hit.clone(false);
        for (const auto contributionIndex : contributionIndices) {
          outputHit.addToContributions(contributions[contributionIndex]);
        }
        calorimeterHits.emplace(hit.getCellID(), std::move(outputHit));
      }
    }
  }

  // Entry selection and explicit entry reads are unique to this component.
  int physicsBX = m_physBX.value();
  if (m_randomBX.value()) {
    physicsBX = std::uniform_int_distribution<int>(1, m_numberOfBunches.value())(engine);
    debug() << "Physics event was placed in bunch crossing " << physicsBX << endmsg;
  }

  for (size_t group = 0; group < m_backgroundEvents->numberOfGroups(); ++group) {
    std::vector<int> bunchCrossings;
    for (int bx = -(physicsBX - 1); bx < m_numberOfBunches.value() - (physicsBX - 1); ++bx) {
      bunchCrossings.push_back(bx);
    }
    std::shuffle(bunchCrossings.begin(), bunchCrossings.end(), engine);

    for (int bxInTrain = 0; bxInTrain < m_numberOfBunches.value(); ++bxInTrain) {
      const int bunchCrossing = bunchCrossings.at(bxInTrain);
      const size_t numberToOverlay = m_poisson.value()[group]
                                         ? std::poisson_distribution<size_t>(m_numberBackground.value()[group])(engine)
                                         : static_cast<size_t>(m_numberBackground.value()[group]);

      const size_t availableEntries = m_backgroundEvents->numberOfEntries(group);
      if (!m_allowReusingEntries.value()[group] && numberToOverlay > availableEntries) {
        throw GaudiException("Background group " + std::to_string(group) + " requested " +
                                 std::to_string(numberToOverlay) + " entries, but only " +
                                 std::to_string(availableEntries) + " are available without replacement",
                             name(), StatusCode::FAILURE);
      }
      const auto selectedEntries =
          selectEntries(availableEntries, numberToOverlay, m_allowReusingEntries.value()[group], engine);

      debug() << "Overlaying " << numberToOverlay << " entries from group " << group << " in bunch crossing "
              << bunchCrossing + physicsBX << endmsg;

      for (const size_t entry : selectedEntries) {
        debug() << "Reading background group " << group << ", entry " << entry << endmsg;
        podio::Frame backgroundEvent = m_backgroundEvents->read(group, entry);
        const auto availableCollections = backgroundEvent.getAvailableCollections();
        const auto timeOffset = bunchCrossing * m_deltaT.value();

        // Background collection handling is the same as OverlayTimingRandomMix.
        if (std::find(availableCollections.begin(), availableCollections.end(), m_MCParticleCollectionName) ==
            availableCollections.end()) {
          warning() << "Collection " << m_MCParticleCollectionName << " not found in background event" << endmsg;
        }

        std::map<int, int> oldToNewParticle;
        std::map<int, std::pair<std::vector<int>, std::vector<int>>> particleRelations;

        if (m_mergeMCParticles) {
          const auto& backgroundParticles =
              backgroundEvent.get<edm4hep::MCParticleCollection>(m_MCParticleCollectionName);
          int outputIndex = outputParticles.size();
          for (size_t i = 0; i < backgroundParticles.size(); ++i) {
            auto outputParticle = backgroundParticles[i].clone(false);
            outputParticle.setTime(backgroundParticles[i].getTime() + timeOffset);
            outputParticle.setOverlay(true);
            outputParticles.push_back(outputParticle);
            for (const auto& parent : backgroundParticles[i].getParents()) {
              particleRelations[outputIndex].first.push_back(parent.getObjectID().index);
            }
            for (const auto& daughter : backgroundParticles[i].getDaughters()) {
              particleRelations[outputIndex].second.push_back(daughter.getObjectID().index);
            }
            oldToNewParticle[i] = outputIndex++;
          }

          for (const auto& [index, relations] : particleRelations) {
            const auto& [parents, daughters] = relations;
            for (const auto parent : parents) {
              if (particleRelations.find(oldToNewParticle[parent]) != particleRelations.end()) {
                outputParticles[index].addToParents(outputParticles[oldToNewParticle[parent]]);
              }
            }
            for (const auto daughter : daughters) {
              if (particleRelations.find(oldToNewParticle[daughter]) != particleRelations.end()) {
                outputParticles[index].addToDaughters(outputParticles[oldToNewParticle[daughter]]);
              }
            }
          }
        }

        for (size_t i = 0; i < simTrackerHits.size(); ++i) {
          const auto collectionName = inputLocations(SIMTRACKERHIT_INDEX_POSITION)[i];
          if (std::find(availableCollections.begin(), availableCollections.end(), collectionName) ==
              availableCollections.end()) {
            warning() << "Collection " << collectionName << " not found in background event" << endmsg;
            continue;
          }

          const auto [windowStart, windowStop] = defineTimeWindow(collectionName);
          if (windowStop <= (bunchCrossing - physicsBX) * m_deltaT.value()) {
            debug() << "Skipping collection " << collectionName << " outside its integration window" << endmsg;
            continue;
          }

          auto& outputCollection = outputTrackerHits[i];
          for (const auto hit : backgroundEvent.get<edm4hep::SimTrackerHitCollection>(collectionName)) {
            const float flightTime = timeOfFlight(hit.getPosition());
            if (hit.getTime() + timeOffset <= windowStart + flightTime ||
                hit.getTime() + timeOffset >= windowStop + flightTime) {
              continue;
            }

            auto outputHit = hit.clone(false);
            outputHit.setOverlay(true);
            outputHit.setTime(hit.getTime() + timeOffset);
            if (m_mergeMCParticles) {
              outputHit.setParticle(outputParticles[oldToNewParticle[hit.getParticle().getObjectID().index]]);
            } else {
              const auto particle = hit.getParticle();
              if (particle.isAvailable()) {
                const auto momentum = particle.getMomentum();
                outputHit.setMomentum(
                    {static_cast<float>(momentum.x), static_cast<float>(momentum.y), static_cast<float>(momentum.z)});
              }
            }
            outputCollection.push_back(outputHit);
          }
        }

        for (size_t i = 0; i < simCaloHits.size(); ++i) {
          const auto collectionName = inputLocations(SIMCALOHIT_INDEX_POSITION)[i];
          if (std::find(availableCollections.begin(), availableCollections.end(), collectionName) ==
              availableCollections.end()) {
            warning() << "Collection " << collectionName << " not found in background event" << endmsg;
            continue;
          }

          const auto [windowStart, windowStop] = defineTimeWindow(collectionName);
          if (windowStop <= (bunchCrossing - physicsBX) * m_deltaT.value()) {
            debug() << "Skipping collection " << collectionName << " outside its integration window" << endmsg;
            continue;
          }

          auto& calorimeterHits = cellIDs[i];
          auto& contributions = outputCaloContributions[i];
          for (const auto hit : backgroundEvent.get<edm4hep::SimCalorimeterHitCollection>(collectionName)) {
            if (calorimeterHits.find(hit.getCellID()) == calorimeterHits.end()) {
              auto outputHit = edm4hep::MutableSimCalorimeterHit();
              bool addHit = false;
              for (const auto& contribution : hit.getContributions()) {
                if (contribution.getTime() + timeOffset > windowStart &&
                    contribution.getTime() + timeOffset < windowStop) {
                  addHit = true;
                  auto outputContribution = contribution.clone(false);
                  if (m_mergeMCParticles) {
                    outputContribution.setParticle(
                        outputParticles[oldToNewParticle[contribution.getParticle().getObjectID().index]]);
                  } else {
                    outputContribution.setParticle(edm4hep::MCParticle());
                  }
                  outputContribution.setTime(contribution.getTime() + timeOffset);
                  outputHit.addToContributions(outputContribution);
                  contributions.push_back(outputContribution);
                }
              }
              if (addHit) {
                outputHit.setCellID(hit.getCellID());
                outputHit.setEnergy(hit.getEnergy());
                outputHit.setPosition(hit.getPosition());
                calorimeterHits[outputHit.getCellID()] = outputHit;
              }
            } else {
              auto& outputHit = calorimeterHits[hit.getCellID()];
              for (const auto& contribution : hit.getContributions()) {
                if (contribution.getTime() + timeOffset > windowStart &&
                    contribution.getTime() + timeOffset < windowStop) {
                  auto outputContribution = contribution.clone(false);
                  if (m_mergeMCParticles) {
                    outputContribution.setParticle(
                        outputParticles[oldToNewParticle[contribution.getParticle().getObjectID().index]]);
                  } else {
                    outputContribution.setParticle(edm4hep::MCParticle());
                  }
                  outputContribution.setTime(contribution.getTime() + timeOffset);
                  outputHit.addToContributions(outputContribution);
                  contributions.push_back(outputContribution);
                }
              }
            }
          }
        }
      }
    }
  }

  // Output assembly is the same as OverlayTimingRandomMix.
  for (const auto& [index, calorimeterHits] : cellIDs) {
    auto outputCollection = edm4hep::SimCalorimeterHitCollection();
    for (const auto& [cellID, hit] : calorimeterHits) {
      outputCollection.push_back(std::move(hit));
    }
    outputCaloHits.emplace_back(std::move(outputCollection));
  }

  debug() << "\n\t\tCollection\t\t|\t\tPre BIB\t\t|\t\tPost "
             "BIB\t\t\n--------------------------------------------------------------------------\n";
  for (size_t i = 0; i < simTrackerHits.size(); ++i) {
    debug() << "\tTrackerHits " << i << "\t|\t\t" << simTrackerHits[i]->size() << "\t\t|\t\t"
            << outputTrackerHits[i].size() << "\n";
  }
  for (size_t i = 0; i < simCaloHits.size(); ++i) {
    debug() << "\tCaloHits " << i << "\t|\t\t" << simCaloHits[i]->size() << "\t\t|\t\t" << outputCaloHits[i].size()
            << "\n";
  }
  debug() << endmsg;

  return std::make_tuple(std::move(outputParticles), std::move(outputTrackerHits), std::move(outputCaloHits),
                         std::move(outputCaloContributions));
}

// the same as OverlayTimingRandomMix::finalize.
StatusCode OverlayTimingRandomEntryMix::finalize() {
  if (m_copyCellIDMetadata) {
    for (const auto& [input, output] :
         {std::make_pair(inputLocations("SimTrackerHits"), outputLocations("OutputSimTrackerHits")),
          std::make_pair(inputLocations("SimCalorimeterHits"), outputLocations("OutputSimCalorimeterHits"))}) {
      for (size_t i = 0; i < input.size(); ++i) {
        const auto value = k4FWCore::getParameter<std::string>(
            podio::collMetadataParamName(input[i], edm4hep::labels::CellIDEncoding), this);
        if (value.has_value()) {
          k4FWCore::putParameter(podio::collMetadataParamName(output[i], edm4hep::labels::CellIDEncoding),
                                 value.value(), this);
        } else {
          warning() << "No metadata found for " << input[i] << " when copying CellID metadata was requested" << endmsg;
        }
      }
    }
  }

  return Gaudi::Algorithm::finalize();
}

DECLARE_COMPONENT(OverlayTimingRandomEntryMix)
