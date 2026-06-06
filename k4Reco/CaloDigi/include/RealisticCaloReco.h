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
#ifndef DIGITIZER_REALISTICCALORECO_H
#define DIGITIZER_REALISTICCALORECO_H 1

#include <k4FWCore/Transformer.h>

#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/CalorimeterHit.h>
#include <edm4hep/CalorimeterHitCollection.h>

#include "CalorimeterHitType.h"
#include "k4Interface/IGeoSvc.h"

#include <string>
#include <vector>

/** === RealisticCaloReco Processor === <br>
    realistic reconstruction of calorimeter hits
    e.g. apply sampling fraction correction
    virtual class, technology indenpendent
    D.Jeans 02/2016.

    24 March 2016: removed gap corrections - to be put into separate processor
    changed relations: now keep relation between reconstructed and simulated hits.
*/

struct RealisticCaloReco : k4FWCore::MultiTransformer<
                               std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>(
                                   const edm4hep::CaloHitSimCaloHitLinkCollection&)> {

public:
  RealisticCaloReco(const std::string& name, ISvcLocator* svcLoc);
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  StatusCode initialize();

  /** Called for every run.
   */
  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
  operator()(const edm4hep::CaloHitSimCaloHitLinkCollection& inputLinks) const;

  /** Called after data processing for clean up.
   */
  StatusCode finalize();

protected:
  float getLayerCalib(int ilayer) const;
  virtual float reconstructEnergy(const edm4hep::CalorimeterHit* hit,
                                  int layer) const = 0; // to be overloaded, technology-specific

  // parameters
  // Grouping of calo layers
  Gaudi::Property<std::vector<int>> m_calLayers{this, "calibration_layergroups", {}, "Grouping of calo layers"};
  // Calibration coefficients for layers groups
  Gaudi::Property<std::vector<float>> m_calibrCoeff{
      this, "calibration_factorsMipGev", {}, "Calibration coefficients (MIP->shower GeV) of layers groups"};
  // Cell ID layer string
  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalCalorimeterReadoutID",
      "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};

  SmartIF<IGeoSvc> m_geoSvc;
};

#endif
