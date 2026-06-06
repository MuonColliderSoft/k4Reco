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

// Calorimeter digitiser for the IDC ECAL and HCAL
// For other detectors/models SimpleCaloDigi should be used
#include "RealisticCaloDigiSilicon.h"

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string>

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

using namespace std;

DECLARE_COMPONENT(RealisticCaloDigiSilicon)

RealisticCaloDigiSilicon::RealisticCaloDigiSilicon(const std::string& name, ISvcLocator* svcLoc)
    : RealisticCaloDigi(name, svcLoc) {}

float RealisticCaloDigiSilicon::convertEnergy(float energy,
                                              int inUnit) const { // convert energy from input to output scale (MIP)
  // converts input energy to MIP scale
  if (inUnit == MIP)
    return energy;
  else if (inUnit == GEVDEP)
    return energy / m_calib_mip;

  throw std::runtime_error("RealisticCaloDigiSilicon::convertEnergy - unknown unit " + std::to_string(inUnit));
}

float RealisticCaloDigiSilicon::digitiseDetectorEnergy(float energy) const {
  // applies extra digitisation to silicon hits
  //  input energy in deposited GeV
  //  output is MIP scale
  float smeared_energy(energy);
  if (m_ehEnergy > 0) {
    // calculate #e-h pairs
    float nehpairs = 1e9 * energy / m_ehEnergy; // check units of energy! _ehEnergy is in eV, energy in GeV
    // fluctuate it by Poisson (actually an overestimate: Fano factor actually makes it smaller, however even this
    // overstimated effect is tiny for our purposes)
    smeared_energy *= m_engine.Poisson(nehpairs) / nehpairs;
  }

  return smeared_energy / m_calib_mip; // convert to MIP units
}
