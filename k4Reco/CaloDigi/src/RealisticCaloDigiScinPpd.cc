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

#include "RealisticCaloDigiScinPpd.h"

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string>

#include "CLHEP/Random/RandBinomial.h"
#include "CLHEP/Random/RandGauss.h"

using namespace std;

DECLARE_COMPONENT(RealisticCaloDigiScinPpd)

RealisticCaloDigiScinPpd::RealisticCaloDigiScinPpd(const std::string& name, ISvcLocator* svcLoc)
    : RealisticCaloDigi(name, svcLoc) {}

float RealisticCaloDigiScinPpd::convertEnergy(float energy,
                                              int inUnit) const { // convert energy from input to output scale (NPE)
  if (inUnit == NPE)
    return energy;
  else if (inUnit == MIP)
    return m_PPD_pe_per_mip * energy;
  else if (inUnit == GEVDEP)
    return m_PPD_pe_per_mip * energy / m_calib_mip;

  throw std::runtime_error("RealisticCaloDigiScinPpd::convertEnergy - unknown unit " + std::to_string(inUnit));
}

float RealisticCaloDigiScinPpd::digitiseDetectorEnergy(float energy) const {
  // input energy in deposited GeV
  // output in npe
  float npe = energy * m_PPD_pe_per_mip / m_calib_mip; // convert to pe scale

  if (m_PPD_n_pixels > 0) {
    // apply average sipm saturation behaviour
    npe = m_PPD_n_pixels * (1.0 - exp(-npe / m_PPD_n_pixels));
    // apply binomial smearing
    float p = npe / m_PPD_n_pixels;             // fraction of hit pixels on SiPM
    npe = m_engine.Binomial(m_PPD_n_pixels, p); // npe now quantised to integer pixels

    if (m_pixSpread > 0) {
      // variations in pixel capacitance
      npe *= m_engine.Gaus(1, m_pixSpread / sqrt(npe));
    }
  }

  return npe;
}
