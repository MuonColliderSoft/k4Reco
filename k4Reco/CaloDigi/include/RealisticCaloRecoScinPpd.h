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
#ifndef REALISTICCALORECOSCINPPD_H
#define REALISTICCALORECOSCINPPD_H 1

#include "RealisticCaloReco.h"

/** === RealisticCaloRecoSilicon Processor === <br>
    realistic reconstruction of scint+PPD calorimeter hits
    D.Jeans 02/2016.
*/

struct RealisticCaloRecoScinPpd final : RealisticCaloReco {
public:
  RealisticCaloRecoScinPpd(const std::string& name, ISvcLocator* svcLoc);

protected:
  float reconstructEnergy(const edm4hep::CalorimeterHit* hit, int layer) const override;

  Gaudi::Property<float> m_PPD_pe_per_mip{this, "ppd_mipPe", 10.0f,
                                          "# Photo-electrons per MIP (scintillator): used to Poisson smear #PEs if >0"};
  Gaudi::Property<int> m_PPD_n_pixels{this, "ppd_npix", 10000,
                                      "Total number of MPPC/SiPM pixels for implementation of saturation effect"};
};

#endif
