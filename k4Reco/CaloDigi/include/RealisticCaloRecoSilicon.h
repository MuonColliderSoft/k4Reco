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
#ifndef REALISTICCALORECOSILICON_H
#define REALISTICCALORECOSILICON_H 1

#include "RealisticCaloReco.h"

/** === RealisticCaloRecoSilicon Processor === <br>
    realistic reconstruction of silicon calorimeter hits
    D.Jeans 02/2016.

    24 March 2016: removed gap corrections - to be put into separate processor

*/

struct RealisticCaloRecoSilicon final : RealisticCaloReco {

public:
  RealisticCaloRecoSilicon(const std::string& name, ISvcLocator* svcLoc);

protected:
  float reconstructEnergy(const edm4hep::CalorimeterHit* hit, int layer) const override;
};

#endif
