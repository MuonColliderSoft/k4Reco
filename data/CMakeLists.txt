#[[
Copyright (c) 2020-2024 Key4hep-Project.

This file is part of Key4hep.
See https://key4hep.github.io/key4hep-doc/ for further info.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
]]

# ECAL thresholds
INSTALL(FILES ECAL_Thresholds_3TeV.root ECAL_Thresholds_10TeV.root
  DESTINATION ${CMAKE_INSTALL_DATADIR}/k4Reco/data)

# HCAL thresholds
INSTALL(FILES HCAL_Thresholds_3TeV.root HCAL_Thresholds_10TeV.root
  DESTINATION ${CMAKE_INSTALL_DATADIR}/k4Reco/data)
