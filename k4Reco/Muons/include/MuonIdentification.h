#pragma once

#include "Gaudi/Property.h"
#include "GaudiKernel/ITHistSvc.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"

#include <edm4hep/TrackCollection.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/ReconstructedParticleCollection.h>

#include "TFormula.h"
#include "TH1F.h"
#include "TH2F.h"

using TransformType = edm4hep::ReconstructedParticleCollection(
                        const edm4hep::TrackCollection&,
                        const edm4hep::CalorimeterHitCollection&
                      );

struct MuonIdentification : public k4FWCore::Transformer<TransformType>
{
public:

  MuonIdentification(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  edm4hep::ReconstructedParticleCollection operator()(
    const edm4hep::TrackCollection& trackCollection,
    const edm4hep::CalorimeterHitCollection& muonHitCollection
  ) const override;
  StatusCode finalize() override;

protected:
  Gaudi::Property<float> m_pt_min { this, "TrackPtMin", 1.0, "Minimum track transverse momentum in GeV" };
  Gaudi::Property<float> m_d0_max { this, "TrackD0Max", 0.1, "Maximum track transverse impact parameter" };
  Gaudi::Property<float> m_z0_max { this, "TrackZ0Max", 0.1, "Maximum track longitudinal impact parameter" };
  Gaudi::Property<float> m_timeResolution { this, "MuonHitsTimeResolution", 0.1,
                                    "Time resolution of the muon detector hits" };
  Gaudi::Property<std::vector<float>> m_xyzResolution { this, "MuonHitsSpatialResolution", {10., 10., 10.},
                                    "Time resolution of the muon detector hits" };
  Gaudi::Property<std::vector<float>> m_deltaRMatch { this, "DeltaRMatch", { 0.2, 0.3 },
                                   "DeltaR for matching tracks with muon detector hits in barrel and endcaps" };
  Gaudi::Property<std::vector<float>> m_timeWindowB { this, "BarrelHitsTimeWindow", { -0.3, 0.3 }, 
                                   "Time window for hits in the muon detecor barrel in ns" };
  Gaudi::Property<std::vector<float>> _timeWindowE { this, "EndcapHitsTimeWindow", { -0.3, 0.3 },
                                   "Time window for hits in the muon detecor endcaps in ns" };
  Gaudi::Property<int> m_nHitsMatch { this, "NHitsMatch", 4, "Minumum number of matching hits in the muon detector" };
  Gaudi::Property<bool> m_fillHistos { this, "FillHistograms", false, "Fill the diagnostic histograms" };

  Gaudi::Property<std::string> m_encodingStringVariable { this, "EncodingStringParameterName",
                            "GlobalTrackerReadoutID",
                            "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  Gaudi::Property<std::string> m_formulaStr { this, "TrackTofCorrFunction",
                                             "(x < 0.5735 || x > 2.57) ? 0.015058 : "
                                             "(x >= 0.5735 && x < 0.626) ? (-9.63278 + 16.8561*x) : "
                                             "(x >= 0.626 && x < 2.516) ? (2.63371 - 4.14427*x + "
                                             "3.42931*x*x - 1.3498*x*x*x + 0.216549*x*x*x*x) : "
                                             "(x >= 2.516 && x < 2.57) ? (43.5031 - 16.9276*x) : 0",
                                             "Correction function for the track time of flight" };

private:
  SmartIF<IGeoSvc>   m_geoSvc;
  SmartIF<ITHistSvc> m_histSvc;

  unsigned int m_muonDetBarrel;
  unsigned int m_muonDetEndcap;

  float m_ecalB_inner_r;
  float m_ecalE_min_z;

  float m_bField = 5.;

  TFormula* m_tof_correction = nullptr;
};

DECLARE_COMPONENT(MuonIdentification)
