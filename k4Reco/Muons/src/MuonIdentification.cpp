#include "MuonIdentification.h"
#include "DD4hep/Detector.h"

#include "Objects/Helix.h"
#include "Pandora/Pandora.h"

#include "math.h"

MuonIdentification::MuonIdentification(const std::string& name, ISvcLocator* svcLoc) :
  Transformer(name, svcLoc,
    {
      KeyValues("InputTrackCollection", { "TrackerHits" }),
      KeyValues("InputMuonHitCollection", { "MuonHits" })
    },
    { KeyValues("OutputMuonCollection", { "RecoMuons" }) }
  )
{}

StatusCode MuonIdentification::initialize()
{
  // --- Get the detector geometry
  m_geoSvc = serviceLocator()->service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  dd4hep::Detector& theDetector = *(m_geoSvc->getDetector());

  // --- Get the muon detector IDs
  m_muonDetBarrel = theDetector.constant<unsigned int>("DetID_Yoke_Barrel");
  m_muonDetEndcap = theDetector.constant<unsigned int>("DetID_Yoke_Endcap");

  // --- Get the ECAL barrel inner radius and the ECAL endcap minimum z
  m_ecalB_inner_r = theDetector.constant<float>("ECalBarrel_inner_radius") / dd4hep::mm;
  m_ecalE_min_z = theDetector.constant<float>("ECalEndcap_min_z") / dd4hep::mm;

  // Get Histogram and Data Services
  m_histSvc = serviceLocator()->service("THistSvc");

  // --- Get the track time-of-flight correction
  m_tof_correction = new TFormula("tof_correction", m_formulaStr.value().c_str());
  
  // --- Get the magnetic field value
  const double pos[3] = {0., 0., 0.}; 
  double bFieldVec[3] = {0., 0., 0.}; 
  theDetector.field().magneticField(pos, bFieldVec);
  m_bField = bFieldVec[2] / dd4hep::tesla;

  return StatusCode::SUCCESS;
}

StatusCode MuonIdentification::finalize()
{
  delete m_tof_correction;
  return StatusCode::SUCCESS;
}

edm4hep::ReconstructedParticleCollection MuonIdentification::operator()(
  const edm4hep::TrackCollection& trackCollection,
  const edm4hep::CalorimeterHitCollection& muonHitCollection
) const
{
  edm4hep::ReconstructedParticleCollection result;

  std::string initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!

  const size_t nHit = trackCollection.size(); 
  for (size_t idx1 = 0; idx1 < nHit ; idx1++)
  {
    edm4hep::Track track = trackCollection.at(idx1);
    edm4hep::TrackState ts_atIP = track.getTrackStates(edm4hep::TrackState::AtIP);
    edm4hep::TrackState ts_atCAL = track.getTrackStates(edm4hep::TrackState::AtCalorimeter);

    const float trk_pt = 0.000299792458 * m_bField / std::fabs(ts_atIP.omega);
    const float trk_cotTheta = ts_atIP.tanLambda;
    const float trk_phi = ts_atIP.phi;
    const float trk_theta = M_PI_2 - std::atan(trk_cotTheta);
    const float trk_d0 = ts_atIP.D0;
    const float trk_z0 = ts_atIP.Z0;

    // --- Track selection
    if (trk_pt < m_pt_min) continue; 
    if (std::fabs(trk_d0) > m_d0_max || std::fabs(trk_z0) > m_z0_max) continue;

    // --- Calulate the track length
    float deltaPhi = fabs(ts_atCAL.phi - trk_phi);
    if (deltaPhi > M_PI) deltaPhi -= 2. * M_PI;
    if (deltaPhi < -M_PI) deltaPhi += 2. * M_PI;
    const float trk_length = deltaPhi / std::fabs(ts_atIP.omega) *
      std::sqrt(1. + trk_cotTheta*trk_cotTheta);

    // --- Get the track intersection point with the ECAL inner surface
    float x_ecal = ts_atCAL.referencePoint[0];
    float y_ecal = ts_atCAL.referencePoint[1];
    float z_ecal = ts_atCAL.referencePoint[2];

    // If the reference point at the calorimeter is not available, extrapolate the track
    // to the ECAL inner surface using a helix
    float tof_correction = 0.;
    if (x_ecal == 0. || y_ecal == 0. || z_ecal == 0.)
    {
      if (std::fabs(ts_atIP.omega) < std::numeric_limits<float>::epsilon()) continue;

      pandora::Helix helix(trk_phi, trk_d0, trk_z0, ts_atIP.omega, trk_cotTheta, m_bField);

      const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());
      pandora::CartesianVector intersectionPoint(0.f, 0.f, 0.f);

      // If the track reaches the ECAL endcap:
      if (fabs(m_ecalB_inner_r * trk_cotTheta) > m_ecalE_min_z )
      {
        const float signPz((helix.GetMomentum().GetZ() > 0.f) ? 1.f : -1.f);
        helix.GetPointInZ(signPz * m_ecalE_min_z, referencePoint, intersectionPoint);
      }
      else // If the track reaches the ECAL barrel:
      {
        helix.GetPointOnCircle(m_ecalB_inner_r, referencePoint, intersectionPoint);
      }

      x_ecal = intersectionPoint.GetX();
      y_ecal = intersectionPoint.GetY();
      z_ecal = intersectionPoint.GetZ();

      // A correction to the particle time of flight is required to account for
      // the approximate analytical extrapolation with the helix
      tof_correction = m_tof_correction->Eval(trk_theta);
    }

    // --- Loop over the muon detector hits
    for(size_t idx2 = 0; idx2 < muonHitCollection.size(); idx2++)
    {
      edm4hep::CalorimeterHit cHit = muonHitCollection.at(idx2);
      unsigned int systemID = bitFieldCoder.get(cHit.getCellID(), "system");
      unsigned int layerID  = bitFieldCoder.get(cHit.getCellID(), "layer");
    }
  }

  return result;
}


