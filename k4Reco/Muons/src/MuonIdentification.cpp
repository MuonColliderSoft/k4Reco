#include "MuonIdentification.h"
#include "DD4hep/Detector.h"

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
  _muonDetBarrel = theDetector.constant<unsigned int>("DetID_Yoke_Barrel");
  _muonDetEndcap = theDetector.constant<unsigned int>("DetID_Yoke_Endcap");

  // Get Histogram and Data Services
	m_histSvc = serviceLocator()->service("THistSvc");

  // --- Get the track time-of-flight correction
  _tof_correction = new TFormula("tof_correction", m_formulaStr.value().c_str());
  
  // --- Get the magnetic field value
  const double pos[3] = {0., 0., 0.}; 
  double bFieldVec[3] = {0., 0., 0.}; 
  theDetector.field().magneticField(pos,bFieldVec);
  _bField = bFieldVec[2]/dd4hep::tesla;

  return StatusCode::SUCCESS;
}

StatusCode MuonIdentification::finalize()
{
  delete _tof_correction;
  return StatusCode::SUCCESS;
}

edm4hep::ReconstructedParticleCollection MuonIdentification::operator()(
  const edm4hep::TrackerHitPlaneCollection& trackerHitCollection,
  const edm4hep::CalorimeterHitCollection& muonHitCollection
) const
{
  edm4hep::ReconstructedParticleCollection result;

  const size_t nHit = trackerHitCollection.size(); 
  for (size_t idx1 = 0; idx1 < nHit ; idx1++)
  {
    edm4hep::TrackerHitPlane tHit = trackerHitCollection.at(idx1);

    for(size_t idx2 = 0; idx2 < muonHitCollection.size(); idx2++)
    {
      edm4hep::CalorimeterHit cHit = muonHitCollection.at(idx2);
    }
  }

  return result;
}


