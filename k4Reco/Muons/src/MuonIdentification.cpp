#include "MuonIdentification.h"
#include "DD4hep/Detector.h"

#include "edm4hep/MutableReconstructedParticle.h"
#include "Objects/Helix.h"
#include "Pandora/Pandora.h"

#include "math.h"
#include <format>

#include "GaudiKernel/RndmGenerators.h"

/* ****************************************************************************
 *  Gaudi algorithm for Muon Identification
 * **************************************************************************** */

MuonIdentification::MuonIdentification(const std::string& name, ISvcLocator* svcLoc) :
  Transformer(name, svcLoc,
    {
      KeyValues("InputTrackCollection", { "TrackerHits" }),
      KeyValues("InputMuonHitCollection", { "MuonHits" })
    },
    { KeyValues("OutputMuonCollection", { "RecoMuons" }) }
  ),
  m_histos(*this)
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

  // --- Get the track time-of-flight correction
  m_tof_correction = new TFormula("tof_correction", m_formulaStr.value().c_str());
  
  // --- Get the magnetic field value
  const double pos[3] = {0., 0., 0.}; 
  double bFieldVec[3] = {0., 0., 0.}; 
  theDetector.field().magneticField(pos, bFieldVec);
  m_bField = bFieldVec[2] / dd4hep::tesla;

  m_histos.initialize();

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

  Rndm::Numbers gauss { randSvc(), Rndm::Gauss(0.0, m_timeResolution) };

  std::string initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!

  const size_t nHit = trackCollection.size(); 
  for (size_t idx1 = 0; idx1 < nHit ; idx1++)
  {
    edm4hep::Track track = trackCollection.at(idx1);
    edm4hep::TrackState ts_atIP = track.getTrackStates(edm4hep::TrackState::AtIP);
    edm4hep::TrackState ts_atCAL = track.getTrackStates(edm4hep::TrackState::AtCalorimeter);

    const float trk_pt = (m_lightSpeed / 1000000) * m_bField / std::fabs(ts_atIP.omega);
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
    int n_matchedHits = 0;
    for(size_t idx2 = 0; idx2 < muonHitCollection.size(); idx2++)
    {
      edm4hep::CalorimeterHit cHit = muonHitCollection.at(idx2);
      unsigned int systemID = bitFieldCoder.get(cHit.getCellID(), "system");
      unsigned int layerID  = bitFieldCoder.get(cHit.getCellID(), "layer");

      float x_centered = cHit.getPosition()[0] - x_ecal;
      float y_centered = cHit.getPosition()[1] - y_ecal;
      float z_centered = cHit.getPosition()[2] - z_ecal;

      float hit_r = std::sqrt(std::pow(x_centered, 2) + std::pow(y_centered, 2));
      float hit_phi = std::atan2(y_centered, x_centered); 
      float hit_theta = std::atan2(hit_r, z_centered);

      float dphi = hit_phi - trk_phi;
      if (dphi > M_PI) dphi -= 2. * M_PI;
      if (dphi < -M_PI) dphi += 2. * M_PI;
      float dtheta = hit_theta - trk_theta;

      float deltaR = std::sqrt(std::pow(dphi, 2) + std::pow(dtheta, 2));

      m_histos.fillDeltaR(systemID, layerID, deltaR, trk_pt, trk_theta);

      if (deltaR < m_deltaRMatch.value()[systemID - m_muonDetBarrel])
      {
        // --- Calculate the particle total time of flight to the muon detector hit
        float dist = std::sqrt(std::pow(x_centered, 2) + std::pow(y_centered, 2) + std::pow(z_centered, 2));
        float trk_p = trk_pt / sin(trk_theta);
        float trk_E = sqrt(std::pow(trk_p, 2) + std::pow(m_muonMass, 2));

        float beta = trk_p / trk_E;
        float tof = (trk_length + dist) / (beta * m_lightSpeed) + tof_correction;

        float deltaT = cHit.getTime() + gauss() - tof;

        m_histos.fillDeltaT(systemID, layerID, deltaT, trk_pt, trk_theta);

        // Time window selection: use barrel or endcap window depending on system
        if (systemID == m_muonDetBarrel)
        {
          if (deltaT < m_timeWindowB[0] || deltaT > m_timeWindowB[1]) continue;
        }
        else if (systemID == m_muonDetEndcap)
        {
          if (deltaT < m_timeWindowE[0] || deltaT > m_timeWindowE[1]) continue;
        }

        n_matchedHits++;
      }
    }

    m_histos.fillMatched(n_matchedHits);

    if (n_matchedHits >= m_nHitsMatch)
    {
      const float charge = ts_atIP.omega > 0.0 ? 1.0 : -1.0;
      const int pdg = charge < 0. ? m_muonPDG : -1 * m_muonPDG; 
      const float mom[3] = { 
        trk_pt * std::cos(ts_atIP.phi),
	      trk_pt * std::sin(ts_atIP.phi),
	      trk_pt * trk_cotTheta
      };
      float energy = std::pow(mom[0], 2) + std::pow(mom[1], 2) + std::pow(mom[2], 2);
      energy = std::sqrt(energy + std::pow(m_muonMass, 2));

      edm4hep::MutableReconstructedParticle muon;
      muon.setMomentum(mom);
      muon.setEnergy(energy);
      muon.setMass(m_muonMass);
      muon.setCharge(charge);
      muon.setPDG(pdg);
      muon.addToTracks(track);

      result.push_back(muon);
    }
  }

  return result;
}

/* ****************************************************************************
 *  Friend class for handling histograms
 * **************************************************************************** */

void MuonIDHistograms::initialize()
{
  if (m_algorithm.m_fillHistos.value())
  {
    m_hDeltaR_Barrel = m_algorithm.histoSvc()->book("MuonSystem", "hDeltaR_Barrel", 
        "#DeltaR between tracks and muon detector hits (barrel);#DeltaR [rad]", 1000, 0.0, 0.5);
    m_hDeltaR_Endcap = m_algorithm.histoSvc()->book("MuonSystem", "hDeltaR_Endcap", 
        "#DeltaR between tracks and muon detector hits (endcap);#DeltaR [rad]", 1000, 0.0, 0.5);
    m_hDeltaT_Barrel = m_algorithm.histoSvc()->book("MuonSystem", "hDeltaT_Barrel",
        "#DeltaT between tracks and muon detector hits (barrel);#DeltaT [ns]", 1000, -10., 10.);
    m_hDeltaT_Endcap = m_algorithm.histoSvc()->book("MuonSystem", "hDeltaT_Endcap",
        "#DeltaT between tracks and muon detector hits (endcap);#DeltaT [ns]", 1000, -10., 10.);

    for (unsigned int ih = 0; ih < m_algorithm.m_muonDetBarrelLayers; ++ih)
    {
      m_hDeltaR_B.push_back(m_algorithm.histoSvc()->book("MuonSystem",
          std::format("hDeltaR_B{}", ih).c_str(),
          std::format("#DeltaR between tracks and muon detector hits (layer B{});#DeltaT [ns]", ih).c_str(),
          1000, 0., 0.5));

      m_hDeltaT_B.push_back(m_algorithm.histoSvc()->book("MuonSystem",
          std::format("hDeltaT_B{}", ih).c_str(),
          std::format("#DeltaT between tracks and muon detector hits (layer B{});#DeltaT [ns]", ih).c_str(),
          1000, -10., 10.));
    }

    for (unsigned int ih = 0; ih < m_algorithm.m_muonDetEndcapLayers; ++ih)
    {
      m_hDeltaR_E.push_back(m_algorithm.histoSvc()->book("MuonSystem",
          std::format("hDeltaR_E{}", ih).c_str(),
          std::format("#DeltaR between tracks and muon detector hits (layer E{});#DeltaR [rad]", ih).c_str(),
          1000, 0., 0.5));

      m_hDeltaT_E.push_back(m_algorithm.histoSvc()->book("MuonSystem",
          std::format("hDeltaT_E{}", ih).c_str(),
          std::format("#DeltaT between tracks and muon detector hits (layer E{});#DeltaT [ns]", ih).c_str(),
          1000, -10., 10.));
    }

    m_hDeltaR_vs_Pt = m_algorithm.histoSvc()->book("MuonSystem", "hDeltaR_vs_Pt",
        "#DeltaR vs p_{T}^{trk};p_{T}^{trk} [GeV];#DeltaR [rad]",
			  200, 0., 100., 200, 0., 0.5);
    m_hDeltaT_vs_Pt = m_algorithm.histoSvc()->book("MuonSystem", "hDeltaT_vs_Pt",
        "#DeltaT vs p_{T}^{trk};p_{T}^{trk} [GeV];#DeltaT [ns]",
	      200, 0., 100., 200, -0.5, 3.5);
    m_hDeltaR_vs_Theta = m_algorithm.histoSvc()->book("MuonSystem", "hDeltaR_vs_Theta",
        "#DeltaR vs #theta_{trk};#theta_{trk} [#circ];#DeltaR [rad]",
				 200, 0., M_PI, 200, 0., 0.5);
    m_hDeltaT_vs_Theta = m_algorithm.histoSvc()->book("MuonSystem", "hDeltaT_vs_Theta",
        "#DeltaT vs #theta_{trk};#theta_{trk} [#circ];#DeltaT [ns]",
				200, 0., M_PI, 200, -0.5, 3.5);

    m_hNhits = m_algorithm.histoSvc()->book("MuonSystem", "hNhits",
        "Number of muon detector hits matched to tracks;N_{hits}", 20, 0., 20.);
  }
}

void MuonIDHistograms::fillDeltaR(unsigned int system, unsigned int layer, float deltaR, float trk_pt, float trk_theta)
{
  if (m_algorithm.m_fillHistos.value())
  {
    if (system == m_algorithm.m_muonDetBarrel)
    {
	    m_hDeltaR_Barrel->fill(deltaR);
	    m_hDeltaR_B[layer]->fill(deltaR);
    }
    else if (system == m_algorithm.m_muonDetEndcap)
    {
	    m_hDeltaR_Endcap->fill(deltaR);
	    m_hDeltaR_E[layer]->fill(deltaR);
    }
    else return;
	  m_hDeltaR_vs_Pt->fill(trk_pt, deltaR);
	  m_hDeltaR_vs_Theta->fill(trk_theta, deltaR);
  }
}

void MuonIDHistograms::fillDeltaT(unsigned int system, unsigned int layer, float deltaT, float trk_pt, float trk_theta)
{
  if (m_algorithm.m_fillHistos.value())
  {
    if (system == m_algorithm.m_muonDetBarrel)
    {
	    m_hDeltaT_Barrel->fill(deltaT);
	    m_hDeltaT_B[layer]->fill(deltaT);
    }
    else if (system == m_algorithm.m_muonDetEndcap)
    {
	    m_hDeltaT_Endcap->fill(deltaT);
	    m_hDeltaT_E[layer]->fill(deltaT);
    }
    else return;
    m_hDeltaT_vs_Pt->fill(trk_pt ,deltaT);
	  m_hDeltaT_vs_Theta->fill(trk_theta, deltaT);
  }
}

void MuonIDHistograms::fillMatched(int num)
{
  if (m_algorithm.m_fillHistos.value())
  {
    m_hNhits->fill(num);
  }
}

