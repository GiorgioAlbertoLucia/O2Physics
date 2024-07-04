// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// TableProducer to generate Trees for pure protons ans pions (from Lambda), Kaons (from Omegas), deuterons (identified with TOF) and He3 (identified with TPC).
// The output trees contain the ITS cluster size information.
//
// Author: Giorgio Alberto Lucia

#include <vector>
#include <utility>
#include <random>
#include <cstdint>
#include <cstring>
#include <unordered_map>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/PIDResponse.h"
#include "DCAFitter/DCAFitterN.h"

#include "PWGLF/DataModel/LFClusterStudiesTable.h"

#include "TDatabasePDG.h"

using namespace ::o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using Track = o2::track::TrackParCov;
using TracksFullIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr>;

using Vec3D = ROOT::Math::SVector<double, 3>;
using momentum = std::array<float, 3>;
using position = std::array<float, 3>;

namespace physics
{
  namespace pdg
  {
    constexpr int Electron = 11;
    constexpr int Pion = 211;
    constexpr int KaonCharged = 321;
    constexpr int Proton = 2212;
    constexpr int Deuteron = 1000010020;
    constexpr int He3 = 1000020030;
    constexpr int Lambda = 3122;
    constexpr int AntiLambda = -3122;
    constexpr int K0s = 310;
    constexpr int XiMinus = 3312;
    constexpr int OmegaMinus = 3334;
  }
}

/**
 * Decay2Bodies class to compute the properties of a 2-body decay
*/
template <class T>
class Decay2Bodies
{
  public:

    Decay2Bodies(T& daughterTrackA, T& daughterTrackB, const position& PV, o2::vertexing::DCAFitterN<2>& fitter) : 
      m_fitter(fitter), m_fitSuccess(true)
    {
      
      m_PV = {PV[0], PV[1], PV[2]};
      auto daughterTrackCovarianceA = getTrackParCov(daughterTrackA);
      auto daughterTrackCovarianceB = getTrackParCov(daughterTrackB);

      int nCand = 0;
        try {
          nCand = m_fitter.process(daughterTrackCovarianceA, daughterTrackCovarianceB);
        } catch (...) {
          LOG(error) << "Exception caught in DCA fitter process call!";
          m_fitSuccess = false;
        }
        if (nCand == 0) {
          m_fitSuccess = false;
        }

      m_trackA = m_fitter.getTrack(0);
      m_trackB = m_fitter.getTrack(1);

      m_trackA.getPxPyPzGlo(m_momA);
      m_trackB.getPxPyPzGlo(m_momB);
      computeMomMother();

      m_decayVtx = m_fitter.getPCACandidate();

      gpu::gpustd::array<float, 2> dcaInfo;
      o2::base::Propagator::Instance()->propagateToDCABxByBz({m_PV[0], m_PV[1], m_PV[2]}, daughterTrackCovarianceA, 2.f, m_fitter.getMatCorrType(), &dcaInfo);
      m_AdcaToPV = std::hypot(dcaInfo[0], dcaInfo[1]);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({m_PV[0], m_PV[1], m_PV[2]}, daughterTrackCovarianceB, 2.f, m_fitter.getMatCorrType(), &dcaInfo);
      m_BdcaToPV = std::hypot(dcaInfo[0], dcaInfo[1]);
    }

    Track getTrackA() const { return m_trackA; }
    Track getTrackB() const { return m_trackB; }

    float p(const char * particle = "mother") const 
    {
      if (strcmp(particle, "mother") == 0) {
        return std::hypot(m_momMother[0], m_momMother[1], m_momMother[2]);
      } else if (strcmp(particle, "A") == 0) {
        return std::hypot(m_momA[0], m_momA[1], m_momA[2]);
      } else if (strcmp(particle, "B") == 0) {
        return std::hypot(m_momB[0], m_momB[1], m_momB[2]);
      } else {
        return -999.f;
      }
    }

    float pt(const char * particle = "mother") const
    {
      if (strcmp(particle, "mother") == 0) {
        return std::hypot(m_momMother[0], m_momMother[1]);
      } else if (strcmp(particle, "A") == 0) {
        return std::hypot(m_momA[0], m_momA[1]);
      } else if (strcmp(particle, "B") == 0) {
        return std::hypot(m_momB[0], m_momB[1]);
      } else {
        return -999.f;
      }
    }

    float eta(const char * particle = "mother") const 
    {
      float mom = 0.f;
      float momZ = 0.f;
      if (strcmp(particle, "mother") == 0) {
        mom = std::hypot(m_momMother[0], m_momMother[1], m_momMother[2]);
        momZ = m_momMother[2];
      } else if (strcmp(particle, "A") == 0) {
        mom = std::hypot(m_momA[0], m_momA[1], m_momA[2]);
        momZ = m_momA[2];
      } else if (strcmp(particle, "B") == 0) {
        mom = std::hypot(m_momB[0], m_momB[1], m_momB[2]);
        momZ = m_momB[2];
      } else {
        return -999.f;
      }
      float proposedEta = 0.5 * std::log((mom + momZ) / (mom - momZ));
      if (std::abs(proposedEta) < static_cast<float>(1e-7)) {
        if (momZ < 0) {
          return -100.f;
        } else {
          return 100.f;
        }
      }
      return proposedEta;
    }

    float phi(const char * particle = "mother") const 
    {
      if (strcmp(particle, "mother") == 0) {
        return std::atan2(m_momMother[1], m_momMother[0]);
      } else if (strcmp(particle, "A") == 0) {
        return std::atan2(m_momA[1], m_momA[0]);
      } else if (strcmp(particle, "B") == 0) {
        return std::atan2(m_momB[1], m_momB[0]);
      } else {
        return -999.f;
      }
    }

    /**
     * Compute the mass of the mother, given two mass hypotheses for the daughters
     *
     * @param massA mass hypothesis for the first daughter
     * @param massB mass hypothesis for the second daughter
     */
    float massMother(const float massA, const float massB) const
    {
      float momA = std::hypot(m_momA[0], m_momA[1], m_momA[2]);
      float momB = std::hypot(m_momB[0], m_momB[1], m_momB[2]);
      float eA = std::hypot(massA, momA);
      float eB = std::hypot(massB, momB);
      float momMother = std::hypot(m_momMother[0], m_momMother[1], m_momMother[2]);
      float eMother = eA + eB;
      return std::sqrt(eMother * eMother - momMother * momMother);
    }

    float radiusMother() const
    {
      return std::hypot(m_decayVtx[0], m_decayVtx[1]);
    }

    /**
     * Compute the DCA of the daughters
    */
    float dcaDaughters() const
    {
      return std::sqrt(m_fitter.getChi2AtPCACandidate());
    }

    /**
     * Compute the DCA of the mother from the primary vertex
     *
     * @param PV primary vertex position
     */
    float dcaMotherToPV() const
    {
      std::array<float, 3> relPos = {m_decayVtx[0] - m_PV[0], m_decayVtx[1] - m_PV[1], m_decayVtx[2] - m_PV[2]};
      float momMother = std::hypot(m_momMother[0], m_momMother[1], m_momMother[2]);
      return std::sqrt((std::pow(relPos[1] * m_momMother[2] - relPos[2] * m_momMother[1], 2) + std::pow(relPos[2] * m_momMother[0] - relPos[0] * m_momMother[2], 2) + std::pow(relPos[0] * m_momMother[1] - relPos[1] * m_momMother[0], 2))) / momMother;
    }

    /**
     * Compute the DCA of the daughters from the primary vertex
     * 
     * @param whichDaughter "A" or "B"
    */
    float dcaDaughtersToPV(const char * whichDaughter = "A") const
    {
      if (strcmp(whichDaughter, "A") == 0) {
        return m_AdcaToPV;
      } else if (strcmp(whichDaughter, "B") == 0) {
        return m_BdcaToPV;
      } else {
        return -999.f;
      }
    }

    /**
     * Return the cosine of the pointing angle of the mother
     * 
    */
    float cosPA() const { return RecoDecay::cpa(m_PV, m_decayVtx, m_momMother); }

    /**
     * Compute the alpha for the Armenteros-Podolanski plot
     */
    float alphaAP() const
    {
      float momTot = std::hypot(m_momMother[0], m_momMother[1], m_momMother[2]);
      float lQlA = (m_momA[0] * m_momMother[0] + m_momA[1] * m_momMother[1] + m_momA[2] * m_momMother[2]) / momTot;
      float lQlNeg = (m_momB[0] * m_momMother[0] + m_momB[1] * m_momMother[1] + m_momB[2] * m_momMother[2]) / momTot;
      return (lQlA - lQlNeg) / (lQlA + lQlNeg);
    }

    /**
     * Compute the qt for the Armenteros-Podolanski plot
     */
    float qtAP() const
    {
      float dp = m_momMother[0] * m_momA[0] + m_momMother[1] * m_momA[1] + m_momMother[2] * m_momA[2];
      float p2V0 = std::hypot(m_momMother[0], m_momMother[1], m_momMother[2]);
      float p2A = std::hypot(m_momA[0], m_momA[1], m_momA[2]);
      return std::sqrt(p2A - dp * dp / p2V0);
    }


    bool fitSuccess() const { return m_fitSuccess; }

  protected:
    void computeMomMother()
    {
      m_momMother[0] = m_momA[0] + m_momB[0];
      m_momMother[1] = m_momA[1] + m_momB[1];
      m_momMother[2] = m_momA[2] + m_momB[2];
    }

  private:
    o2::vertexing::DCAFitterN<2> m_fitter;

    std::array<float, 3> m_momA;      // momentum of the first daughter
    std::array<float, 3> m_momB;      // momentum of the second daughter
    std::array<float, 3> m_momMother; // momentum of the mother

    std::array<float, 3> m_PV;        // primary vertex
    Vec3D m_decayVtx;

    Track m_trackA;
    Track m_trackB;

    float m_AdcaToPV;
    float m_BdcaToPV;

    bool m_fitSuccess;

};

class BitMap
{
  public:
    BitMap() = default;

    BitMap(const std::vector<std::string>& variables) 
    {
      for (const auto& var : variables) {
        m_bits[var] = false;
      }
    }

    void set(const std::string& var)
    {
      m_bits[var] = 1;
    }

    void clear(const std::string& var)
    {
      m_bits[var] = 0;
    }

    bool get(const std::string& var) const
    {
      return m_bits.at(var);
    }

    /**
     * Check if the candidate is ambiguous (two candidates are true)
    */
    bool isAmbiguous() const
    {
      int count = 0;
        for (const auto& pair : m_bits) {
            if (pair.second) {
                count++;
                if (count > 1) {
                    return true;
                }
            }
        }
        return false;
    }

  private:
    std::unordered_map<std::string, bool> m_bits;
};

/**
 * Struct with the information of a candidate V0 to be stored in the tree
*/
struct CandidateV0 {

  // mother information
  int64_t globalIndex_V0 = -999;

  float p_V0 = -999.f;
  float pt_V0 = -999.f;
  float eta_V0 = -999.f;
  float phi_V0 = -999.f;
  
  float mass_V0 = -999.f;
  int pdgCode_V0 = -999;
  
  float radius_V0 = -999.f;
  float dcaPV_V0 = -999.f;
  float cosPA_V0 = -999.f;
  float alphaAP_V0 = -999.f;
  float qtAP_V0 = -999.f;

  // pos daughter information
  int64_t globalIndex_pos = -999;
  
  float p_pos = -999.f;
  float pt_pos = -999.f;
  float eta_pos = -999.f;
  float phi_pos = -999.f;
  float pTPC_pos = -999.f;
  
  int pdgCode_pos = -999;
  
  float dcaPV_pos = -999.f;
  uint32_t itsClSize_pos = 0u;
  uint16_t tpcSignal_pos = 0u;
  uint16_t nTPCcls_pos = 0u;

  float nSigmaTPCe_pos = -999.f;
  float nSigmaTPCpi_pos = -999.f;
  float nSigmaTPCk_pos = -999.f;
  float nSigmaTPCp_pos = -999.f;
  float nSigmaTPCdeu_pos = -999.f;
  float nSigmaTPChe3_pos = -999.f;

  float nSigmaTOFe_pos = -999.f;
  float nSigmaTOFpi_pos = -999.f;
  float nSigmaTOFk_pos = -999.f;
  float nSigmaTOFp_pos = -999.f;
  float nSigmaTOFdeu_pos = -999.f;
  float nSigmaTOFhe3_pos = -999.f;


  float chi2its_pos = -999.f;
  float chi2tpc_pos = -999.f;
  float chi2matchingITSTPC_pos = -999.f;


  // neg daughter information
  int64_t globalIndex_neg = -999;

  float p_neg = -999.f;
  float pt_neg = -999.f;
  float eta_neg = -999.f;
  float phi_neg = -999.f;
  float pTPC_neg = -999.f;
  
  int pdgCode_neg = -999;
  
  float dcaPV_neg = -999.f;
  uint32_t itsClSize_neg = 0u;
  uint16_t tpcSignal_neg = 0u;
  uint16_t nTPCcls_neg = 0u;

  float nSigmaTPCe_neg = -999.f;
  float nSigmaTPCpi_neg = -999.f;
  float nSigmaTPCk_neg = -999.f;
  float nSigmaTPCp_neg = -999.f;
  float nSigmaTPCdeu_neg = -999.f;
  float nSigmaTPChe3_neg = -999.f;

  float nSigmaTOFe_neg = -999.f;
  float nSigmaTOFpi_neg = -999.f;
  float nSigmaTOFk_neg = -999.f;
  float nSigmaTOFp_neg = -999.f;
  float nSigmaTOFdeu_neg = -999.f;
  float nSigmaTOFhe3_neg = -999.f;

  float chi2its_neg = -999.f;
  float chi2tpc_neg = -999.f;
  float chi2matchingITSTPC_neg = -999.f;
};

struct CandidateCascade {
  // mother information
  int64_t globalIndex_mother = -999;

  float p_mother = -999.f;
  float pt_mother = -999.f;
  float eta_mother = -999.f;
  float phi_mother = -999.f;
  
  float mass_mother = -999.f;
  int pdgCode_mother = -999;
  
  float radius_mother = -999.f;
  float dcaPV_mother = -999.f;
  float cosPA_mother = -999.f;
  float alphaAP_mother = -999.f;
  float qtAP_mother = -999.f;

  // K daughter information
  int64_t globalIndex_K = -999;
  
  float p_K = -999.f;
  float pt_K = -999.f;
  float eta_K = -999.f;
  float phi_K = -999.f;
  float pTPC_K = -999.f;
  
  int pdgCode_K = -999;

  float dcaPV_K = -999.f;
  uint32_t itsClSize_K = 0u;
  uint16_t tpcSignal_K = 0u;
  uint16_t nTPCcls_K = 0u;

  float nSigmaTPCe_K = -999.f;
  float nSigmaTPCpi_K = -999.f;
  float nSigmaTPCk_K = -999.f;
  float nSigmaTPCp_K = -999.f;
  float nSigmaTPCdeu_K = -999.f;
  float nSigmaTPChe3_K = -999.f;

  float nSigmaTOFe_K = -999.f;
  float nSigmaTOFpi_K = -999.f;
  float nSigmaTOFk_K = -999.f;
  float nSigmaTOFp_K = -999.f;
  float nSigmaTOFdeu_K = -999.f;
  float nSigmaTOFhe3_K = -999.f;

  float chi2its_K = -999.f;
  float chi2tpc_K = -999.f;
  float chi2matchingITSTPC_K = -999.f;
};

struct LFTreeCreator {

  Produces<o2::aod::ClStV0Table> m_v0Table;
  Produces<o2::aod::ClStCascTable> m_KTable;

  Service<o2::ccdb::BasicCCDBManager> m_ccdb;
  int m_runNumber;
  float m_d_bz;
  uint32_t m_randomSeed = 0.;

  Configurable<bool> setting_fillV0{"fillV0", true, "Fill the V0 tree"};
  Configurable<bool> setting_fillK{"fillK", true, "Fill the K tree"};

  Configurable<int> setting_materialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

  Configurable<float> setting_zVtxMax{"zVtxMax", 10.f, "Maximum z vertex position"};

  Configurable<float> setting_downscaleFactor{"downscaleFactor", 1.f, "Downscale factor for the V0 candidates"};
  Configurable<bool> setting_applyAdditionalEvSel{"applyAdditionalEvSel", false, "Apply additional event selection"};

  Configurable<float> v0track_nClsItsMin{"v0track_NclsItsMin", 0.f, "Minimum number of ITS clusters for the V0 daughters"};
  Configurable<float> v0track_nClsTpcMin{"v0track_NclsTpcMin", 100.f, "Minimum number of TPC clusters for the V0 daughters"};
  Configurable<float> v0track_nClsTpcMaxShared{"v0track_NclsTpcMaxShared", 5.f, "Maximum number of shared TPC clusters for the V0 daughters"};

  //Configurable<float> v0setting_etaMaxV0{"etaMaxV0", 0.8f, "Maximum eta for the V0 daughters"};
  Configurable<float> v0setting_etaMaxV0dau{"etaMaxV0dau", 0.8f, "Maximum eta for the V0 daughters"};
  Configurable<float> v0setting_dcaV0daughters{"v0setting_dcaV0daughters", 0.5f, "DCA between the V0 daughters"};
  Configurable<float> v0setting_dcaV0toPV{"v0setting_dcaV0fromPV", 1.f, "DCA of the V0 to the primary vertex"};
  Configurable<float> v0setting_dcaDaughtersToPV{"v0setting_dcaDaughtersToPV", 1.f, "DCA of the daughters to the primary vertex"};
  Configurable<float> v0setting_radiusMax{"v0setting_radiusMax", 100.f, "Maximum radius of the V0 accepted"};
  Configurable<float> v0setting_radiusMin{"v0setting_radiusMin", 5.f, "Minimum radius of the V0 accepted"};
  Configurable<float> v0setting_cosPA{"v0setting_cosPA", 0.99f, "Cosine of the pointing angle of the V0"};
  Configurable<float> v0setting_nsigmatpc{"v0setting_nsigmaTPC", 4.f, "Number of sigmas for the TPC PID"};
  Configurable<float> v0setting_massWindowLambda{"v0setting_massWindowLambda", 0.02f, "Mass window for the Lambda"};
  Configurable<float> v0setting_massWindowK0s{"v0setting_massWindowK0s", 0.02f, "Mass window for the K0s"};
  
  Configurable<float> cascsetting_dcaCascDaughters{"v0setting_dcaV0daughters", 0.1f, "DCA between the V0 daughters"};
  Configurable<float> cascsetting_cosPA{"v0setting_cosPA", 0.99f, "Cosine of the pointing angle of the V0"};
  Configurable<float> cascsetting_massWindowOmega{"v0setting_massWindowOmega", 0.01f, "Mass window for the Omega"};
  Configurable<float> cascsetting_massWindowXi{"v0setting_massWindowXi", 0.01f, "Mass window for the Xi"};

  Preslice<aod::V0s> m_perCollisionV0 = o2::aod::v0::collisionId;
  Preslice<aod::Cascades> m_perCollisionCascade = o2::aod::cascade::collisionId;

  HistogramRegistry m_hCheck{
    "LFTreeCreator",
    {
      {"armenteros_plot_before_selections", "Armenteros-Podolanski plot; #alpha; q_{T} (GeV/c)", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 0.2f}}}},
      {"armenteros_plot", "Armenteros-Podolanski plot; #alpha; q_{T} (GeV/c)", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 0.2f}}}},
      {"Xi_vs_Omega", "Mass Xi vs Omega; mass Omega (GeV/c^{2}); mass Xi (GeV/c^{2})", {HistType::kTH2F, {{100, 1.f, 2.f}, {100, 1.f, 2.f}}}},
      {"massOmega", "Mass Omega; mass Omega (GeV/c^{2}); counts", {HistType::kTH1F, {{100, 1.f, 2.f}}}},
      {"massOmegaWithBkg", "Mass Omega with Background; mass Omega (GeV/c^{2}); counts", {HistType::kTH1F, {{100, 1.f, 2.f}}}},
      {"massLambda", "Mass Lambda; mass Lambda (GeV/c^{2}); counts", {HistType::kTH1F, {{100, 1.f, 2.f}}}},
      {"zVtx", "Binning for the vertex z in cm", {HistType::kTH1F, {{100, -20.f, 20.f}}}}
    },
    OutputObjHandlingPolicy::AnalysisObject, 
    false, 
    true};  // check histograms

  std::vector<CandidateV0> m_candidateV0s;
  std::vector<CandidateCascade> m_candidateCascades;
  o2::vertexing::DCAFitterN<2> m_fitter;

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (m_runNumber == bc.runNumber()) {
      return;
    }

    auto timestamp = bc.timestamp();
    o2::parameters::GRPMagField* grpmag = 0x0;

    auto grpmagPath{"GLO/Config/GRPMagField"};
    grpmag = m_ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
    if (!grpmag) {
      LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField for timestamp " << timestamp;
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);

    // Fetch magnetic field from ccdb for current collision
    m_d_bz = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << m_d_bz << " kG";
    m_runNumber = bc.runNumber();
    m_fitter.setBz(m_d_bz);

    // o2::base::Propagator::Instance()->setMatLUT(lut);
  }

  void init(o2::framework::InitContext& ic)
  {
    m_runNumber = 0;
    m_d_bz = 0;

    m_ccdb->setURL("http://alice-ccdb.cern.ch");
    m_ccdb->setCaching(true);
    m_ccdb->setLocalObjectValidityChecking();
    m_ccdb->setFatalWhenNull(false);
    // lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    m_fitter.setPropagateToPCA(true);
    m_fitter.setMaxR(200.);
    m_fitter.setMinParamChange(1e-3);
    m_fitter.setMinRelChi2Change(0.9);
    m_fitter.setMaxDZIni(4);
    m_fitter.setMaxDXYIni(4);
    m_fitter.setMaxChi2(1e9);
    m_fitter.setUseAbsDCA(true);
    m_fitter.setWeightedFinalPCA(false);
    int mat{static_cast<int>(setting_materialCorrection)};
    m_fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

  }

  /**
   * Select the V0 daughters based on the quality cuts
  */
  template <typename T>
  bool qualitySelectionV0Daughter(const T& track)
  {
    if (std::abs(track.eta()) > v0setting_etaMaxV0dau) {
      return false;
    }
    if (track.itsNCls() < v0track_nClsItsMin ||
        track.tpcNClsFound() < v0track_nClsTpcMin ||
        track.tpcNClsCrossedRows() < v0track_nClsTpcMin ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcNClsShared() > v0track_nClsTpcMaxShared) {
      return false;
    }
    return true;
  }

  template<class C, class T>
  void fillV0Tree(const C& collision, const T& /*tracks*/, const aod::V0s& v0s) 
  {
    m_candidateV0s.clear();

    // V0s
    for (const auto& v0: v0s) {
      auto posTrack = v0.posTrack_as<T>();
      auto negTrack = v0.negTrack_as<T>();

      bool qualityTestPos = qualitySelectionV0Daughter(posTrack);
      bool qualityTestNeg = qualitySelectionV0Daughter(negTrack);

      if (!qualityTestPos || !qualityTestNeg) {
        continue;
      }

      std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};
      Decay2Bodies<T> decayV0(posTrack, negTrack, PV, m_fitter);
      if (!decayV0.fitSuccess()) {
        continue;
      }

      m_hCheck.fill(HIST("armenteros_plot_before_selections"), decayV0.alphaAP(), decayV0.qtAP());

      // V0 selections
      float ptV0 = decayV0.pt();
      /*
      if (ptV0 < v0setting_ptMin || ptV0 > v0setting_ptMax) {
        continue;
      }
      */
      /*
      if (std::abs(decayV0.eta()) > v0setting_etaMaxV0) {
        continue;
      }
      */
      if (decayV0.dcaDaughters() > v0setting_dcaV0daughters) {
        continue;
      }
      float radiusV0 = decayV0.radiusMother();
      if (radiusV0 > v0setting_radiusMax || radiusV0 < v0setting_radiusMin) {
        continue;
      }
      if (std::abs(decayV0.dcaMotherToPV()) > v0setting_dcaV0toPV) {
        continue;
      }
      if (decayV0.cosPA() < v0setting_cosPA) {
        continue;
      }

      // mass hypothesis
      float massLambdaV0 = decayV0.massMother(o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged);
      float massK0sV0 = decayV0.massMother(o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged);

      BitMap selectionMap({"Lambda", "AntiLambda", "K0s"});
      if (std::abs(massK0sV0 - o2::constants::physics::MassK0Short) > v0setting_massWindowK0s) {
        selectionMap.clear("K0s");
      }
      if ((std::abs(massLambdaV0 - o2::constants::physics::MassLambda0) > v0setting_massWindowLambda) && (decayV0.alphaAP() > 0)) { 
        selectionMap.clear("Lambda");
      }
      if ((std::abs(massLambdaV0 - o2::constants::physics::MassLambda0) > v0setting_massWindowLambda) && (decayV0.alphaAP() < 0)) {
        selectionMap.clear("AntiLambda");
      }
      /*
      if (TESTBIT(selectionMap, Lambda) && ( ((posTrack.hasTPC()) && (std::abs(posTrack.tpcNSigmaPr()) > v0setting_nsigmatpc)) || ((negTrack.hasTPC()) && (std::abs(negTrack.tpcNSigmaPi()) > v0setting_nsigmatpc)))) {
        CLRBIT(selectionMap, Lambda);
      }
      */
      /*
      if (TESTBIT(selectionMap, AntiLambda) && ( ((posTrack.hasTPC()) && (std::abs(posTrack.tpcNSigmaPi()) > v0setting_nsigmatpc)) || ((negTrack.hasTPC()) && (std::abs(negTrack.tpcNSigmaPr()) > v0setting_nsigmatpc)))) {
        CLRBIT(selectionMap, AntiLambda);
      }
      */

      float massV0 = -999.f;
      int pdgCodeV0 = -999;
      int pdgCodePos = -999;
      int pdgCodeNeg = -999;
      if (selectionMap.get("Lambda")) {
        massV0 = massLambdaV0;
        pdgCodeV0 = physics::pdg::Lambda;
        pdgCodePos = physics::pdg::Proton;
        pdgCodeNeg = -physics::pdg::Pion;
      } else if (selectionMap.get("AntiLambda")) {
        massV0 = massLambdaV0;
        pdgCodeV0 = -physics::pdg::Lambda;
        pdgCodePos = -physics::pdg::Proton;
        pdgCodeNeg = physics::pdg::Pion;
      } else if (selectionMap.get("K0s")) {
        massV0 = massK0sV0;
        pdgCodeV0 = physics::pdg::K0s;
        pdgCodePos = physics::pdg::Pion;
        pdgCodeNeg = -physics::pdg::Pion;
      } else {
        continue;
      }

      // skip ambiguous V0s
      if (selectionMap.isAmbiguous()) {
        continue;
      }

      auto dcaToPVpos = decayV0.dcaDaughtersToPV("A");
      if (dcaToPVpos < v0setting_dcaDaughtersToPV /*&& std::abs(dcaInfo[0]) < v0setting_dcaDaughtersToPV*/) {
        continue;
      }      
      float dcaToPVneg = decayV0.dcaDaughtersToPV("B");
      if (dcaToPVneg < v0setting_dcaDaughtersToPV /*&& std::abs(dcaInfo[0]) < v0setting_dcaDaughtersToPV*/) {
        continue;
      }

      // fill the candidate and add it to the Armenteros plot
      m_hCheck.fill(HIST("armenteros_plot"), decayV0.alphaAP(), decayV0.qtAP());

      auto posPropagatedTrack = decayV0.getTrackA();
      auto negPropagatedTrack = decayV0.getTrackB();

      CandidateV0 candidate;
      candidate.globalIndex_V0 = v0.globalIndex();
      candidate.p_V0 = decayV0.p();
      candidate.pt_V0 = decayV0.pt();
      candidate.eta_V0 = decayV0.eta();
      candidate.phi_V0 = decayV0.phi();
      candidate.mass_V0 = massV0;
      candidate.pdgCode_V0 = pdgCodeV0;
      candidate.radius_V0 = radiusV0;
      candidate.dcaPV_V0 = decayV0.dcaV0toPV();
      candidate.cosPA_V0 = decayV0.cosPA();
      candidate.alphaAP_V0 = decayV0.alphaAP();
      candidate.qtAP_V0 = decayV0.qtAP();
      
      candidate.globalIndex_pos = posTrack.globalIndex();
      candidate.p_pos = posTrack.p();
      candidate.pt_pos = posTrack.pt();
      candidate.eta_pos = posTrack.eta();
      candidate.phi_pos = posTrack.phi();
      candidate.pTPC_pos = posTrack.tpcInnerParam();
      candidate.pdgCode_pos = pdgCodePos;
      candidate.dcaPV_pos = decayV0.dcaDaughtersToPV("A");
      candidate.itsClSize_pos = posTrack.itsClusterSizes();
      candidate.tpcSignal_pos = posTrack.tpcSignal();
      candidate.nTPCcls_pos = posTrack.tpcNClsFound();
      candidate.nSigmaTPCe_pos = posTrack.tpcNSigmaEl();
      candidate.nSigmaTPCpi_pos = posTrack.tpcNSigmaPi();
      candidate.nSigmaTPCk_pos = posTrack.tpcNSigmaKa();
      candidate.nSigmaTPCp_pos = posTrack.tpcNSigmaPr();
      candidate.nSigmaTPCdeu_pos = posTrack.tpcNSigmaDe();
      candidate.nSigmaTPChe3_pos = posTrack.tpcNSigmaHe();
      candidate.nSigmaTPCe_pos = posTrack.tofNSigmaEl();
      candidate.nSigmaTPCpi_pos = posTrack.tofNSigmaPi();
      candidate.nSigmaTPCk_pos = posTrack.tofNSigmaKa();
      candidate.nSigmaTPCp_pos = posTrack.tofNSigmaPr();
      candidate.nSigmaTPCdeu_pos = posTrack.tofNSigmaDe();
      candidate.nSigmaTPChe3_pos = posTrack.tofNSigmaHe();
      candidate.chi2its_pos = posTrack.itsChi2NCl();
      candidate.chi2tpc_pos = posTrack.tpcChi2NCl();
      candidate.chi2matchingITSTPC_pos = posTrack.getChi2Match(); // to add

      candidate.globalIndex_neg = negTrack.globalIndex();
      candidate.p_neg = negTrack.p();
      candidate.pt_neg = negTrack.pt();
      candidate.eta_neg = negTrack.eta();
      candidate.phi_neg = negTrack.phi();
      candidate.pTPC_neg = negTrack.tpcInnerParam();
      candidate.pdgCode_neg = pdgCodeNeg;
      candidate.dcaPV_neg = decayV0.dcaDaughtersToPV("B");
      candidate.itsClSize_neg = negTrack.itsClusterSizes();
      candidate.tpcSignal_neg = negTrack.tpcSignal();
      candidate.nTPCcls_neg = negTrack.tpcNClsFound();
      candidate.nSigmaTPCe_neg = negTrack.tpcNSigmaEl();
      candidate.nSigmaTPCpi_neg = negTrack.tpcNSigmaPi();
      candidate.nSigmaTPCk_neg = negTrack.tpcNSigmaKa();
      candidate.nSigmaTPCp_neg = negTrack.tpcNSigmaPr();
      candidate.nSigmaTPCdeu_neg = negTrack.tpcNSigmaDe();
      candidate.nSigmaTPChe3_neg = negTrack.tpcNSigmaHe();
      candidate.nSigmaTPCe_neg = negTrack.tofNSigmaEl();
      candidate.nSigmaTPCpi_neg = negTrack.tofNSigmaPi();
      candidate.nSigmaTPCk_neg = negTrack.tofNSigmaKa();
      candidate.nSigmaTPCp_neg = negTrack.tofNSigmaPr();
      candidate.nSigmaTPCdeu_neg = negTrack.tofNSigmaDe();
      candidate.nSigmaTPChe3_neg = negTrack.tofNSigmaHe();
      candidate.chi2its_neg = negTrack.itsChi2NCl();
      candidate.chi2tpc_neg = negTrack.tpcChi2NCl();
      candidate.chi2matchingITSTPC_neg = negTrack.getChi2Match(); // to add
      
      m_candidateV0s.push_back(candidate);
    }
  }

  template<class C, class T>
  void fillKTree(const C& collision, const T& /*tracks*/, const aod::Cascades& cascades) 
  {
    m_candidateCascades.clear();
    
    for (auto& cascade: cascades) {

      auto v0track = cascade.template v0_as<aod::V0s>();
      auto bachelorTrack = cascade.template bachelor_as<T>(); 
      auto bachelorTrackPar = getTrackPar(bachelorTrack);

      std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};

      Decay2Bodies<T> decay(v0track, bachelorTrack, PV, m_fitter);

      if (!decay.fitSuccess()) {
        continue;
      }
      
      if (decay.dcaDaughters() > cascsetting_dcaCascDaughters) {
        continue;
      }
      if (decay.cosPA() < cascsetting_cosPA) {
        continue;
      }

      float massXi = decay.massMother(o2::constants::physics::MassLambda0, o2::constants::physics::MassPionCharged);
      float massOmega = decay.massMother(o2::constants::physics::MassLambda0, o2::constants::physics::MassKaonCharged);

      m_hCheck.fill(HIST("Xi_vs_Omega"), massOmega, massXi);
      if (std::abs(massOmega - o2::constants::physics::MassOmegaMinus) > cascsetting_massWindowOmega) {
        continue;
      }
      m_hCheck.fill(HIST("massOmegaWithBkg"), massOmega);

      if (std::abs(massXi - o2::constants::physics::MassXiMinus) < cascsetting_massWindowXi) {
        continue;
      } // enhance purity by rejecting Xi background
      m_hCheck.fill(HIST("massOmega"), massOmega);

      auto KPropagatedTrack = decay.getTrackB();

      CandidateCascade candidate;
      candidate.globalIndex_mother = cascade.globalIndex();
      candidate.p_mother = decay.p("mother");
      candidate.pt_mother = decay.pt("mother");
      candidate.eta_mother = decay.eta("mother");
      candidate.phi_mother = decay.phi("mother");
      candidate.mass_mother = massXi;
      candidate.pdgCode_mother = physics::pdg::XiMinus;
      candidate.radius_mother = decay.radiusMother();
      candidate.dcaPV_mother = decay.dcaMotherToPV();
      candidate.cosPA_mother = decay.cosPA();

      candidate.globalIndex_K = bachelorTrack.globalIndex();
      candidate.p_K = bachelorTrack.p();
      candidate.pt_K = bachelorTrack.pt();
      candidate.eta_K = bachelorTrack.eta();
      candidate.phi_K = bachelorTrack.phi();
      candidate.pdgCode_K = physics::pdg::KaonCharged;
      candidate.dcaPV_K = decay.dcaDaughtersToPV("A");
      candidate.itsClSize_K = bachelorTrack.itsClusterSizes();
      candidate.tpcSignal_K = bachelorTrack.tpcSignal();
      candidate.nTPCcls_K = bachelorTrack.tpcNClsFound();
      candidate.nSigmaTPCe_K = bachelorTrack.tpcNSigmaEl();
      candidate.nSigmaTPCpi_K = bachelorTrack.tpcNSigmaPi();
      candidate.nSigmaTPCk_K = bachelorTrack.tpcNSigmaKa();
      candidate.nSigmaTPCp_K = bachelorTrack.tpcNSigmaPr();
      candidate.nSigmaTPCdeu_K = bachelorTrack.tpcNSigmaDe();
      candidate.nSigmaTPChe3_K = bachelorTrack.tpcNSigmaHe();
      candidate.chi2its_K = bachelorTrack.itsChi2NCl();
      candidate.chi2tpc_K = bachelorTrack.tpcChi2NCl();
      candidate.chi2matchingITSTPC_K = bachelorTrack.getChi2Match(); // to add

      m_candidateCascades.push_back(candidate);

    }
  }

  void processRun3(soa::Join<aod::Collisions, aod::EvSels> const& collisions, TracksFullIU const& tracks, aod::V0s const& v0s, aod::Cascades const& cascades) 
  {
    for (auto& collision: collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8()) {
        continue;
      }

      if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
        continue;
      }

      if (std::abs(collision.posZ()) > setting_zVtxMax) {
        continue;
      }

      m_hCheck.fill(HIST("zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto v0Table_thisCollision = v0s.sliceBy(m_perCollisionV0, collIdx);
      auto cascTable_thisCollision = cascades.sliceBy(m_perCollisionCascade, collIdx);
      v0Table_thisCollision.bindExternalIndices(&tracks);
      cascTable_thisCollision.bindExternalIndices(&tracks);
      cascTable_thisCollision.bindExternalIndices(&v0s);

      if (setting_fillV0)
        fillV0Tree(collision, tracks, v0Table_thisCollision);
      if (setting_fillK)
        fillKTree(collision, tracks, cascTable_thisCollision);

      for (auto& v0: m_candidateV0s) {
        m_v0Table(
          v0.p_V0,
          //v0.pt_V0,
          v0.eta_V0,
          v0.phi_V0,
          //v0.mass_V0,
          v0.pdgCode_V0,
          v0.radius_V0,
          v0.dcaPV_V0,
          v0.cosPA_V0,
          v0.alphaAP_V0,
          v0.qtAP_V0,
          v0.p_pos,
          //v0.pt_pos,
          v0.eta_pos,
          v0.phi_pos,
          v0.pTPC_pos,
          //v0.pdgCode_pos,
          v0.dcaPV_pos,
          v0.itsClSize_pos,
          v0.tpcSignal_pos,
          v0.nTPCcls_pos,
          v0.nSigmaTPCe_pos,
          v0.nSigmaTPCpi_pos,
          v0.nSigmaTPCk_pos,
          v0.nSigmaTPCp_pos,
          v0.nSigmaTPCdeu_pos,
          v0.nSigmaTPChe3_pos,
          v0.nSigmaTOFe_pos,
          v0.nSigmaTOFpi_pos,
          v0.nSigmaTOFk_pos,
          v0.nSigmaTOFp_pos,
          v0.nSigmaTOFdeu_pos,
          v0.nSigmaTOFhe3_pos,
          v0.chi2its_pos,
          v0.chi2tpc_pos,
          v0.chi2matchingITSTPC_pos,
          v0.p_neg,
          //v0.pt_neg,
          v0.eta_neg,
          v0.phi_neg,
          v0.pTPC_neg,
          //v0.pdgCode_neg,
          v0.dcaPV_neg,
          v0.itsClSize_neg,
          v0.tpcSignal_neg,
          v0.nTPCcls_neg,
          v0.nSigmaTPCe_neg,
          v0.nSigmaTPCpi_neg,
          v0.nSigmaTPCk_neg,
          v0.nSigmaTPCp_neg,
          v0.nSigmaTPCdeu_neg,
          v0.nSigmaTPChe3_neg,
          v0.nSigmaTOFe_neg,
          v0.nSigmaTOFpi_neg,
          v0.nSigmaTOFk_neg,
          v0.nSigmaTOFp_neg,
          v0.nSigmaTOFdeu_neg,
          v0.nSigmaTOFhe3_neg,
          v0.chi2its_neg,
          v0.chi2tpc_neg,
          v0.chi2matchingITSTPC_neg
        );
      }

      for (auto& cascade: m_candidateCascades) {
        m_KTable(
          cascade.p_mother,
          //cascade.pt_mother,
          cascade.eta_mother,
          cascade.phi_mother,
          cascade.mass_mother,
          //cascade.pdgCode_mother,
          cascade.radius_mother,
          cascade.dcaPV_mother,
          cascade.cosPA_mother,
          cascade.alphaAP_mother,
          cascade.qtAP_mother,
          cascade.p_K,
          //cascade.pt_K,
          cascade.eta_K,
          cascade.phi_K,
          cascade.pTPC_K,
          cascade.pdgCode_K,
          cascade.dcaPV_K,
          cascade.itsClSize_K,
          cascade.tpcSignal_K,
          cascade.nTPCcls_K,
          cascade.nSigmaTPCe_K,
          cascade.nSigmaTPCpi_K,
          cascade.nSigmaTPCk_K,
          cascade.nSigmaTPCp_K,
          cascade.nSigmaTPCdeu_K,
          cascade.nSigmaTPChe3_K,
          cascade.nSigmaTOFe_K,
          cascade.nSigmaTOFpi_K,
          cascade.nSigmaTOFk_K,
          cascade.nSigmaTOFp_K,
          cascade.nSigmaTOFdeu_K,
          cascade.nSigmaTOFhe3_K,
          cascade.chi2its_K,
          cascade.chi2tpc_K,
          cascade.chi2matchingITSTPC_K
        );
      }
    }
  }


};