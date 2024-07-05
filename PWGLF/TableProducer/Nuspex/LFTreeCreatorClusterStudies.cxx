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
#include "TPDGCode.h"

using namespace ::o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using Track = o2::track::TrackParCov;
using TracksFullIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCEl, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe, aod::pidTPCHe, aod::pidTOFEl, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe, aod::pidTOFHe>;
using CollisionsCustom = soa::Join<aod::Collisions, aod::EvSels>;

// =========================================================================================================

namespace physics
{
  namespace pdg 
  {
    constexpr int Deuteron = 1000010020;
    constexpr int He3 = 1000020030;
  }
}

namespace BetheBloch
{

  constexpr double defaultParams[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
  static const std::vector<std::string> parNames{"p0", "p1", "p2", "p3", "p4", "resolution"};

}

// =========================================================================================================

enum V0Type : uint8_t {
  K0s = 0,
  Lambda,
  AntiLambda
};

enum CascadeType : uint8_t {
  XiMinus = 0,
  OmegaMinus
};

// =========================================================================================================

/**
 * Struct with the information of a candidate De to be stored in the tree
*/
struct CandidateDe
{
  int64_t globalIndex_De = -999;
  
  float p_De = -999.f;
  float pt_De = -999.f;
  float eta_De = -999.f;
  float phi_De = -999.f;
  float pTPC_De = -999.f;
  
  int pdgCode_De = -999;

  float dcaPV_De = -999.f;
  uint32_t itsClSize_De = 0u;
  uint16_t tpcSignal_De = 0u;
  uint16_t nTPCcls_De = 0u;

  float nSigmaTPCe_De = -999.f;
  float nSigmaTPCpi_De = -999.f;
  float nSigmaTPCk_De = -999.f;
  float nSigmaTPCp_De = -999.f;
  float nSigmaTPCdeu_De = -999.f;
  float nSigmaTPChe3_De = -999.f;

  float nSigmaTOFe_De = -999.f;
  float nSigmaTOFpi_De = -999.f;
  float nSigmaTOFk_De = -999.f;
  float nSigmaTOFp_De = -999.f;
  float nSigmaTOFdeu_De = -999.f;
  float nSigmaTOFhe3_De = -999.f;

  float chi2its_De = -999.f;
  float chi2tpc_De = -999.f;
  bool hasTPC_De = true;
};

/**
 * Struct with the information of a candidate He3 to be stored in the tree
*/
struct CandidateHe3 
{
  int64_t globalIndex_He3 = -999;
  
  float p_He3 = -999.f;
  float pt_He3 = -999.f;
  float eta_He3 = -999.f;
  float phi_He3 = -999.f;
  float pTPC_He3 = -999.f;
  
  int pdgCode_He3 = -999;

  float dcaPV_He3 = -999.f;
  uint32_t itsClSize_He3 = 0u;
  uint16_t tpcSignal_He3 = 0u;
  uint16_t nTPCcls_He3 = 0u;

  float nSigmaTPCe_He3 = -999.f;
  float nSigmaTPCpi_He3 = -999.f;
  float nSigmaTPCHe3_He3 = -999.f;
  float nSigmaTPCp_He3 = -999.f;
  float nSigmaTPCdeu_He3 = -999.f;
  float nSigmaTPChe3_He3 = -999.f;

  float nSigmaTOFe_He3 = -999.f;
  float nSigmaTOFpi_He3 = -999.f;
  float nSigmaTOFHe3_He3 = -999.f;
  float nSigmaTOFp_He3 = -999.f;
  float nSigmaTOFdeu_He3 = -999.f;
  float nSigmaTOFhe3_He3 = -999.f;

  float chi2its_He3 = -999.f;
  float chi2tpc_He3 = -999.f;
  bool hasTPC_He3 = true;
};

// =========================================================================================================

struct LFTreeCreatorClusterStudies {

  Service<o2::ccdb::BasicCCDBManager> m_ccdb;
  int m_runNumber;
  float m_d_bz;
  uint32_t m_randomSeed = 0.;

  Configurable<bool> setting_fillV0{"fillV0", true, "Fill the V0 tree"};
  Configurable<bool> setting_fillK{"fillK", true, "Fill the K tree"};
  Configurable<bool> setting_fillDe{"fillDe", true, "Fill the De tree"};
  Configurable<bool> setting_fillHe3{"fillHe3", true, "Fill the He3 tree"};

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

  Configurable<float> desetting_nsigmatpc{"desetting_nsigmaCutTPC", 4.f, "Number of sigmas for the TPC PID"};
  Configurable<bool> he3setting_compensatePIDinTracking{"he3setting_compensatePIDinTracking", true, "Compensate PID in tracking"};
  Configurable<float> he3setting_nsigmatpc{"he3setting_nsigmaCutTPC", 4.f, "Number of sigmas for the TPC PID"};

  // Bethe Bloch parameters
  std::array<float, 6> m_BBparamsDe, m_BBparamsHe;
  Configurable<LabeledArray<double>> setting_BetheBlochParamsDe{"cfgBetheBlochParamsDe", {BetheBloch::defaultParams[0], 1, 6, {"De"}, BetheBloch::parNames}, "TPC Bethe-Bloch parameterisation for Deuteron"};
  Configurable<LabeledArray<double>> setting_BetheBlochParams{"cfgBetheBlochParams", {BetheBloch::defaultParams[0], 1, 6, {"He3"}, BetheBloch::parNames}, "TPC Bethe-Bloch parameterisation for He3"};

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

  Produces<o2::aod::ClStV0Table> m_v0Table;
  Produces<o2::aod::ClStCascTable> m_KTable;
  Produces<o2::aod::ClStNucTable> m_nucleiTable;

  struct V0TrackParCov {
    int64_t globalIndex;
    Track trackParCov;
  };
  std::vector<V0TrackParCov> m_v0TrackParCovs;

  o2::vertexing::DCAFitterN<2> m_fitter;

  // =========================================================================================================

  template <typename T>
  bool initializeFitter(const T& trackParCovA, const T& trackParCovB) 
  {
    int nCand = 0;
        try {
          nCand = m_fitter.process(trackParCovA, trackParCovB);
        } catch (...) {
          LOG(error) << "Exception caught in DCA fitter process call!";
          return false;
        }
        if (nCand == 0) {
          return false;
        }

      return true;
  }

  void computeTrackMomentum(const int itrack, std::array<float, 3>& mom)
  {
    auto fittedTrack = m_fitter.getTrack(itrack);
    fittedTrack.getPxPyPzGlo(mom);
  }

  void computeMotherMomentum(const std::array<float, 3>& momA, const std::array<float, 3>& momB, std::array<float, 3>& momMother)
  {
    momMother[0] = momA[0] + momB[0];
    momMother[1] = momA[1] + momB[1];
    momMother[2] = momA[2] + momB[2];
  }
  
  /**
   * Compute the alpha for the Armenteros-Podolanski plot
  */
  float computeAlphaAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momA, const std::array<float, 3>& momB)
  {
    float momTot = std::hypot(momMother[0], momMother[1], momMother[2]);
    float lQlA = (momA[0] * momMother[0] + momA[1] * momMother[1] + momA[2] * momMother[2]) / momTot;
    float lQlNeg = (momB[0] * momMother[0] + momB[1] * momMother[1] + momB[2] * momMother[2]) / momTot;
    return (lQlA - lQlNeg) / (lQlA + lQlNeg);
  }

  /**
   * Compute the qt for the Armenteros-Podolanski plot
  */
  float computeQtAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momA)
  {
    float dp = momMother[0] * momA[0] + momMother[1] * momA[1] + momMother[2] * momA[2];
    float p2V0 = std::hypot(momMother[0], momMother[1], momMother[2]);
    float p2A = std::hypot(momA[0], momA[1], momA[2]);
    return std::sqrt(p2A - dp * dp / p2V0);
  }

  float dcaMotherToPV(const std::array<float, 3>& decayVtx, const std::array<float, 3>& PV, std::array<float, 3> momMother) const
  {
    std::array<float, 3> relPos = {decayVtx[0] - PV[0], decayVtx[1] - PV[1], decayVtx[2] - PV[2]};
    float lmomMotherl = std::hypot(momMother[0], momMother[1], momMother[2]);
    return std::sqrt((std::pow(relPos[1] * momMother[2] - relPos[2] * momMother[1], 2) + std::pow(relPos[2] * momMother[0] - relPos[0] * momMother[2], 2) + std::pow(relPos[0] * momMother[1] - relPos[1] * momMother[0], 2))) / lmomMotherl;
  }
  
  template <typename T>
  float dcaDaughterToPV(const std::array<float, 3>& PV, T& trackParCov, gpu::gpustd::array<float, 2>& dcaInfo)
  {
      o2::base::Propagator::Instance()->propagateToDCABxByBz({PV[0], PV[1], PV[2]}, trackParCov, 2.f, m_fitter.getMatCorrType(), &dcaInfo);
      return std::hypot(dcaInfo[0], dcaInfo[1]);
  }
  
  float computeMassMother(const float massA, const float massB, const std::array<float, 3>& momA, const std::array<float, 3>& momB, const std::array<float, 3>& momMother) const
  {
    float eA = std::hypot(massA, std::hypot(momA[0], momA[1], momA[2]));
    float eB = std::hypot(massB, std::hypot(momB[0], momB[1], momB[2]));
    float lmomMotherl = std::hypot(momMother[0], momMother[1], momMother[2]);
    float eMother = eA + eB;
    return std::sqrt(eMother * eMother - lmomMotherl * lmomMotherl);
  }

  // =========================================================================================================

  bool collisionSelection(const CollisionsCustom::iterator& collision)
  {
    if (!collision.sel8()) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (std::abs(collision.posZ()) > setting_zVtxMax) {
      return false;
    }
    return true;
  }

  // =========================================================================================================

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

  bool qualitySelectionV0(const double dcaV0toPV, const double dcaV0daughters, const double radiusV0, const double cosPA)
  {
    if (dcaV0daughters > v0setting_dcaV0daughters) {
        return false;
      }
      if (radiusV0 > v0setting_radiusMax || radiusV0 < v0setting_radiusMin) {
        return false;
      }
      if (dcaV0toPV > v0setting_dcaV0toPV) {
        return false;
      }
      if (cosPA < v0setting_cosPA) {
        return false;
      }
      return true;
  }

  bool qualitySelectionCascade(const double dcaCascDaughters, const double cosPA)
  {
    if (dcaCascDaughters > cascsetting_dcaCascDaughters) {
      return false;
    }
    if (cosPA < cascsetting_cosPA) {
      return false;
    }
    return true;
  }

  // =========================================================================================================

  template <typename T>
  float computeNSigmaDe(const T& candidate)
  {
    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(candidate.tpcInnerParam() * 2.f / constants::physics::MassDeuteron), m_BBparamsDe[0], m_BBparamsDe[1], m_BBparamsDe[2], m_BBparamsDe[3], m_BBparamsDe[4]);
    double resoTPC{expTPCSignal * m_BBparamsDe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename T>
  bool selectionPIDDe(const T& candidate)
  {
    auto nSigmaDe = computeNSigmaDe(candidate);
    if (std::abs(nSigmaDe) < desetting_nsigmatpc) {
      return true;
    }
    return false;
  }

  // =========================================================================================================

  template <typename T>
  float computeNSigmaHe3(const T& candidate)
  {
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParam = (heliumPID && he3setting_compensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(correctedTPCinnerParam * 2.f / constants::physics::MassHelium3), m_BBparamsHe[0], m_BBparamsHe[1], m_BBparamsHe[2], m_BBparamsHe[3], m_BBparamsHe[4]);
    double resoTPC{expTPCSignal * m_BBparamsHe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename T>
  bool selectionPIDHe3(const T& candidate)
  {
    auto nSigmaHe3 = computeNSigmaHe3(candidate);
    if (std::abs(nSigmaHe3) < he3setting_nsigmatpc) {
      return true;
    }
    return false;
  }

  // =========================================================================================================

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
  
  void fillV0Tree(const std::array<float, 3>& PV, const aod::V0s& v0s) 
  {
    m_v0TrackParCovs.clear();

    for (const auto& v0: v0s) {

      auto posTrack = v0.posTrack_as<TracksFullIU>();
      auto negTrack = v0.negTrack_as<TracksFullIU>();

      if (!qualitySelectionV0Daughter(posTrack) || !qualitySelectionV0Daughter(negTrack)) {
        continue;
      }

      auto daughterTrackCovarianceA = getTrackParCov(posTrack);
      auto daughterTrackCovarianceB = getTrackParCov(negTrack);

      if (!initializeFitter(daughterTrackCovarianceA, daughterTrackCovarianceB)) {
        continue;
      }

      std::array<float, 3> momPos, momNeg, momMother;
      computeTrackMomentum(0, momPos);
      computeTrackMomentum(1, momNeg);
      computeMotherMomentum(momPos, momNeg, momMother);
      ROOT::Math::SVector<double, 3> vec_decayVtx = m_fitter.getPCACandidate();
      std::array<float, 3> decayVtx = {static_cast<float>(vec_decayVtx[0]), static_cast<float>(vec_decayVtx[1]), static_cast<float>(vec_decayVtx[2])};

      float alphaAP = computeAlphaAP(momMother, momPos, momNeg);
      float qtAP = computeQtAP(momMother, momPos);

      m_hCheck.fill(HIST("armenteros_plot_before_selections"), alphaAP, qtAP);

      // V0 selections
      //float ptV0 = decayV0.pt();
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
      float dcaV0daughters = std::sqrt(m_fitter.getChi2AtPCACandidate());
      float radiusV0 = std::hypot(decayVtx[0], decayVtx[1]);
      float dcaV0toPV = std::hypot(decayVtx[0] - PV[0], decayVtx[1] - PV[1]);
      float cosPA = RecoDecay::cpa(PV, decayVtx, momMother);

      if (!qualitySelectionV0(dcaV0toPV, dcaV0daughters, radiusV0, cosPA)) {
        continue;
      }
      
      // mass hypothesis
      float massLambdaV0 = computeMassMother(o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, momPos, momNeg, momMother);
      float massK0sV0 = computeMassMother(o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged, momPos, momNeg, momMother);

      uint8_t selectionMap{BIT(Lambda) | BIT(AntiLambda) | BIT(K0s)};
      if (std::abs(massK0sV0 - o2::constants::physics::MassK0Short) > v0setting_massWindowK0s) {
        CLRBIT(selectionMap, K0s);
      }
      if ((std::abs(massLambdaV0 - o2::constants::physics::MassLambda0) > v0setting_massWindowLambda) && (alphaAP > 0)) { 
        CLRBIT(selectionMap, Lambda);
      }
      if ((std::abs(massLambdaV0 - o2::constants::physics::MassLambda0) > v0setting_massWindowLambda) && (alphaAP < 0)) {
        CLRBIT(selectionMap, AntiLambda);
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

      if (selectionMap == 0 || (selectionMap & (selectionMap - 1)) != 0) {
        continue;
      }

      int pdgCodeV0 = -999;
      if (TESTBIT(selectionMap, Lambda)) {
        pdgCodeV0 = kLambda0;
      } else if (TESTBIT(selectionMap, AntiLambda)) {
        pdgCodeV0 = -kLambda0;
      } else if (TESTBIT(selectionMap, K0s)) {
        pdgCodeV0 = kK0Short;
      } else {
        continue;
      }

      gpu::gpustd::array<float, 2> dcaInfo;
      float dcaToPVpos = dcaDaughterToPV(PV, daughterTrackCovarianceA, dcaInfo);
      if (dcaToPVpos < v0setting_dcaDaughtersToPV /*&& std::abs(dcaInfo[0]) < v0setting_dcaDaughtersToPV*/) {
        continue;
      }      
      float dcaToPVneg = dcaDaughterToPV(PV, daughterTrackCovarianceB, dcaInfo);
      if (dcaToPVneg < v0setting_dcaDaughtersToPV /*&& std::abs(dcaInfo[0]) < v0setting_dcaDaughtersToPV*/) {
        continue;
      }

      // fill the candidate and add it to the Armenteros plot
      m_hCheck.fill(HIST("armenteros_plot"), alphaAP, qtAP);

      V0TrackParCov v0TrackParCov{v0.globalIndex(), m_fitter.createParentTrackParCov()};
      m_v0TrackParCovs.push_back(v0TrackParCov);

      m_v0Table(
          std::hypot(momMother[0], momMother[1], momMother[2]),           // p_V0
          //v0.pt_V0,
          RecoDecay::eta(momMother),                                      // eta_V0
          RecoDecay::phi(momMother),                                      // phi_V0
          //v0.mass_V0,                             
          pdgCodeV0,                                                      // pdgCode_V0
          radiusV0,                                                       // radius_V0  
          dcaV0toPV,                                                      // dcaPV_V0
          cosPA,                                                          // cosPA_V0
          alphaAP,                                                        // alphaAP_V0
          qtAP,                                                           // qtAP_V0
          std::hypot(momPos[0], momPos[1], momPos[2]) * posTrack.sign(),  // p_pos
          //v0.pt_pos,
          RecoDecay::eta(momPos),                                         // eta_pos  
          RecoDecay::phi(momPos),                                         // phi_pos
          posTrack.tpcInnerParam() * posTrack.sign(),                     // pTPC_pos
          //v0.pdgCode_pos,
          dcaToPVpos,                                                     // dcaPV_pos
          posTrack.itsClusterSizes(),                                     // itsClSize_pos
          posTrack.tpcSignal(),                                           // tpcSignal_pos
          posTrack.tpcNClsFound(),                                        // nTPCcls_pos
          posTrack.tpcNSigmaEl(),                                         // nSigmaTPCe_pos
          posTrack.tpcNSigmaPi(),                                         // nSigmaTPCpi_pos
          posTrack.tpcNSigmaKa(),                                         // nSigmaTPCk_pos
          posTrack.tpcNSigmaPr(),                                         // nSigmaTPCp_pos
          posTrack.tpcNSigmaDe(),                                         // nSigmaTPCdeu_pos
          posTrack.tpcNSigmaHe(),                                         // nSigmaTPChe3_pos
          posTrack.tofNSigmaEl(),                                         // nSigmaTOFe_pos
          posTrack.tofNSigmaPi(),                                         // nSigmaTOFpi_pos
          posTrack.tofNSigmaKa(),                                         // nSigmaTOFk_pos
          posTrack.tofNSigmaPr(),                                         // nSigmaTOFp_pos
          posTrack.tofNSigmaDe(),                                         // nSigmaTOFdeu_pos
          posTrack.tofNSigmaHe(),                                         // nSigmaTOFhe3_pos
          posTrack.itsChi2NCl(),                                          // chi2its_pos
          posTrack.tpcChi2NCl(),                                          // chi2tpc_pos
          posTrack.hasTPC(),                                              // hasTPC_pos
          std::hypot(momNeg[0], momNeg[1], momNeg[2]) * negTrack.sign(),  // p_neg
          //v0.pt_neg,
          RecoDecay::eta(momNeg),                                         // eta_neg
          RecoDecay::phi(momNeg),                                         // phi_neg
          negTrack.tpcInnerParam() * negTrack.sign(),                     // pTPC_neg
          //v0.pdgCode_neg,
          dcaToPVneg,                                                     // dcaPV_neg
          negTrack.itsClusterSizes(),                                     // itsClSize_neg
          negTrack.tpcSignal(),                                           // tpcSignal_neg
          negTrack.tpcNClsFound(),                                        // nTPCcls_neg
          negTrack.tpcNSigmaEl(),                                         // nSigmaTPCe_neg
          negTrack.tpcNSigmaPi(),                                         // nSigmaTPCpi_neg
          negTrack.tpcNSigmaKa(),                                         // nSigmaTPCk_neg 
          negTrack.tpcNSigmaPr(),                                         // nSigmaTPCp_neg
          negTrack.tpcNSigmaDe(),                                         // nSigmaTPCdeu_neg
          negTrack.tpcNSigmaHe(),                                         // nSigmaTPChe3_neg
          negTrack.tofNSigmaEl(),                                         // nSigmaTOFe_neg
          negTrack.tofNSigmaPi(),                                         // nSigmaTOFpi_neg
          negTrack.tofNSigmaKa(),                                         // nSigmaTOFk_neg
          negTrack.tofNSigmaPr(),                                         // nSigmaTOFp_neg
          negTrack.tofNSigmaDe(),                                         // nSigmaTOFdeu_neg
          negTrack.tofNSigmaHe(),                                         // nSigmaTOFhe3_neg
          negTrack.itsChi2NCl(),                                          // chi2its_neg
          negTrack.tpcChi2NCl(),                                          // chi2tpc_neg
          negTrack.hasTPC()                                               // hasTPC_neg
          );
    }
  }

  void fillKTree(const std::array<float, 3>& PV, const aod::Cascades& cascades) 
  {
    
    for (auto& cascade: cascades) {

      auto v0Track = cascade.template v0_as<aod::V0s>();
      auto bachelorTrack = cascade.template bachelor_as<TracksFullIU>(); 
      auto bachelorTrackPar = getTrackPar(bachelorTrack);

      auto itv0 = std::find_if(m_v0TrackParCovs.begin(), m_v0TrackParCovs.end(), [&](const V0TrackParCov& v0) { return v0.globalIndex == v0Track.globalIndex(); });
      if (itv0 == m_v0TrackParCovs.end()) {
        continue;
      }
      
      auto v0TrackCovariance = itv0->trackParCov;
      auto bachelorTrackCovariance = getTrackParCov(bachelorTrack);

      if (!initializeFitter(v0TrackCovariance, bachelorTrackCovariance)) {
        continue;
      }

      std::array<float, 3> momV0, momBachelor, momMother;
      computeTrackMomentum(0, momV0);
      computeTrackMomentum(1, momBachelor);
      computeMotherMomentum(momV0, momBachelor, momMother);
      ROOT::Math::SVector<double, 3> vec_decayVtx = m_fitter.getPCACandidate();
      std::array<float, 3> decayVtx = {static_cast<float>(vec_decayVtx[0]), static_cast<float>(vec_decayVtx[1]), static_cast<float>(vec_decayVtx[2])};

      float dcaV0daughters = std::sqrt(m_fitter.getChi2AtPCACandidate());
      float radius = std::hypot(decayVtx[0], decayVtx[1]);
      float cosPA = RecoDecay::cpa(PV, decayVtx, momMother);
      float dcaCascToPV = dcaMotherToPV(decayVtx, PV, momMother);

      if (!qualitySelectionCascade(dcaV0daughters, cosPA)) {
        continue;
      }

      gpu::gpustd::array<float, 2> dcaInfo;
      float dcaToPVbachelor = dcaDaughterToPV(PV, bachelorTrackCovariance, dcaInfo);
      
      float massXi = computeMassMother(o2::constants::physics::MassLambda0, o2::constants::physics::MassPionCharged, momV0, momBachelor, momMother);
      float massOmega = computeMassMother(o2::constants::physics::MassLambda0, o2::constants::physics::MassKaonCharged, momV0, momBachelor, momMother);

      m_hCheck.fill(HIST("Xi_vs_Omega"), massOmega, massXi);
      if (std::abs(massOmega - o2::constants::physics::MassOmegaMinus) > cascsetting_massWindowOmega) {
        continue;
      }
      m_hCheck.fill(HIST("massOmegaWithBkg"), massOmega);

      if (std::abs(massXi - o2::constants::physics::MassXiMinus) < cascsetting_massWindowXi) {
        continue;
      } // enhance purity by rejecting Xi background
      m_hCheck.fill(HIST("massOmega"), massOmega);

      m_KTable(
        std::hypot(momMother[0], momMother[1], momMother[2]),                               // p_mother
        //cascade.pt_mother,                    
        RecoDecay::eta(momMother),                                                          // eta_mother
        RecoDecay::phi(momMother),                                                          // phi_mother
        massXi,                                                                             // mass_mother
        //kXiMinus * int(bachelorTrack.sign()),                                               // pdgCode_mother
        radius,                                                                             // radius_mother
        dcaCascToPV,                                                                        // dcaPV_mother
        cosPA,                                                                              // cosPA_mother
        std::hypot(momBachelor[0], momBachelor[1], momBachelor[2]) * bachelorTrack.sign(),  // p_K
        //cascade.pt_K,                   
        RecoDecay::eta(momBachelor),                                                        // eta_K
        RecoDecay::phi(momBachelor),                                                        // phi_K
        bachelorTrack.tpcInnerParam() * bachelorTrack.sign(),                               // pTPC_K
        kKPlus * int(bachelorTrack.sign()),                                                 // pdgCode_K
        dcaToPVbachelor,                                                                    // dcaPV_K
        bachelorTrack.itsClusterSizes(),                                                    // itsClSize_K
        bachelorTrack.tpcSignal(),                                                          // tpcSignal_K
        bachelorTrack.tpcNClsFound(),                                                       // nTPCcls_K
        bachelorTrack.tpcNSigmaEl(),                                                        // nSigmaTPCe_K
        bachelorTrack.tpcNSigmaPi(),                                                        // nSigmaTPCpi_K
        bachelorTrack.tpcNSigmaKa(),                                                        // nSigmaTPCk_K
        bachelorTrack.tpcNSigmaPr(),                                                        // nSigmaTPCp_K
        bachelorTrack.tpcNSigmaDe(),                                                        // nSigmaTPCdeu_K
        bachelorTrack.tpcNSigmaHe(),                                                        // nSigmaTPChe3_K
        bachelorTrack.tofNSigmaEl(),                                                        // nSigmaTOFe_K
        bachelorTrack.tofNSigmaPi(),                                                        // nSigmaTOFpi_K
        bachelorTrack.tofNSigmaKa(),                                                        // nSigmaTOFk_K
        bachelorTrack.tofNSigmaPr(),                                                        // nSigmaTOFp_K
        bachelorTrack.tofNSigmaDe(),                                                        // nSigmaTOFdeu_K
        bachelorTrack.tofNSigmaHe(),                                                        // nSigmaTOFhe3_K
        bachelorTrack.itsChi2NCl(),                                                         // chi2its_K
        bachelorTrack.tpcChi2NCl(),                                                         // chi2tpc_K
        bachelorTrack.hasTPC()                                                              // hasTPC_K
      );                    

    }
  }

  void fillDeTree(const TracksFullIU& tracks)
  {
    for (auto& track: tracks) {
      if (track.pidForTracking() != o2::track::PID::Deuteron) {
        continue;
      }
      if (!selectionPIDDe(track)) {
        continue;
      }

      m_nucleiTable(
        track.p() * track.sign(),               // p_De,
        //track.pt() * track.sign(),            // pt_De,
        track.eta(),                            // eta_De,
        track.phi(),                            // phi_De,
        track.tpcInnerParam() * track.sign(),   // pTPC_De,
        physics::pdg::Deuteron,                 // pdgCode_De,
        0.,                                     // dcaPV_De,
        track.itsClusterSizes(),                // itsClSize_De,
        track.tpcSignal(),                      // tpcSignal_De,
        track.tpcNClsFound(),                   // nTPCcls_De,
        track.tpcNSigmaEl(),                    // nSigmaTPCe_De,
        track.tpcNSigmaPi(),                    // nSigmaTPCpi_De,
        track.tpcNSigmaKa(),                    // nSigmaTPCk_De,
        track.tpcNSigmaPr(),                    // nSigmaTPCp_De,
        track.tpcNSigmaDe(),                    // nSigmaTPCdeu_De,
        track.tpcNSigmaHe(),                    // nSigmaTPChe3_De,
        track.tofNSigmaEl(),                    // nSigmaTOFe_De,
        track.tofNSigmaPi(),                    // nSigmaTOFpi_De,
        track.tofNSigmaKa(),                    // nSigmaTOFk_De,
        track.tofNSigmaPr(),                    // nSigmaTOFp_De,
        track.tofNSigmaDe(),                    // nSigmaTOFdeu_De,
        track.tofNSigmaHe(),                    // nSigmaTOFhe3_De,
        track.itsChi2NCl(),                     // chi2its_De,
        track.tpcChi2NCl(),                     // chi2tpc_De,
        track.hasTPC()                          // hasTPC_De
      );
    }
  }

  void fillHe3Tree(const TracksFullIU& tracks)
  {
    for (auto& track: tracks) {
      if (track.pidForTracking() != o2::track::PID::Helium3) {
        continue;
      }
      if (!selectionPIDHe3(track)) {
        continue;
      }
      
      m_nucleiTable(
        track.p() * track.sign(),               // p_He3,
        //track.pt() * track.sign(),            // pt_He3,
        track.eta(),                            // eta_He3,
        track.phi(),                            // phi_He3,
        track.tpcInnerParam() * track.sign(),   // pTPC_He3,
        physics::pdg::Deuteron,                 // pdgCode_He3,
        0.,                                     // dcaPV_He3,
        track.itsClusterSizes(),                // itsClSize_He3,
        track.tpcSignal(),                      // tpcSignal_He3,
        track.tpcNClsFound(),                   // nTPCcls_He3,
        track.tpcNSigmaEl(),                    // nSigmaTPCe_He3,
        track.tpcNSigmaPi(),                    // nSigmaTPCpi_He3,
        track.tpcNSigmaKa(),                    // nSigmaTPCk_He3,
        track.tpcNSigmaPr(),                    // nSigmaTPCp_He3,
        track.tpcNSigmaDe(),                    // nSigmaTPCdeu_He3,
        track.tpcNSigmaHe(),                    // nSigmaTPChe3_He3,
        track.tofNSigmaEl(),                    // nSigmaTOFe_He3,
        track.tofNSigmaPi(),                    // nSigmaTOFpi_He3,
        track.tofNSigmaKa(),                    // nSigmaTOFk_He3,
        track.tofNSigmaPr(),                    // nSigmaTOFp_He3,
        track.tofNSigmaDe(),                    // nSigmaTOFdeu_He3,
        track.tofNSigmaHe(),                    // nSigmaTOFhe3_He3,
        track.itsChi2NCl(),                     // chi2its_He3,
        track.tpcChi2NCl(),                     // chi2tpc_He3,
        track.hasTPC()                          // hasTPC_He3,
      );
    }

  }

  void processRun3(CollisionsCustom::iterator const& collision, TracksFullIU const& tracks, aod::V0s const& v0s, aod::Cascades const& cascades) 
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!collisionSelection(collision)) {
      return;
    }

    m_hCheck.fill(HIST("zVtx"), collision.posZ());
    std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};

    const uint64_t collIdx = collision.globalIndex();
    auto v0Table_thisCollision = v0s.sliceBy(m_perCollisionV0, collIdx);
    auto cascTable_thisCollision = cascades.sliceBy(m_perCollisionCascade, collIdx);
    v0Table_thisCollision.bindExternalIndices(&tracks);
    cascTable_thisCollision.bindExternalIndices(&tracks);
    cascTable_thisCollision.bindExternalIndices(&v0s);

    if (setting_fillV0)
      fillV0Tree(PV, v0Table_thisCollision);
    if (setting_fillK && setting_fillV0) // the v0 loops are needed for the Ks
      fillKTree(PV, cascTable_thisCollision);
    if (setting_fillDe)
      fillDeTree(tracks);
    if (setting_fillHe3)
      fillHe3Tree(tracks);

  }
  PROCESS_SWITCH(LFTreeCreatorClusterStudies, processRun3, "process Run 3", false);

};