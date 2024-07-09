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

namespace physics
{
namespace pdg
{
constexpr int Deuteron = 1000010020;
constexpr int He3 = 1000020030;
} // namespace pdg
} // namespace physics

namespace BetheBloch
{

constexpr double defaultParams[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> parNames{"p0", "p1", "p2", "p3", "p4", "resolution"};

} // namespace BetheBloch

enum V0Type : uint8_t {
  K0s = 0,
  Lambda,
  AntiLambda,
  Photon
};

enum CascadeType : uint8_t {
  XiMinus = 0,
  OmegaMinus
};

enum Selections {
  kNoCut = 0,
  kSel8,
  kVtxZ,
  kAll
};

enum V0Selections {
  kV0NoCut = 0,
  kV0DaughterQuality,
  kV0DCA,
  kV0Radius,
  kV0dcaPV,
  kV0CosPA,
  kV0PID,
  kV0All
};

enum CascSelections {
  kCascDCA = 0,
  kCascCosPA,
  kAcceptedOmega,
  kRejectedXi,
  kCascAll
};

enum DeSelections {
  kDeNoCut = 0,
  kDePIDforTrk,
  kDePID,
  kDeAll
};

enum He3Selections {
  kHe3NoCut = 0,
  kHe3PIDforTrk,
  kHe3PID,
  kHe3All
};

struct CandidateV0 {
  float p_V0 = -999.f;
  float eta_V0 = -999.f;
  float phi_V0 = -999.f;
  int pdgCode_V0 = -999.f;
  float radius_V0 = -999.f;
  float dcaPV_V0 = -999.f;
  float cosPA_V0 = -999.f;
  float alphaAP_V0 = -999.f;
  float qtAP_V0 = -999.f;
  int mc_pdgCode_V0 = -999.f;

  float p_pos = -999.f;
  float eta_pos = -999.f;
  float phi_pos = -999.f;
  float pTPC_pos = -999.f;
  float dcaPV_pos = -999.f;
  uint32_t itsClsize_pos = 0;
  uint16_t tpcSignal_pos = 0;
  uint8_t tpcNcls_pos = 0;
  float nSigmaTPCel_pos = -999.f;
  float nSigmaTPCpi_pos = -999.f;
  float nSigmaTPCka_pos = -999.f;
  float nSigmaTPCpr_pos = -999.f;
  float nSigmaTPCde_pos = -999.f;
  float nSigmaTPChe_pos = -999.f;
  float nSigmaTOFel_pos = -999.f;
  float nSigmaTOFpi_pos = -999.f;
  float nSigmaTOFka_pos = -999.f;
  float nSigmaTOFpr_pos = -999.f;
  float nSigmaTOFde_pos = -999.f;
  float nSigmaTOFhe_pos = -999.f;
  float chi2its_pos = -999.f;
  float chi2tpc_pos = -999.f;
  bool hasTPC_pos = false;
  int mc_pdgCode_pos = -999.f;

  float p_neg = -999.f;
  float eta_neg = -999.f;
  float phi_neg = -999.f;
  float pTPC_neg = -999.f;
  float dcaPV_neg = -999.f;
  uint32_t itsClsize_neg = 0;
  uint16_t tpcSignal_neg = 0;
  uint8_t tpcNcls_neg = 0;
  float nSigmaTPCel_neg = -999.f;
  float nSigmaTPCpi_neg = -999.f;
  float nSigmaTPCka_neg = -999.f;
  float nSigmaTPCpr_neg = -999.f;
  float nSigmaTPCde_neg = -999.f;
  float nSigmaTPChe_neg = -999.f;
  float nSigmaTOFel_neg = -999.f;
  float nSigmaTOFpi_neg = -999.f;
  float nSigmaTOFka_neg = -999.f;
  float nSigmaTOFpr_neg = -999.f;
  float nSigmaTOFde_neg = -999.f;
  float nSigmaTOFhe_neg = -999.f;
  float chi2its_neg = -999.f;
  float chi2tpc_neg = -999.f;
  bool hasTPC_neg = false;
  int mc_pdgCode_neg = -999.f;
};

struct CandidateK {
  float p_casc = -999.f;
  float eta_casc = -999.f;
  float phi_casc = -999.f;
  int pdgCode_casc = -999.f;
  float radius_casc = -999.f;
  float dcaPV_casc = -999.f;
  float cosPA_casc = -999.f;
  int mc_pdgCode_casc = -999.f;

  float p_K = -999.f;
  float eta_K = -999.f;
  float phi_K = -999.f;
  float pTPC_K = -999.f;
  float dcaPV_K = -999.f;
  uint32_t itsClsize_K = 0;
  uint16_t tpcSignal_K = 0;
  uint8_t tpcNcls_K = 0;
  float nSigmaTPCel_K = -999.f;
  float nSigmaTPCpi_K = -999.f;
  float nSigmaTPCka_K = -999.f;
  float nSigmaTPCpr_K = -999.f;
  float nSigmaTPCde_K = -999.f;
  float nSigmaTPChe_K = -999.f;
  float nSigmaTOFel_K = -999.f;
  float nSigmaTOFpi_K = -999.f;
  float nSigmaTOFka_K = -999.f;
  float nSigmaTOFpr_K = -999.f;
  float nSigmaTOFde_K = -999.f;
  float nSigmaTOFhe_K = -999.f;
  float chi2its_K = -999.f;
  float chi2tpc_K = -999.f;
  bool hasTPC_K = false;
  int mc_pdgCode_K = -999.f;
};

struct candidateDe {
  float p_de = -999.f;
  float eta_de = -999.f;
  float phi_de = -999.f;
  float pTPC_de = -999.f;
  float dcaToPV_de = -999.f;
  uint32_t itsClsize_de = 0;
  uint16_t tpcSignal_de = 0;
  uint8_t tpcNcls_de = 0;
  float nSigmaTPCel_de = -999.f;
  float nSigmaTPCpi_de = -999.f;
  float nSigmaTPCka_de = -999.f;
  float nSigmaTPCpr_de = -999.f;
  float nSigmaTPCde_de = -999.f;
  float nSigmaTPChe_de = -999.f;
  float nSigmaTOFel_de = -999.f;
  float nSigmaTOFpi_de = -999.f;
  float nSigmaTOFka_de = -999.f;
  float nSigmaTOFpr_de = -999.f;
  float nSigmaTOFde_de = -999.f;
  float nSigmaTOFhe_de = -999.f;
  float chi2its_de = -999.f;
  float chi2tpc_de = -999.f;
  bool hasTPC_de = false;
  int mc_pdgCode_de = -999.f;
};

struct candidateHe {
  float p_he = -999.f;
  float eta_he = -999.f;
  float phi_he = -999.f;
  float pTPC_he = -999.f;
  float dcaToPV_he = -999.f;
  uint32_t itsClsize_he = 0;
  uint16_t tpcSignal_he = 0;
  uint8_t tpcNcls_he = 0;
  float nSigmaTPCel_he = -999.f;
  float nSigmaTPCpi_he = -999.f;
  float nSigmaTPCka_he = -999.f;
  float nSigmaTPCpr_he = -999.f;
  float nSigmaTPCde_he = -999.f;
  float nSigmaTPChe_he = -999.f;
  float nSigmaTOFel_he = -999.f;
  float nSigmaTOFpi_he = -999.f;
  float nSigmaTOFka_he = -999.f;
  float nSigmaTOFpr_he = -999.f;
  float nSigmaTOFde_he = -999.f;
  float nSigmaTOFhe_he = -999.f;
  float chi2its_he = -999.f;
  float chi2tpc_he = -999.f;
  bool hasTPC_he = false;
  int mc_pdgCode_he = -999.f;
};

struct LfTreeCreatorClusterStudies {

  Service<o2::ccdb::BasicCCDBManager> m_ccdb;
  int m_runNumber;
  int m_collisionCounter = 0;
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

  // Configurable<float> v0setting_etaMaxV0{"etaMaxV0", 0.8f, "Maximum eta for the V0 daughters"};
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

  Configurable<float> cascsetting_dcaCascDaughters{"casc_setting_dcaV0daughters", 0.1f, "DCA between the V0 daughters"};
  Configurable<float> cascsetting_cosPA{"casc_setting_cosPA", 0.99f, "Cosine of the pointing angle of the V0"};
  Configurable<float> cascsetting_massWindowOmega{"casc_setting_massWindowOmega", 0.01f, "Mass window for the Omega"};
  Configurable<float> cascsetting_massWindowXi{"casc_setting_massWindowXi", 0.01f, "Mass window for the Xi"};

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
    {{"collision_selections", "Collision selection; selection; counts", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}},
     {"v0_selections", "V0 selection; selection; counts", {HistType::kTH1F, {{V0Selections::kV0All, -0.5, static_cast<double>(V0Selections::kV0All) - 0.5}}}},
     {"casc_selections", "Cascade selection; selection; counts", {HistType::kTH1F, {{CascSelections::kCascAll, -0.5, static_cast<double>(CascSelections::kCascAll) - 0.5}}}},
     {"de_selections", "Deuteron track selection; selection; counts", {HistType::kTH1F, {{DeSelections::kDeAll, -0.5, static_cast<double>(DeSelections::kDeAll) - 0.5}}}},
     {"he3_selections", "He3 track selection; selection; counts", {HistType::kTH1F, {{He3Selections::kHe3All, -0.5, static_cast<double>(He3Selections::kHe3All) - 0.5}}}},
     {"invmass_Lambda", "#Lambda invariant mass; m (GeV/#it{c}^{2}); counts", {HistType::kTH1F, {{200, 0.7f, 1.5f}}}},
     {"invmass_K0s", "K0s invariant mass; m (GeV/#it{c}^{2}); counts", {HistType::kTH1F, {{200, 0.4f, 0.6f}}}},
     {"armenteros_plot_before_selections", "Armenteros-Podolanski plot; #alpha; q_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 0.3f}}}},
     {"armenteros_plot", "Armenteros-Podolanski plot; #alpha; q_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 3.f}}}},
     {"Xi_vs_Omega", "Mass Xi vs Omega; mass Omega (GeV/#it{c}^{2}); mass Xi (GeV/#it{c}^{2})", {HistType::kTH2F, {{100, 1.f, 2.f}, {100, 1.f, 2.f}}}},
     {"massOmega", "Mass Omega; mass Omega (GeV/#it{c}^{2}); #it{c}ounts", {HistType::kTH1F, {{100, 1.f, 2.f}}}},
     {"massOmegaWithBkg", "Mass Omega with Background; mass Omega (GeV/#it{c}^{2}); #it{c}ounts", {HistType::kTH1F, {{100, 1.f, 2.f}}}},
     {"massLambda", "Mass Lambda; mass Lambda (GeV/#it{c}^{2}); #it{c}ounts", {HistType::kTH1F, {{100, 1.f, 2.f}}}},
     {"BetheBlochDe", "Bethe-Bloch Deuteron; p (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
     {"BetheBlochHe3", "Bethe-Bloch He3; p (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
     {"zVtx", "Binning for the vertex z in cm", {HistType::kTH1F, {{100, -20.f, 20.f}}}}},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true}; // check histograms

  Produces<o2::aod::ClStV0Table> m_v0Table;
  Produces<o2::aod::ClStCascTable> m_KTable;
  Produces<o2::aod::ClStNucTable> m_nucleiTable;

  struct V0TrackParCov {
    int64_t globalIndex;
    Track trackParCov;
  };
  std::vector<V0TrackParCov> m_v0TrackParCovs;

  o2::vertexing::DCAFitterN<2> m_fitter;

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

  bool collisionSelection(const CollisionsCustom::iterator& collision)
  {
    m_hCheck.fill(HIST("collision_selections"), Selections::kNoCut);
    // if (!collision.sel8()) {
    //   return false;
    // }
    m_hCheck.fill(HIST("collision_selections"), Selections::kSel8);
    // if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
    //   return false;
    // }
    if (std::abs(collision.posZ()) > setting_zVtxMax) {
      return false;
    }
    m_hCheck.fill(HIST("collision_selections"), Selections::kVtxZ);
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
    m_hCheck.fill(HIST("v0_selections"), V0Selections::kV0DCA);
    if (radiusV0 > v0setting_radiusMax || radiusV0 < v0setting_radiusMin) {
      return false;
    }
    m_hCheck.fill(HIST("v0_selections"), V0Selections::kV0Radius);
    if (dcaV0toPV > v0setting_dcaV0toPV) {
      return false;
    }
    m_hCheck.fill(HIST("v0_selections"), V0Selections::kV0dcaPV);
    if (cosPA < v0setting_cosPA) {
      return false;
    }
    m_hCheck.fill(HIST("v0_selections"), V0Selections::kV0CosPA);
    return true;
  }

  bool qualitySelectionCascade(const double dcaCascDaughters, const double cosPA)
  {
    if (dcaCascDaughters > cascsetting_dcaCascDaughters) {
      return false;
    }
    m_hCheck.fill(HIST("casc_selections"), CascSelections::kCascDCA);
    if (cosPA < cascsetting_cosPA) {
      return false;
    }
    m_hCheck.fill(HIST("casc_selections"), CascSelections::kCascCosPA);
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

  void fillV0Cand(const std::array<float, 3>& PV, const aod::V0s::iterator& v0, CandidateV0& candV0)
  {
    m_hCheck.fill(HIST("v0_selections"), V0Selections::kV0NoCut);

    auto posTrack = v0.posTrack_as<TracksFullIU>();
    auto negTrack = v0.negTrack_as<TracksFullIU>();
    if (!qualitySelectionV0Daughter(posTrack) || !qualitySelectionV0Daughter(negTrack)) {
      return;
    }
    m_hCheck.fill(HIST("v0_selections"), V0Selections::kV0DaughterQuality);

    auto daughterTrackCovarianceA = getTrackParCov(posTrack);
    auto daughterTrackCovarianceB = getTrackParCov(negTrack);
    if (!initializeFitter(daughterTrackCovarianceA, daughterTrackCovarianceB)) {
      return;
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
    // float ptV0 = decayV0.pt();
    /*
    if (ptV0 < v0setting_ptMin || ptV0 > v0setting_ptMax) {
      return;
    }
    */
    /*
    if (std::abs(decayV0.eta()) > v0setting_etaMaxV0) {
      return;
    }
    */

    // float dcaV0daughters = std::sqrt(m_fitter.getChi2AtPCACandidate());
    float radiusV0 = std::hypot(decayVtx[0], decayVtx[1]);
    float dcaV0toPV = std::hypot(decayVtx[0] - PV[0], decayVtx[1] - PV[1]);
    float cosPA = RecoDecay::cpa(PV, decayVtx, momMother);
    // if (!qualitySelectionV0(dcaV0toPV, dcaV0daughters, radiusV0, cosPA)) {
    //   return;
    // }

    // mass hypothesis
    float massLambdaV0 = computeMassMother(o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, momPos, momNeg, momMother);
    float massK0sV0 = computeMassMother(o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged, momPos, momNeg, momMother);
    // float massPhotonV0 = computeMassMother(o2::constants::physics::MassElectron, o2::constants::physics::MassElectron, momPos, momNeg, momMother);
    m_hCheck.fill(HIST("invmass_Lambda"), massLambdaV0);
    m_hCheck.fill(HIST("invmass_K0s"), massK0sV0);

    uint8_t selectionMap{BIT(Lambda) | BIT(AntiLambda) | BIT(K0s) | BIT(Photon)};
    if (!v0.isPhotonV0()) {
      CLRBIT(selectionMap, Photon);
    }
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
      return;
    }
    m_hCheck.fill(HIST("v0_selections"), V0Selections::kV0PID);

    int pdgCodeV0 = -999;
    if (TESTBIT(selectionMap, Lambda)) {
      pdgCodeV0 = kLambda0;
    } else if (TESTBIT(selectionMap, AntiLambda)) {
      pdgCodeV0 = -kLambda0;
    } else if (TESTBIT(selectionMap, K0s)) {
      pdgCodeV0 = kK0Short;
    } else if (TESTBIT(selectionMap, Photon)) {
      pdgCodeV0 = kGamma;
    } else {
      return;
    }

    gpu::gpustd::array<float, 2> dcaInfo;
    float dcaToPVpos = dcaDaughterToPV(PV, daughterTrackCovarianceA, dcaInfo);
    if (dcaToPVpos < v0setting_dcaDaughtersToPV /*&& std::abs(dcaInfo[0]) < v0setting_dcaDaughtersToPV*/) {
      return;
    }
    float dcaToPVneg = dcaDaughterToPV(PV, daughterTrackCovarianceB, dcaInfo);
    if (dcaToPVneg < v0setting_dcaDaughtersToPV /*&& std::abs(dcaInfo[0]) < v0setting_dcaDaughtersToPV*/) {
      return;
    }

    // fill the candidate and add it to the Armenteros plot
    m_hCheck.fill(HIST("armenteros_plot"), alphaAP, qtAP);

    V0TrackParCov v0TrackParCov{v0.globalIndex(), m_fitter.createParentTrackParCov()};
    m_v0TrackParCovs.push_back(v0TrackParCov);

    candV0.p_V0 = std::hypot(momMother[0], momMother[1], momMother[2]);
    candV0.eta_V0 = RecoDecay::eta(momMother);
    candV0.phi_V0 = RecoDecay::phi(momMother);
    candV0.pdgCode_V0 = pdgCodeV0;
    candV0.radius_V0 = radiusV0;
    candV0.dcaPV_V0 = dcaV0toPV;
    candV0.cosPA_V0 = cosPA;
    candV0.alphaAP_V0 = alphaAP;
    candV0.qtAP_V0 = qtAP;
    candV0.p_pos = std::hypot(momPos[0], momPos[1], momPos[2]) * posTrack.sign();
    candV0.eta_pos = RecoDecay::eta(momPos);
    candV0.phi_pos = RecoDecay::phi(momPos);
    candV0.pTPC_pos = posTrack.tpcInnerParam() * posTrack.sign();
    candV0.dcaPV_pos = dcaToPVpos;
    candV0.itsClsize_pos = posTrack.itsClusterSizes();
    candV0.tpcSignal_pos = posTrack.tpcSignal();
    candV0.tpcNcls_pos = posTrack.tpcNClsFound();
    candV0.nSigmaTPCel_pos = posTrack.tpcNSigmaEl();
    candV0.nSigmaTPCpi_pos = posTrack.tpcNSigmaPi();
    candV0.nSigmaTPCka_pos = posTrack.tpcNSigmaKa();
    candV0.nSigmaTPCpr_pos = posTrack.tpcNSigmaPr();
    candV0.nSigmaTPCde_pos = posTrack.tpcNSigmaDe();
    candV0.nSigmaTPChe_pos = posTrack.tpcNSigmaHe();
    candV0.nSigmaTOFel_pos = posTrack.tofNSigmaEl();
    candV0.nSigmaTOFpi_pos = posTrack.tofNSigmaPi();
    candV0.nSigmaTOFka_pos = posTrack.tofNSigmaKa();
    candV0.nSigmaTOFpr_pos = posTrack.tofNSigmaPr();
    candV0.nSigmaTOFde_pos = posTrack.tofNSigmaDe();
    candV0.nSigmaTOFhe_pos = posTrack.tofNSigmaHe();
    candV0.chi2its_pos = posTrack.itsChi2NCl();
    candV0.chi2tpc_pos = posTrack.tpcChi2NCl();
    candV0.hasTPC_pos = posTrack.hasTPC();
    candV0.p_neg = std::hypot(momNeg[0], momNeg[1], momNeg[2]) * negTrack.sign();
    candV0.eta_neg = RecoDecay::eta(momNeg);
    candV0.phi_neg = RecoDecay::phi(momNeg);
    candV0.pTPC_neg = negTrack.tpcInnerParam() * negTrack.sign();
    candV0.dcaPV_neg = dcaToPVneg;
    candV0.itsClsize_neg = negTrack.itsClusterSizes();
    candV0.tpcSignal_neg = negTrack.tpcSignal();
    candV0.tpcNcls_neg = negTrack.tpcNClsFound();
    candV0.nSigmaTPCel_neg = negTrack.tpcNSigmaEl();
    candV0.nSigmaTPCpi_neg = negTrack.tpcNSigmaPi();
    candV0.nSigmaTPCka_neg = negTrack.tpcNSigmaKa();
    candV0.nSigmaTPCpr_neg = negTrack.tpcNSigmaPr();
    candV0.nSigmaTPCde_neg = negTrack.tpcNSigmaDe();
    candV0.nSigmaTPChe_neg = negTrack.tpcNSigmaHe();
    candV0.nSigmaTOFel_neg = negTrack.tofNSigmaEl();
    candV0.nSigmaTOFpi_neg = negTrack.tofNSigmaPi();
    candV0.nSigmaTOFka_neg = negTrack.tofNSigmaKa();
    candV0.nSigmaTOFpr_neg = negTrack.tofNSigmaPr();
    candV0.nSigmaTOFde_neg = negTrack.tofNSigmaDe();
    candV0.nSigmaTOFhe_neg = negTrack.tofNSigmaHe();
    candV0.chi2its_neg = negTrack.itsChi2NCl();
    candV0.chi2tpc_neg = negTrack.tpcChi2NCl();
    candV0.hasTPC_neg = negTrack.hasTPC();
  }

  void fillV0Table(const CandidateV0& candV0)
  {
    m_v0Table(
      candV0.p_V0,
      candV0.eta_V0,
      candV0.phi_V0,
      candV0.pdgCode_V0,
      candV0.radius_V0,
      candV0.dcaPV_V0,
      candV0.cosPA_V0,
      candV0.alphaAP_V0,
      candV0.qtAP_V0,
      candV0.mc_pdgCode_V0,
      candV0.p_pos,
      candV0.eta_pos,
      candV0.phi_pos,
      candV0.pTPC_pos,
      candV0.dcaPV_pos,
      candV0.itsClsize_pos,
      candV0.tpcSignal_pos,
      candV0.tpcNcls_pos,
      candV0.nSigmaTPCel_pos,
      candV0.nSigmaTPCpi_pos,
      candV0.nSigmaTPCka_pos,
      candV0.nSigmaTPCpr_pos,
      candV0.nSigmaTPCde_pos,
      candV0.nSigmaTPChe_pos,
      candV0.nSigmaTOFel_pos,
      candV0.nSigmaTOFpi_pos,
      candV0.nSigmaTOFka_pos,
      candV0.nSigmaTOFpr_pos,
      candV0.nSigmaTOFde_pos,
      candV0.nSigmaTOFhe_pos,
      candV0.chi2its_pos,
      candV0.chi2tpc_pos,
      candV0.hasTPC_pos,
      candV0.mc_pdgCode_pos,
      candV0.p_neg,
      candV0.eta_neg,
      candV0.phi_neg,
      candV0.pTPC_neg,
      candV0.dcaPV_neg,
      candV0.itsClsize_neg,
      candV0.tpcSignal_neg,
      candV0.tpcNcls_neg,
      candV0.nSigmaTPCel_neg,
      candV0.nSigmaTPCpi_neg,
      candV0.nSigmaTPCka_neg,
      candV0.nSigmaTPCpr_neg,
      candV0.nSigmaTPCde_neg,
      candV0.nSigmaTPChe_neg,
      candV0.nSigmaTOFel_neg,
      candV0.nSigmaTOFpi_neg,
      candV0.nSigmaTOFka_neg,
      candV0.nSigmaTOFpr_neg,
      candV0.nSigmaTOFde_neg,
      candV0.nSigmaTOFhe_neg,
      candV0.chi2its_neg,
      candV0.chi2tpc_neg,
      candV0.hasTPC_neg,
      candV0.mc_pdgCode_neg);
  }

  void fillKCand(const std::array<float, 3>& PV, const aod::Cascades::iterator& cascade, CandidateK& candK)
  {
    m_hCheck.fill(HIST("casc_selections"), Selections::kNoCut);

    auto v0Track = cascade.template v0_as<aod::V0s>();
    auto bachelorTrack = cascade.template bachelor_as<TracksFullIU>();
    auto bachelorTrackPar = getTrackPar(bachelorTrack);

    auto itv0 = std::find_if(m_v0TrackParCovs.begin(), m_v0TrackParCovs.end(), [&](const V0TrackParCov& v0) { return v0.globalIndex == v0Track.globalIndex(); });
    if (itv0 == m_v0TrackParCovs.end()) {
      return;
    }

    auto v0TrackCovariance = itv0->trackParCov;
    auto bachelorTrackCovariance = getTrackParCov(bachelorTrack);
    if (!initializeFitter(v0TrackCovariance, bachelorTrackCovariance)) {
      return;
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
      return;
    }
    gpu::gpustd::array<float, 2> dcaInfo;
    float dcaToPVbachelor = dcaDaughterToPV(PV, bachelorTrackCovariance, dcaInfo);

    float massXi = computeMassMother(o2::constants::physics::MassLambda0, o2::constants::physics::MassPionCharged, momV0, momBachelor, momMother);
    float massOmega = computeMassMother(o2::constants::physics::MassLambda0, o2::constants::physics::MassKaonCharged, momV0, momBachelor, momMother);
    m_hCheck.fill(HIST("Xi_vs_Omega"), massOmega, massXi);
    if (std::abs(massOmega - o2::constants::physics::MassOmegaMinus) > cascsetting_massWindowOmega) {
      return;
    }
    m_hCheck.fill(HIST("massOmegaWithBkg"), massOmega);
    m_hCheck.fill(HIST("casc_selections"), CascSelections::kAcceptedOmega);
    if (std::abs(massXi - o2::constants::physics::MassXiMinus) < cascsetting_massWindowXi) {
      return;
    } // enhance purity by rejecting Xi background
    m_hCheck.fill(HIST("massOmega"), massOmega);
    m_hCheck.fill(HIST("casc_selections"), CascSelections::kRejectedXi);

    candK.p_casc = std::hypot(momMother[0], momMother[1], momMother[2]);
    candK.eta_casc = RecoDecay::eta(momMother);
    candK.phi_casc = RecoDecay::phi(momMother);
    candK.pdgCode_casc = kOmegaMinus * int(bachelorTrack.sign());
    candK.radius_casc = radius;
    candK.dcaPV_casc = dcaCascToPV;
    candK.cosPA_casc = cosPA;
    candK.p_K = std::hypot(momBachelor[0], momBachelor[1], momBachelor[2]) * bachelorTrack.sign();
    candK.eta_K = RecoDecay::eta(momBachelor);
    candK.phi_K = RecoDecay::phi(momBachelor);
    candK.pTPC_K = bachelorTrack.tpcInnerParam() * bachelorTrack.sign();
    candK.dcaPV_K = dcaToPVbachelor;
    candK.itsClsize_K = bachelorTrack.itsClusterSizes();
    candK.tpcSignal_K = bachelorTrack.tpcSignal();
    candK.tpcNcls_K = bachelorTrack.tpcNClsFound();
    candK.nSigmaTPCel_K = bachelorTrack.tpcNSigmaEl();
    candK.nSigmaTPCpi_K = bachelorTrack.tpcNSigmaPi();
    candK.nSigmaTPCka_K = bachelorTrack.tpcNSigmaKa();
    candK.nSigmaTPCpr_K = bachelorTrack.tpcNSigmaPr();
    candK.nSigmaTPCde_K = bachelorTrack.tpcNSigmaDe();
    candK.nSigmaTPChe_K = bachelorTrack.tpcNSigmaHe();
    candK.nSigmaTOFel_K = bachelorTrack.tofNSigmaEl();
    candK.nSigmaTOFpi_K = bachelorTrack.tofNSigmaPi();
    candK.nSigmaTOFka_K = bachelorTrack.tofNSigmaKa();
    candK.nSigmaTOFpr_K = bachelorTrack.tofNSigmaPr();
    candK.nSigmaTOFde_K = bachelorTrack.tofNSigmaDe();
    candK.nSigmaTOFhe_K = bachelorTrack.tofNSigmaHe();
    candK.chi2its_K = bachelorTrack.itsChi2NCl();
    candK.chi2tpc_K = bachelorTrack.tpcChi2NCl();
    candK.hasTPC_K = bachelorTrack.hasTPC();
  }

  void fillKTable(const CandidateK& candK)
  {
    m_KTable(
      candK.p_casc,          // p_casc
      candK.eta_casc,        // eta_casc
      candK.phi_casc,        // phi_casc
      candK.pdgCode_casc,    // pdgCode_casc
      candK.radius_casc,     // radius_casc
      candK.dcaPV_casc,      // dcaPV_casc
      candK.cosPA_casc,      // cosPA_casc
      candK.mc_pdgCode_casc, // mc_pdgCode_casc
      candK.p_K,             // p_K
      candK.eta_K,           // eta_K
      candK.phi_K,           // phi_K
      candK.pTPC_K,          // pTPC_K
      candK.dcaPV_K,         // dcaPV_K
      candK.itsClsize_K,     // itsClSize_K
      candK.tpcSignal_K,     // tpcSignal_K
      candK.tpcNcls_K,       // nTPCcls_K
      candK.nSigmaTPCel_K,   // nSigmaTPCe_K
      candK.nSigmaTPCpi_K,   // nSigmaTPCpi_K
      candK.nSigmaTPCka_K,   // nSigmaTPCk_K
      candK.nSigmaTPCpr_K,   // nSigmaTPCp_K
      candK.nSigmaTPCde_K,   // nSigmaTPCdeu_K
      candK.nSigmaTPChe_K,   // nSigmaTPChe3_K
      candK.nSigmaTOFel_K,   // nSigmaTOFe_K
      candK.nSigmaTOFpi_K,   // nSigmaTOFpi_K
      candK.nSigmaTOFka_K,   // nSigmaTOFk_K
      candK.nSigmaTOFpr_K,   // nSigmaTOFp_K
      candK.nSigmaTOFde_K,   // nSigmaTOFdeu_K
      candK.nSigmaTOFhe_K,   // nSigmaTOFhe3_K
      candK.chi2its_K,       // chi2its_K
      candK.chi2tpc_K,       // chi2tpc_K
      candK.hasTPC_K,        // hasTPC_K
      candK.mc_pdgCode_K     // mc_pdgCode_K
    );
  }

  void fillDeTree(const TracksFullIU& tracks)
  {
    for (auto& track : tracks) {

      m_hCheck.fill(HIST("de_selections"), DeSelections::kDeNoCut);
      if (track.pidForTracking() != o2::track::PID::Deuteron) {
        continue;
      }
      m_hCheck.fill(HIST("de_selections"), DeSelections::kDePIDforTrk);
      if (!selectionPIDDe(track)) {
        continue;
      }
      m_hCheck.fill(HIST("de_selections"), DeSelections::kDePID);
      m_hCheck.fill(HIST("BetheBlochDe"), track.p() * track.sign(), track.tpcSignal());

      m_nucleiTable(
        track.p() * track.sign(), // p_De,
        // track.pt() * track.sign(),            // pt_De,
        track.eta(),                          // eta_De,
        track.phi(),                          // phi_De,
        track.tpcInnerParam() * track.sign(), // pTPC_De,
        physics::pdg::Deuteron,               // pdgCode_De,
        0.,                                   // dcaPV_De,
        track.itsClusterSizes(),              // itsClSize_De,
        track.tpcSignal(),                    // tpcSignal_De,
        track.tpcNClsFound(),                 // nTPCcls_De,
        track.tpcNSigmaEl(),                  // nSigmaTPCe_De,
        track.tpcNSigmaPi(),                  // nSigmaTPCpi_De,
        track.tpcNSigmaKa(),                  // nSigmaTPCk_De,
        track.tpcNSigmaPr(),                  // nSigmaTPCp_De,
        track.tpcNSigmaDe(),                  // nSigmaTPCdeu_De,
        track.tpcNSigmaHe(),                  // nSigmaTPChe3_De,
        track.tofNSigmaEl(),                  // nSigmaTOFe_De,
        track.tofNSigmaPi(),                  // nSigmaTOFpi_De,
        track.tofNSigmaKa(),                  // nSigmaTOFk_De,
        track.tofNSigmaPr(),                  // nSigmaTOFp_De,
        track.tofNSigmaDe(),                  // nSigmaTOFdeu_De,
        track.tofNSigmaHe(),                  // nSigmaTOFhe3_De,
        track.itsChi2NCl(),                   // chi2its_De,
        track.tpcChi2NCl(),                   // chi2tpc_De,
        track.hasTPC()                        // hasTPC_De
      );
    }
  }

  void fillHe3Tree(const TracksFullIU& tracks)
  {
    for (auto& track : tracks) {

      m_hCheck.fill(HIST("he3_selections"), He3Selections::kHe3NoCut);

      if (track.pidForTracking() != o2::track::PID::Helium3) {
        continue;
      }
      m_hCheck.fill(HIST("he3_selections"), He3Selections::kHe3PIDforTrk);
      m_hCheck.fill(HIST("BetheBlochHe3"), track.p() * track.sign(), track.tpcSignal());
      if (!selectionPIDHe3(track)) {
        continue;
      }
      m_hCheck.fill(HIST("he3_selections"), He3Selections::kHe3PID);
      // m_hCheck.fill(HIST("BetheBlochHe3"), track.p() * track.sign(), track.tpcSignal());

      m_nucleiTable(
        track.p() * track.sign(), // p_He3,
        // track.pt() * track.sign(),            // pt_He3,
        track.eta(),                          // eta_He3,
        track.phi(),                          // phi_He3,
        track.tpcInnerParam() * track.sign(), // pTPC_He3,
        physics::pdg::Deuteron,               // pdgCode_He3,
        0.,                                   // dcaPV_He3,
        track.itsClusterSizes(),              // itsClSize_He3,
        track.tpcSignal(),                    // tpcSignal_He3,
        track.tpcNClsFound(),                 // nTPCcls_He3,
        track.tpcNSigmaEl(),                  // nSigmaTPCe_He3,
        track.tpcNSigmaPi(),                  // nSigmaTPCpi_He3,
        track.tpcNSigmaKa(),                  // nSigmaTPCk_He3,
        track.tpcNSigmaPr(),                  // nSigmaTPCp_He3,
        track.tpcNSigmaDe(),                  // nSigmaTPCdeu_He3,
        track.tpcNSigmaHe(),                  // nSigmaTPChe3_He3,
        track.tofNSigmaEl(),                  // nSigmaTOFe_He3,
        track.tofNSigmaPi(),                  // nSigmaTOFpi_He3,
        track.tofNSigmaKa(),                  // nSigmaTOFk_He3,
        track.tofNSigmaPr(),                  // nSigmaTOFp_He3,
        track.tofNSigmaDe(),                  // nSigmaTOFdeu_He3,
        track.tofNSigmaHe(),                  // nSigmaTOFhe3_He3,
        track.itsChi2NCl(),                   // chi2its_He3,
        track.tpcChi2NCl(),                   // chi2tpc_He3,
        track.hasTPC()                        // hasTPC_He3,
      );
    }
  }

  void processRun3(CollisionsCustom const& collisions, TracksFullIU const& tracks, aod::V0s const& v0s, aod::Cascades const& cascades, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      m_collisionCounter++;
      // LOG(info) << "Processing collision " << m_collisionCounter << " with zVtx = " << collision.posZ();

      if (!collisionSelection(collision)) {
        continue;
      }

      m_hCheck.fill(HIST("zVtx"), collision.posZ());
      std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};

      const uint64_t collIdx = collision.globalIndex();
      auto v0Table_thisCollision = v0s.sliceBy(m_perCollisionV0, collIdx);
      auto cascTable_thisCollision = cascades.sliceBy(m_perCollisionCascade, collIdx);
      v0Table_thisCollision.bindExternalIndices(&tracks);
      cascTable_thisCollision.bindExternalIndices(&tracks);
      cascTable_thisCollision.bindExternalIndices(&v0s);

      if (setting_fillV0) {
        m_v0TrackParCovs.clear();
        for (auto& v0 : v0Table_thisCollision) {
          CandidateV0 candV0;
          fillV0Cand(PV, v0, candV0);
          fillV0Table(candV0);
        }
      }
      // fillV0Tree(PV, v0Table_thisCollision);
      if (setting_fillK && setting_fillV0) { // the v0 loops are needed for the Ks
        for (auto& cascade : cascTable_thisCollision) {
          CandidateK candK;
          fillKCand(PV, cascade, candK);
          fillKTable(candK);
        }
      }
      if (setting_fillDe)
        fillDeTree(tracks);
      if (setting_fillHe3)
        fillHe3Tree(tracks);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processRun3, "process Run 3", false);

}; // LfTreeCreatorClusterStudies

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LfTreeCreatorClusterStudies>(cfgc)};
}
