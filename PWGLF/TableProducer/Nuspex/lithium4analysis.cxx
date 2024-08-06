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
// Analysis task for anti-lithium4 analysis

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include <cmath>
#include <array>
#include <cstdlib>

#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "CCDB/BasicCCDBManager.h"

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "DCAFitter/DCAFitterN.h"
#include "PWGLF/DataModel/LFLithium4Tables.h"
#include "PWGLF/Utils/svPoolCreator.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using McIter = aod::McParticles::iterator;
using CollBracket = o2::math_utils::Bracket<int>;
using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;
using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TOFSignal, aod::TOFEvTime>;
using TrackCandidatesMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TOFSignal, aod::TOFEvTime, aod::McTrackLabels>;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};

constexpr float he3Mass = o2::constants::physics::MassHelium3;
constexpr float protonMass = o2::constants::physics::MassProton;
constexpr int li4PDG = 1000030040;
constexpr int prPDG = 2212;
constexpr int hePDG = 1000020030;

enum Selections {
  kNoCuts = 0,
  kGlobalTrack,
  kTrackCuts,
  kPID,
  kAll
};

} // namespace

struct Lithium4Candidate {

  float sign = 0.f;
  std::array<float, 3> momHe3 = {99.f, 99.f, 99.f};
  std::array<float, 3> momPr = {99.f, 99.f, 99.f};

  float recoPtHe3() const { return sign * std::hypot(momHe3[0], momHe3[1]); }
  float recoPhiHe3() const { return std::atan2(momHe3[1], momHe3[0]); }
  float recoEtaHe3() const { return std::asinh(momHe3[2] / recoPtHe3()); }
  float recoPtPr() const { return sign * std::hypot(momPr[0], momPr[1]); }
  float recoPhiPr() const { return std::atan2(momPr[1], momPr[0]); }
  float recoEtaPr() const { return std::asinh(momPr[2] / recoPtPr()); }

  float DCAxyHe3 = -10.f;
  float DCAzHe3 = -10.f;
  float DCAxyPr = -10.f;
  float DCAzPr = -10.f;

  uint16_t tpcSignalHe3 = 0u;
  float momHe3TPC = -99.f;
  uint16_t tpcSignalPr = 0u;
  float momPrTPC = -99.f;

  uint8_t nTPCClustersHe3 = 0u;
  float nSigmaHe3 = -10.f;
  float nSigmaPr = -10.f;

  float massTOFHe3 = -10;
  float massTOFPr = -10;

  uint32_t PIDtrkHe3 = 0xFFFFF; // PID in tracking
  uint32_t PIDtrkPr = 0xFFFFF;

  uint32_t itsClSizeHe3 = 0u;
  uint32_t itsClSizePr = 0u;

  uint8_t sharedClustersHe3 = 0u;
  uint8_t sharedClustersPr = 0u;

  bool isBkgUS = false;
  bool isBkgEM = false;

  float momHe3MC = -99.f;
  float momPrMC = -99.f;

  float l4PtMC = -99.f;
  float l4MassMC = -10.f;

  float invMass = -10.f;
};

struct lithium4analysis {

  Produces<aod::Lithium4Table> m_outputDataTable;
  Produces<aod::Lithium4TableMC> m_outputMCTable;
  Service<o2::ccdb::BasicCCDBManager> m_ccdb;

  // Selections
  Configurable<float> setting_cutVertex{"setting_cutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> setting_cutPT{"setting_cutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> setting_cutMaxPrPT{"setting_cutMaxPrPT", 1.8, "Max PT cut on proton"};
  Configurable<float> setting_cutEta{"setting_cutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> setting_cutDCAxy{"setting_cutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> setting_cutDCAz{"setting_cutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> setting_nsigmaCutTPC{"setting_nsigmaCutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> setting_nsigmaCutTOF{"setting_nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> setting_noMixedEvents{"setting_noMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> setting_enableBkgUS{"setting_enableBkgUS", false, "Enable US background"};
  Configurable<bool> setting_isMC{"setting_isMC", false, "Run MC"};
  Configurable<bool> setting_skipAmbiTracks{"setting_skipAmbiTracks", false, "Skip ambiguous tracks"};
  Configurable<float> setting_customVertexerTimeMargin{"setting_customVertexerTimeMargin", 800, "Time margin for custom vertexer (ns)"};

  Configurable<LabeledArray<double>> setting_BetheBlochParams{"setting_BetheBlochParams", {betheBlochDefault[0], 1, 6, {"He3"}, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
  Configurable<bool> setting_compensatePIDinTracking{"setting_compensatePIDinTracking", false, "If true, divide tpcInnerParam by the electric charge"};
  Configurable<int> setting_materialCorrection{"setting_materialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Material correction type"};

  // CCDB options
  Configurable<double> ccdbsetting_d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdbsetting_ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbsetting_grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> ccdbsetting_grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> ccdbsetting_lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> ccdbsetting_geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> ccdbsetting_pidPath{"pidPath", "", "Path to the PID response object"};

  o2::pid::tof::Beta<TrackCandidates::iterator> m_responseBeta;
  o2::pid::tof::Beta<TrackCandidatesMC::iterator> m_responseBetaMC;

  Preslice<TrackCandidates> perCol = aod::track::collisionId;
  Preslice<TrackCandidatesMC> perColMC = aod::track::collisionId;

  // binning for EM background
  ConfigurableAxis axisVertex{"axisVertex", {30, -10, 10}, "vertex axis for bin"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;
  BinningType binningOnPositions{{axisVertex}, true};
  SliceCache cache;
  SameKindPair<EventCandidates, TrackCandidates, BinningType> pair{binningOnPositions, setting_noMixedEvents, -1, &cache};

  // Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<2> m_fitter;
  svPoolCreator m_svPoolCreator{hePDG, prPDG};

  int m_runNumber = 0;
  float m_d_bz = 0.0f;
  std::array<float, 6> m_BBparamsHe;

  std::vector<int> m_recoCollisionIDs;
  std::vector<bool> m_goodCollisions;
  std::vector<SVCand> m_trackPairs;

  HistogramRegistry m_hAnalysis{
    "histos",
    {{"hCentrality", "Centrality distribution", {HistType::kTH1F, {{2001, -0.5, 2000.5}}}},
     {"hVtxZ", "Vertex distribution in Z;Z (cm)", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
     {"hNcontributor", "Number of primary vertex contributor", {HistType::kTH1F, {{2000, 0.0f, 2000.0f}}}},
     {"hDCAxyHe3", ";DCA_{xy} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hDCAzHe3", ";DCA_{z} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hLitInvMass", "; M(^{3}He + p) (GeV/#it{c}^{2})", {HistType::kTH1F, {{50, 3.74f, 3.85f}}}},
     {"hHe3Pt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{200, 0.0f, 6.0f}}}},
     {"hProtonPt", "Pt distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{200, 0.0f, 3.0f}}}},
     {"h2dEdxHe3candidates", "dEdx distribution; Signed #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
     {"h2NsigmaHe3TPC", "NsigmaHe3 TPC distribution; Signed #it{p}/#it{z} (GeV/#it{c}); n#sigma_{TPC}({}^{3}He)", {HistType::kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaProtonTPC", "NsigmaProton TPC distribution; Signed #it{p}/#it{z} (GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaProtonTOF", "NsigmaProton TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}}}},
     {"hTrackSel", "Accepted tracks", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}}},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  void init(o2::framework::InitContext&)
  {
    m_runNumber = 0;
    m_d_bz = 0;

    m_ccdb->setURL(ccdbsetting_ccdburl);
    m_ccdb->setCaching(true);
    m_ccdb->setLocalObjectValidityChecking();
    m_ccdb->setFatalWhenNull(false);

    m_fitter.setPropagateToPCA(true);
    m_fitter.setMaxR(200.);
    m_fitter.setMinParamChange(1e-3);
    m_fitter.setMinRelChi2Change(0.9);
    m_fitter.setMaxDZIni(1e9);
    m_fitter.setMaxChi2(1e9);
    m_fitter.setUseAbsDCA(true);
    int mat{static_cast<int>(setting_materialCorrection)};
    m_fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

    m_svPoolCreator.setTimeMargin(setting_customVertexerTimeMargin);
    if (setting_skipAmbiTracks) {
      m_svPoolCreator.setSkipAmbiTracks();
    }

    for (int i = 0; i < 5; i++) {
      m_BBparamsHe[i] = setting_BetheBlochParams->get("He3", Form("p%i", i));
    }
    m_BBparamsHe[5] = setting_BetheBlochParams->get("He3", "resolution");

    std::vector<std::string> selection_labels = {"All", "Global track", "Track selection", "PID {}^{3}He"};
    for (int i = 0; i < Selections::kAll; i++) {
      m_hAnalysis.get<TH1>(HIST("hTrackSel"))->GetXaxis()->SetBinLabel(i + 1, selection_labels[i].c_str());
    }
  }

  void initCCDB(const aod::BCsWithTimestamps::iterator& bc)
  {
    if (m_runNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();

    /* Francesco */
    o2::parameters::GRPObject* grpo = m_ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbsetting_grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (ccdbsetting_d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        m_d_bz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << m_d_bz << " kZG";
      } else {
        m_d_bz = ccdbsetting_d_bz_input;
      }
    } else {
      grpmag = m_ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbsetting_grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdbsetting_grpmagPath << " of object GRPMagField and " << ccdbsetting_grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (ccdbsetting_d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        m_d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << m_d_bz << " kZG";
      } else {
        m_d_bz = ccdbsetting_d_bz_input;
      }
    }

    if (!ccdbsetting_pidPath.value.empty()) {
      auto he3pid = m_ccdb->getForTimeStamp<std::array<float, 6>>(ccdbsetting_pidPath.value + "_He3", run3grp_timestamp);
      std::copy(he3pid->begin(), he3pid->end(), m_BBparamsHe.begin());
    } else {
      for (int i = 0; i < 5; i++) {
        m_BBparamsHe[i] = setting_BetheBlochParams->get("He3", Form("p%i", i));
      }
      m_BBparamsHe[5] = setting_BetheBlochParams->get("He3", "resolution");
    }
    m_fitter.setBz(m_d_bz);
    m_runNumber = bc.runNumber();
  }

  // ==================================================================================================================

  template <typename C>
  void selectCollisions(const C& collisions)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!collision.sel8() || std::abs(collision.posZ()) > setting_cutVertex) {
        continue;
      }
      m_goodCollisions[collision.globalIndex()] = true;
    }
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {

    if (candidate.itsNCls() < 5 ||
        candidate.tpcNClsFound() < 90 || // candidate.tpcNClsFound() < 70 ||
        candidate.tpcNClsCrossedRows() < 70 ||
        candidate.tpcNClsCrossedRows() < 0.8 * candidate.tpcNClsFindable() ||
        candidate.tpcChi2NCl() > 4.f ||
        candidate.itsChi2NCl() > 36.f) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool selectionPIDProton(const T& candidate)
  {
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPr()) < setting_nsigmaCutTOF && std::abs(candidate.tpcNSigmaPr()) < setting_nsigmaCutTPC) {
        m_hAnalysis.fill(HIST("h2NsigmaProtonTPC"), candidate.tpcInnerParam(), candidate.tpcNSigmaPr());
        m_hAnalysis.fill(HIST("h2NsigmaProtonTOF"), candidate.p(), candidate.tofNSigmaPr());
        return true;
      }
    } else if (std::abs(candidate.tpcNSigmaPr()) < setting_nsigmaCutTPC) {
      m_hAnalysis.fill(HIST("h2NsigmaProtonTPC"), candidate.tpcInnerParam(), candidate.tpcNSigmaPr());
      return true;
    }
    return false;
  }

  template <typename T>
  float computeNSigmaHe3(const T& candidate)
  {
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;

    float correctedTPCinnerParam = (heliumPID && setting_compensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(correctedTPCinnerParam * 2.f / constants::physics::MassHelium3), m_BBparamsHe[0], m_BBparamsHe[1], m_BBparamsHe[2], m_BBparamsHe[3], m_BBparamsHe[4]);

    double resoTPC{expTPCSignal * m_BBparamsHe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename T>
  bool selectionPIDHe3(const T& candidate)
  {
    auto nSigmaHe3 = computeNSigmaHe3(candidate);
    if (std::abs(nSigmaHe3) < setting_nsigmaCutTPC) {
      return true;
    }
    return false;
  }

  // ==================================================================================================================

  template <typename T1, typename T2>
  bool fillCandidateInfo(const T1& trackHe3, const T2& trackPr, Lithium4Candidate& li4cand, bool mix)
  {
    li4cand.momHe3 = std::array{2 * trackHe3.px(), 2 * trackHe3.py(), 2 * trackHe3.pz()};
    li4cand.momPr = std::array{trackPr.px(), trackPr.py(), trackPr.pz()};

    float invMass = RecoDecay::m(array{li4cand.momHe3, li4cand.momPr}, std::array{o2::constants::physics::MassHelium3, o2::constants::physics::MassProton});

    if (invMass < 3.74 || invMass > 3.85 || trackPr.pt() > setting_cutMaxPrPT) {
      return false;
    }

    li4cand.sign = trackHe3.sign();

    li4cand.DCAxyHe3 = trackHe3.dcaXY();
    li4cand.DCAzHe3 = trackHe3.dcaZ();
    li4cand.DCAxyPr = trackPr.dcaXY();
    li4cand.DCAzPr = trackPr.dcaZ();

    li4cand.tpcSignalHe3 = trackHe3.tpcSignal();
    bool heliumPID = trackHe3.pidForTracking() == o2::track::PID::Helium3 || trackHe3.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParamHe3 = (heliumPID && setting_compensatePIDinTracking) ? trackHe3.tpcInnerParam() / 2.f : trackHe3.tpcInnerParam();
    li4cand.momHe3TPC = correctedTPCinnerParamHe3;
    li4cand.tpcSignalPr = trackPr.tpcSignal();
    li4cand.momPrTPC = trackPr.tpcInnerParam();

    li4cand.nTPCClustersHe3 = trackHe3.tpcNClsFound();
    li4cand.nSigmaHe3 = computeNSigmaHe3(trackHe3);
    li4cand.nSigmaPr = trackPr.tpcNSigmaPr();

    li4cand.PIDtrkHe3 = trackHe3.pidForTracking();
    li4cand.PIDtrkPr = trackPr.pidForTracking();

    li4cand.itsClSizeHe3 = trackHe3.itsClusterSizes();
    li4cand.itsClSizePr = trackPr.itsClusterSizes();

    li4cand.sharedClustersHe3 = trackHe3.tpcNClsShared();
    li4cand.sharedClustersPr = trackPr.tpcNClsShared();

    li4cand.isBkgUS = trackHe3.sign() * trackPr.sign() < 0;
    li4cand.isBkgEM = mix;

    li4cand.invMass = invMass;

    return true;
  }

  template <typename T1>
  void fillCandidateBeta(const T1& trackHe3, const T1& trackPr, Lithium4Candidate& li4cand)
  // void fillCandidateBeta(const TrackCandidates::iterator& trackHe3, const TrackCandidates::iterator& trackPr, Lithium4Candidate& li4cand)
  {
    if (trackHe3.hasTOF()) {
      float beta = m_responseBeta.GetBeta(trackHe3);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      bool heliumPID = trackHe3.pidForTracking() == o2::track::PID::Helium3 || trackHe3.pidForTracking() == o2::track::PID::Alpha;
      float correctedTPCinnerParamHe3 = (heliumPID && setting_compensatePIDinTracking) ? trackHe3.tpcInnerParam() / 2.f : trackHe3.tpcInnerParam();
      li4cand.massTOFHe3 = correctedTPCinnerParamHe3 * 2.f * std::sqrt(1.f / (beta * beta) - 1.f);
    }
    if (trackPr.hasTOF()) {
      float beta = m_responseBeta.GetBeta(trackPr);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      li4cand.massTOFPr = trackPr.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
    }
  }

  bool fillMcCandidateInfo(Lithium4Candidate& li4cand, const McIter& mctrackHe3, const McIter& mctrackPr, const McIter& mothertrack)
  {
    li4cand.momHe3MC = mctrackHe3.pt() * (mctrackHe3.pdgCode() > 0 ? 1 : -1);
    li4cand.momPrMC = mctrackPr.pt() * (mctrackPr.pdgCode() > 0 ? 1 : -1);
    li4cand.l4PtMC = mothertrack.pt() * (mothertrack.pdgCode() > 0 ? 1 : -1);
    const double eLit = mctrackHe3.e() + mctrackPr.e();
    li4cand.l4MassMC = std::sqrt(eLit * eLit - mothertrack.p() * mothertrack.p());
    return true;
  }

  template <typename T1>
  void fillMcCandidateBeta(const T1& trackHe3, const T1& trackPr, Lithium4Candidate& li4cand)
  {
    if (trackHe3.hasTOF()) {
      float beta = m_responseBetaMC.GetBeta(trackHe3);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      bool heliumPID = trackHe3.pidForTracking() == o2::track::PID::Helium3 || trackHe3.pidForTracking() == o2::track::PID::Alpha;
      float correctedTPCinnerParamHe3 = (heliumPID && setting_compensatePIDinTracking) ? trackHe3.tpcInnerParam() / 2.f : trackHe3.tpcInnerParam();
      li4cand.massTOFHe3 = correctedTPCinnerParamHe3 * 2.f * std::sqrt(1.f / (beta * beta) - 1.f);
    }
    if (trackPr.hasTOF()) {
      float beta = m_responseBetaMC.GetBeta(trackPr);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      li4cand.massTOFPr = trackPr.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
    }
  }

  template <typename T>
  void pairTracksSameEvent(const T& tracks)
  {
    for (auto track0 : tracks) {

      m_hAnalysis.fill(HIST("hTrackSel"), Selections::kNoCuts);
      bool heliumPID = track0.pidForTracking() == o2::track::PID::Helium3 || track0.pidForTracking() == o2::track::PID::Alpha;

      float correctedTPCinnerParam = (heliumPID && setting_compensatePIDinTracking) ? track0.tpcInnerParam() / 2.f : track0.tpcInnerParam();
      m_hAnalysis.fill(HIST("h2dEdxHe3candidates"), correctedTPCinnerParam * 2.f, track0.tpcSignal());

      if (!track0.isGlobalTrackWoDCA()) {
        continue;
      }
      m_hAnalysis.fill(HIST("hTrackSel"), Selections::kGlobalTrack);

      if (!selectionTrack(track0)) {
        continue;
      }
      m_hAnalysis.fill(HIST("hTrackSel"), Selections::kTrackCuts);

      if (!selectionPIDHe3(track0)) {
        continue;
      }
      m_hAnalysis.fill(HIST("hTrackSel"), Selections::kPID);

      for (auto track1 : tracks) {

        if (track0 == track1) {
          continue;
        }

        if (!setting_enableBkgUS) {
          if (track0.sign() * track1.sign() < 0) {
            continue;
          }
        }

        if (!track1.isGlobalTrackWoDCA() || !selectionTrack(track1) || !selectionPIDProton(track1)) {
          continue;
        }

        SVCand trackPair;
        trackPair.tr0Idx = track0.globalIndex();
        trackPair.tr1Idx = track1.globalIndex();
        m_trackPairs.push_back(trackPair);
      }
    }
  }

  void pairTracksEventMixing()
  {
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (!c1.sel8() || !c2.sel8()) {
        continue;
      }
      m_hAnalysis.fill(HIST("hNcontributor"), c1.numContrib());
      m_hAnalysis.fill(HIST("hVtxZ"), c1.posZ());

      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!t1.isGlobalTrackWoDCA() || !selectionTrack(t1) || !t2.isGlobalTrackWoDCA() || !selectionTrack(t2)) {
          continue;
        }

        TrackCandidates::iterator he3Cand, protonCand;
        if (selectionPIDHe3(t1) && selectionPIDProton(t2)) {
          he3Cand = t1, protonCand = t2;
        } else if (selectionPIDHe3(t2) && selectionPIDProton(t1)) {
          he3Cand = t2, protonCand = t1;
        } else {
          continue;
        }

        bool heliumPID = he3Cand.pidForTracking() == o2::track::PID::Helium3 || he3Cand.pidForTracking() == o2::track::PID::Alpha;
        float correctedTPCinnerParam = (heliumPID && setting_compensatePIDinTracking) ? he3Cand.tpcInnerParam() / 2.f : he3Cand.tpcInnerParam();
        m_hAnalysis.fill(HIST("h2dEdxHe3candidates"), correctedTPCinnerParam * 2.f, he3Cand.tpcSignal());

        SVCand trackPair;
        trackPair.tr0Idx = he3Cand.globalIndex();
        trackPair.tr1Idx = protonCand.globalIndex();
        m_trackPairs.push_back(trackPair);
      }
    }
  }

  void fillTable(const Lithium4Candidate& li4cand, bool isMC = false)
  {
    if (!isMC) {
      m_outputDataTable(
        li4cand.recoPtHe3(),
        li4cand.recoEtaHe3(),
        li4cand.recoPhiHe3(),
        li4cand.recoPtPr(),
        li4cand.recoEtaPr(),
        li4cand.recoPhiPr(),
        li4cand.DCAxyHe3,
        li4cand.DCAzHe3,
        li4cand.DCAxyPr,
        li4cand.DCAzPr,
        li4cand.tpcSignalHe3,
        li4cand.momHe3TPC,
        li4cand.tpcSignalPr,
        li4cand.momPrTPC,
        li4cand.nTPCClustersHe3,
        li4cand.nSigmaHe3,
        li4cand.nSigmaPr,
        li4cand.massTOFHe3,
        li4cand.massTOFPr,
        li4cand.PIDtrkHe3,
        li4cand.PIDtrkPr,
        li4cand.itsClSizeHe3,
        li4cand.itsClSizePr,
        li4cand.sharedClustersHe3,
        li4cand.sharedClustersPr,
        li4cand.isBkgUS,
        li4cand.isBkgEM);
    } else {
      m_outputMCTable(
        li4cand.recoPtHe3(),
        li4cand.recoEtaHe3(),
        li4cand.recoPhiHe3(),
        li4cand.recoPtPr(),
        li4cand.recoEtaPr(),
        li4cand.recoPhiPr(),
        li4cand.DCAxyHe3,
        li4cand.DCAzHe3,
        li4cand.DCAxyPr,
        li4cand.DCAzPr,
        li4cand.tpcSignalHe3,
        li4cand.momHe3TPC,
        li4cand.tpcSignalPr,
        li4cand.momPrTPC,
        li4cand.nTPCClustersHe3,
        li4cand.nSigmaHe3,
        li4cand.nSigmaPr,
        li4cand.massTOFHe3,
        li4cand.massTOFPr,
        li4cand.PIDtrkHe3,
        li4cand.PIDtrkPr,
        li4cand.itsClSizeHe3,
        li4cand.itsClSizePr,
        li4cand.sharedClustersHe3,
        li4cand.sharedClustersPr,
        li4cand.isBkgUS,
        li4cand.isBkgEM,
        li4cand.momHe3MC,
        li4cand.momPrMC,
        li4cand.l4PtMC,
        li4cand.l4MassMC);
    }
  }

  void fillHistograms(const Lithium4Candidate& l4cand)
  {
    int candSign = l4cand.sign;
    m_hAnalysis.fill(HIST("hHe3Pt"), l4cand.recoPtHe3());
    m_hAnalysis.fill(HIST("hProtonPt"), l4cand.recoPtPr());
    m_hAnalysis.fill(HIST("hLitInvMass"), l4cand.invMass);
    m_hAnalysis.fill(HIST("hDCAxyHe3"), l4cand.DCAxyHe3);
    m_hAnalysis.fill(HIST("hDCAzHe3"), l4cand.DCAzHe3);
    m_hAnalysis.fill(HIST("h2NsigmaHe3TPC"), candSign * l4cand.momHe3TPC, l4cand.nSigmaHe3);
    m_hAnalysis.fill(HIST("h2NsigmaProtonTPC"), candSign * l4cand.momPrTPC, l4cand.nSigmaPr);
    m_hAnalysis.fill(HIST("h2NsigmaProtonTOF"), l4cand.recoPtPr(), l4cand.nSigmaPr);
  }

  // ==================================================================================================================

  void processSameEvent(const soa::Join<aod::Collisions, aod::EvSels>& collisions, const TrackCandidates& tracks, const aod::BCs&)
  {
    for (auto& collision : collisions) {

      m_trackPairs.clear();

      if (!collision.sel8() || std::abs(collision.posZ()) > setting_cutVertex) {
        continue;
      }
      m_hAnalysis.fill(HIST("hNcontributor"), collision.numContrib());
      m_hAnalysis.fill(HIST("hVtxZ"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto TrackTable_thisCollision = tracks.sliceBy(perCol, collIdx);
      TrackTable_thisCollision.bindExternalIndices(&tracks);

      pairTracksSameEvent(TrackTable_thisCollision);

      for (auto& trackPair : m_trackPairs) {

        auto track0 = TrackTable_thisCollision.rawIteratorAt(trackPair.tr0Idx);
        auto track1 = TrackTable_thisCollision.rawIteratorAt(trackPair.tr1Idx);

        Lithium4Candidate li4cand;
        if (!fillCandidateInfo(track0, track1, li4cand, /* mix */ false)) {
          continue;
        }
        fillCandidateBeta(track0, track1, li4cand);
        fillHistograms(li4cand);
        fillTable(li4cand, false);
      }
    }
  }
  PROCESS_SWITCH(lithium4analysis, processSameEvent, "Process Same event", false);

  void processMixedEvent(EventCandidates& /*collisions*/, const TrackCandidates& tracks)
  {
    m_trackPairs.clear();
    pairTracksEventMixing();

    for (auto& trackPair : m_trackPairs) {

      auto track0 = tracks.rawIteratorAt(trackPair.tr0Idx);
      auto track1 = tracks.rawIteratorAt(trackPair.tr1Idx);

      Lithium4Candidate li4cand;
      if (!fillCandidateInfo(track0, track1, li4cand, /* mix */ true)) {
        continue;
      }
      fillCandidateBeta(track0, track1, li4cand);
      fillHistograms(li4cand);
      fillTable(li4cand, false);
    }
  }
  PROCESS_SWITCH(lithium4analysis, processMixedEvent, "Process Mixed event", false);

  void processMC(const soa::Join<aod::Collisions, aod::EvSels>& collisions, const aod::BCs&, const TrackCandidatesMC& tracks, const aod::McParticles& mcParticles)
  {
    std::vector<unsigned int> filledMothers;

    for (auto& collision : collisions) {

      m_trackPairs.clear();

      if (!collision.sel8() || std::abs(collision.posZ()) > setting_cutVertex) {
        continue;
      }

      m_hAnalysis.fill(HIST("hNcontributor"), collision.numContrib());
      m_hAnalysis.fill(HIST("hVtxZ"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto TrackTable_thisCollision = tracks.sliceBy(perColMC, collIdx);
      TrackTable_thisCollision.bindExternalIndices(&tracks);

      pairTracksSameEvent(TrackTable_thisCollision);

      for (auto& trackPair : m_trackPairs) {

        auto track0 = TrackTable_thisCollision.rawIteratorAt(trackPair.tr0Idx);
        auto track1 = TrackTable_thisCollision.rawIteratorAt(trackPair.tr1Idx);

        if (!track0.has_mcParticle() || !track1.has_mcParticle()) {
          continue;
        }

        auto mctrackHe3 = track0.mcParticle();
        auto mctrackPr = track1.mcParticle();

        if (std::abs(mctrackHe3.pdgCode()) != hePDG || std::abs(mctrackPr.pdgCode()) != prPDG) {
          continue;
        }

        for (auto& mothertrack : mctrackHe3.mothers_as<aod::McParticles>()) {
          for (auto& mothertrackPr : mctrackPr.mothers_as<aod::McParticles>()) {

            if (mothertrack != mothertrackPr || std::abs(mothertrack.pdgCode()) != li4PDG || std::abs(mothertrack.y()) > 1) {
              continue;
            }

            Lithium4Candidate li4cand;
            if (!fillCandidateInfo(track0, track1, li4cand, /* mix */ false)) {
              continue;
            }
            fillMcCandidateBeta(track0, track1, li4cand);
            fillMcCandidateInfo(li4cand, mctrackHe3, mctrackPr, mothertrack);
            fillHistograms(li4cand);
            fillTable(li4cand, true);
            filledMothers.push_back(mothertrack.globalIndex());
          }
        }
      }
    }

    for (auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.pdgCode()) != li4PDG || std::abs(mcParticle.y()) > 1 || mcParticle.isPhysicalPrimary() == false) {
        continue;
      }

      if (std::find(filledMothers.begin(), filledMothers.end(), mcParticle.globalIndex()) != filledMothers.end()) {
        continue;
      }

      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      bool daughtHe3(false), daughtPr(false);
      McIter mcHe3, mcPr;
      for (auto kCurrentDaughter : kDaughters) {
        if (std::abs(kCurrentDaughter.pdgCode()) == hePDG) {
          daughtHe3 = true;
          mcHe3 = kCurrentDaughter;
        } else if (std::abs(kCurrentDaughter.pdgCode()) == prPDG) {
          daughtPr = true;
          mcPr = kCurrentDaughter;
        }
      }
      if (daughtHe3 && daughtPr) {
        Lithium4Candidate li4cand;
        fillMcCandidateInfo(li4cand, mcHe3, mcPr, mcParticle);
        fillTable(li4cand, true);
      }
    }
  }
  PROCESS_SWITCH(lithium4analysis, processMC, "Process MC", false);

  void processSameEventPools(const soa::Join<aod::Collisions, aod::EvSels>& collisions, const TrackCandidates& tracks, const aod::AmbiguousTracks& ambiguousTracks, const aod::BCsWithTimestamps& bcs)
  {
    m_goodCollisions.clear();
    m_goodCollisions.resize(collisions.size(), false);
    selectCollisions(collisions);

    m_svPoolCreator.clearPools();
    m_svPoolCreator.fillBC2Coll(collisions, bcs);

    for (auto& track : tracks) {

      if (!selectionTrack(track))
        continue;

      bool selPr = selectionPIDProton(track);
      bool selHe = selectionPIDHe3(track);
      if ((!selPr && !selHe) || (selPr && selHe))
        continue;

      int pdgHypo = selHe ? hePDG : prPDG;

      m_svPoolCreator.appendTrackCand(track, collisions, pdgHypo, ambiguousTracks, bcs);
    }
    auto& svPool = m_svPoolCreator.getSVCandPool(collisions);

    for (auto& svCand : svPool) {
      auto heTrack = tracks.rawIteratorAt(svCand.tr0Idx);
      auto prTrack = tracks.rawIteratorAt(svCand.tr1Idx);
      auto collIdxs = svCand.collBracket;
      Lithium4Candidate li4cand;

      if (!fillCandidateInfo(heTrack, prTrack, li4cand, /* mix */ false))
        continue;
      fillCandidateBeta(heTrack, prTrack, li4cand);
      fillHistograms(li4cand);
      fillTable(li4cand, false);
    }
  }
  PROCESS_SWITCH(lithium4analysis, processSameEventPools, "Process Same event pools", false);

  void processMcPools(const soa::Join<aod::Collisions, aod::EvSels>& collisions, const TrackCandidatesMC& tracks, const aod::AmbiguousTracks& ambiguousTracks, const aod::BCsWithTimestamps& bcs)
  {
    m_goodCollisions.clear();
    m_goodCollisions.resize(collisions.size(), false);
    selectCollisions(collisions);

    m_svPoolCreator.clearPools();
    m_svPoolCreator.fillBC2Coll(collisions, bcs);

    std::vector<unsigned int> filledMothers;

    for (auto& track : tracks) {

      if (!selectionTrack(track))
        continue;

      bool selPr = selectionPIDProton(track);
      bool selHe = selectionPIDHe3(track);
      if ((!selPr && !selHe) || (selPr && selHe))
        continue;

      int pdgHypo = selHe ? hePDG : prPDG;

      m_svPoolCreator.appendTrackCand(track, collisions, pdgHypo, ambiguousTracks, bcs);
    }
    auto& svPool = m_svPoolCreator.getSVCandPool(collisions);

    for (auto& svCand : svPool) {
      auto heTrack = tracks.rawIteratorAt(svCand.tr0Idx);
      auto prTrack = tracks.rawIteratorAt(svCand.tr1Idx);
      auto collIdxs = svCand.collBracket;
      Lithium4Candidate li4cand;

      if (!fillCandidateInfo(heTrack, prTrack, li4cand, /* mix */ false))
        continue;
      fillMcCandidateBeta(heTrack, prTrack, li4cand);
      fillHistograms(li4cand);
      fillTable(li4cand, true);
      // cfr hyperRecoTask.cxx:561
      // filledMothers.push_back(mothertrack.globalIndex());
    }

    // for (auto& mcParticle : mcParticles) {
    //
    //  if (std::abs(mcParticle.pdgCode()) != li4PDG || std::abs(mcParticle.y()) > 1 || mcParticle.isPhysicalPrimary() == false) {
    //    continue;
    //  }
    //
    //  if (std::find(filledMothers.begin(), filledMothers.end(), mcParticle.globalIndex()) != filledMothers.end()) {
    //    continue;
    //  }
    //
    //  auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
    //  bool daughtHe3(false), daughtPr(false);
    //  McIter mcHe3, mcPr;
    //  for (auto kCurrentDaughter : kDaughters) {
    //    if (std::abs(kCurrentDaughter.pdgCode()) == hePDG) {
    //      daughtHe3 = true;
    //      mcHe3 = kCurrentDaughter;
    //    } else if (std::abs(kCurrentDaughter.pdgCode()) == prPDG) {
    //      daughtPr = true;
    //      mcPr = kCurrentDaughter;
    //    }
    //  }
    //  if (daughtHe3 && daughtPr) {
    //    Lithium4Candidate li4cand;
    //    fillMcCandidateInfo(li4cand, mcHe3, mcPr, mcParticle);
    //    fillTable(li4cand, true);
    //  }
    //}
  }
  PROCESS_SWITCH(lithium4analysis, processMcPools, "Process MC pools", false);
};

WorkflowSpec defineDataProcessing(const ConfigContext& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lithium4analysis>(cfgc, TaskName{"lithium4analysis"})};
}
