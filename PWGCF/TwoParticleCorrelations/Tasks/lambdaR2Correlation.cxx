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

/// \file lambdaCorrelationAnalysis.cxx
/// \brief R2 correlation of Lambda baryons.
/// \author Yash Patley <yash.patley@cern.ch>

#include <TLorentzVector.h>
#include <TRandom.h>

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace o2::aod
{
namespace lambdacollision
{
DECLARE_SOA_COLUMN(Cent, cent, float);
}
DECLARE_SOA_TABLE(LambdaCollisions, "AOD", "LAMBDACOLS", o2::soa::Index<>,
                  lambdacollision::Cent,
                  aod::collision::PosX,
                  aod::collision::PosY,
                  aod::collision::PosZ);
using LambdaCollision = LambdaCollisions::iterator;

namespace lambdatrack
{
DECLARE_SOA_INDEX_COLUMN(LambdaCollision, lambdaCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(PosTrackId, postrackid, int64_t);
DECLARE_SOA_COLUMN(NegTrackId, negtrackid, int64_t);
DECLARE_SOA_COLUMN(V0Type, v0type, int8_t);
DECLARE_SOA_COLUMN(MassWindow, masswindow, int8_t);
} // namespace lambdatrack
DECLARE_SOA_TABLE(LambdaTracks, "AOD", "LAMBDATRACKS", o2::soa::Index<>,
                  lambdatrack::LambdaCollisionId,
                  lambdatrack::Pt,
                  lambdatrack::Rap,
                  lambdatrack::Phi,
                  lambdatrack::Mass,
                  lambdatrack::PosTrackId,
                  lambdatrack::NegTrackId,
                  lambdatrack::V0Type,
                  lambdatrack::MassWindow);
using LambdaTrack = LambdaTracks::iterator;
} // namespace o2::aod

enum PidType {
  kPion = 0,
  kProton
};

enum ParticleType {
  kLambda = 0,
  kAntiLambda
};

enum ParticlePairType {
  kLambdaAntiLambda = 0,
  kLambdaLambda,
  kAntiLambdaAntiLambda
};

enum MassWindowType {
  kCentralWindow = 0,
  kLeftWindow,
  kRightWindow
};

struct lambdaCorrTableProducer {

  Produces<aod::LambdaCollisions> lambdaCollisionTable;
  Produces<aod::LambdaTracks> lambdaTrackTable;

  // Collisions
  Configurable<float> cfg_z_vtx{"cfg_z_vtx", 10.0, "z vertex cut"};
  Configurable<bool> cfg_trigger_tvx_sel{"cfg_trigger_tvx_sel", false, "Trigger Time and Vertex Selection"};
  Configurable<bool> cfg_tf_border{"cfg_tf_border", false, "Timeframe Border Selection"};
  Configurable<bool> cfg_noitsro_border{"cfg_noitsro_border", false, "No ITSRO Border Cut"};
  Configurable<bool> cfg_sel8_sel{"cfg_sel8_sel", true, "Sel8 (T0A + T0C) Selection"};
  Configurable<bool> cfg_itstpc_vtx{"cfg_itstpc_vtx", false, "ITS+TPC Vertex Selection"};
  Configurable<bool> cfg_pileup_reject{"cfg_pileup_reject", false, "Pileup rejection"};
  Configurable<bool> cfg_zvtx_time_diff{"cfg_zvtx_time_diff", false, "z-vtx time diff selection"};

  // Tracks
  Configurable<float> cfg_eta_cut{"cfg_eta_cut", 0.8, "Pseudorapidity cut"};
  Configurable<int> cfg_min_crossed_rows{"cfg_min_crossed_rows", 70, "min crossed rows"};
  Configurable<double> cfg_tpc_nsigma{"cfg_tpc_nsigma", 3.0, "TPC NSigma Selection Cut"};
  Configurable<double> cfg_tof_nsigma{"cfg_tof_nsigma", 3.0, "TOF NSigma Selection Cut"};

  // V0s
  Configurable<double> cfg_min_dca_V0_daughters{"cfg_min_dca_V0_daughters", 1.0, "min DCA between V0 daughters"};
  Configurable<double> cfg_min_dca_pos_to_PV{"cfg_min_dca_pos_to_PV", 0.1, "Minimum V0 Positive Track DCAr cut to PV"};
  Configurable<double> cfg_min_dca_neg_to_PV{"cfg_min_dca_neg_to_PV", 0.1, "Minimum V0 Negative Track DCAr cut to PV"};
  Configurable<double> cfg_min_dca_V0_to_PV{"cfg_min_dca_V0_to_PV", 0.6, "Minimum DCA V0 to PV"};
  Configurable<double> cfg_min_V0_radius{"cfg_min_V0_radius", 0.0, "Minimum V0 radius from PV"};
  Configurable<double> cfg_max_V0_radius{"cfg_max_V0_radius", 50.0, "Maximum V0 radius from PV"};
  Configurable<double> cfg_min_ctau{"cfg_min_ctau", 0.0, "Minimum ctau"};
  Configurable<double> cfg_max_ctau{"cfg_max_ctau", 50.0, "Maximum ctau"};
  Configurable<double> cfg_min_V0_cosPA{"cfg_min_V0_cosPA", 0.995, "Minimum V0 CosPA to PV"};
  Configurable<double> cfg_lambda_mass_window{"cfg_lambda_mass_window", 0.01, "Mass Window to select Lambda"};
  Configurable<double> cfg_kshort_rej{"cfg_kshort_rej", 0.005, "Reject K0Short Candidates"};

  // V0s kinmatic acceptance
  Configurable<float> cfg_v0_pt_min{"cfg_v0_pt_min", 0.5, "Minimum V0 pT"};
  Configurable<float> cfg_v0_pt_max{"cfg_v0_pt_max", 2.5, "Minimum V0 pT"};
  Configurable<float> cfg_v0_rap_max{"cfg_v0_rap_max", 0.8, "|rap| cut"};

  // bool eta/rapidity
  Configurable<bool> cfg_do_eta_analysis{"cfg_do_eta_analysis", false, "Eta Analysis"};

  // lambda mass windows
  Configurable<float> cfg_lambda_mass_min{"cfg_lambda_mass_min", 1.11, "Minimum Central Window"};
  Configurable<float> cfg_lambda_mass_max{"cfg_lambda_mass_max", 1.12, "Maximum Central Window"};
  Configurable<float> cfg_lambda_left_min{"cfg_lambda_left_min", 1.08, "Minimum Left Window"};
  Configurable<float> cfg_lambda_left_max{"cfg_lambda_left_max", 1.10, "Maximum Left Window"};
  Configurable<float> cfg_lambda_right_min{"cfg_lambda_right_min", 1.13, "Minimum Right Window"};
  Configurable<float> cfg_lambda_right_max{"cfg_lambda_right_max", 1.15, "Maximum Right Window"};

  // global variable declaration
  std::map<MassWindowType, float> mass_map_min;
  std::map<MassWindowType, float> mass_map_max;

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // global variable
    mass_map_min = {{kCentralWindow, cfg_lambda_mass_min}, {kLeftWindow, cfg_lambda_left_min}, {kRightWindow, cfg_lambda_right_min}};
    mass_map_max = {{kCentralWindow, cfg_lambda_mass_max}, {kLeftWindow, cfg_lambda_left_max}, {kRightWindow, cfg_lambda_right_max}};

    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");

    const AxisSpec axisV0Mass(100, 1.06, 1.16, "Inv Mass (GeV/#it{c}^{2})");
    const AxisSpec axisV0Pt(200, 0., 5., "p_{T} (GeV/#it{c})");
    const AxisSpec axisV0Rap(16, -1., 1., "rap");
    const AxisSpec axisV0Phi(36, 0., 2. * TMath::Pi(), "#phi (rad)");

    const AxisSpec axisRadius(500, 0, 100, "r(cm)");
    const AxisSpec axisCosPA(120, 0.97, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaV0PV(200, 0, 2., "dca (cm)");
    const AxisSpec axisDcaProngPV(200, 0, 20., "dca (cm)");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daug DCA (cm^{2})");
    const AxisSpec axisCTau(500, 0, 100, "c#tau (cm/#it{c})");

    const AxisSpec axisMomPID(80, 0, 4, "p (GeV/#it{c})");
    const AxisSpec axisNsigma(401, -10.025, 10.025, {"n#sigma"});
    const AxisSpec axisdEdx(360, 20, 200, "#frac{dE}{dx}");

    // Create Histograms.
    // QA
    histos.add("QA_Checks/h1d_lambda_mass", "M_{#Lambda}", kTH1F, {axisV0Mass});

    // QA Lambda
    histos.add("QA_Sel_Lambda/h1d_V0_inv_mass", "V_{0} mass", kTH1F, {axisV0Mass});
    histos.add("QA_Sel_Lambda/h1d_V0_pt", "V_{0} p_{T}", kTH1F, {axisV0Pt});
    histos.add("QA_Sel_Lambda/h1d_V0_eta", "#eta-distribution", kTH1F, {axisV0Rap});
    histos.add("QA_Sel_Lambda/h1d_V0_rap", "y-distribution", kTH1F, {axisV0Rap});
    histos.add("QA_Sel_Lambda/h1d_V0_phi", "#phi-distribution", kTH1F, {axisV0Phi});

    histos.add("QA_Sel_Lambda/h1d_dca_V0_daughters", "DCA between V0 daughters", kTH1F, {axisDcaDau});
    histos.add("QA_Sel_Lambda/h1d_dca_pos_to_PV", "DCA positive prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA_Sel_Lambda/h1d_dca_neg_to_PV", "DCA negative prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA_Sel_Lambda/h1d_dca_V0_to_PV", "DCA V0 to PV", kTH1F, {axisDcaV0PV});
    histos.add("QA_Sel_Lambda/h1d_V0_cospa", "cos(#theta_{PA})", kTH1F, {axisCosPA});
    histos.add("QA_Sel_Lambda/h1d_V0_radius", "V_{0} Decay Radius in XY plane", kTH1F, {axisRadius});
    histos.add("QA_Sel_Lambda/h1d_V0_ctau", "V_{0} c#tau", kTH1F, {axisCTau});

    histos.add("QA_Sel_Lambda/h2d_pr_dEdx_vs_p", "TPC Signal Proton", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA_Sel_Lambda/h2d_pi_dEdx_vs_p", "TPC Signal Pion", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA_Sel_Lambda/h2d_pr_nsigma_tpc", "TPC Proton", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_pi_nsigma_tpc", "TPC Pion", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_pr_nsigma_tof", "TOF Proton", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_pi_nsigma_tof", "TOF Pion", kTH2F, {axisMomPID, axisNsigma});

    // QA Anti-Lambda
    histos.addClone("QA_Sel_Lambda/", "QA_Sel_AntiLambda/");
  }

  template <typename C>
  bool selCol(C const& col)
  {

    if (fabs(col.posZ()) > cfg_z_vtx) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kIsTriggerTVX) && cfg_trigger_tvx_sel) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kNoTimeFrameBorder) && cfg_tf_border) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kNoITSROFrameBorder) && cfg_noitsro_border) {
      return false;
    }

    if (!col.sel8() && cfg_sel8_sel) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kIsVertexITSTPC) && cfg_itstpc_vtx) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kNoSameBunchPileup) && cfg_pileup_reject) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && cfg_zvtx_time_diff) {
      return false;
    }

    return true;
  }

  template <typename V, typename T>
  bool topologicalCutsV0(V const& v0, T const&)
  {

    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    if (fabs(postrack.eta()) > cfg_eta_cut) {
      return false;
    }

    if (fabs(negtrack.eta()) > cfg_eta_cut) {
      return false;
    }

    if (postrack.tpcNClsCrossedRows() < cfg_min_crossed_rows) {
      return false;
    }

    if (negtrack.tpcNClsCrossedRows() < cfg_min_crossed_rows) {
      return false;
    }

    if (v0.dcaV0daughters() > cfg_min_dca_V0_daughters) {
      return false;
    }

    if (fabs(v0.dcapostopv()) < cfg_min_dca_pos_to_PV) {
      return false;
    }

    if (fabs(v0.dcanegtopv()) < cfg_min_dca_neg_to_PV) {
      return false;
    }

    if (v0.dcav0topv() > cfg_min_dca_V0_to_PV) {
      return false;
    }

    if ((v0.v0radius() > cfg_max_V0_radius) || (v0.v0radius() < cfg_min_V0_radius)) {
      return false;
    }

    if (v0.v0cosPA() < cfg_min_V0_cosPA) {
      return false;
    }

    return true;
  }

  template <PidType pid, typename T>
  bool selPIDTrack(T const& track)
  {

    bool selTPCv0type = false, selTOFv0type = false;
    float tpcNSigma = 0., tofNSigma = 0.;

    switch (pid) {

      case kPion:
        tpcNSigma = track.tpcNSigmaPi();
        tofNSigma = track.tofNSigmaPi();
        break;

      case kProton:
        tpcNSigma = track.tpcNSigmaPr();
        tofNSigma = track.tofNSigmaPr();
        break;
    }

    if (track.hasTOF()) {
      if (fabs(tofNSigma) < cfg_tof_nsigma) {
        selTOFv0type = true;
      }
      if (fabs(tpcNSigma) < cfg_tpc_nsigma) {
        selTPCv0type = true;
      }
    } else {
      selTOFv0type = true;
      if (fabs(tpcNSigma) < cfg_tpc_nsigma) {
        selTPCv0type = true;
      }
    }

    if (selTPCv0type && selTOFv0type) {
      return true;
    }

    return false;
  }

  template <ParticleType part, typename C, typename V, typename T>
  void fillQALambda(C const& col, V const& v0, T const&)
  {

    static constexpr std::string_view sub_dir[] = {"QA_Sel_Lambda/", "QA_Sel_AntiLambda/"};

    // ctau
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;

    // daugthers
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_V0_daughters"), v0.dcaV0daughters());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_pos_to_PV"), fabs(v0.dcapostopv()));
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_neg_to_PV"), fabs(v0.dcanegtopv()));
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_V0_to_PV"), fabs(v0.dcav0topv()));
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_cospa"), v0.v0cosPA());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_radius"), v0.v0radius());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_ctau"), ctau);
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_pt"), v0.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_eta"), v0.eta());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_rap"), v0.yLambda());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_phi"), v0.phi());
    if constexpr (part == kLambda) {
      histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_inv_mass"), v0.mLambda());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_pr_dEdx_vs_p"), postrack.tpcInnerParam(), postrack.tpcSignal());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_pi_dEdx_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcSignal());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_pr_nsigma_tpc"), postrack.tpcInnerParam(), postrack.tpcNSigmaPr());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_pi_nsigma_tpc"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPi());
      if (postrack.hasTOF()) {
        histos.fill(HIST(sub_dir[part]) + HIST("h2d_pr_nsigma_tof"), postrack.tofExpMom(), postrack.tofNSigmaPr());
      }
      if (negtrack.hasTOF()) {
        histos.fill(HIST(sub_dir[part]) + HIST("h2d_pi_nsigma_tof"), negtrack.tofExpMom(), negtrack.tofNSigmaPi());
      }
    } else {
      histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_inv_mass"), v0.mAntiLambda());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_pi_dEdx_vs_p"), postrack.tpcInnerParam(), postrack.tpcSignal());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_pr_dEdx_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcSignal());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_pi_nsigma_tpc"), postrack.tpcInnerParam(), postrack.tpcNSigmaPi());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_pr_nsigma_tpc"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPr());
      if (postrack.hasTOF()) {
        histos.fill(HIST(sub_dir[part]) + HIST("h2d_pi_nsigma_tof"), postrack.tofExpMom(), postrack.tofNSigmaPi());
      }
      if (negtrack.hasTOF()) {
        histos.fill(HIST(sub_dir[part]) + HIST("h2d_pr_nsigma_tof"), negtrack.tofExpMom(), negtrack.tofNSigmaPr());
      }
    }
  }

  template <ParticleType v0part, PidType pos_prong, PidType neg_prong, MassWindowType masswin, typename C, typename V, typename T>
  void selV0Particle(C const& collision, V const& v0track, T const& tracks)
  {

    // initialize variables
    auto postrack = v0track.template posTrack_as<T>();
    auto negtrack = v0track.template negTrack_as<T>();
    float mass = 0., rap = 0.;

    if (v0part == kLambda) {
      mass = v0track.mLambda();
    }
    if (v0part == kAntiLambda) {
      mass = v0track.mAntiLambda();
    }

    if (cfg_do_eta_analysis) {
      rap = v0track.eta();
    } else {
      rap = v0track.yLambda();
    }

    // apply mass window selection [global]
    if ((fabs(mass - MassLambda0) >= cfg_lambda_mass_window) || (fabs(v0track.mK0Short() - MassK0Short) <= cfg_kshort_rej)) {
      return;
    }

    // apply daughter particle id
    if (!selPIDTrack<pos_prong>(postrack) || !selPIDTrack<neg_prong>(negtrack)) {
      return;
    }

    // apply kinematic acceptance
    if (v0track.pt() <= cfg_v0_pt_min || v0track.pt() >= cfg_v0_pt_max) {
      return;
    }

    if (cfg_do_eta_analysis && (fabs(v0track.eta()) >= cfg_v0_rap_max)) {
      return;
    } else if (fabs(v0track.yLambda()) >= cfg_v0_rap_max) {
      return;
    }

    if (v0part == kLambda) {
      histos.fill(HIST("QA_Checks/h1d_lambda_mass"), v0track.mLambda());
    }

    if (v0part == kAntiLambda) {
      histos.fill(HIST("QA_Checks/h1d_lambda_mass"), v0track.mAntiLambda());
    }

    // apply further mass window selection [central, left, right]
    if (mass > mass_map_min[masswin] && mass < mass_map_max[masswin]) {
      if (masswin == kCentralWindow) {
        fillQALambda<v0part>(collision, v0track, tracks);
      }
      lambdaTrackTable(lambdaCollisionTable.lastIndex(), v0track.pt(), rap, v0track.phi(), mass, postrack.index(), negtrack.index(), v0part, masswin);
    }
  }

  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;

  void process(Collisions::iterator const& collision, aod::V0Datas const& V0s, Tracks const& tracks)
  {

    // select collision
    if (!selCol(collision)) {
      return;
    }

    lambdaCollisionTable(collision.centFT0M(), collision.posX(), collision.posY(), collision.posZ());

    for (auto const& v0 : V0s) {

      // apply topological cuts on v0 candidates
      if (!topologicalCutsV0(v0, tracks)) {
        continue;
      }

      selV0Particle<kLambda, kProton, kPion, kCentralWindow>(collision, v0, tracks);
      selV0Particle<kLambda, kProton, kPion, kLeftWindow>(collision, v0, tracks);
      selV0Particle<kLambda, kProton, kPion, kRightWindow>(collision, v0, tracks);
      selV0Particle<kAntiLambda, kPion, kProton, kCentralWindow>(collision, v0, tracks);
      selV0Particle<kAntiLambda, kPion, kProton, kLeftWindow>(collision, v0, tracks);
      selV0Particle<kAntiLambda, kPion, kProton, kRightWindow>(collision, v0, tracks);
    }
  }
};

struct lambdaCorrelationAnalysis {

  // Global Configurables
  Configurable<int> cfg_nRapBins{"cfg_nRapBins", 16, "N Rapidity Bins"};
  Configurable<float> cfg_Rap_Min{"cfg_Rap_Min", -0.8, "Minimum Rapidity"};
  Configurable<float> cfg_Rap_Max{"cfg_Rap_Max", 0.8, "Maximum Rapidity"};
  Configurable<int> cfg_nPhiBins{"cfg_nPhiBins", 64, "N Phi Bins"};
  Configurable<float> cfg_Phi_Min{"cfg_Phi_Min", 0, "Minimum Phi"};
  Configurable<float> cfg_Phi_Max{"cfg_Phi_Max", 2 * TMath::Pi(), "Maximum Phi"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Initialize global variables
  float nrapbins = 0.;
  float kminrap = 0.;
  float kmaxrap = 0.;
  float nphibins = 0.;
  float kminphi = 0.;
  float kmaxphi = 0.;
  float rapbinwidth = 0.;
  float phibinwidth = 0.;

  void init(InitContext const&)
  {
    nrapbins = static_cast<float>(cfg_nRapBins);
    kminrap = static_cast<float>(cfg_Rap_Min);
    kmaxrap = static_cast<float>(cfg_Rap_Max);
    nphibins = static_cast<float>(cfg_nPhiBins);
    kminphi = static_cast<float>(cfg_Phi_Min);
    kmaxphi = static_cast<float>(cfg_Phi_Max);

    rapbinwidth = (kmaxrap - kminrap) / nrapbins;
    phibinwidth = (kmaxphi - kminphi) / nphibins;

    int knrapphibins = static_cast<int>(cfg_nRapBins) * static_cast<int>(cfg_nPhiBins);
    float kminrapphi = 0.;
    float kmaxrapphi = knrapphibins;

    const AxisSpec axisPosZ(220, -11, 11, "V_{z} (cm)");
    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisMass(100, 1.06, 1.16, "Inv Mass (GeV/#it{c}^{2})");
    const AxisSpec axisRap(cfg_nRapBins, cfg_Rap_Min, cfg_Rap_Max, "rap");
    const AxisSpec axisPhi(cfg_nPhiBins, cfg_Phi_Min, cfg_Phi_Max, "#phi (rad)");
    const AxisSpec axisRapPhi(knrapphibins, kminrapphi, kmaxrapphi, "rap #phi");

    // Create Histograms.
    // Event
    histos.add("Event/h1d_posz", "V_{Z} Distribution", kTH1F, {axisPosZ});
    histos.add("Event/h1d_ft0m_mult_percentile", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/h1d_lambda_tot_mult", "#Lambda/#bar{#Lambda} - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/h1d_lambda_multiplicity", "#Lambda - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/h1d_antilambda_multiplicity", "#bar{#Lambda} - Multiplicity", kTH1I, {axisMult});

    // Lambda
    histos.add("Lambda/h1d_inv_mass", "M_{p#pi}", kTH1F, {axisMass});

    // Anti-Lambda
    histos.addClone("Lambda/", "AntiLambda/");

    // single and two particle densities
    histos.add("Lambda_Mass/h2d_n1_LaP", "#rho_{1}^{#Lambda}", kTH2D, {axisRap, axisPhi});
    histos.add("Lambda_Mass/h2d_n1_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2D, {axisRap, axisPhi});
    histos.add("Lambda_Mass/h2d_n2_LaP_LaM", "#rho_{2}^{#Lambda - #bar{#Lambda}}", kTH2D, {axisRapPhi, axisRapPhi});
    histos.add("Lambda_Mass/h2d_n2_LaP_LaP", "#rho_{2}^{#Lambda - #Lambda}", kTH2D, {axisRapPhi, axisRapPhi});
    histos.add("Lambda_Mass/h2d_n2_LaM_LaM", "#rho_{2}^{#bar{#Lambda} - #bar{#Lambda}}", kTH2D, {axisRapPhi, axisRapPhi});

    histos.addClone("Lambda_Mass/", "Lambda_Right/");
    histos.addClone("Lambda_Mass/", "Lambda_Left/");
  }

  template <ParticlePairType part_pair, MassWindowType mass_win, typename U>
  void fillPairHistos(U& p1, U& p2)
  {

    static constexpr std::string_view sub_dir_type[] = {"Lambda_Mass/", "Lambda_Left/", "Lambda_Right/"};
    static constexpr std::string_view sub_dir_hist[] = {"h2d_n2_LaP_LaM", "h2d_n2_LaP_LaP", "h2d_n2_LaM_LaM"};

    int rapbin1 = static_cast<int>((p1.rap() - kminrap) / rapbinwidth);
    int rapbin2 = static_cast<int>((p2.rap() - kminrap) / rapbinwidth);

    int phibin1 = static_cast<int>(p1.phi() / phibinwidth);
    int phibin2 = static_cast<int>(p2.phi() / phibinwidth);

    if (rapbin1 >= 0 && rapbin2 >= 0 && phibin1 >= 0 && phibin2 >= 0 && rapbin1 < nrapbins && rapbin2 < nrapbins && phibin1 < nphibins && phibin2 < nphibins) {

      int rapphix = rapbin1 * nphibins + phibin1;
      int rapphiy = rapbin2 * nphibins + phibin2;

      histos.fill(HIST(sub_dir_type[mass_win]) + HIST(sub_dir_hist[part_pair]), rapphix + 0.5, rapphiy + 0.5);
    }
  }

  template <ParticleType part, MassWindowType masswin, typename T>
  void analyzeSingles(T const& tracks)
  {

    static constexpr std::string_view sub_dir_part[] = {"Lambda/", "AntiLambda/"};
    static constexpr std::string_view sub_dir_mass_win[] = {"Lambda_Mass/", "Lambda_Left/", "Lambda_Right/"};
    static constexpr std::string_view sub_dir_hist[] = {"h2d_n1_LaP", "h2d_n1_LaM"};

    for (auto const& track : tracks) {
      histos.fill(HIST(sub_dir_part[part]) + HIST("h1d_inv_mass"), track.mass());
      histos.fill(HIST(sub_dir_mass_win[masswin]) + HIST(sub_dir_hist[part]), track.rap(), track.phi());
    }
  }

  template <ParticlePairType partpair, MassWindowType masswin, bool samelambda, typename T>
  void analyzePairs(T const& trks_1, T const& trks_2)
  {

    for (auto const& trk_1 : trks_1) {
      for (auto const& trk_2 : trks_2) {
        if (samelambda && (trk_1.index() == trk_2.index()) && (trk_1.postrackid() == trk_2.postrackid()) && (trk_1.negtrackid() == trk_2.negtrackid())) {
          continue;
        }
        fillPairHistos<partpair, masswin>(trk_1, trk_2);
      }
    }
  }

  using Lambda_Collisions = aod::LambdaCollisions;
  using Lambda_Tracks = aod::LambdaTracks;

  SliceCache cache;

  Partition<Lambda_Tracks> part_lambda_tracks = (aod::lambdatrack::v0type == (int8_t)kLambda && aod::lambdatrack::masswindow == (int8_t)kCentralWindow);
  Partition<Lambda_Tracks> part_anti_lambda_tracks = (aod::lambdatrack::v0type == (int8_t)kAntiLambda && aod::lambdatrack::masswindow == (int8_t)kCentralWindow);

  Partition<Lambda_Tracks> part_lambda_tracks_left_masswin = (aod::lambdatrack::v0type == (int8_t)kLambda && aod::lambdatrack::masswindow == (int8_t)kLeftWindow);
  Partition<Lambda_Tracks> part_anti_lambda_tracks_left_masswin = (aod::lambdatrack::v0type == (int8_t)kAntiLambda && aod::lambdatrack::masswindow == (int8_t)kLeftWindow);

  Partition<Lambda_Tracks> part_lambda_tracks_right_masswin = (aod::lambdatrack::v0type == (int8_t)kLambda && aod::lambdatrack::masswindow == (int8_t)kRightWindow);
  Partition<Lambda_Tracks> part_anti_lambda_tracks_right_masswin = (aod::lambdatrack::v0type == (int8_t)kAntiLambda && aod::lambdatrack::masswindow == (int8_t)kRightWindow);

  void process(Lambda_Collisions::iterator const& collision, Lambda_Tracks const&)
  {

    histos.fill(HIST("Event/h1d_posz"), collision.posZ());
    histos.fill(HIST("Event/h1d_ft0m_mult_percentile"), collision.cent());

    auto lambda_tracks = part_lambda_tracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto anti_lambda_tracks = part_anti_lambda_tracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);

    auto lambda_tracks_left = part_lambda_tracks_left_masswin->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto anti_lambda_tracks_left = part_anti_lambda_tracks_left_masswin->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);

    auto lambda_tracks_right = part_lambda_tracks_right_masswin->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto anti_lambda_tracks_right = part_anti_lambda_tracks_right_masswin->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);

    analyzeSingles<kLambda, kCentralWindow>(lambda_tracks);
    analyzeSingles<kLambda, kLeftWindow>(lambda_tracks_left);
    analyzeSingles<kLambda, kRightWindow>(lambda_tracks_right);

    analyzeSingles<kAntiLambda, kCentralWindow>(anti_lambda_tracks);
    analyzeSingles<kAntiLambda, kLeftWindow>(anti_lambda_tracks_left);
    analyzeSingles<kAntiLambda, kRightWindow>(anti_lambda_tracks_right);

    analyzePairs<kLambdaAntiLambda, kCentralWindow, false>(lambda_tracks, anti_lambda_tracks);
    analyzePairs<kLambdaAntiLambda, kLeftWindow, false>(lambda_tracks_left, anti_lambda_tracks_left);
    analyzePairs<kLambdaAntiLambda, kRightWindow, false>(lambda_tracks_right, anti_lambda_tracks_right);

    analyzePairs<kLambdaLambda, kCentralWindow, true>(lambda_tracks, lambda_tracks);
    analyzePairs<kLambdaLambda, kLeftWindow, true>(lambda_tracks_left, lambda_tracks_left);
    analyzePairs<kLambdaLambda, kRightWindow, true>(lambda_tracks_right, lambda_tracks_right);

    analyzePairs<kAntiLambdaAntiLambda, kCentralWindow, true>(anti_lambda_tracks, anti_lambda_tracks);
    analyzePairs<kAntiLambdaAntiLambda, kLeftWindow, true>(anti_lambda_tracks_left, anti_lambda_tracks_left);
    analyzePairs<kAntiLambdaAntiLambda, kRightWindow, true>(anti_lambda_tracks_right, anti_lambda_tracks_right);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdaCorrTableProducer>(cfgc),
    adaptAnalysisTask<lambdaCorrelationAnalysis>(cfgc)};
}
