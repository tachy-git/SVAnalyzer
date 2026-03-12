// -*- C++ -*-
// Package:    MyAnalysis/staMuon
// Class:      staMuon

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "FWCore/Utilities/interface/ESInputTag.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <vector>
#include <cmath>
#include <array>
#include <limits>
#include <string>

static constexpr float kInvalid = -999.f;
static constexpr int   kNFits   = 1;

class staMuon : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit staMuon(const edm::ParameterSet&);
  ~staMuon() override {}

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  // == Tokens ================================================================
  edm::EDGetTokenT<std::vector<pat::Muon>>         muonToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>>         disMuonToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>>      pvToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  TTree* tree_;

  // == event-level: gen J/psi (pdgId == 443) =================================
  std::vector<float> gen_jpsi_pt_;
  std::vector<float> gen_jpsi_eta_;
  std::vector<float> gen_jpsi_phi_;

  // == event-level: gen dimuon system ========================================
  std::vector<float> gen_dimu_pt_;
  std::vector<float> gen_dimu_eta_;
  std::vector<float> gen_dimu_phi_;
  std::vector<float> gen_dimu_mass_;
  std::vector<float> gen_dimu_dR_;

  // == event-level: gen muons ================================================
  std::vector<float> gen_mu_pt_;
  std::vector<float> gen_mu_eta_;
  std::vector<float> gen_mu_phi_;
  std::vector<int>   gen_mu_pid_;
  std::vector<int>   gen_mu_charge_;
  std::vector<float> gen_mu_pt_m_;
  std::vector<float> gen_mu_pt_p_;
  std::vector<float> gen_mu_px_;
  std::vector<float> gen_mu_py_;
  std::vector<float> gen_mu_pz_;

  // == pair-level scalars ====================================================
  std::vector<int>   coll_;
  std::vector<float> invm_;
  std::vector<float> invm_bef_;

  // DCA between the two muon tracks — one scalar per OS pair
  std::vector<float> dca2d_;
  std::vector<float> dca3d_;
  std::vector<float> dca3d_cp_rho_;

  // == gen-reco muon matching ================================================
  std::vector<std::vector<float>> genReco_dEta_;
  std::vector<std::vector<float>> genReco_dPhi_;
  std::vector<std::vector<float>> genReco_dR_;
  std::vector<std::vector<float>> genReco_ptRatio_;
  std::vector<std::vector<float>> genReco_pxRatio_;
  std::vector<std::vector<float>> genReco_pyRatio_;

  // == muon-level: original PAT muon =========================================
  std::vector<std::vector<float>> pt_mu_;
  std::vector<std::vector<float>> eta_mu_;
  std::vector<std::vector<int>>   charge_mu_;
  std::vector<std::vector<float>> ip_mu_;
  std::vector<std::vector<float>> time_mu_;
  std::vector<std::vector<int>>   cscHit_mu_;
  std::vector<std::vector<int>>   dtHit_mu_;
  std::vector<std::vector<int>>   dir_mu_;
  std::vector<std::vector<float>> normChi_mu_;
  std::vector<std::vector<float>> px_mu_;
  std::vector<std::vector<float>> py_mu_;
  std::vector<std::vector<float>> pz_mu_;
  std::vector<std::vector<float>> phi_mu_;

  // == reco muon pT split by charge ==========================================
  std::vector<float> reco_mu_pt_m_;
  std::vector<float> reco_mu_pt_p_;

  // == per-fit: pair-level ===================================================
  std::array<std::vector<int>,   kNFits> isValid_vtx_;
  std::array<std::vector<float>, kNFits> invm_aft_;
  std::array<std::vector<float>, kNFits> Lxy_;
  std::array<std::vector<float>, kNFits> LxyErr_;
  std::array<std::vector<float>, kNFits> vx_;
  std::array<std::vector<float>, kNFits> vy_;
  std::array<std::vector<float>, kNFits> normChi_vtx_;
  std::array<std::vector<float>, kNFits> prob_vtx_;

  // == per-fit: muon-level kinematics ========================================
  std::array<std::vector<std::vector<float>>, kNFits> pt_refit_;
  std::array<std::vector<std::vector<float>>, kNFits> eta_refit_;
  std::array<std::vector<std::vector<float>>, kNFits> relPtUnc_refit_;

  // == helpers ===============================================================
  static bool passSASelection(const pat::Muon& mu) {
    if (!mu.isStandAloneMuon())    return false;
    if ( mu.isTrackerMuon())       return false;
    if (mu.pt() <= 5.)             return false;
    if (std::abs(mu.eta()) >= 2.4) return false;
    return true;
  }

  static float myDeltaR(float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = phi1 - phi2;
    while (dphi >  M_PI) dphi -= 2.f * M_PI;
    while (dphi < -M_PI) dphi += 2.f * M_PI;
    return std::sqrt(deta * deta + dphi * dphi);
  }

  static float deltaPhi(float phi1, float phi2) {
    float dphi = phi1 - phi2;
    while (dphi >  M_PI) dphi -= 2.f * M_PI;
    while (dphi < -M_PI) dphi += 2.f * M_PI;
    return dphi;
  }

  void clearVectors();

  struct MuonVars {
    float pt, eta, phi, ip, time, normChi;
    float px, py, pz;
    int   charge, cscHit, dtHit, dir;

    static MuonVars from(const pat::Muon& mu) {
      MuonVars v;
      v.pt     = mu.pt();
      v.eta    = mu.eta();
      v.phi    = mu.phi();
      v.charge = mu.charge();
      v.ip     = mu.dB();
      v.time   = mu.time().timeAtIpInOut;
      v.px     = mu.px();
      v.py     = mu.py();
      v.pz     = mu.pz();
      const auto& hp = mu.bestTrack()->hitPattern();
      v.cscHit  = hp.numberOfValidMuonCSCHits();
      v.dtHit   = hp.numberOfValidMuonDTHits();
      v.dir     = static_cast<int>(mu.bestTrack()->seedDirection());
      v.normChi = mu.bestTrack()->normalizedChi2();
      return v;
    }
  };

  struct RefitVars {
    float pt, eta, phi, relPtUnc;

    static RefitVars from(const reco::TransientTrack& trk,
                          const GlobalVector& gMom) {
      const reco::Track& t = trk.track();
      RefitVars v;
      v.pt       = t.pt();
      v.eta      = t.eta();
      v.phi      = std::atan2(gMom.y(), gMom.x());
      v.relPtUnc = (t.pt() > 0.f) ? t.ptError() / t.pt() : -1.f;
      return v;
    }
    static RefitVars invalid() { return {kInvalid, kInvalid, kInvalid, kInvalid}; }
  };

  void pushMuonPair(const MuonVars& v1, const MuonVars& v2);
};


// == Constructor ===============================================================
staMuon::staMuon(const edm::ParameterSet& iConfig) {
  usesResource("TFileService");
  muonToken_    = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
  disMuonToken_ = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("disMuons"));
  pvToken_      = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"));
  genToken_     = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));
  ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(
                       edm::ESInputTag("", "TransientTrackBuilder"));
}


// == beginJob ==================================================================
void staMuon::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "staMuon tree");

  // event-level gen
  tree_->Branch("gen_jpsi_pt",    &gen_jpsi_pt_);
  tree_->Branch("gen_jpsi_eta",   &gen_jpsi_eta_);
  tree_->Branch("gen_jpsi_phi",   &gen_jpsi_phi_);
  tree_->Branch("gen_dimu_pt",    &gen_dimu_pt_);
  tree_->Branch("gen_dimu_eta",   &gen_dimu_eta_);
  tree_->Branch("gen_dimu_phi",   &gen_dimu_phi_);
  tree_->Branch("gen_dimu_mass",  &gen_dimu_mass_);
  tree_->Branch("gen_dimu_dR",    &gen_dimu_dR_);
  tree_->Branch("gen_mu_pt",      &gen_mu_pt_);
  tree_->Branch("gen_mu_eta",     &gen_mu_eta_);
  tree_->Branch("gen_mu_phi",     &gen_mu_phi_);
  tree_->Branch("gen_mu_pid",     &gen_mu_pid_);
  tree_->Branch("gen_mu_charge",  &gen_mu_charge_);
  tree_->Branch("gen_mu_pt_m",    &gen_mu_pt_m_);
  tree_->Branch("gen_mu_pt_p",    &gen_mu_pt_p_);
  tree_->Branch("gen_mu_px",      &gen_mu_px_);
  tree_->Branch("gen_mu_py",      &gen_mu_py_);
  tree_->Branch("gen_mu_pz",      &gen_mu_pz_);

  // pair-level scalars
  tree_->Branch("coll",     &coll_);
  tree_->Branch("invm",     &invm_);
  tree_->Branch("invm_bef", &invm_bef_);

  // DCA between the two muon tracks
  tree_->Branch("dca2d",        &dca2d_);
  tree_->Branch("dca3d",        &dca3d_);
  tree_->Branch("dca3d_cp_rho", &dca3d_cp_rho_);

  // gen-reco muon matching per pair
  tree_->Branch("genReco_dEta",    &genReco_dEta_);
  tree_->Branch("genReco_dPhi",    &genReco_dPhi_);
  tree_->Branch("genReco_dR",      &genReco_dR_);
  tree_->Branch("genReco_ptRatio", &genReco_ptRatio_);
  tree_->Branch("genReco_pxRatio", &genReco_pxRatio_);
  tree_->Branch("genReco_pyRatio", &genReco_pyRatio_);

  // muon-level original
  tree_->Branch("pt_mu",      &pt_mu_);
  tree_->Branch("eta_mu",     &eta_mu_);
  tree_->Branch("charge_mu",  &charge_mu_);
  tree_->Branch("ip_mu",      &ip_mu_);
  tree_->Branch("time_mu",    &time_mu_);
  tree_->Branch("cscHit_mu",  &cscHit_mu_);
  tree_->Branch("dtHit_mu",   &dtHit_mu_);
  tree_->Branch("dir_mu",     &dir_mu_);
  tree_->Branch("normChi_mu", &normChi_mu_);
  tree_->Branch("px_mu",      &px_mu_);
  tree_->Branch("py_mu",      &py_mu_);
  tree_->Branch("pz_mu",      &pz_mu_);
  tree_->Branch("phi_mu",     &phi_mu_);

  // reco muon pT split by charge
  tree_->Branch("reco_mu_pt_m", &reco_mu_pt_m_);
  tree_->Branch("reco_mu_pt_p", &reco_mu_pt_p_);

  // per-fit branches
  for (int f = 0; f < kNFits; ++f) {
    const std::string s = std::to_string(f);
    tree_->Branch(("isValid_vtx_"    + s).c_str(), &isValid_vtx_[f]);
    tree_->Branch(("invm_aft_"       + s).c_str(), &invm_aft_[f]);
    tree_->Branch(("Lxy_"            + s).c_str(), &Lxy_[f]);
    tree_->Branch(("LxyErr_"         + s).c_str(), &LxyErr_[f]);
    tree_->Branch(("vx_"             + s).c_str(), &vx_[f]);
    tree_->Branch(("vy_"             + s).c_str(), &vy_[f]);
    tree_->Branch(("normChi_vtx_"    + s).c_str(), &normChi_vtx_[f]);
    tree_->Branch(("prob_vtx_"       + s).c_str(), &prob_vtx_[f]);
    tree_->Branch(("pt_refit_"       + s).c_str(), &pt_refit_[f]);
    tree_->Branch(("eta_refit_"      + s).c_str(), &eta_refit_[f]);
    tree_->Branch(("relPtUnc_refit_" + s).c_str(), &relPtUnc_refit_[f]);
  }
}


// == clearVectors ==============================================================
void staMuon::clearVectors() {
  gen_jpsi_pt_.clear();  gen_jpsi_eta_.clear();  gen_jpsi_phi_.clear();
  gen_dimu_pt_.clear();  gen_dimu_eta_.clear();  gen_dimu_phi_.clear();
  gen_dimu_mass_.clear(); gen_dimu_dR_.clear();
  gen_mu_pt_.clear();    gen_mu_eta_.clear();
  gen_mu_phi_.clear();   gen_mu_pid_.clear();    gen_mu_charge_.clear();
  gen_mu_pt_m_.clear();  gen_mu_pt_p_.clear();
  gen_mu_px_.clear();    gen_mu_py_.clear();     gen_mu_pz_.clear();

  coll_.clear();
  invm_.clear();
  invm_bef_.clear();

  dca2d_.clear();
  dca3d_.clear();
  dca3d_cp_rho_.clear();

  genReco_dEta_.clear();
  genReco_dPhi_.clear();
  genReco_dR_.clear();
  genReco_ptRatio_.clear();
  genReco_pxRatio_.clear();
  genReco_pyRatio_.clear();

  pt_mu_.clear();     eta_mu_.clear();    charge_mu_.clear();
  ip_mu_.clear();     time_mu_.clear();
  cscHit_mu_.clear(); dtHit_mu_.clear();
  dir_mu_.clear();    normChi_mu_.clear();
  px_mu_.clear();     py_mu_.clear();
  pz_mu_.clear();     phi_mu_.clear();

  reco_mu_pt_m_.clear();
  reco_mu_pt_p_.clear();

  for (int f = 0; f < kNFits; ++f) {
    isValid_vtx_[f].clear();
    invm_aft_[f].clear();
    Lxy_[f].clear();
    LxyErr_[f].clear();
    vx_[f].clear();
    vy_[f].clear();
    normChi_vtx_[f].clear();
    prob_vtx_[f].clear();
    pt_refit_[f].clear();
    eta_refit_[f].clear();
    relPtUnc_refit_[f].clear();
  }
}


// == pushMuonPair ==============================================================
void staMuon::pushMuonPair(const MuonVars& v1, const MuonVars& v2) {
  pt_mu_.push_back     ({v1.pt,      v2.pt});
  eta_mu_.push_back    ({v1.eta,     v2.eta});
  charge_mu_.push_back ({v1.charge,  v2.charge});
  ip_mu_.push_back     ({v1.ip,      v2.ip});
  time_mu_.push_back   ({v1.time,    v2.time});
  cscHit_mu_.push_back ({v1.cscHit,  v2.cscHit});
  dtHit_mu_.push_back  ({v1.dtHit,   v2.dtHit});
  dir_mu_.push_back    ({v1.dir,     v2.dir});
  normChi_mu_.push_back({v1.normChi, v2.normChi});
  px_mu_.push_back     ({v1.px,      v2.px});
  py_mu_.push_back     ({v1.py,      v2.py});
  pz_mu_.push_back     ({v1.pz,      v2.pz});
  phi_mu_.push_back    ({v1.phi,     v2.phi});
}


// == analyze ===================================================================
void staMuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  clearVectors();

  edm::Handle<std::vector<pat::Muon>>         muons;
  edm::Handle<std::vector<pat::Muon>>         disMuons;
  edm::Handle<std::vector<reco::Vertex>>      vertices;
  edm::Handle<std::vector<reco::GenParticle>> genParticles;

  iEvent.getByToken(muonToken_,    muons);
  iEvent.getByToken(disMuonToken_, disMuons);
  iEvent.getByToken(pvToken_,      vertices);
  iEvent.getByToken(genToken_,     genParticles);

  const auto& ttBuilder = iSetup.getData(ttBuilderToken_);

  // == fill event-level gen branches =========================================
  if (genParticles.isValid()) {
    const reco::GenParticle* genMuPlus  = nullptr;
    const reco::GenParticle* genMuMinus = nullptr;

    for (const auto& gp : *genParticles) {
      const int pid = gp.pdgId();

      if (pid == 443) {
        gen_jpsi_pt_.push_back(gp.pt());
        gen_jpsi_eta_.push_back(gp.eta());
        gen_jpsi_phi_.push_back(gp.phi());
      }

      if (std::abs(pid) == 13) {
        gen_mu_pt_.push_back(gp.pt());
        gen_mu_eta_.push_back(gp.eta());
        gen_mu_phi_.push_back(gp.phi());
        gen_mu_pid_.push_back(pid);
        gen_mu_charge_.push_back(gp.charge());
        gen_mu_px_.push_back(static_cast<float>(gp.px()));
        gen_mu_py_.push_back(static_cast<float>(gp.py()));
        gen_mu_pz_.push_back(static_cast<float>(gp.pz()));

        if (pid == 13)  gen_mu_pt_m_.push_back(gp.pt());
        if (pid == -13) gen_mu_pt_p_.push_back(gp.pt());

        if (pid == -13) {
          if (!genMuPlus  || gp.pt() > genMuPlus->pt()) genMuPlus  = &gp;
        } else {
          if (!genMuMinus || gp.pt() > genMuMinus->pt()) genMuMinus = &gp;
        }
      }
    }

    if (genMuPlus && genMuMinus) {
      TLorentzVector gp1, gp2;
      gp1.SetPtEtaPhiM(genMuPlus->pt(),  genMuPlus->eta(),  genMuPlus->phi(),  0.105658f);
      gp2.SetPtEtaPhiM(genMuMinus->pt(), genMuMinus->eta(), genMuMinus->phi(), 0.105658f);
      TLorentzVector gDimu = gp1 + gp2;
      gen_dimu_pt_.push_back(gDimu.Pt());
      gen_dimu_eta_.push_back(gDimu.Eta());
      gen_dimu_phi_.push_back(gDimu.Phi());
      gen_dimu_mass_.push_back(gDimu.M());
      gen_dimu_dR_.push_back(myDeltaR(genMuPlus->eta(),  genMuPlus->phi(),
                                      genMuMinus->eta(), genMuMinus->phi()));
    }
  }

  if (!vertices.isValid() || vertices->empty()) { tree_->Fill(); return; }
  const reco::Vertex& pv = vertices->at(0);
  const reco::Vertex::Error& ve = pv.error();
  const GlobalError pvErr(ve.At(0,0), ve.At(1,0), ve.At(1,1),
                          ve.At(2,0), ve.At(2,1), ve.At(2,2));

  // == process one muon collection ===========================================
  auto processMuonCollection = [&](const std::vector<pat::Muon>& coll, int collId) {

    for (const auto& mu : coll) {
      const reco::TrackRef ot = mu.outerTrack();
      if (ot.isNull()) continue;
      const float pt = ot->pt();
      const int q = mu.charge();
      if (q == -1)      reco_mu_pt_m_.push_back(pt);
      else if (q == +1) reco_mu_pt_p_.push_back(pt);
    }

    std::vector<const pat::Muon*> muonsPtr;
    muonsPtr.reserve(coll.size());
    for (const auto& mu : coll) muonsPtr.push_back(&mu);

    coll_.push_back(static_cast<int>(muonsPtr.size()));

    constexpr float mMu = 0.105658f;

    // OuterTrack-based invariant mass for all OS combinations with outerTrack pt > 10
    for (size_t i = 0; i < muonsPtr.size(); ++i) {
      for (size_t j = i + 1; j < muonsPtr.size(); ++j) {
        const pat::Muon* mu1 = muonsPtr[i];
        const pat::Muon* mu2 = muonsPtr[j];

        if (mu1->charge() * mu2->charge() >= 0) continue;

        const reco::TrackRef ot1 = mu1->outerTrack();
        const reco::TrackRef ot2 = mu2->outerTrack();
        if (ot1.isNull() || ot2.isNull()) continue;

        if (ot1->pt() <= 10.0f || ot2->pt() <= 10.0f) continue;

        TLorentzVector p1, p2;
        p1.SetPtEtaPhiM(ot1->pt(), ot1->eta(), ot1->phi(), mMu);
        p2.SetPtEtaPhiM(ot2->pt(), ot2->eta(), ot2->phi(), mMu);

        invm_.push_back((p1 + p2).M());
      }
    }

    // Pair loop
    std::vector<const pat::Muon*> sa;
    sa.reserve(coll.size());
    for (const auto& mu : coll) sa.push_back(&mu);

    for (size_t i = 0; i < sa.size(); ++i) {
      for (size_t j = i + 1; j < sa.size(); ++j) {
        const pat::Muon* mu1 = sa[i];
        const pat::Muon* mu2 = sa[j];
        if (mu1->charge() * mu2->charge() >= 0) continue;

        TLorentzVector p1, p2;
        p1.SetPtEtaPhiM(mu1->pt(), mu1->eta(), mu1->phi(), mMu);
        p2.SetPtEtaPhiM(mu2->pt(), mu2->eta(), mu2->phi(), mMu);

        const auto mv1 = MuonVars::from(*mu1);
        const auto mv2 = MuonVars::from(*mu2);

        invm_bef_.push_back((p1 + p2).M());
        pushMuonPair(mv1, mv2);

        // ---- gen-reco matching (charge-matched) ----------------------------
        {
          std::array<float, 2> matchDEta    = {{kInvalid, kInvalid}};
          std::array<float, 2> matchDPhi    = {{kInvalid, kInvalid}};
          std::array<float, 2> matchDR      = {{kInvalid, kInvalid}};
          std::array<float, 2> matchPtRatio = {{kInvalid, kInvalid}};
          std::array<float, 2> matchPxRatio = {{kInvalid, kInvalid}};
          std::array<float, 2> matchPyRatio = {{kInvalid, kInvalid}};

          const std::array<const pat::Muon*, 2> recoMuArr = {{mu1, mu2}};

          if (genParticles.isValid()) {
            for (int k = 0; k < 2; ++k) {
              const float recoEta    = recoMuArr[k]->eta();
              const float recoPhi    = recoMuArr[k]->phi();
              const float recoPt     = recoMuArr[k]->pt();
              const float recoPx     = recoMuArr[k]->px();
              const float recoPy     = recoMuArr[k]->py();
              const int   recoCharge = recoMuArr[k]->charge();
              float bestDR = std::numeric_limits<float>::max();

              for (const auto& gp : *genParticles) {
                if (std::abs(gp.pdgId()) != 13)  continue;
                if (gp.charge() != recoCharge)   continue;  // charge match

                const float dR = myDeltaR(recoEta, recoPhi,
                                          static_cast<float>(gp.eta()),
                                          static_cast<float>(gp.phi()));
                if (dR < bestDR) {
                  bestDR           = dR;
                  matchDEta[k]     = recoEta - static_cast<float>(gp.eta());
                  matchDPhi[k]     = deltaPhi(recoPhi, static_cast<float>(gp.phi()));
                  matchDR[k]       = dR;
                  matchPtRatio[k]  = (gp.pt() > 0.f)
                                       ? recoPt / static_cast<float>(gp.pt())
                                       : kInvalid;
                  matchPxRatio[k]  = (std::abs(gp.px()) > 1e-3)
                                       ? recoPx / static_cast<float>(gp.px())
                                       : kInvalid;
                  matchPyRatio[k]  = (std::abs(gp.py()) > 1e-3)
                                       ? recoPy / static_cast<float>(gp.py())
                                       : kInvalid;
                }
              }
            }
          }

          genReco_dEta_.push_back   ({matchDEta[0],    matchDEta[1]});
          genReco_dPhi_.push_back   ({matchDPhi[0],    matchDPhi[1]});
          genReco_dR_.push_back     ({matchDR[0],      matchDR[1]});
          genReco_ptRatio_.push_back({matchPtRatio[0], matchPtRatio[1]});
          genReco_pxRatio_.push_back({matchPxRatio[0], matchPxRatio[1]});
          genReco_pyRatio_.push_back({matchPyRatio[0], matchPyRatio[1]});
        }

        // ---- DCA between the two muon tracks --------------------------------
        reco::TransientTrack tt1 = ttBuilder.build(mu1->bestTrack());
        reco::TransientTrack tt2 = ttBuilder.build(mu2->bestTrack());

        float dca2d = kInvalid;
        {
          ClosestApproachInRPhi cApp;
          cApp.calculate(tt1.initialFreeState(), tt2.initialFreeState());
          if (cApp.status())
            dca2d = static_cast<float>(cApp.distance());
        }

        float dca3d = kInvalid;
        float dca3d_cp_rho = kInvalid;
        {
          TwoTrackMinimumDistance ttmd;
          if (ttmd.calculate(tt1.initialFreeState(), tt2.initialFreeState())) {
            dca3d = static_cast<float>(ttmd.distance());
            const GlobalPoint cp = ttmd.crossingPoint();
            dca3d_cp_rho = static_cast<float>(std::sqrt(cp.x() * cp.x() + cp.y() * cp.y()));
          }
        }

        dca2d_.push_back(dca2d);
        dca3d_.push_back(dca3d);
        dca3d_cp_rho_.push_back(dca3d_cp_rho);

        // == Kalman vertex fitting ============================================
        KalmanVertexFitter kvf(true, true);
        std::vector<reco::TransientTrack> currentTracks = {tt1, tt2};

        for (int ifit = 0; ifit < kNFits; ++ifit) {
          const GlobalPoint seedPoint(300.f, 0.f, 0.f);
          AlgebraicSymMatrix33 mat;
          mat[0][0] = 1.f;
          mat[1][1] = 1.f;
          mat[2][2] = 1.f;
          const GlobalError seedError(mat);

          TransientVertex tv = kvf.vertex(currentTracks);
          const bool fitOk = tv.isValid() && tv.refittedTracks().size() == 2;

          if (!fitOk) {
            for (int jfit = ifit; jfit < kNFits; ++jfit) {
              isValid_vtx_[jfit].push_back(0);
              invm_aft_[jfit].push_back(kInvalid);
              Lxy_[jfit].push_back(kInvalid);
              LxyErr_[jfit].push_back(kInvalid);
              vx_[jfit].push_back(kInvalid);
              vy_[jfit].push_back(kInvalid);
              normChi_vtx_[jfit].push_back(kInvalid);
              prob_vtx_[jfit].push_back(kInvalid);
              pt_refit_[jfit].push_back      ({kInvalid, kInvalid});
              eta_refit_[jfit].push_back     ({kInvalid, kInvalid});
              relPtUnc_refit_[jfit].push_back({kInvalid, kInvalid});
            }
            break;
          }

          const auto& refitTrks = tv.refittedTracks();

          const auto gm1 = refitTrks[0].impactPointState().globalMomentum();
          const auto gm2 = refitTrks[1].impactPointState().globalMomentum();
          TLorentzVector r1, r2;
          r1.SetXYZM(gm1.x(), gm1.y(), gm1.z(), mMu);
          r2.SetXYZM(gm2.x(), gm2.y(), gm2.z(), mMu);

          const GlobalPoint svPos = tv.position();
          const GlobalError svErr = tv.positionError();
          const GlobalPoint disp(svPos.x(), svPos.y(), 0.f);
          const GlobalError totErr = svErr + pvErr;

          const float normChi = (tv.degreesOfFreedom() > 0)
                                  ? tv.totalChiSquared() / tv.degreesOfFreedom()
                                  : kInvalid;
          const float prob = TMath::Prob(tv.totalChiSquared(),
                                         static_cast<int>(tv.degreesOfFreedom()));

          isValid_vtx_[ifit].push_back(1);
          invm_aft_[ifit].push_back((r1 + r2).M());
          Lxy_[ifit].push_back(disp.perp());
          LxyErr_[ifit].push_back(std::sqrt(totErr.rerr(disp)));
          vx_[ifit].push_back(svPos.x());
          vy_[ifit].push_back(svPos.y());
          normChi_vtx_[ifit].push_back(normChi);
          prob_vtx_[ifit].push_back(prob);

          const auto rv1 = RefitVars::from(refitTrks[0], gm1);
          const auto rv2 = RefitVars::from(refitTrks[1], gm2);
          pt_refit_[ifit].push_back      ({rv1.pt,       rv2.pt});
          eta_refit_[ifit].push_back     ({rv1.eta,      rv2.eta});
          relPtUnc_refit_[ifit].push_back({rv1.relPtUnc, rv2.relPtUnc});

          if (ifit + 1 < kNFits)
            currentTracks = {refitTrks[0], refitTrks[1]};
        }
      }
    }

    (void)collId;
  };

  if (disMuons.isValid()) processMuonCollection(*disMuons, 1);

  tree_->Fill();
}

DEFINE_FWK_MODULE(staMuon);
