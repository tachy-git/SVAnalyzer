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

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "FWCore/Utilities/interface/ESInputTag.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <vector>
#include <cmath>

// Sentinel value for quantities undefined when vertex fit fails
static constexpr float kInvalid = -999.f;

class staMuon : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit staMuon(const edm::ParameterSet&);
  ~staMuon() override {}

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  // Tokens
  edm::EDGetTokenT<std::vector<pat::Muon>>    muonToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>>    disMuonToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  // Tree
  TTree* tree_;

  // ==================================================================
  // Branch variables
  //
  // pair-level:  vector<T>            -- one entry per OS dimuon pair
  // muon-level:  vector<vector<T>>    -- one inner vector per pair,
  //                                      always size 2: [0]=mu1, [1]=mu2
  //
  // When vertex fit fails:
  //   isValid_vtx_ = 0, fit-derived pair branches = kInvalid
  //   refit muon branches = {kInvalid, kInvalid}
  // ==================================================================

  // --- pair-level ---
  std::vector<int>   coll_;           // 0=slimmedMuons, 1=slimmedDisplacedMuons
  std::vector<int>   isValid_vtx_;    // 1=fit succeeded, 0=fit failed
  std::vector<float> invm_bef_;       // dimuon inv. mass before vertex fit
  std::vector<float> invm_aft_;       // dimuon inv. mass after vertex fit  (kInvalid if failed)
  std::vector<float> Lxy_;            // transverse displacement |SV-PV|    (kInvalid if failed)
  std::vector<float> LxyErr_;         // uncertainty on Lxy                 (kInvalid if failed)
  std::vector<float> normChi_vtx_;    // vertex chi2/ndof                   (kInvalid if failed)
  std::vector<float> prob_vtx_;       // vertex fit probability              (kInvalid if failed)

  // --- muon-level: PAT muon / original track (always filled) ---
  std::vector<std::vector<float>> pt_mu_;
  std::vector<std::vector<float>> eta_mu_;
  std::vector<std::vector<float>> ip_mu_;          // dB() -- 2D IP wrt beam spot
  std::vector<std::vector<float>> time_mu_;        // timeAtIpInOut
  std::vector<std::vector<int>>   cscHit_mu_;      // number of valid CSC hits
  std::vector<std::vector<int>>   dtHit_mu_;       // number of valid DT hits
  std::vector<std::vector<int>>   dir_mu_;
  std::vector<std::vector<float>> normChi_mu_;

  // --- muon-level: refitted track (kInvalid if fit failed) ---
  std::vector<std::vector<float>> pt_refit_;
  std::vector<std::vector<float>> eta_refit_;
  std::vector<std::vector<float>> relPtUnc_refit_;

  // ==================================================================

  static bool passSASelection(const pat::Muon& mu) {
    if (!mu.isStandAloneMuon())    return false;
    if ( mu.isTrackerMuon())       return false;
    if (mu.pt() <= 5.)             return false;
    if (std::abs(mu.eta()) >= 2.4) return false;
    return true;
  }

  void clearVectors();

  // Structs to collect variables before pushing to branches
  struct MuonVars {
    float pt, eta, ip, time, normChi;
    int   cscHit, dtHit, dir;

    static MuonVars from(const pat::Muon& mu) {
      MuonVars v;
      v.pt    = mu.pt();
      v.eta   = mu.eta();
      v.ip    = mu.dB();
      v.time  = mu.time().timeAtIpInOut;
      const auto& hp = mu.bestTrack()->hitPattern();
      v.cscHit   = hp.numberOfValidMuonCSCHits();
      v.dtHit    = hp.numberOfValidMuonDTHits();
      v.dir      = static_cast<int>(mu.bestTrack()->seedDirection());
      v.normChi  = mu.bestTrack()->normalizedChi2();
      return v;
    }
  };

  struct RefitVars {
    float pt, eta, relPtUnc;

    static RefitVars from(const reco::TransientTrack& trk) {
      const reco::Track& t = trk.track();
      RefitVars v;
      v.pt       = t.pt();
      v.eta      = t.eta();
      v.relPtUnc = (t.pt() > 0.f) ? t.ptError() / t.pt() : -1.f;
      return v;
    }

    static RefitVars invalid() {
      RefitVars v;
      v.pt = v.eta = v.relPtUnc = kInvalid;
      return v;
    }
  };

  void pushMuonPair(const MuonVars& v1, const MuonVars& v2);
  void pushRefitPair(const RefitVars& v1, const RefitVars& v2);
};


// -------------------- Constructor --------------------
staMuon::staMuon(const edm::ParameterSet& iConfig) {
  usesResource("TFileService");
  muonToken_      = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
  disMuonToken_   = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("disMuons"));
  pvToken_        = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"));
  ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(
                       edm::ESInputTag("", "TransientTrackBuilder"));
}


// -------------------- beginJob --------------------
void staMuon::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "staMuon tree");

  // pair-level
  tree_->Branch("coll",        &coll_);
  tree_->Branch("isValid_vtx", &isValid_vtx_);
  tree_->Branch("invm_bef",    &invm_bef_);
  tree_->Branch("invm_aft",    &invm_aft_);
  tree_->Branch("Lxy",         &Lxy_);
  tree_->Branch("LxyErr",      &LxyErr_);
  tree_->Branch("normChi_vtx", &normChi_vtx_);
  tree_->Branch("prob_vtx",    &prob_vtx_);

  // muon-level (original) -- vector<vector<T>>
  tree_->Branch("pt_mu",       &pt_mu_);
  tree_->Branch("eta_mu",      &eta_mu_);
  tree_->Branch("ip_mu",       &ip_mu_);
  tree_->Branch("time_mu",     &time_mu_);
  tree_->Branch("cscHit_mu",   &cscHit_mu_);
  tree_->Branch("dtHit_mu",    &dtHit_mu_);
  tree_->Branch("dir_mu",      &dir_mu_);
  tree_->Branch("normChi_mu",  &normChi_mu_);

  // muon-level (refit) -- vector<vector<T>>
  tree_->Branch("pt_refit",       &pt_refit_);
  tree_->Branch("eta_refit",      &eta_refit_);
  tree_->Branch("relPtUnc_refit", &relPtUnc_refit_);
}


// -------------------- clearVectors --------------------
void staMuon::clearVectors() {
  coll_.clear();        isValid_vtx_.clear();
  invm_bef_.clear();    invm_aft_.clear();
  Lxy_.clear();         LxyErr_.clear();
  normChi_vtx_.clear(); prob_vtx_.clear();

  pt_mu_.clear();       eta_mu_.clear();
  ip_mu_.clear();       time_mu_.clear();
  cscHit_mu_.clear();   dtHit_mu_.clear(); dir_mu_.clear(); normChi_mu_.clear();

  pt_refit_.clear();    eta_refit_.clear();
  relPtUnc_refit_.clear();
}


// -------------------- pushMuonPair --------------------
void staMuon::pushMuonPair(const MuonVars& v1, const MuonVars& v2) {
  pt_mu_.push_back({v1.pt,       v2.pt});
  eta_mu_.push_back({v1.eta,     v2.eta});
  ip_mu_.push_back({v1.ip,       v2.ip});
  time_mu_.push_back({v1.time,   v2.time});
  cscHit_mu_.push_back({v1.cscHit, v2.cscHit});
  dtHit_mu_.push_back({v1.dtHit,   v2.dtHit});
  dir_mu_.push_back({v1.dir, v2.dir});
  normChi_mu_.push_back({v1.normChi,   v2.normChi});
}


// -------------------- pushRefitPair --------------------
void staMuon::pushRefitPair(const RefitVars& v1, const RefitVars& v2) {
  pt_refit_.push_back({v1.pt,         v2.pt});
  eta_refit_.push_back({v1.eta,       v2.eta});
  relPtUnc_refit_.push_back({v1.relPtUnc, v2.relPtUnc});
}


// -------------------- analyze --------------------
void staMuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  clearVectors();

  edm::Handle<std::vector<pat::Muon>>    muons;
  edm::Handle<std::vector<pat::Muon>>    disMuons;
  edm::Handle<std::vector<reco::Vertex>> vertices;

  iEvent.getByToken(muonToken_,    muons);
  iEvent.getByToken(disMuonToken_, disMuons);
  iEvent.getByToken(pvToken_,      vertices);

  const auto& ttBuilder = iSetup.getData(ttBuilderToken_);

  if (vertices->empty()) return;
  const reco::Vertex& pv = vertices->at(0);
  const GlobalPoint pvPos(pv.x(), pv.y(), pv.z());
  const reco::Vertex::Error& ve = pv.error();
  const GlobalError pvErr(ve.At(0,0), ve.At(1,0), ve.At(1,1),
                           ve.At(2,0), ve.At(2,1), ve.At(2,2));

  // -----------------------------------------------------------------------
  auto processMuonCollection = [&](const std::vector<pat::Muon>& coll, int collId) {

    std::vector<const pat::Muon*> sa;
    sa.reserve(coll.size());
    for (const auto& mu : coll)
      if (passSASelection(mu)) sa.push_back(&mu);

    for (size_t i = 0; i < sa.size(); ++i) {
      for (size_t j = i + 1; j < sa.size(); ++j) {
        const pat::Muon* mu1 = sa[i];
        const pat::Muon* mu2 = sa[j];
        if (mu1->charge() * mu2->charge() >= 0) continue;

        // ---- invariant mass before fit ----
        TLorentzVector p1, p2;
        p1.SetPtEtaPhiM(mu1->pt(), mu1->eta(), mu1->phi(), 0.105658f);
        p2.SetPtEtaPhiM(mu2->pt(), mu2->eta(), mu2->phi(), 0.105658f);
        float invm_bef = (p1 + p2).M();

        // ---- PAT muon info (always available) ----
        const auto mv1 = MuonVars::from(*mu1);
        const auto mv2 = MuonVars::from(*mu2);

        // ---- Kalman Vertex Fit ----
        KalmanVertexFitter kvf(true, true);
        std::vector<reco::TransientTrack> tTracks;
        tTracks.emplace_back(ttBuilder.build(mu1->bestTrack()));
        tTracks.emplace_back(ttBuilder.build(mu2->bestTrack()));
        TransientVertex tv = kvf.vertex(tTracks);

        const bool fitOk = tv.isValid() && tv.refittedTracks().size() == 2;

        // ---- always fill: pair header + PAT muon branches ----
        coll_.push_back(collId);
        invm_bef_.push_back(invm_bef);
        pushMuonPair(mv1, mv2);

        if (!fitOk) {
          isValid_vtx_.push_back(0);
          invm_aft_.push_back(kInvalid);
          Lxy_.push_back(kInvalid);
          LxyErr_.push_back(kInvalid);
          normChi_vtx_.push_back(kInvalid);
          prob_vtx_.push_back(kInvalid);
          pushRefitPair(RefitVars::invalid(), RefitVars::invalid());
          continue;
        }

        // ---- fit succeeded ----
        const auto& refitTrks = tv.refittedTracks();

        // invariant mass after fit
        constexpr float mMu = 0.105658f;
        const auto gm1 = refitTrks[0].impactPointState().globalMomentum();
        const auto gm2 = refitTrks[1].impactPointState().globalMomentum();
        TLorentzVector r1, r2;
        r1.SetXYZM(gm1.x(), gm1.y(), gm1.z(), mMu);
        r2.SetXYZM(gm2.x(), gm2.y(), gm2.z(), mMu);
        float invm_aft = (r1 + r2).M();

        // Lxy
        const GlobalPoint svPos = tv.position();
        const GlobalError svErr = tv.positionError();
        const GlobalPoint disp(svPos.x() - pvPos.x(),
                               svPos.y() - pvPos.y(),
                               0.f);
        const GlobalError totErr = svErr + pvErr;
        float lxy    = disp.perp();
        float lxyErr = std::sqrt(totErr.rerr(disp));

        // vertex chi2 / prob
        float normChi_vtx = (tv.degreesOfFreedom() > 0)
                              ? tv.totalChiSquared() / tv.degreesOfFreedom()
                              : kInvalid;
        float prob_vtx = TMath::Prob(tv.totalChiSquared(),
                                     static_cast<int>(tv.degreesOfFreedom()));

        // fill pair-level fit branches
        isValid_vtx_.push_back(1);
        invm_aft_.push_back(invm_aft);
        Lxy_.push_back(lxy);
        LxyErr_.push_back(lxyErr);
        normChi_vtx_.push_back(normChi_vtx);
        prob_vtx_.push_back(prob_vtx);

        // fill refit branches
        pushRefitPair(RefitVars::from(refitTrks[0]),
                      RefitVars::from(refitTrks[1]));
      }
    }
  };

  if (muons.isValid())    processMuonCollection(*muons,    0);
  if (disMuons.isValid()) processMuonCollection(*disMuons, 1);

  tree_->Fill();
}


// -------------------- Register Module --------------------
DEFINE_FWK_MODULE(staMuon);
