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
#include "DataFormats/PatCandidates/interface/Jet.h"
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
#include <limits>
#include <string>

static constexpr float kInvalid = -999.f;

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
  edm::EDGetTokenT<std::vector<pat::Jet>>          jetToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  TTree* tree_;

  // == event-level: gen muons ================================================
  std::vector<float> gen_mu_pt_;
  std::vector<float> gen_mu_eta_;
  std::vector<float> gen_mu_phi_;
  std::vector<int>   gen_mu_charge_;

  // == pair-level scalars ====================================================
  std::vector<float> invm_;
  std::vector<float> invm_bef_;

  // DCA between the two muon tracks
  std::vector<float> dca3d_;
  std::vector<float> dca3d_cp_rho_;

  // == muon-level: original PAT muon =========================================
  std::vector<std::vector<float>> pt_mu_;
  std::vector<std::vector<float>> eta_mu_;
  std::vector<std::vector<float>> phi_mu_;
  std::vector<std::vector<int>>   charge_mu_;
  std::vector<std::vector<float>> ip_mu_;
  std::vector<std::vector<float>> time_mu_;
  std::vector<std::vector<int>>   cscHit_mu_;
  std::vector<std::vector<int>>   dtHit_mu_;
  std::vector<std::vector<float>> normChi_mu_;

  // == single vertex fit =====================================================
  std::vector<int>   isValid_vtx_;
  std::vector<float> invm_aft_;
  std::vector<float> Lxy_;
  std::vector<float> LxyErr_;
  std::vector<float> vx_;
  std::vector<float> vy_;
  std::vector<float> normChi_vtx_;
  std::vector<float> prob_vtx_;
  std::vector<std::vector<float>> pt_refit_;
  std::vector<std::vector<float>> eta_refit_;
  std::vector<std::vector<float>> relPtUnc_refit_;

  // == jet information =======================================================
  std::vector<float> jet_pt_;
  std::vector<float> jet_eta_;
  std::vector<float> jet_phi_;
  std::vector<int>   jet_partonFlavour_;

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
    int   charge, cscHit, dtHit;

    static MuonVars from(const pat::Muon& mu) {
      MuonVars v;
      v.pt      = mu.pt();
      v.eta     = mu.eta();
      v.phi     = mu.phi();
      v.charge  = mu.charge();
      v.ip      = mu.dB();
      v.time    = mu.time().timeAtIpInOut;
      const auto& hp = mu.bestTrack()->hitPattern();
      v.cscHit  = hp.numberOfValidMuonCSCHits();
      v.dtHit   = hp.numberOfValidMuonDTHits();
      v.normChi = mu.bestTrack()->normalizedChi2();
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
    static RefitVars invalid() { return {kInvalid, kInvalid, kInvalid}; }
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
  jetToken_     = consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"));
  ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(
                       edm::ESInputTag("", "TransientTrackBuilder"));
}


// == beginJob ==================================================================
void staMuon::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "staMuon tree");

  // event-level gen
  tree_->Branch("gen_mu_pt",      &gen_mu_pt_);
  tree_->Branch("gen_mu_eta",     &gen_mu_eta_);
  tree_->Branch("gen_mu_phi",     &gen_mu_phi_);
  tree_->Branch("gen_mu_charge",  &gen_mu_charge_);

  // pair-level scalars
  tree_->Branch("invm",     &invm_);
  tree_->Branch("invm_bef", &invm_bef_);

  // DCA
  tree_->Branch("dca3d",        &dca3d_);
  tree_->Branch("dca3d_cp_rho", &dca3d_cp_rho_);

  // muon-level
  tree_->Branch("pt_mu",      &pt_mu_);
  tree_->Branch("eta_mu",     &eta_mu_);
  tree_->Branch("phi_mu",     &phi_mu_);
  tree_->Branch("charge_mu",  &charge_mu_);
  tree_->Branch("ip_mu",      &ip_mu_);
  tree_->Branch("time_mu",    &time_mu_);
  tree_->Branch("cscHit_mu",  &cscHit_mu_);
  tree_->Branch("dtHit_mu",   &dtHit_mu_);
  tree_->Branch("normChi_mu", &normChi_mu_);

  // vertex fit (single)
  tree_->Branch("isValid_vtx",    &isValid_vtx_);
  tree_->Branch("invm_aft",       &invm_aft_);
  tree_->Branch("Lxy",            &Lxy_);
  tree_->Branch("LxyErr",         &LxyErr_);
  tree_->Branch("vx",             &vx_);
  tree_->Branch("vy",             &vy_);
  tree_->Branch("normChi_vtx",    &normChi_vtx_);
  tree_->Branch("prob_vtx",       &prob_vtx_);
  tree_->Branch("pt_refit",       &pt_refit_);
  tree_->Branch("eta_refit",      &eta_refit_);
  tree_->Branch("relPtUnc_refit", &relPtUnc_refit_);

  // jet
  tree_->Branch("jet_pt",            &jet_pt_);
  tree_->Branch("jet_eta",           &jet_eta_);
  tree_->Branch("jet_phi",           &jet_phi_);
  tree_->Branch("jet_partonFlavour", &jet_partonFlavour_);
}


// == clearVectors ==============================================================
void staMuon::clearVectors() {
  gen_mu_pt_.clear();    gen_mu_eta_.clear();
  gen_mu_phi_.clear();   gen_mu_charge_.clear();

  invm_.clear();
  invm_bef_.clear();

  dca3d_.clear();
  dca3d_cp_rho_.clear();

  pt_mu_.clear();     eta_mu_.clear();    phi_mu_.clear();
  charge_mu_.clear(); ip_mu_.clear();     time_mu_.clear();
  cscHit_mu_.clear(); dtHit_mu_.clear();  normChi_mu_.clear();

  isValid_vtx_.clear();
  invm_aft_.clear();
  Lxy_.clear();
  LxyErr_.clear();
  vx_.clear();
  vy_.clear();
  normChi_vtx_.clear();
  prob_vtx_.clear();
  pt_refit_.clear();
  eta_refit_.clear();
  relPtUnc_refit_.clear();

  jet_pt_.clear();
  jet_eta_.clear();
  jet_phi_.clear();
  jet_partonFlavour_.clear();
}


// == pushMuonPair ==============================================================
void staMuon::pushMuonPair(const MuonVars& v1, const MuonVars& v2) {
  pt_mu_.push_back     ({v1.pt,      v2.pt});
  eta_mu_.push_back    ({v1.eta,     v2.eta});
  phi_mu_.push_back    ({v1.phi,     v2.phi});
  charge_mu_.push_back ({v1.charge,  v2.charge});
  ip_mu_.push_back     ({v1.ip,      v2.ip});
  time_mu_.push_back   ({v1.time,    v2.time});
  cscHit_mu_.push_back ({v1.cscHit,  v2.cscHit});
  dtHit_mu_.push_back  ({v1.dtHit,   v2.dtHit});
  normChi_mu_.push_back({v1.normChi, v2.normChi});
}


// == analyze ===================================================================
void staMuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  clearVectors();

  edm::Handle<std::vector<pat::Muon>>         muons;
  edm::Handle<std::vector<pat::Muon>>         disMuons;
  edm::Handle<std::vector<reco::Vertex>>      vertices;
  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  edm::Handle<std::vector<pat::Jet>>          jets;

  iEvent.getByToken(muonToken_,    muons);
  iEvent.getByToken(disMuonToken_, disMuons);
  iEvent.getByToken(pvToken_,      vertices);
  iEvent.getByToken(genToken_,     genParticles);
  iEvent.getByToken(jetToken_,     jets);

  const auto& ttBuilder = iSetup.getData(ttBuilderToken_);

  // == fill event-level gen branches =========================================
  if (genParticles.isValid()) {
    // collect all gen muons
    std::vector<const reco::GenParticle*> genMuons;
    for (const auto& gp : *genParticles) {
      if (std::abs(gp.pdgId()) != 13) continue;
      gen_mu_pt_.push_back(gp.pt());
      gen_mu_eta_.push_back(gp.eta());
      gen_mu_phi_.push_back(gp.phi());
      gen_mu_charge_.push_back(gp.charge());
      genMuons.push_back(&gp);
    }
  }

  // == fill jet branches =====================================================
  if (jets.isValid()) {
    for (const auto& jet : *jets) {
      jet_pt_.push_back(jet.pt());
      jet_eta_.push_back(jet.eta());
      jet_phi_.push_back(jet.phi());
      jet_partonFlavour_.push_back(jet.partonFlavour());
    }
  }

  if (!vertices.isValid() || vertices->empty()) { tree_->Fill(); return; }
  const reco::Vertex& pv = vertices->at(0);
  const reco::Vertex::Error& ve = pv.error();
  const GlobalError pvErr(ve.At(0,0), ve.At(1,0), ve.At(1,1),
                          ve.At(2,0), ve.At(2,1), ve.At(2,2));

  // == process muon collection ===============================================
  auto processMuonCollection = [&](const std::vector<pat::Muon>& coll) {

    constexpr float mMu = 0.105658f;

    std::vector<const pat::Muon*> muonsPtr;
    muonsPtr.reserve(coll.size());
    for (const auto& mu : coll){
      muonsPtr.push_back(&mu);
    }

    // OuterTrack-based invariant mass for OS combinations with outerTrack pt > 10
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
    std::vector<const pat::Muon*> SAmuonsPtr;
    SAmuonsPtr.reserve(coll.size());
    for (const auto& mu : coll){
      if( !passSASelection(mu) ) continue;
      SAmuonsPtr.push_back(&mu);
    }


    // Pair loop
    for (size_t i = 0; i < SAmuonsPtr.size(); ++i) {
      for (size_t j = i + 1; j < SAmuonsPtr.size(); ++j) {
        const pat::Muon* mu1 = SAmuonsPtr[i];
        const pat::Muon* mu2 = SAmuonsPtr[j];
        if (mu1->charge() * mu2->charge() >= 0) continue;

        TLorentzVector p1, p2;
        p1.SetPtEtaPhiM(mu1->pt(), mu1->eta(), mu1->phi(), mMu);
        p2.SetPtEtaPhiM(mu2->pt(), mu2->eta(), mu2->phi(), mMu);

        const auto mv1 = MuonVars::from(*mu1);
        const auto mv2 = MuonVars::from(*mu2);

        invm_bef_.push_back((p1 + p2).M());
        pushMuonPair(mv1, mv2);

        // ---- DCA between the two muon tracks --------------------------------
        reco::TransientTrack tt1 = ttBuilder.build(mu1->bestTrack());
        reco::TransientTrack tt2 = ttBuilder.build(mu2->bestTrack());

        float dca3d       = kInvalid;
        float dca3d_cp_rho = kInvalid;
        {
          TwoTrackMinimumDistance ttmd;
          if (ttmd.calculate(tt1.initialFreeState(), tt2.initialFreeState())) {
            dca3d = static_cast<float>(ttmd.distance());
            const GlobalPoint cp = ttmd.crossingPoint();
            dca3d_cp_rho = static_cast<float>(
                std::sqrt(cp.x() * cp.x() + cp.y() * cp.y()));
          }
        }
        dca3d_.push_back(dca3d);
        dca3d_cp_rho_.push_back(dca3d_cp_rho);

        // == Kalman vertex fitting (single) ===================================
        KalmanVertexFitter kvf(true, true);
        const std::vector<reco::TransientTrack> tracks = {tt1, tt2};
        TransientVertex tv = kvf.vertex(tracks);

        const bool fitOk = tv.isValid() && tv.refittedTracks().size() == 2;
        if (!fitOk) {
          isValid_vtx_.push_back(0);
          invm_aft_.push_back(kInvalid);
          Lxy_.push_back(kInvalid);
          LxyErr_.push_back(kInvalid);
          vx_.push_back(kInvalid);
          vy_.push_back(kInvalid);
          normChi_vtx_.push_back(kInvalid);
          prob_vtx_.push_back(kInvalid);
          pt_refit_.push_back      ({kInvalid, kInvalid});
          eta_refit_.push_back     ({kInvalid, kInvalid});
          relPtUnc_refit_.push_back({kInvalid, kInvalid});
          continue;
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

        isValid_vtx_.push_back(1);
        invm_aft_.push_back((r1 + r2).M());
        Lxy_.push_back(disp.perp());
        LxyErr_.push_back(std::sqrt(totErr.rerr(disp)));
        vx_.push_back(svPos.x());
        vy_.push_back(svPos.y());
        normChi_vtx_.push_back(normChi);
        prob_vtx_.push_back(prob);

        const auto rv1 = RefitVars::from(refitTrks[0]);
        const auto rv2 = RefitVars::from(refitTrks[1]);
        pt_refit_.push_back      ({rv1.pt,       rv2.pt});
        eta_refit_.push_back     ({rv1.eta,      rv2.eta});
        relPtUnc_refit_.push_back({rv1.relPtUnc, rv2.relPtUnc});
      }
    }
  };

  if (disMuons.isValid()) processMuonCollection(*disMuons);

  tree_->Fill();
}

DEFINE_FWK_MODULE(staMuon);
