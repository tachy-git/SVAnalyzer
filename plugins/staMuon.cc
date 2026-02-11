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

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/ESInputTag.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TMath.h"

#include <algorithm>
#include <vector>
#include <cmath>

class staMuon : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit staMuon(const edm::ParameterSet&);
  ~staMuon() override {}

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  // Tokens
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> disMuonToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  // Tree
  TTree* tree_;

  // Per-muon (flat) branches: concatenate slimmedMuons then slimmedDisplacedMuons
  std::vector<int>   mu_coll_;            // 0=slimmedMuons, 1=slimmedDisplacedMuons
  std::vector<int>   mu_charge_;

  // Track info (store -1 when missing)
  std::vector<float> inner_outerPt_;
  std::vector<float> inner_outerEta_;
  std::vector<float> inner_outerPhi_;

  std::vector<float> outer_innerPt_;
  std::vector<float> outer_innerEta_;
  std::vector<float> outer_innerPhi_;

  std::vector<float> dr_;
  std::vector<float> ptRatio_;

  // Standalone-pair info (within each collection separately)
  std::vector<float> sa_mll_;
  std::vector<float> sa_mll2_;
  std::vector<float> sa_mll3_;
  std::vector<float> sa_mu1_pt_;
  std::vector<float> sa_mu2_pt_;
  std::vector<float> sa_mu1_ip_;
  std::vector<float> sa_mu2_ip_;
  std::vector<int>   sa_coll_;            // 0=slimmedMuons, 1=slimmedDisplacedMuons

  std::vector<float> sa_Lxy_;
  std::vector<float> sa_prob_;

  // Helpers
  static bool passSASelection(const pat::Muon& mu) {
    if (!mu.isStandAloneMuon()) return false;
    if (mu.isTrackerMuon()) return false;                 // SA and NOT TrackerMuon
    if (mu.pt() <= 10.0) return false;
    if (std::abs(mu.eta()) >= 2.4) return false;
    return true;
  }

  void clearVectors();

};



// -------------------- Constructor --------------------
staMuon::staMuon(const edm::ParameterSet& iConfig) {
  usesResource("TFileService");

  muonToken_     = consumes<std::vector<pat::Muon>>( iConfig.getParameter<edm::InputTag>("muons") );
  disMuonToken_     = consumes<std::vector<pat::Muon>>( iConfig.getParameter<edm::InputTag>("disMuons") );
  ttBuilderToken_ =
      esConsumes<TransientTrackBuilder, TransientTrackRecord>(
          edm::ESInputTag("", "TransientTrackBuilder"));
}



// -------------------- beginJob --------------------
void staMuon::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "");

  tree_->Branch("mu_coll",        &mu_coll_);
  tree_->Branch("mu_charge",      &mu_charge_);

  tree_->Branch("inner_outerPt",  &inner_outerPt_);
  tree_->Branch("inner_outerEta", &inner_outerEta_);
  tree_->Branch("inner_outerPhi", &inner_outerPhi_);

  tree_->Branch("outer_innerPt",  &outer_innerPt_);
  tree_->Branch("outer_innerEta", &outer_innerEta_);
  tree_->Branch("outer_innerPhi", &outer_innerPhi_);

  tree_->Branch("dr", &dr_);
  tree_->Branch("ptRatio", &ptRatio_);

  tree_->Branch("sa_coll",   &sa_coll_);
  tree_->Branch("sa_mll",    &sa_mll_);
  tree_->Branch("sa_mll2",    &sa_mll2_);
  tree_->Branch("sa_mll3",    &sa_mll3_);
  tree_->Branch("sa_mu1_pt", &sa_mu1_pt_);
  tree_->Branch("sa_mu2_pt", &sa_mu2_pt_);
  tree_->Branch("sa_mu1_ip", &sa_mu1_ip_);
  tree_->Branch("sa_mu2_ip", &sa_mu2_ip_);

  tree_->Branch("sa_Lxy",    &sa_Lxy_);
  tree_->Branch("sa_prob",    &sa_prob_);
}



// -------------------- clearVectors --------------------
void staMuon::clearVectors() {
  mu_coll_.clear();
  mu_charge_.clear();

  inner_outerPt_.clear();
  inner_outerEta_.clear();
  inner_outerPhi_.clear();

  outer_innerPt_.clear();
  outer_innerEta_.clear();
  outer_innerPhi_.clear();

  dr_.clear();
  ptRatio_.clear();

  sa_coll_.clear();
  sa_mll_.clear();
  sa_mll2_.clear();
  sa_mll3_.clear();
  sa_mu1_pt_.clear();
  sa_mu2_pt_.clear();
  sa_mu1_ip_.clear();
  sa_mu2_ip_.clear();

  sa_Lxy_.clear();
  sa_prob_.clear();
}

// -------------------- analyze --------------------
void staMuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  clearVectors();

  edm::Handle<std::vector<pat::Muon>> muons;
  edm::Handle<std::vector<pat::Muon>> disMuons;

  iEvent.getByToken(muonToken_, muons);
  iEvent.getByToken(disMuonToken_, disMuons);

  //edm::ESHandle<TransientTrackBuilder> ttBuilder;
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);
  const auto& ttBuilder = iSetup.getData(ttBuilderToken_);

  // Helper lambda to process a collection
  auto processMuonCollection = [&](const std::vector<pat::Muon>& coll, int collId) {
    // 1) Check the global muons
    for (const auto& mu : coll) {
      if( ! mu.isGlobalMuon() ) continue;
      mu_coll_.push_back(collId);
      mu_charge_.push_back(mu.charge());


      float i_pt = -1.f, i_eta = -99.f, i_phi = -99.f;
      float o_pt = -1.f, o_eta = -99.f, o_phi = -99.f;

      // inner track
      if (mu.innerTrack().isNonnull()) {
        const reco::Track* intrk = mu.innerTrack().get();
        /*if( intrk->outerOk() ){
          i_pt  = intrk->outerMomentum().rho();
          i_eta = intrk->outerPosition().eta();
          i_phi = intrk->outerPosition().phi();
        }*/
        i_pt = intrk->pt();
        i_eta = intrk->eta();
        i_phi = intrk->phi();
      }
      // outer track
      if (mu.outerTrack().isNonnull()) {
        const reco::Track* outtrk = mu.outerTrack().get();
        /*if( outtrk->innerOk() ){
          o_pt  = outtrk->innerMomentum().rho();
          o_eta = outtrk->innerPosition().eta();
          o_phi = outtrk->innerPosition().phi();
        }*/
        o_pt = outtrk->pt();
        o_eta = outtrk->eta();
        o_phi = outtrk->phi();
      }

      inner_outerPt_.push_back(i_pt);
      inner_outerEta_.push_back(i_eta);
      inner_outerPhi_.push_back(i_phi);
      outer_innerPt_.push_back(o_pt);
      outer_innerEta_.push_back(o_eta);
      outer_innerPhi_.push_back(o_phi);

      float deta = i_eta - o_eta;
      float dphi = TVector2::Phi_mpi_pi(i_phi-o_phi);
      double dr = std::sqrt(deta*deta + dphi*dphi);
      dr_.push_back(dr);

      if( i_pt<0 || o_pt<0 ) continue;
      double ptRatio = ( i_pt - o_pt ) / i_pt ;
      ptRatio_.push_back(ptRatio);

    }

    // 2) Invariant mass of dimuon (standalone)
    std::vector<const pat::Muon*> sa;
    sa.reserve(coll.size());
    for (const auto& mu : coll) {
      if (!passSASelection(mu)) continue;
      sa.push_back(&mu);
    }

    // OS pairs -> invariant mass
    for (size_t i = 0; i < sa.size(); ++i) {
      for (size_t j = i + 1; j < sa.size(); ++j) {
        if (sa[i]->charge() * sa[j]->charge() >= 0) continue;

        TLorentzVector p1, p2;
        p1.SetPtEtaPhiM(sa[i]->pt(), sa[i]->eta(), sa[i]->phi(), 0.105658);
        p2.SetPtEtaPhiM(sa[j]->pt(), sa[j]->eta(), sa[j]->phi(), 0.105658);

        float mll = (p1 + p2).M();

        sa_coll_.push_back(collId);
        sa_mll_.push_back(mll);

        // Save the two muon pT (ordered)
        float pt1 = sa[i]->pt();
        float pt2 = sa[j]->pt();
        float ip1 = sa[i]->dB();
        float ip2 = sa[j]->dB();
        if (pt2 > pt1) std::swap(pt1, pt2);
        sa_mu1_pt_.push_back(pt1);
        sa_mu2_pt_.push_back(pt2);
        sa_mu1_ip_.push_back(ip1);
        sa_mu2_ip_.push_back(ip2);

        if( ip1 < 1. || ip2 < 1.) continue;
        sa_mll2_.push_back(mll);

        // Kalman Vertex Fitter
        KalmanVertexFitter kvf(true,true);
        std::vector<reco::TransientTrack> tTracks;

        tTracks.emplace_back(ttBuilder.build( sa[i]->bestTrack() ));
        tTracks.emplace_back(ttBuilder.build( sa[j]->bestTrack() ));

        TransientVertex tv = kvf.vertex(tTracks);
        if (!tv.isValid()) continue;

        //GlobalPoint v = tv.position();
        auto v = tv.position();
        float Lxy = std::hypot(v.x(), v.y());
        float prob = TMath::Prob(tv.totalChiSquared(), tv.degreesOfFreedom());

        sa_Lxy_.push_back(Lxy);
        sa_prob_.push_back(prob);

        // Refitting
        const auto& refittedTracks = tv.refittedTracks();
        if (refittedTracks.size() < 2) return;
        /*
        for (auto const& rt : refittedTracks) {
          GlobalVector p = rt.impactPointState().globalMomentum();
        }
        */
        const double mMu = 0.105658;
        auto m1 = refittedTracks[0].impactPointState().globalMomentum();
        auto m2 = refittedTracks[1].impactPointState().globalMomentum();
        TLorentzVector mu1, mu2;
        mu1.SetXYZM(m1.x(), m1.y(), m1.z(), mMu);
        mu2.SetXYZM(m2.x(), m2.y(), m2.z(), mMu);
        double mass = (mu1 + mu2).M();
        sa_mll3_.push_back(mass);

      }
    }
  };

  if (muons.isValid()) {
    processMuonCollection(*muons, 0);
  }
  if (disMuons.isValid()) {
    processMuonCollection(*disMuons, 1);
  }

  tree_->Fill();
}



// -------------------- Register Module --------------------
DEFINE_FWK_MODULE(staMuon);
