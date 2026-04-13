// -*- C++ -*-

#include <memory>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <cstdlib>

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/ESInputTag.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"

class staMuon : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit staMuon(const edm::ParameterSet&);
  ~staMuon() override {}

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void clearVectors();

  static constexpr float kInvalid = -999.f;

  struct MuonVars {
    float pt, eta, phi, ip, time, normChi;
    int charge, dir, cscHit, dtHit;

    static MuonVars from(const pat::Muon& mu) {
      MuonVars v;
      v.pt = mu.pt();
      v.eta = mu.eta();
      v.phi = mu.phi();
      v.charge = mu.charge();
      v.ip = mu.dB();
      v.time = mu.time().timeAtIpInOut;
      v.dir = mu.time().direction();

      if (mu.bestTrack()) {
        const auto& hp = mu.bestTrack()->hitPattern();
        v.cscHit = hp.numberOfValidMuonCSCHits();
        v.dtHit = hp.numberOfValidMuonDTHits();
        v.normChi = mu.bestTrack()->normalizedChi2();
      } else {
        v.cscHit = static_cast<int>(kInvalid);
        v.dtHit = static_cast<int>(kInvalid);
        v.normChi = kInvalid;
      }
      return v;
    }
  };

  struct RefitVars {
    float pt, eta, relPtUnc;

    static RefitVars from(const reco::TransientTrack& trk) {
      const reco::Track& t = trk.track();
      RefitVars v;
      v.pt = t.pt();
      v.eta = t.eta();
      v.relPtUnc = (t.pt() > 0.f ? t.ptError() / t.pt() : kInvalid);
      return v;
    }
  };

  struct GenMuonVars {
    float pt, eta, phi;
    int charge;

    static GenMuonVars from(const reco::GenParticle& gp) {
      GenMuonVars v;
      v.pt = gp.pt();
      v.eta = gp.eta();
      v.phi = gp.phi();
      v.charge = gp.charge();
      return v;
    }
  };

  static bool passSAMuonSelection(const pat::Muon& mu) {
    if (!mu.isStandAloneMuon()) return false;
    if (mu.isTrackerMuon()) return false;
    if (!mu.bestTrack()) return false;
    if (mu.pt() < 8.) return false;
    return true;
  }

  void pushMuonPair(const MuonVars& v1, const MuonVars& v2);
  void pushGenMuonPair(const GenMuonVars& v1, const GenMuonVars& v2);

  edm::EDGetTokenT<std::vector<pat::Muon>> disMuonToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  TTree* tree_;

  // event-level multiplicities
  int nDisMuon_;
  int nGenMuon_;

  // displacedMuon collection, standalone muons only
  std::vector<float> invm_mu_bef_;
  std::vector<float> invm_mu_aft_;
  std::vector<float> Lxy_mu_;
  std::vector<float> LxyErr_mu_;
  std::vector<float> normChi_mu_vtx_;
  std::vector<float> vtxProb_mu_;
  std::vector<int>   isValid_mu_vtx_;
  std::vector<float> dR_mu_;
  std::vector<float> dca3d_mu_;

  // reco muon kinematics
  std::vector<std::vector<float>> pt_mu_;
  std::vector<std::vector<float>> eta_mu_;
  std::vector<std::vector<float>> phi_mu_;
  std::vector<std::vector<int>>   charge_mu_;
  std::vector<std::vector<float>> ip_mu_;
  std::vector<std::vector<float>> time_mu_;
  std::vector<std::vector<int>>   dir_mu_;
  std::vector<std::vector<int>>   cscHit_mu_;
  std::vector<std::vector<int>>   dtHit_mu_;
  std::vector<std::vector<float>> normChi_mu_;

  // refit kinematics
  std::vector<std::vector<float>> pt_refit_mu_;
  std::vector<std::vector<float>> eta_refit_mu_;
  std::vector<std::vector<float>> relPtUnc_refit_mu_;

  // gen muon pair kinematics
  std::vector<std::vector<float>> gen_pt_mu_;
  std::vector<std::vector<float>> gen_eta_mu_;
  std::vector<std::vector<float>> gen_phi_mu_;
  std::vector<std::vector<int>>   gen_charge_mu_;
  std::vector<float> dR_gen_mu_;
};

constexpr float staMuon::kInvalid;

staMuon::staMuon(const edm::ParameterSet& iConfig) {
  usesResource("TFileService");

  disMuonToken_ =
      consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("disMuons"));
  pvToken_ =
      consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"));
  genParticleToken_ =
      consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));

  ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(
      edm::ESInputTag("", "TransientTrackBuilder"));
}

void staMuon::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "staMuon tree");

  tree_->Branch("nDisMuon", &nDisMuon_, "nDisMuon/I");

  tree_->Branch("invm_mu_bef", &invm_mu_bef_);
  tree_->Branch("invm_mu_aft", &invm_mu_aft_);
  tree_->Branch("Lxy_mu", &Lxy_mu_);
  tree_->Branch("LxyErr_mu", &LxyErr_mu_);
  tree_->Branch("normChi_mu_vtx", &normChi_mu_vtx_);
  tree_->Branch("vtxProb_mu", &vtxProb_mu_);
  tree_->Branch("isValid_mu_vtx", &isValid_mu_vtx_);
  tree_->Branch("dR_mu", &dR_mu_);
  tree_->Branch("dca3d_mu", &dca3d_mu_);

  tree_->Branch("pt_mu", &pt_mu_);
  tree_->Branch("eta_mu", &eta_mu_);
  tree_->Branch("phi_mu", &phi_mu_);
  tree_->Branch("charge_mu", &charge_mu_);
  tree_->Branch("ip_mu", &ip_mu_);
  tree_->Branch("time_mu", &time_mu_);
  tree_->Branch("dir_mu", &dir_mu_);
  tree_->Branch("cscHit_mu", &cscHit_mu_);
  tree_->Branch("dtHit_mu", &dtHit_mu_);
  tree_->Branch("normChi_mu", &normChi_mu_);

  tree_->Branch("pt_refit_mu", &pt_refit_mu_);
  tree_->Branch("eta_refit_mu", &eta_refit_mu_);
  tree_->Branch("relPtUnc_refit_mu", &relPtUnc_refit_mu_);

  tree_->Branch("gen_pt_mu", &gen_pt_mu_);
  tree_->Branch("gen_eta_mu", &gen_eta_mu_);
  tree_->Branch("gen_phi_mu", &gen_phi_mu_);
  tree_->Branch("gen_charge_mu", &gen_charge_mu_);
  tree_->Branch("dR_gen_mu", &dR_gen_mu_);
}

void staMuon::clearVectors() {
  nDisMuon_ = 0;

  invm_mu_bef_.clear();
  invm_mu_aft_.clear();
  Lxy_mu_.clear();
  LxyErr_mu_.clear();
  normChi_mu_vtx_.clear();
  vtxProb_mu_.clear();
  isValid_mu_vtx_.clear();
  dR_mu_.clear();
  dca3d_mu_.clear();

  pt_mu_.clear();
  eta_mu_.clear();
  phi_mu_.clear();
  charge_mu_.clear();
  ip_mu_.clear();
  time_mu_.clear();
  dir_mu_.clear();
  cscHit_mu_.clear();
  dtHit_mu_.clear();
  normChi_mu_.clear();

  pt_refit_mu_.clear();
  eta_refit_mu_.clear();
  relPtUnc_refit_mu_.clear();

  gen_pt_mu_.clear();
  gen_eta_mu_.clear();
  gen_phi_mu_.clear();
  gen_charge_mu_.clear();
  dR_gen_mu_.clear();
}

void staMuon::pushMuonPair(const MuonVars& v1, const MuonVars& v2) {
  pt_mu_.push_back({v1.pt, v2.pt});
  eta_mu_.push_back({v1.eta, v2.eta});
  phi_mu_.push_back({v1.phi, v2.phi});
  charge_mu_.push_back({v1.charge, v2.charge});
  ip_mu_.push_back({v1.ip, v2.ip});
  time_mu_.push_back({v1.time, v2.time});
  dir_mu_.push_back({v1.dir, v2.dir});
  cscHit_mu_.push_back({v1.cscHit, v2.cscHit});
  dtHit_mu_.push_back({v1.dtHit, v2.dtHit});
  normChi_mu_.push_back({v1.normChi, v2.normChi});
}

void staMuon::pushGenMuonPair(const GenMuonVars& v1, const GenMuonVars& v2) {
  gen_pt_mu_.push_back({v1.pt, v2.pt});
  gen_eta_mu_.push_back({v1.eta, v2.eta});
  gen_phi_mu_.push_back({v1.phi, v2.phi});
  gen_charge_mu_.push_back({v1.charge, v2.charge});
}

void staMuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  clearVectors();

  edm::Handle<std::vector<pat::Muon>> disMuons;
  edm::Handle<std::vector<reco::Vertex>> vertices;
  edm::Handle<reco::GenParticleCollection> genParticles;

  iEvent.getByToken(disMuonToken_, disMuons);
  iEvent.getByToken(pvToken_, vertices);
  iEvent.getByToken(genParticleToken_, genParticles);

  constexpr float mMu = 0.105658f;

  // ============================================================
  // 1. GEN muons first
  // ============================================================
  std::vector<const reco::GenParticle*> genMuPtrs;
  if (genParticles.isValid()) {
    for (const auto& gp : *genParticles) {
      if (std::abs(gp.pdgId()) != 13) continue;
      if (gp.status() != 1) continue;
      if (!gp.isLastCopy()) continue;
      genMuPtrs.push_back(&gp);
    }
  }

  // ============================================================
  // 2. Need PV for vertex-related quantities
  // ============================================================
  if (!vertices.isValid() || vertices->empty()) {
    return;
  }

  const auto& ttBuilder = iSetup.getData(ttBuilderToken_);

  const reco::Vertex& pv = vertices->at(0);
  const reco::Vertex::Error& ve = pv.error();
  const GlobalError pvErr(ve.At(0,0), ve.At(1,0), ve.At(1,1),
                          ve.At(2,0), ve.At(2,1), ve.At(2,2));

  // ============================================================
  // 4. displacedMuon collection only
  //    standaloneMuon && !trackerMuon
  // ============================================================
  if (disMuons.isValid()) {
    nDisMuon_ = static_cast<int>(disMuons->size());
    std::vector<const pat::Muon*> saMuPtrs;
    saMuPtrs.reserve(disMuons->size());

    for (const auto& mu : *disMuons) {
      if (!passSAMuonSelection(mu)) continue;
      saMuPtrs.push_back(&mu);
    }

    for (size_t i = 0; i < saMuPtrs.size(); ++i) {
      for (size_t j = i + 1; j < saMuPtrs.size(); ++j) {
        const pat::Muon* mu1 = saMuPtrs[i];
        const pat::Muon* mu2 = saMuPtrs[j];

        if (mu1->charge() * mu2->charge() >= 0) continue;

        TLorentzVector p1_bef, p2_bef;
        p1_bef.SetPtEtaPhiM(mu1->pt(), mu1->eta(), mu1->phi(), mMu);
        p2_bef.SetPtEtaPhiM(mu2->pt(), mu2->eta(), mu2->phi(), mMu);
        invm_mu_bef_.push_back((p1_bef + p2_bef).M());

        dR_mu_.push_back(reco::deltaR(
            mu1->eta(), mu1->phi(), mu2->eta(), mu2->phi()));

        const auto mv1 = MuonVars::from(*mu1);
        const auto mv2 = MuonVars::from(*mu2);
        pushMuonPair(mv1, mv2);

        reco::TransientTrack tt1 = ttBuilder.build(mu1->bestTrack());
        reco::TransientTrack tt2 = ttBuilder.build(mu2->bestTrack());

        {
          float dca3d = kInvalid;

          TwoTrackMinimumDistance ttmd;
          if (ttmd.calculate(tt1.initialFreeState(), tt2.initialFreeState())) {
            dca3d = static_cast<float>(ttmd.distance());
          }
          dca3d_mu_.push_back(dca3d);
        }

        KalmanVertexFitter kvf(true, true);
        std::vector<reco::TransientTrack> tracks = {tt1, tt2};
        TransientVertex tv = kvf.vertex(tracks);

        if (!tv.isValid() || tv.refittedTracks().size() != 2) {
          isValid_mu_vtx_.push_back(0);
          invm_mu_aft_.push_back(kInvalid);
          Lxy_mu_.push_back(kInvalid);
          LxyErr_mu_.push_back(kInvalid);
          normChi_mu_vtx_.push_back(kInvalid);
          vtxProb_mu_.push_back(kInvalid);

          pt_refit_mu_.push_back({kInvalid, kInvalid});
          eta_refit_mu_.push_back({kInvalid, kInvalid});
          relPtUnc_refit_mu_.push_back({kInvalid, kInvalid});
          continue;
        }

        const auto& refitTrks = tv.refittedTracks();
        const auto gm1 = refitTrks[0].impactPointState().globalMomentum();
        const auto gm2 = refitTrks[1].impactPointState().globalMomentum();

        TLorentzVector p1_aft, p2_aft;
        p1_aft.SetXYZM(gm1.x(), gm1.y(), gm1.z(), mMu);
        p2_aft.SetXYZM(gm2.x(), gm2.y(), gm2.z(), mMu);

        const GlobalPoint svPos = tv.position();
        const GlobalError svErr = tv.positionError();
        const GlobalPoint disp(svPos.x(), svPos.y(), 0.f);
        const GlobalError totErr = svErr + pvErr;

        const float normChi =
            (tv.degreesOfFreedom() > 0 ? tv.totalChiSquared() / tv.degreesOfFreedom() : kInvalid);

        float vtxProb = kInvalid;
        if (tv.degreesOfFreedom() > 0) {
          vtxProb = static_cast<float>(
              TMath::Prob(tv.totalChiSquared(), static_cast<int>(tv.degreesOfFreedom())));
        }

        isValid_mu_vtx_.push_back(1);
        invm_mu_aft_.push_back((p1_aft + p2_aft).M());
        Lxy_mu_.push_back(disp.perp());
        LxyErr_mu_.push_back(std::sqrt(totErr.rerr(disp)));
        normChi_mu_vtx_.push_back(normChi);
        vtxProb_mu_.push_back(vtxProb);

        const auto rv1 = RefitVars::from(refitTrks[0]);
        const auto rv2 = RefitVars::from(refitTrks[1]);

        pt_refit_mu_.push_back({rv1.pt, rv2.pt});
        eta_refit_mu_.push_back({rv1.eta, rv2.eta});
        relPtUnc_refit_mu_.push_back({rv1.relPtUnc, rv2.relPtUnc});
      }
    }
  }

  tree_->Fill();
}

DEFINE_FWK_MODULE(staMuon);
