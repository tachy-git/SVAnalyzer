// -*- C++ -*-

#include <memory>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <cstdlib>
#include <utility>

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
    int dir, cscHit, dtHit;

    static MuonVars from(const pat::Muon& mu) {
      MuonVars v;
      v.pt = mu.pt();
      v.eta = mu.eta();
      v.phi = mu.phi();
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
    float pt, eta, phi, relPtUnc;

    static RefitVars from(const reco::TransientTrack& trk) {
      const reco::Track& t = trk.track();
      RefitVars v;
      v.pt = t.pt();
      v.eta = t.eta();
      v.phi = t.phi();
      v.relPtUnc = (t.pt() > 0.f ? t.ptError() / t.pt() : kInvalid);
      return v;
    }
  };

  struct GenMuonVars {
    float pt, eta, phi;

    static GenMuonVars from(const reco::GenParticle& gp) {
      GenMuonVars v;
      v.pt = gp.pt();
      v.eta = gp.eta();
      v.phi = gp.phi();
      return v;
    }
  };

  template <typename T>
  static std::pair<const T*, const T*> orderedByCharge(const T* a, const T* b) {
    if (!a || !b) return {nullptr, nullptr};
    const int qa = a->charge();
    const int qb = b->charge();
    if (qa * qb >= 0) return {nullptr, nullptr};
    return (qa > 0) ? std::make_pair(a, b) : std::make_pair(b, a);
  }

  static bool passSAMuonSelection(const pat::Muon& mu) {
    if (!mu.isStandAloneMuon()) return false;
    if (mu.isTrackerMuon()) return false;
    if (!mu.bestTrack()) return false;
    if (mu.pt() < 8.) return false;
    return true;
  }

  void pushMuonPair(const MuonVars& vPos, const MuonVars& vNeg);
  void pushGenMuonPair(const GenMuonVars& vPos, const GenMuonVars& vNeg);

  edm::EDGetTokenT<std::vector<pat::Muon>> disMuonToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  TTree* tree_;

  int nDisMuon_;
  int nGenMuon_;

  float invm_gen_mu_;

  std::vector<float> invm_mu_bef_;
  std::vector<float> invm_mu_aft_;
  std::vector<float> Lxy_mu_;
  std::vector<float> LxyErr_mu_;
  std::vector<float> normChi_mu_vtx_;
  std::vector<float> vtxProb_mu_;
  std::vector<int>   isValid_mu_vtx_;
  std::vector<float> dR_mu_;
  std::vector<float> dca3d_mu_;
  std::vector<float> dR_gen_;
  std::vector<float> dR_refit_;

  std::vector<std::vector<float>> pt_mu_;
  std::vector<std::vector<float>> eta_mu_;
  std::vector<std::vector<float>> phi_mu_;
  std::vector<std::vector<float>> ip_mu_;
  std::vector<std::vector<float>> time_mu_;
  std::vector<std::vector<int>>   dir_mu_;
  std::vector<std::vector<int>>   cscHit_mu_;
  std::vector<std::vector<int>>   dtHit_mu_;
  std::vector<std::vector<float>> normChi_mu_;

  std::vector<std::vector<float>> pt_refit_mu_;
  std::vector<std::vector<float>> eta_refit_mu_;
  std::vector<std::vector<float>> phi_refit_mu_;
  std::vector<std::vector<float>> relPtUnc_refit_mu_;

  std::vector<std::vector<float>> gen_pt_mu_;
  std::vector<std::vector<float>> gen_eta_mu_;
  std::vector<std::vector<float>> gen_phi_mu_;
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
  tree_->Branch("nGenMuon", &nGenMuon_, "nGenMuon/I");
  tree_->Branch("invm_gen_mu", &invm_gen_mu_, "invm_gen_mu/F");

  tree_->Branch("invm_mu_bef", &invm_mu_bef_);
  tree_->Branch("invm_mu_aft", &invm_mu_aft_);
  tree_->Branch("Lxy_mu", &Lxy_mu_);
  tree_->Branch("LxyErr_mu", &LxyErr_mu_);
  tree_->Branch("normChi_mu_vtx", &normChi_mu_vtx_);
  tree_->Branch("vtxProb_mu", &vtxProb_mu_);
  tree_->Branch("isValid_mu_vtx", &isValid_mu_vtx_);
  tree_->Branch("dR_mu", &dR_mu_);
  tree_->Branch("dca3d_mu", &dca3d_mu_);
  tree_->Branch("dR_gen", &dR_gen_);
  tree_->Branch("dR_refit", &dR_refit_);

  tree_->Branch("pt_mu", &pt_mu_);
  tree_->Branch("eta_mu", &eta_mu_);
  tree_->Branch("phi_mu", &phi_mu_);
  tree_->Branch("ip_mu", &ip_mu_);
  tree_->Branch("time_mu", &time_mu_);
  tree_->Branch("dir_mu", &dir_mu_);
  tree_->Branch("cscHit_mu", &cscHit_mu_);
  tree_->Branch("dtHit_mu", &dtHit_mu_);
  tree_->Branch("normChi_mu", &normChi_mu_);

  tree_->Branch("pt_refit_mu", &pt_refit_mu_);
  tree_->Branch("eta_refit_mu", &eta_refit_mu_);
  tree_->Branch("phi_refit_mu", &phi_refit_mu_);
  tree_->Branch("relPtUnc_refit_mu", &relPtUnc_refit_mu_);

  tree_->Branch("gen_pt_mu", &gen_pt_mu_);
  tree_->Branch("gen_eta_mu", &gen_eta_mu_);
  tree_->Branch("gen_phi_mu", &gen_phi_mu_);
}

void staMuon::clearVectors() {
  nDisMuon_ = 0;
  nGenMuon_ = 0;
  invm_gen_mu_ = kInvalid;

  invm_mu_bef_.clear();
  invm_mu_aft_.clear();
  Lxy_mu_.clear();
  LxyErr_mu_.clear();
  normChi_mu_vtx_.clear();
  vtxProb_mu_.clear();
  isValid_mu_vtx_.clear();
  dR_mu_.clear();
  dca3d_mu_.clear();
  dR_gen_.clear();
  dR_refit_.clear();

  pt_mu_.clear();
  eta_mu_.clear();
  phi_mu_.clear();
  ip_mu_.clear();
  time_mu_.clear();
  dir_mu_.clear();
  cscHit_mu_.clear();
  dtHit_mu_.clear();
  normChi_mu_.clear();

  pt_refit_mu_.clear();
  eta_refit_mu_.clear();
  phi_refit_mu_.clear();
  relPtUnc_refit_mu_.clear();

  gen_pt_mu_.clear();
  gen_eta_mu_.clear();
  gen_phi_mu_.clear();
}

void staMuon::pushMuonPair(const MuonVars& vPos, const MuonVars& vNeg) {
  pt_mu_.push_back({vPos.pt, vNeg.pt});
  eta_mu_.push_back({vPos.eta, vNeg.eta});
  phi_mu_.push_back({vPos.phi, vNeg.phi});
  ip_mu_.push_back({vPos.ip, vNeg.ip});
  time_mu_.push_back({vPos.time, vNeg.time});
  dir_mu_.push_back({vPos.dir, vNeg.dir});
  cscHit_mu_.push_back({vPos.cscHit, vNeg.cscHit});
  dtHit_mu_.push_back({vPos.dtHit, vNeg.dtHit});
  normChi_mu_.push_back({vPos.normChi, vNeg.normChi});
}

void staMuon::pushGenMuonPair(const GenMuonVars& vPos, const GenMuonVars& vNeg) {
  gen_pt_mu_.push_back({vPos.pt, vNeg.pt});
  gen_eta_mu_.push_back({vPos.eta, vNeg.eta});
  gen_phi_mu_.push_back({vPos.phi, vNeg.phi});
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

  std::vector<const reco::GenParticle*> genMuPtrs;
  if (genParticles.isValid()) {
    for (const auto& gp : *genParticles) {
      if (std::abs(gp.pdgId()) != 13) continue;
      if (gp.status() != 1) continue;
      if (!gp.isLastCopy()) continue;
      genMuPtrs.push_back(&gp);
    }
  }
  nGenMuon_ = static_cast<int>(genMuPtrs.size());

  if (!vertices.isValid() || vertices->empty()) {
    return;
  }

  if (!disMuons.isValid()) {
    return;
  }

  nDisMuon_ = static_cast<int>(disMuons->size());
  if (nDisMuon_ < 2) {
    return;
  }

  const auto& ttBuilder = iSetup.getData(ttBuilderToken_);

  const reco::Vertex& pv = vertices->at(0);
  const reco::Vertex::Error& ve = pv.error();
  const GlobalError pvErr(ve.At(0,0), ve.At(1,0), ve.At(1,1),
                          ve.At(2,0), ve.At(2,1), ve.At(2,2));

  const reco::GenParticle* genPos = nullptr;
  const reco::GenParticle* genNeg = nullptr;
  for (size_t i = 0; i < genMuPtrs.size(); ++i) {
    for (size_t j = i + 1; j < genMuPtrs.size(); ++j) {
      const auto ordered = orderedByCharge(genMuPtrs[i], genMuPtrs[j]);
      if (ordered.first && ordered.second) {
        genPos = ordered.first;
        genNeg = ordered.second;
        break;
      }
    }
    if (genPos && genNeg) break;
  }

  if (genPos && genNeg) {
    TLorentzVector pGenPos, pGenNeg;
    pGenPos.SetPtEtaPhiM(genPos->pt(), genPos->eta(), genPos->phi(), mMu);
    pGenNeg.SetPtEtaPhiM(genNeg->pt(), genNeg->eta(), genNeg->phi(), mMu);
    invm_gen_mu_ = (pGenPos + pGenNeg).M();
  } else {
    invm_gen_mu_ = kInvalid;
  }

  std::vector<const pat::Muon*> saMuPtrs;
  saMuPtrs.reserve(disMuons->size());

  for (const auto& mu : *disMuons) {
    if (!passSAMuonSelection(mu)) continue;
    saMuPtrs.push_back(&mu);
  }

  for (size_t i = 0; i < saMuPtrs.size(); ++i) {
    for (size_t j = i + 1; j < saMuPtrs.size(); ++j) {
      const auto orderedReco = orderedByCharge(saMuPtrs[i], saMuPtrs[j]);
      const pat::Muon* muPos = orderedReco.first;
      const pat::Muon* muNeg = orderedReco.second;
      if (!muPos || !muNeg) continue;

      TLorentzVector pPos_bef, pNeg_bef;
      pPos_bef.SetPtEtaPhiM(muPos->pt(), muPos->eta(), muPos->phi(), mMu);
      pNeg_bef.SetPtEtaPhiM(muNeg->pt(), muNeg->eta(), muNeg->phi(), mMu);
      invm_mu_bef_.push_back((pPos_bef + pNeg_bef).M());

      dR_mu_.push_back(static_cast<float>(
          reco::deltaR(muPos->eta(), muPos->phi(), muNeg->eta(), muNeg->phi())));

      const auto mvPos = MuonVars::from(*muPos);
      const auto mvNeg = MuonVars::from(*muNeg);
      pushMuonPair(mvPos, mvNeg);

      reco::TransientTrack ttPos = ttBuilder.build(muPos->bestTrack());
      reco::TransientTrack ttNeg = ttBuilder.build(muNeg->bestTrack());

      {
        float dca3d = kInvalid;
        TwoTrackMinimumDistance ttmd;
        if (ttmd.calculate(ttPos.initialFreeState(), ttNeg.initialFreeState())) {
          dca3d = static_cast<float>(ttmd.distance());
        }
        dca3d_mu_.push_back(dca3d);
      }

      if (genPos && genNeg) {
        pushGenMuonPair(GenMuonVars::from(*genPos), GenMuonVars::from(*genNeg));

        dR_gen_.push_back(static_cast<float>(
            reco::deltaR(genPos->eta(), genPos->phi(), genNeg->eta(), genNeg->phi())));
      } else {
        gen_pt_mu_.push_back({kInvalid, kInvalid});
        gen_eta_mu_.push_back({kInvalid, kInvalid});
        gen_phi_mu_.push_back({kInvalid, kInvalid});
        dR_gen_.push_back(kInvalid);
      }

      KalmanVertexFitter kvf(true, true);
      std::vector<reco::TransientTrack> tracks = {ttPos, ttNeg};
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
        phi_refit_mu_.push_back({kInvalid, kInvalid});
        relPtUnc_refit_mu_.push_back({kInvalid, kInvalid});
        dR_refit_.push_back(kInvalid);
        continue;
      }

      const auto& refitTrks = tv.refittedTracks();
      const auto orderedRefit = orderedByCharge(&refitTrks[0], &refitTrks[1]);
      const reco::TransientTrack* refitPos = orderedRefit.first;
      const reco::TransientTrack* refitNeg = orderedRefit.second;

      if (!refitPos || !refitNeg) {
        isValid_mu_vtx_.push_back(0);
        invm_mu_aft_.push_back(kInvalid);
        Lxy_mu_.push_back(kInvalid);
        LxyErr_mu_.push_back(kInvalid);
        normChi_mu_vtx_.push_back(kInvalid);
        vtxProb_mu_.push_back(kInvalid);

        pt_refit_mu_.push_back({kInvalid, kInvalid});
        eta_refit_mu_.push_back({kInvalid, kInvalid});
        phi_refit_mu_.push_back({kInvalid, kInvalid});
        relPtUnc_refit_mu_.push_back({kInvalid, kInvalid});
        dR_refit_.push_back(kInvalid);
        continue;
      }

      const auto gmPos = refitPos->impactPointState().globalMomentum();
      const auto gmNeg = refitNeg->impactPointState().globalMomentum();

      const float refitEtaPos = static_cast<float>(gmPos.eta());
      const float refitPhiPos = static_cast<float>(gmPos.phi());
      const float refitEtaNeg = static_cast<float>(gmNeg.eta());
      const float refitPhiNeg = static_cast<float>(gmNeg.phi());

      TLorentzVector pPos_aft, pNeg_aft;
      pPos_aft.SetXYZM(gmPos.x(), gmPos.y(), gmPos.z(), mMu);
      pNeg_aft.SetXYZM(gmNeg.x(), gmNeg.y(), gmNeg.z(), mMu);

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
      invm_mu_aft_.push_back((pPos_aft + pNeg_aft).M());
      Lxy_mu_.push_back(disp.perp());
      LxyErr_mu_.push_back(std::sqrt(totErr.rerr(disp)));
      normChi_mu_vtx_.push_back(normChi);
      vtxProb_mu_.push_back(vtxProb);

      const auto rvPos = RefitVars::from(*refitPos);
      const auto rvNeg = RefitVars::from(*refitNeg);

      pt_refit_mu_.push_back({rvPos.pt, rvNeg.pt});
      eta_refit_mu_.push_back({rvPos.eta, rvNeg.eta});
      phi_refit_mu_.push_back({rvPos.phi, rvNeg.phi});
      relPtUnc_refit_mu_.push_back({rvPos.relPtUnc, rvNeg.relPtUnc});

      dR_refit_.push_back(static_cast<float>(
          reco::deltaR(refitEtaPos, refitPhiPos, refitEtaNeg, refitPhiNeg)));
    }
  }

  tree_->Fill();
}

DEFINE_FWK_MODULE(staMuon);
