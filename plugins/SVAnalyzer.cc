// -*- C++ -*-
// Package:    MyAnalysis/SVAnalyzer
// Class:      SVAnalyzer

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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "TTree.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TMath.h"

#include <vector>
#include <cmath>

class SVAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SVAnalyzer(const edm::ParameterSet&);
  ~SVAnalyzer() override {}

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  // Tokens
  edm::EDGetTokenT<std::vector<pat::Jet>> jetToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate>> svToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;

  // Tree
  TTree* tree_;

  // Branch variables
  std::vector<float> lxy_SV_;
  std::vector<float> dR_ptl_;
  std::vector<std::vector<int>> mom_ptl_;
  std::vector<float> invm_ptl_;
  std::vector<float> mass_SV_;
  std::vector<std::vector<float>> pt_trk_;
  std::vector<float> pt_SV_;
  std::vector<std::vector<float>> eta_trk_;
  std::vector<float> normChi2_SV_;
  std::vector<std::vector<float>> normChi2_trk_;
  std::vector<std::vector<float>> dxySig_trk_;
  std::vector<std::vector<int>> hits_trk_;
  std::vector<std::vector<float>> ip2d_mu_;
  std::vector<float> vtxProb_SV_;
  std::vector<float> dotProd_ptl_;
  std::vector<float> dotProd_SV_;
  std::vector<int> partFlav_jet_;

  // Helpers
  float deltaR(const reco::Candidate::LorentzVector& v1,
               const reco::Candidate::LorentzVector& v2);

  void clearVectors();

  // PV selection struct + function
  struct BestPVInfo {
    const reco::Vertex* pv;
    float angle;
    float lxy;
  };

  BestPVInfo findBestPV(const std::vector<reco::Vertex>& pvs,
                        const reco::VertexCompositePtrCandidate& sv);
};



// -------------------- Constructor --------------------
SVAnalyzer::SVAnalyzer(const edm::ParameterSet& iConfig) {
  usesResource("TFileService");

  jetToken_  = consumes<std::vector<pat::Jet>>( iConfig.getParameter<edm::InputTag>("jets") );
  pvToken_   = consumes<std::vector<reco::Vertex>>( iConfig.getParameter<edm::InputTag>("vertices") );
  svToken_   = consumes<std::vector<reco::VertexCompositePtrCandidate>>( iConfig.getParameter<edm::InputTag>("secondaryVertices") );
  muonToken_ = consumes<std::vector<pat::Muon>>( iConfig.getParameter<edm::InputTag>("muons") );
}



// -------------------- beginJob --------------------
void SVAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "Secondary Vertex Analysis");

  tree_->Branch("lxy_SV",       &lxy_SV_);
  tree_->Branch("dR_ptl",       &dR_ptl_);
  tree_->Branch("mom_ptl",      &mom_ptl_);
  tree_->Branch("invm_ptl",     &invm_ptl_);
  tree_->Branch("mass_SV",      &mass_SV_);
  tree_->Branch("pt_trk",       &pt_trk_);
  tree_->Branch("pt_SV",        &pt_SV_);
  tree_->Branch("eta_trk",      &eta_trk_);
  tree_->Branch("normChi2_SV",  &normChi2_SV_);
  tree_->Branch("normChi2_trk", &normChi2_trk_);
  tree_->Branch("dxySig_trk",   &dxySig_trk_);
  tree_->Branch("hits_trk",     &hits_trk_);
  tree_->Branch("ip2d_mu",      &ip2d_mu_);
  tree_->Branch("vtxProb_SV",   &vtxProb_SV_);
  tree_->Branch("dotProd_ptl",  &dotProd_ptl_);
  tree_->Branch("dotProd_SV",   &dotProd_SV_);
  tree_->Branch("partFlav_jet", &partFlav_jet_);
}



// -------------------- clearVectors --------------------
void SVAnalyzer::clearVectors() {
  lxy_SV_.clear();
  dR_ptl_.clear();
  mom_ptl_.clear();
  invm_ptl_.clear();
  mass_SV_.clear();
  pt_trk_.clear();
  pt_SV_.clear();
  eta_trk_.clear();
  normChi2_SV_.clear();
  normChi2_trk_.clear();
  dxySig_trk_.clear();
  hits_trk_.clear();
  ip2d_mu_.clear();
  vtxProb_SV_.clear();
  dotProd_ptl_.clear();
  dotProd_SV_.clear();
  partFlav_jet_.clear();
}



// -------------------- deltaR --------------------
float SVAnalyzer::deltaR(const reco::Candidate::LorentzVector& v1,
                         const reco::Candidate::LorentzVector& v2) {
  float deta = v1.eta() - v2.eta();
  float dphi = TVector2::Phi_mpi_pi(v1.phi() - v2.phi());
  return std::sqrt(deta*deta + dphi*dphi);
}



// -------------------- findBestPV --------------------
SVAnalyzer::BestPVInfo
SVAnalyzer::findBestPV(const std::vector<reco::Vertex>& pvs,
                       const reco::VertexCompositePtrCandidate& sv)
{
  TVector3 svMom(sv.px(), sv.py(), sv.pz());
  TVector3 svPos(sv.vx(), sv.vy(), sv.vz());

  float bestAngle = 9999.f;
  float bestLxy   = -1.f;
  const reco::Vertex* bestPV = nullptr;

  for (const auto& pv : pvs) {
    TVector3 pvPos(pv.x(), pv.y(), pv.z());
    TVector3 diff = svPos - pvPos;
    float angle = svMom.Angle(diff);

    if (std::abs(angle) < std::abs(bestAngle)) {
      bestAngle = angle;
      bestLxy   = diff.Perp();
      bestPV    = &pv;
    }
  }

  return {bestPV, bestAngle, bestLxy};
}



// -------------------- analyze --------------------
void SVAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {

  clearVectors();

  edm::Handle<std::vector<pat::Jet>> jets;
  edm::Handle<std::vector<reco::Vertex>> vertices;
  edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> svs;
  edm::Handle<std::vector<pat::Muon>> muons;

  iEvent.getByToken(jetToken_,  jets);
  iEvent.getByToken(pvToken_,   vertices);
  iEvent.getByToken(svToken_,   svs);
  iEvent.getByToken(muonToken_, muons);

  if (!jets.isValid() || !vertices.isValid() || !svs.isValid() || !muons.isValid())
    return;
  if (vertices->empty()) return;


  for (const auto& jet : *jets) {

    for (const auto& sv : *svs) {

      if (deltaR(sv.p4(), jet.p4()) > 0.4) continue;

      size_t nSrcCand = sv.numberOfSourceCandidatePtrs();
      if (nSrcCand != 2) continue;

      std::vector<const reco::Candidate*> ptls;

      for (size_t i = 0; i < nSrcCand; ++i) {
        const reco::Candidate* cand = sv.sourceCandidatePtr(i).get();
        if (std::abs(cand->pdgId()) == 13)
          ptls.push_back(cand);
      }

      if (ptls.size() != 2) continue;
      if (ptls[0]->charge() * ptls[1]->charge() >= 0) continue;

      // Match to PAT muons
      std::vector<const pat::Muon*> matchedMuons(2, nullptr);
      std::vector<const reco::Track*> trks(2, nullptr);

      for (int i = 0; i < 2; ++i) {
        float bestDR = 999.;
        for (const auto& mu : *muons) {
          float dr = deltaR(mu.p4(), ptls[i]->p4());
          if (dr < bestDR) {
            bestDR = dr;
            matchedMuons[i] = &mu;
          }
        }
        if (matchedMuons[i] && matchedMuons[i]->bestTrack())
          trks[i] = matchedMuons[i]->bestTrack();
      }

      if (!trks[0] || !trks[1]) continue;

      // Choose best PV for this SV
      auto info = findBestPV(*vertices, sv);
      float anglePV = info.angle;
      float lxyPV   = info.lxy;

      // Dimuon kinematics
      reco::Candidate::LorentzVector dimuon = ptls[0]->p4() + ptls[1]->p4();
      float invm = dimuon.M();
      TVector3 dimuVec(dimuon.px(), dimuon.py(), dimuon.pz());

      float massSV = sv.mass();
      float pSV    = sv.p();
      float corrMassSV =
        std::sqrt(massSV*massSV + pSV*pSV*std::sin(anglePV)*std::sin(anglePV))
        + pSV*std::sin(anglePV);

      // Fill branches
      lxy_SV_.push_back(lxyPV);
      dR_ptl_.push_back(deltaR(ptls[0]->p4(), ptls[1]->p4()));
      pt_SV_.push_back(sv.pt());
      mass_SV_.push_back(corrMassSV);
      invm_ptl_.push_back(invm);
      normChi2_SV_.push_back(sv.vertexNormalizedChi2());
      vtxProb_SV_.push_back(TMath::Prob(sv.vertexChi2(), (int)sv.vertexNdof()));
      dotProd_ptl_.push_back(std::cos(dimuVec.Angle(TVector3(sv.vx(), sv.vy(), sv.vz()))));
      dotProd_SV_.push_back(std::cos(anglePV));
      partFlav_jet_.push_back(jet.partonFlavour());

      // Track-level info
      std::vector<int> moms;
      std::vector<float> pts, etas, chi2, dxysig, ip2d;
      std::vector<int> hits;

      for (int i = 0; i < 2; ++i) {
        moms.push_back(matchedMuons[i]->simMotherPdgId());
        pts.push_back(trks[i]->pt());
        etas.push_back(trks[i]->eta());
        chi2.push_back(trks[i]->normalizedChi2());
        dxysig.push_back(trks[i]->dxy(info.pv->position()) / trks[i]->dxyError());
        hits.push_back(trks[i]->numberOfValidHits());
        ip2d.push_back(matchedMuons[i]->dB(pat::Muon::PV2D));
      }

      mom_ptl_.push_back(moms);
      pt_trk_.push_back(pts);
      eta_trk_.push_back(etas);
      normChi2_trk_.push_back(chi2);
      dxySig_trk_.push_back(dxysig);
      hits_trk_.push_back(hits);
      ip2d_mu_.push_back(ip2d);
    }
  }

  tree_->Fill();
}



// -------------------- Register Module --------------------
DEFINE_FWK_MODULE(SVAnalyzer);
