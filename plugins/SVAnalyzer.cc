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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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
  edm::EDGetTokenT<std::vector<pat::Electron>> electronToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genToken_;

  // Tree
  TTree* tree_;

  // Jet branches
  std::vector<float> jet_pt_;
  std::vector<float> jet_eta_;
  std::vector<float> jet_phi_;
  std::vector<int> jet_partonFlavour_;

  // SV branches
  std::vector<float> sv_mass_;
  std::vector<float> sv_pt_;
  std::vector<float> sv_eta_;
  std::vector<float> sv_phi_;
  std::vector<float> sv_vx_;
  std::vector<float> sv_vy_;
  std::vector<float> sv_lxy_;
  std::vector<float> sv_cosAngle_;
  std::vector<float> sv_chi2_;
  std::vector<float> sv_ndof_;

  // Lepton branches (nested: one vector per SV, each containing 2 leptons)
  std::vector<std::vector<float>> lep_pt_;
  std::vector<std::vector<float>> lep_eta_;
  std::vector<std::vector<float>> lep_phi_;
  std::vector<std::vector<float>> lep_ip2d_;
  
  // Muon-specific branches (only filled if lepton is muon, otherwise -1)
  std::vector<std::vector<int>> lep_isGlobalMuon_;
  std::vector<std::vector<int>> lep_isStandAloneMuon_;
  std::vector<std::vector<int>> lep_simMotherPdgId_;
  
  // GEN matching branches (nested: one vector per SV, each containing 2 leptons)
  std::vector<std::vector<float>> gen_pt_;
  std::vector<std::vector<float>> gen_eta_;
  std::vector<std::vector<float>> gen_phi_;
  std::vector<std::vector<int>> gen_motherPdgId_;

  // Helpers
  float deltaR(const reco::Candidate::LorentzVector& v1,
               const reco::Candidate::LorentzVector& v2);

  void clearVectors();

  // Matching functions
  template<typename T>
  const T* findClosestLepton(const std::vector<T>& leptons, 
                             const reco::Candidate* cand);

  const reco::GenParticle* findClosestGen(const std::vector<reco::GenParticle>& genParticles,
                                          const reco::Candidate* lepton);

  // PV selection struct + function
  struct BestPVInfo {
    //const reco::Vertex* pv;
    float angle;
    float lxy;
  };

  BestPVInfo findBestPV(const std::vector<reco::Vertex>& pvs,
                        const reco::VertexCompositePtrCandidate& sv);
};



// -------------------- Constructor --------------------
SVAnalyzer::SVAnalyzer(const edm::ParameterSet& iConfig) {
  usesResource("TFileService");

  jetToken_      = consumes<std::vector<pat::Jet>>( iConfig.getParameter<edm::InputTag>("jets") );
  pvToken_       = consumes<std::vector<reco::Vertex>>( iConfig.getParameter<edm::InputTag>("vertices") );
  svToken_       = consumes<std::vector<reco::VertexCompositePtrCandidate>>( iConfig.getParameter<edm::InputTag>("secondaryVertices") );
  muonToken_     = consumes<std::vector<pat::Muon>>( iConfig.getParameter<edm::InputTag>("muons") );
  electronToken_ = consumes<std::vector<pat::Electron>>( iConfig.getParameter<edm::InputTag>("electrons") );
  genToken_      = consumes<std::vector<reco::GenParticle>>( iConfig.getParameter<edm::InputTag>("genParticles") );
}



// -------------------- beginJob --------------------
void SVAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "Secondary Vertex Analysis");

  // Jet branches
  tree_->Branch("jet_pt",            &jet_pt_);
  tree_->Branch("jet_eta",           &jet_eta_);
  tree_->Branch("jet_phi",           &jet_phi_);
  tree_->Branch("jet_partonFlavour", &jet_partonFlavour_);

  // SV branches
  tree_->Branch("sv_mass",     &sv_mass_);
  tree_->Branch("sv_pt",       &sv_pt_);
  tree_->Branch("sv_eta",      &sv_eta_);
  tree_->Branch("sv_phi",      &sv_phi_);
  tree_->Branch("sv_vx",       &sv_vx_);
  tree_->Branch("sv_vy",       &sv_vy_);
  tree_->Branch("sv_lxy",      &sv_lxy_);
  tree_->Branch("sv_cosAngle", &sv_cosAngle_);
  tree_->Branch("sv_chi2",     &sv_chi2_);
  tree_->Branch("sv_ndof",     &sv_ndof_);

  // Lepton branches
  tree_->Branch("lep_pt",   &lep_pt_);
  tree_->Branch("lep_eta",  &lep_eta_);
  tree_->Branch("lep_phi",  &lep_phi_);
  tree_->Branch("lep_ip2d", &lep_ip2d_);
  
  // Muon-specific branches
  tree_->Branch("lep_isGlobalMuon",     &lep_isGlobalMuon_);
  tree_->Branch("lep_isStandAloneMuon", &lep_isStandAloneMuon_);
  tree_->Branch("lep_simMotherPdgId",   &lep_simMotherPdgId_);

  // GEN branches
  tree_->Branch("gen_pt",          &gen_pt_);
  tree_->Branch("gen_eta",         &gen_eta_);
  tree_->Branch("gen_phi",         &gen_phi_);
  tree_->Branch("gen_motherPdgId", &gen_motherPdgId_);
}



// -------------------- clearVectors --------------------
void SVAnalyzer::clearVectors() {
  // Jets
  jet_pt_.clear();
  jet_eta_.clear();
  jet_phi_.clear();
  jet_partonFlavour_.clear();

  // SVs
  sv_mass_.clear();
  sv_pt_.clear();
  sv_eta_.clear();
  sv_phi_.clear();
  sv_vx_.clear();
  sv_vy_.clear();
  sv_lxy_.clear();
  sv_cosAngle_.clear();
  sv_chi2_.clear();
  sv_ndof_.clear();

  // Leptons
  lep_pt_.clear();
  lep_eta_.clear();
  lep_phi_.clear();
  lep_ip2d_.clear();
  lep_isGlobalMuon_.clear();
  lep_isStandAloneMuon_.clear();
  lep_simMotherPdgId_.clear();

  // GEN
  gen_pt_.clear();
  gen_eta_.clear();
  gen_phi_.clear();
  gen_motherPdgId_.clear();
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
  //const reco::Vertex* bestPV = nullptr;

  for (const auto& pv : pvs) {
    TVector3 pvPos(pv.x(), pv.y(), pv.z());
    TVector3 diff = svPos - pvPos;
    float angle = svMom.Angle(diff);

    if (std::abs(angle) < std::abs(bestAngle)) {
      bestAngle = angle;
      bestLxy   = diff.Perp();
      //bestPV    = &pv;
    }
  }
  return {bestAngle, bestLxy};
}



// -------------------- findClosestLepton --------------------
template<typename T>
const T* SVAnalyzer::findClosestLepton(const std::vector<T>& leptons, 
                                       const reco::Candidate* cand)
{
  float bestDR = 999.;
  const T* bestLepton = nullptr;
  
  for (const auto& lep : leptons) {
    float dr = deltaR(lep.p4(), cand->p4());
    if (dr < bestDR) {
      bestDR = dr;
      bestLepton = &lep;
    }
  }
  return bestLepton;
}



// -------------------- findClosestGen --------------------
const reco::GenParticle* SVAnalyzer::findClosestGen(
    const std::vector<reco::GenParticle>& genParticles,
    const reco::Candidate* lepton)
{
  if (!lepton) return nullptr;

  float bestDR = 999.;
  const reco::GenParticle* bestGen = nullptr;

  int pdgId = lepton->pdgId();
  for (const auto& gen : genParticles) {
    if (gen.pdgId() != pdgId) continue;
    float dr = deltaR(gen.p4(), lepton->p4());
    if (dr < bestDR) {
      bestDR = dr;
      bestGen = &gen;
    }
  }
  if (bestDR < 0.1) {
    return bestGen;
  }
  return nullptr;
}



// -------------------- analyze --------------------
void SVAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {

  clearVectors();

  edm::Handle<std::vector<pat::Jet>> jets;
  edm::Handle<std::vector<reco::Vertex>> vertices;
  edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> svs;
  edm::Handle<std::vector<pat::Muon>> muons;
  edm::Handle<std::vector<pat::Electron>> electrons;
  edm::Handle<std::vector<reco::GenParticle>> genParticles;

  iEvent.getByToken(jetToken_,      jets);
  iEvent.getByToken(pvToken_,       vertices);
  iEvent.getByToken(svToken_,       svs);
  iEvent.getByToken(muonToken_,     muons);
  iEvent.getByToken(electronToken_, electrons);
  iEvent.getByToken(genToken_,      genParticles);

  if (!jets.isValid() || !vertices.isValid() || !svs.isValid())
    return;
  if (vertices->empty()) return;

  // ========== 1. Save all jet information ==========
  for (const auto& jet : *jets) {
    if( jet.pt() < 10. || std::abs(jet.eta()) > 2.4 ) continue;
    jet_pt_.push_back(jet.pt());
    jet_eta_.push_back(jet.eta());
    jet_phi_.push_back(jet.phi());
    jet_partonFlavour_.push_back(jet.partonFlavour());
  }


  // ========== 2. Process secondary vertices ==========
  for (const auto& sv : *svs) {

    size_t nSrcCand = sv.numberOfSourceCandidatePtrs();
    if (nSrcCand != 2) continue;

    std::vector<const reco::Candidate*> cands;
    for (size_t i = 0; i < nSrcCand; ++i) {
      const reco::Candidate* cand = sv.sourceCandidatePtr(i).get();
      cands.push_back(cand);
    }

    // Check for opposite sign
    if (cands[0]->charge() * cands[1]->charge() >= 0) continue;

    // Determine if dimuon or dielectron based on pdgId
    bool isDimuon = (std::abs(cands[0]->pdgId()) == 13 && std::abs(cands[1]->pdgId()) == 13);
    bool isDielectron = (std::abs(cands[0]->pdgId()) == 11 && std::abs(cands[1]->pdgId()) == 11);
    if (!isDimuon && !isDielectron) continue;

    // Apply pT > 5 GeV requirement
    if (cands[0]->pt() < 5.0 || cands[1]->pt() < 5.0) continue;

    // Match to PAT muons or electrons using helper function
    std::vector<const pat::Muon*> matchedMuons(2, nullptr);
    std::vector<const pat::Electron*> matchedElectrons(2, nullptr);

    if (isDimuon && muons.isValid()) {
      for (int i = 0; i < 2; ++i) {
        matchedMuons[i] = findClosestLepton(*muons, cands[i]);
      }
    } 
    else if (isDielectron && electrons.isValid()) {
      for (int i = 0; i < 2; ++i) {
        matchedElectrons[i] = findClosestLepton(*electrons, cands[i]);
      }
    }

    // Both leptons must be matched
    if (isDimuon && (!matchedMuons[0] || !matchedMuons[1])) continue;
    if (isDielectron && (!matchedElectrons[0] || !matchedElectrons[1])) continue;

    // Choose best PV for this SV
    auto info = findBestPV(*vertices, sv);

    // Calculate invariant mass
    reco::Candidate::LorentzVector dilepton = cands[0]->p4() + cands[1]->p4();
    float invMass = dilepton.M();

    // Fill SV branches
    sv_mass_.push_back(invMass);
    sv_pt_.push_back(sv.pt());
    sv_eta_.push_back(sv.eta());
    sv_phi_.push_back(sv.phi());
    sv_vx_.push_back(sv.vx());
    sv_vy_.push_back(sv.vy());
    sv_lxy_.push_back(info.lxy);
    sv_cosAngle_.push_back(std::cos(info.angle));
    sv_chi2_.push_back(sv.vertexChi2());
    sv_ndof_.push_back(sv.vertexNdof());

    // ========== 3. Fill lepton information ==========
    std::vector<float> pts, etas, phis, ip2ds;
    std::vector<int> isGlobal, isStandAlone, simMother;

    for (int i = 0; i < 2; ++i) {
      if (isDimuon && matchedMuons[i]) {
        pts.push_back(matchedMuons[i]->pt());
        etas.push_back(matchedMuons[i]->eta());
        phis.push_back(matchedMuons[i]->phi());
        ip2ds.push_back(matchedMuons[i]->dB(pat::Muon::PV2D));
        
        isGlobal.push_back(matchedMuons[i]->isGlobalMuon() ? 1 : 0);
        isStandAlone.push_back(matchedMuons[i]->isStandAloneMuon() ? 1 : 0);
        simMother.push_back(matchedMuons[i]->simMotherPdgId());
      } 
      else if (isDielectron && matchedElectrons[i]) {
        pts.push_back(matchedElectrons[i]->pt());
        etas.push_back(matchedElectrons[i]->eta());
        phis.push_back(matchedElectrons[i]->phi());
        ip2ds.push_back(matchedElectrons[i]->dB(pat::Electron::PV2D));
        
        isGlobal.push_back(-1);
        isStandAlone.push_back(-1);
        simMother.push_back(-999);
      }
    }
    lep_pt_.push_back(pts);
    lep_eta_.push_back(etas);
    lep_phi_.push_back(phis);
    lep_ip2d_.push_back(ip2ds);
    lep_isGlobalMuon_.push_back(isGlobal);
    lep_isStandAloneMuon_.push_back(isStandAlone);
    lep_simMotherPdgId_.push_back(simMother);

    // ========== 4. GEN matching ==========
    std::vector<float> gen_pts, gen_etas, gen_phis;
    std::vector<int> gen_mothers;

    if (genParticles.isValid()) {
      for (int i = 0; i < 2; ++i) {
        const reco::Candidate* lepton = nullptr;
        
        if ( isDimuon ) {
          lepton = matchedMuons[i];
        } else if ( isDielectron ) {
          lepton = matchedElectrons[i];
        }

        const reco::GenParticle* bestGen = findClosestGen(*genParticles, lepton);

        if (bestGen) {
          gen_pts.push_back(bestGen->pt());
          gen_etas.push_back(bestGen->eta());
          gen_phis.push_back(bestGen->phi());
          
          // Get first mother pdgId
          if (bestGen->numberOfMothers() > 0) {
            gen_mothers.push_back(bestGen->mother(0)->pdgId());
          } else {
            gen_mothers.push_back(-1);
          }
        } else {
          gen_pts.push_back(-1.0);
          gen_etas.push_back(-1.0);
          gen_phis.push_back(-1.0);
          gen_mothers.push_back(-1);
        }
      }
    } else {
      // If genParticles not valid, fill with -1
      for (int i = 0; i < 2; ++i) {
        gen_pts.push_back(-1.0);
        gen_etas.push_back(-1.0);
        gen_phis.push_back(-1.0);
        gen_mothers.push_back(-1);
      }
    }

    gen_pt_.push_back(gen_pts);
    gen_eta_.push_back(gen_etas);
    gen_phi_.push_back(gen_phis);
    gen_motherPdgId_.push_back(gen_mothers);
  }

  tree_->Fill();
}



// -------------------- Register Module --------------------
DEFINE_FWK_MODULE(SVAnalyzer);
