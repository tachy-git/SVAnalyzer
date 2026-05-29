// studySimMuon.cc

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TTree.h"

#include <cmath>
#include <set>
#include <vector>

class studySimMuon : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit studySimMuon(const edm::ParameterSet&);
  ~studySimMuon() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void clearBranches();
  void clearNotMatchedBranches();

  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<std::vector<pat::Jet>> jetToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;

  // ======= TTree for gen-matched reco muons (one entry per gen muon) =======
  TTree* tree_;

  float gen_muon_pt_;
  float gen_muon_eta_;
  float gen_muon_phi_;

  std::vector<float> reco_muon_dR_;
  std::vector<float> reco_muon_ptRatio_dR0p1_;

  int n_reco_muon_dR0p1_ptRatio0p9to1p1_;
  std::vector<int> reco_muon_isGlobalMuon_;
  std::vector<int> reco_muon_isTrackerMuon_;
  std::vector<int> reco_muon_isStandaloneMuon_;
  std::vector<int> reco_muon_isCaloMuon_;
  std::vector<int> reco_muon_isPFMuon_;
  std::vector<int> reco_muon_isRPCMuon_;
  std::vector<int> reco_muon_isGEMMuon_;
  std::vector<int> reco_muon_isME0Muon_;
  std::vector<int> reco_muon_simExtType_;
  std::vector<int> reco_muon_simPdgId_;

  int closest_reco_muon_isGlobalMuon_;
  int closest_reco_muon_isTrackerMuon_;
  int closest_reco_muon_isStandaloneMuon_;
  int closest_reco_muon_isCaloMuon_;
  int closest_reco_muon_isPFMuon_;
  int closest_reco_muon_isRPCMuon_;
  int closest_reco_muon_isGEMMuon_;
  int closest_reco_muon_isME0Muon_;
  int closest_reco_muon_simExtType_;
  int closest_reco_muon_simPdgId_;

  // ======= TTree for not-matched reco muons (one entry per reco muon) =======
  TTree* notMatchedTree_;

  float notMatched_reco_muon_pt_;
  float notMatched_reco_muon_eta_;
  float notMatched_reco_muon_phi_;
  int notMatched_reco_muon_isGlobalMuon_;
  int notMatched_reco_muon_isTrackerMuon_;
  int notMatched_reco_muon_isStandaloneMuon_;
  int notMatched_reco_muon_isCaloMuon_;
  int notMatched_reco_muon_isPFMuon_;
  int notMatched_reco_muon_isRPCMuon_;
  int notMatched_reco_muon_isGEMMuon_;
  int notMatched_reco_muon_isME0Muon_;
  int notMatched_reco_muon_simExtType_;
  int notMatched_reco_muon_simPdgId_;
};

studySimMuon::studySimMuon(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      jetToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      genParticleToken_(
          consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
      tree_(nullptr),
      gen_muon_pt_(0.f),
      gen_muon_eta_(0.f),
      gen_muon_phi_(0.f),
      n_reco_muon_dR0p1_ptRatio0p9to1p1_(0),
      closest_reco_muon_isGlobalMuon_(-1),
      closest_reco_muon_isTrackerMuon_(-1),
      closest_reco_muon_isStandaloneMuon_(-1),
      closest_reco_muon_isCaloMuon_(-1),
      closest_reco_muon_isPFMuon_(-1),
      closest_reco_muon_isRPCMuon_(-1),
      closest_reco_muon_isGEMMuon_(-1),
      closest_reco_muon_isME0Muon_(-1),
      closest_reco_muon_simExtType_(-1),
      closest_reco_muon_simPdgId_(0),
      notMatchedTree_(nullptr),
      notMatched_reco_muon_pt_(0.f),
      notMatched_reco_muon_eta_(0.f),
      notMatched_reco_muon_phi_(0.f),
      notMatched_reco_muon_isGlobalMuon_(-1),
      notMatched_reco_muon_isTrackerMuon_(-1),
      notMatched_reco_muon_isStandaloneMuon_(-1),
      notMatched_reco_muon_isCaloMuon_(-1),
      notMatched_reco_muon_isPFMuon_(-1),
      notMatched_reco_muon_isRPCMuon_(-1),
      notMatched_reco_muon_isGEMMuon_(-1),
      notMatched_reco_muon_isME0Muon_(-1),
      notMatched_reco_muon_simExtType_(-1),
      notMatched_reco_muon_simPdgId_(0) {
  usesResource("TFileService");
}

void studySimMuon::beginJob() {
  edm::Service<TFileService> fs;

  // ======= Gen-matched reco muon tree =======
  tree_ = fs->make<TTree>("tree", "Gen muon to reco muon matching");

  tree_->Branch("gen_muon_pt", &gen_muon_pt_, "gen_muon_pt/F");
  tree_->Branch("gen_muon_eta", &gen_muon_eta_, "gen_muon_eta/F");
  tree_->Branch("gen_muon_phi", &gen_muon_phi_, "gen_muon_phi/F");

  tree_->Branch("reco_muon_dR", &reco_muon_dR_);
  tree_->Branch("reco_muon_ptRatio_dR0p1", &reco_muon_ptRatio_dR0p1_);

  tree_->Branch("n_reco_muon_dR0p1_ptRatio0p9to1p1",
                &n_reco_muon_dR0p1_ptRatio0p9to1p1_,
                "n_reco_muon_dR0p1_ptRatio0p9to1p1/I");
  tree_->Branch("reco_muon_isGlobalMuon", &reco_muon_isGlobalMuon_);
  tree_->Branch("reco_muon_isTrackerMuon", &reco_muon_isTrackerMuon_);
  tree_->Branch("reco_muon_isStandaloneMuon", &reco_muon_isStandaloneMuon_);
  tree_->Branch("reco_muon_isCaloMuon", &reco_muon_isCaloMuon_);
  tree_->Branch("reco_muon_isPFMuon", &reco_muon_isPFMuon_);
  tree_->Branch("reco_muon_isRPCMuon", &reco_muon_isRPCMuon_);
  tree_->Branch("reco_muon_isGEMMuon", &reco_muon_isGEMMuon_);
  tree_->Branch("reco_muon_isME0Muon", &reco_muon_isME0Muon_);
  tree_->Branch("reco_muon_simExtType", &reco_muon_simExtType_);
  tree_->Branch("reco_muon_simPdgId", &reco_muon_simPdgId_);

  tree_->Branch("closest_reco_muon_isGlobalMuon",
                &closest_reco_muon_isGlobalMuon_,
                "closest_reco_muon_isGlobalMuon/I");
  tree_->Branch("closest_reco_muon_isTrackerMuon",
                &closest_reco_muon_isTrackerMuon_,
                "closest_reco_muon_isTrackerMuon/I");
  tree_->Branch("closest_reco_muon_isStandaloneMuon",
                &closest_reco_muon_isStandaloneMuon_,
                "closest_reco_muon_isStandaloneMuon/I");
  tree_->Branch("closest_reco_muon_isCaloMuon",
                &closest_reco_muon_isCaloMuon_,
                "closest_reco_muon_isCaloMuon/I");
  tree_->Branch("closest_reco_muon_isPFMuon",
                &closest_reco_muon_isPFMuon_,
                "closest_reco_muon_isPFMuon/I");
  tree_->Branch("closest_reco_muon_isRPCMuon",
                &closest_reco_muon_isRPCMuon_,
                "closest_reco_muon_isRPCMuon/I");
  tree_->Branch("closest_reco_muon_isGEMMuon",
                &closest_reco_muon_isGEMMuon_,
                "closest_reco_muon_isGEMMuon/I");
  tree_->Branch("closest_reco_muon_isME0Muon",
                &closest_reco_muon_isME0Muon_,
                "closest_reco_muon_isME0Muon/I");
  tree_->Branch("closest_reco_muon_simExtType",
                &closest_reco_muon_simExtType_,
                "closest_reco_muon_simExtType/I");
  tree_->Branch("closest_reco_muon_simPdgId",
                &closest_reco_muon_simPdgId_,
                "closest_reco_muon_simPdgId/I");

  // ======= Not-matched reco muon tree (one entry per reco muon) =======
  notMatchedTree_ = fs->make<TTree>("notMatchedTree", "Reco muons not matched to any gen muon");

  notMatchedTree_->Branch("pt", &notMatched_reco_muon_pt_, "pt/F");
  notMatchedTree_->Branch("eta", &notMatched_reco_muon_eta_, "eta/F");
  notMatchedTree_->Branch("phi", &notMatched_reco_muon_phi_, "phi/F");
  notMatchedTree_->Branch("isGlobalMuon", &notMatched_reco_muon_isGlobalMuon_, "isGlobalMuon/I");
  notMatchedTree_->Branch("isTrackerMuon", &notMatched_reco_muon_isTrackerMuon_, "isTrackerMuon/I");
  notMatchedTree_->Branch("isStandaloneMuon", &notMatched_reco_muon_isStandaloneMuon_, "isStandaloneMuon/I");
  notMatchedTree_->Branch("isCaloMuon", &notMatched_reco_muon_isCaloMuon_, "isCaloMuon/I");
  notMatchedTree_->Branch("isPFMuon", &notMatched_reco_muon_isPFMuon_, "isPFMuon/I");
  notMatchedTree_->Branch("isRPCMuon", &notMatched_reco_muon_isRPCMuon_, "isRPCMuon/I");
  notMatchedTree_->Branch("isGEMMuon", &notMatched_reco_muon_isGEMMuon_, "isGEMMuon/I");
  notMatchedTree_->Branch("isME0Muon", &notMatched_reco_muon_isME0Muon_, "isME0Muon/I");
  notMatchedTree_->Branch("simExtType", &notMatched_reco_muon_simExtType_, "simExtType/I");
  notMatchedTree_->Branch("simPdgId", &notMatched_reco_muon_simPdgId_, "simPdgId/I");
}

void studySimMuon::clearBranches() {
  gen_muon_pt_ = 0.f;
  gen_muon_eta_ = 0.f;
  gen_muon_phi_ = 0.f;

  reco_muon_dR_.clear();
  reco_muon_ptRatio_dR0p1_.clear();

  n_reco_muon_dR0p1_ptRatio0p9to1p1_ = 0;
  reco_muon_isGlobalMuon_.clear();
  reco_muon_isTrackerMuon_.clear();
  reco_muon_isStandaloneMuon_.clear();
  reco_muon_isCaloMuon_.clear();
  reco_muon_isPFMuon_.clear();
  reco_muon_isRPCMuon_.clear();
  reco_muon_isGEMMuon_.clear();
  reco_muon_isME0Muon_.clear();
  reco_muon_simExtType_.clear();
  reco_muon_simPdgId_.clear();

  closest_reco_muon_isGlobalMuon_ = -1;
  closest_reco_muon_isTrackerMuon_ = -1;
  closest_reco_muon_isStandaloneMuon_ = -1;
  closest_reco_muon_isCaloMuon_ = -1;
  closest_reco_muon_isPFMuon_ = -1;
  closest_reco_muon_isRPCMuon_ = -1;
  closest_reco_muon_isGEMMuon_ = -1;
  closest_reco_muon_isME0Muon_ = -1;
  closest_reco_muon_simExtType_ = -1;
  closest_reco_muon_simPdgId_ = 0;
}

void studySimMuon::clearNotMatchedBranches() {
  notMatched_reco_muon_pt_ = 0.f;
  notMatched_reco_muon_eta_ = 0.f;
  notMatched_reco_muon_phi_ = 0.f;
  notMatched_reco_muon_isGlobalMuon_ = -1;
  notMatched_reco_muon_isTrackerMuon_ = -1;
  notMatched_reco_muon_isStandaloneMuon_ = -1;
  notMatched_reco_muon_isCaloMuon_ = -1;
  notMatched_reco_muon_isPFMuon_ = -1;
  notMatched_reco_muon_isRPCMuon_ = -1;
  notMatched_reco_muon_isGEMMuon_ = -1;
  notMatched_reco_muon_isME0Muon_ = -1;
  notMatched_reco_muon_simExtType_ = -1;
  notMatched_reco_muon_simPdgId_ = 0;
}

void studySimMuon::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<std::vector<pat::Muon>> muons;
  edm::Handle<std::vector<pat::Jet>> jets;
  edm::Handle<reco::GenParticleCollection> genParticles;

  iEvent.getByToken(muonToken_, muons);
  iEvent.getByToken(jetToken_, jets);
  iEvent.getByToken(genParticleToken_, genParticles);

  if (!muons.isValid() || !jets.isValid() || !genParticles.isValid()) {
    return;
  }

  // ======= Pass 1: per gen muon, fill tree + collect matched reco indices =======
  std::set<size_t> matchedRecoIndices;

  for (const auto& genMuon : *genParticles) {
    if (std::abs(genMuon.pdgId()) != 13) continue;
    if (genMuon.status() != 1) continue;
    if (genMuon.pt() <= 30. || std::abs(genMuon.eta()) > 0.8) continue;

    bool genMuonIsInJet = false;
    for (const auto& jet : *jets) {
      if (reco::deltaR(genMuon.eta(), genMuon.phi(), jet.eta(), jet.phi()) < 0.3) {
        genMuonIsInJet = true;
        break;
      }
    }
    if (!genMuonIsInJet) continue;

    clearBranches();
    gen_muon_pt_ = genMuon.pt();
    gen_muon_eta_ = genMuon.eta();
    gen_muon_phi_ = genMuon.phi();

    float smallestMatchedDR = 999.f;

    for (size_t i = 0; i < muons->size(); ++i) {
      const auto& muon = (*muons)[i];
      const float dR = reco::deltaR(genMuon.eta(), genMuon.phi(), muon.eta(), muon.phi());
      const float ptRatio = muon.pt() / genMuon.pt();

      reco_muon_dR_.push_back(dR);

      if (dR < 0.1) {
        reco_muon_ptRatio_dR0p1_.push_back(ptRatio);
      }

      if (dR < 0.1 && ptRatio > 0.9 && ptRatio < 1.1) {
        matchedRecoIndices.insert(i);
        ++n_reco_muon_dR0p1_ptRatio0p9to1p1_;
        reco_muon_isGlobalMuon_.push_back(static_cast<int>(muon.isGlobalMuon()));
        reco_muon_isTrackerMuon_.push_back(static_cast<int>(muon.isTrackerMuon()));
        reco_muon_isStandaloneMuon_.push_back(static_cast<int>(muon.isStandAloneMuon()));
        reco_muon_isCaloMuon_.push_back(static_cast<int>(muon.isCaloMuon()));
        reco_muon_isPFMuon_.push_back(static_cast<int>(muon.isPFMuon()));
        reco_muon_isRPCMuon_.push_back(static_cast<int>(muon.isRPCMuon()));
        reco_muon_isGEMMuon_.push_back(static_cast<int>(muon.isGEMMuon()));
        reco_muon_isME0Muon_.push_back(static_cast<int>(muon.isME0Muon()));
        reco_muon_simExtType_.push_back(static_cast<int>(muon.simExtType()));
        reco_muon_simPdgId_.push_back(muon.simPdgId());

        if (dR < smallestMatchedDR) {
          smallestMatchedDR = dR;
          closest_reco_muon_isGlobalMuon_ = static_cast<int>(muon.isGlobalMuon());
          closest_reco_muon_isTrackerMuon_ = static_cast<int>(muon.isTrackerMuon());
          closest_reco_muon_isStandaloneMuon_ = static_cast<int>(muon.isStandAloneMuon());
          closest_reco_muon_isCaloMuon_ = static_cast<int>(muon.isCaloMuon());
          closest_reco_muon_isPFMuon_ = static_cast<int>(muon.isPFMuon());
          closest_reco_muon_isRPCMuon_ = static_cast<int>(muon.isRPCMuon());
          closest_reco_muon_isGEMMuon_ = static_cast<int>(muon.isGEMMuon());
          closest_reco_muon_isME0Muon_ = static_cast<int>(muon.isME0Muon());
          closest_reco_muon_simExtType_ = static_cast<int>(muon.simExtType());
          closest_reco_muon_simPdgId_ = muon.simPdgId();
        }
      }
    }

    if (n_reco_muon_dR0p1_ptRatio0p9to1p1_ == 0) {
      edm::LogPrint("studySimMuon")
          << "No matched reco muon: Run=" << iEvent.id().run()
          << " Lumi=" << iEvent.id().luminosityBlock()
          << " Event=" << iEvent.id().event()
          << " genMuon pt=" << genMuon.pt()
          << " eta=" << genMuon.eta()
          << " phi=" << genMuon.phi();
    }

    tree_->Fill();
  }

  // ======= Pass 2: fill notMatchedTree for all unmatched reco muons =======
  for (size_t i = 0; i < muons->size(); ++i) {
    if (matchedRecoIndices.count(i) > 0) continue;

    const auto& muon = (*muons)[i];
    clearNotMatchedBranches();

    notMatched_reco_muon_pt_ = muon.pt();
    notMatched_reco_muon_eta_ = muon.eta();
    notMatched_reco_muon_phi_ = muon.phi();
    notMatched_reco_muon_isGlobalMuon_ = static_cast<int>(muon.isGlobalMuon());
    notMatched_reco_muon_isTrackerMuon_ = static_cast<int>(muon.isTrackerMuon());
    notMatched_reco_muon_isStandaloneMuon_ = static_cast<int>(muon.isStandAloneMuon());
    notMatched_reco_muon_isCaloMuon_ = static_cast<int>(muon.isCaloMuon());
    notMatched_reco_muon_isPFMuon_ = static_cast<int>(muon.isPFMuon());
    notMatched_reco_muon_isRPCMuon_ = static_cast<int>(muon.isRPCMuon());
    notMatched_reco_muon_isGEMMuon_ = static_cast<int>(muon.isGEMMuon());
    notMatched_reco_muon_isME0Muon_ = static_cast<int>(muon.isME0Muon());
    notMatched_reco_muon_simExtType_ = static_cast<int>(muon.simExtType());
    notMatched_reco_muon_simPdgId_ = muon.simPdgId();

    notMatchedTree_->Fill();
  }
}

void studySimMuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("genParticles", edm::InputTag("prunedGenParticles"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(studySimMuon);
