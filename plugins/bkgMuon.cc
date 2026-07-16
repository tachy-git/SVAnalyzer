// bkgMuon.cc

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "TTree.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace {
template <typename Func>
void ignoreMissingProduct(Func&& func) {
  try {
    func();
  } catch (const edm::Exception& ex) {
    if (ex.categoryCode() != edm::errors::ProductNotFound) {
      throw;
    }
  }
}

struct MuonBranches {
  int nMuon = 0;
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<float> dRJet;
  std::vector<int> charge;
  std::vector<int> genMuonMatched;
  std::vector<int> genPdgId;
  std::vector<float> genDR;
  std::vector<int> genMotherPdgId;
  std::vector<int> simExtType;
  std::vector<int> simPdgId;
  std::vector<int> isGlobalMuon;
  std::vector<int> isTrackerMuon;
  std::vector<int> isStandaloneMuon;
  std::vector<int> numberOfMatches;
  std::vector<float> segmentCompatibility;
  std::vector<float> dxy;
  std::vector<float> dz;
  std::vector<int> recHitsSize;
  std::vector<int> sharedSegmentsMax;
  std::vector<int> sharedPartner;
  std::vector<int> hasSharedSegment;
  std::vector<int> numberOfMatchedStations;
  std::vector<int> numberOfChambers;
  std::vector<int> numberOfChambersCSCorDT;
  std::vector<int> numberOfMatchedRPCLayers;
  std::vector<int> numberOfStandAloneMuonHits;
  std::vector<int> numberOfStandAloneMuonMatchedStations;
  std::vector<float> standAloneMuonNormChi2;
  std::vector<int> timeNdof;
  std::vector<float> timeAtIpInOut;
  std::vector<int> rpcTimeNdof;
  std::vector<float> rpcTimeAtIpInOut;
  std::vector<float> inverseBeta;
  std::vector<float> inverseBetaErr;

  void clear() {
    nMuon = 0;
    pt.clear();
    eta.clear();
    phi.clear();
    dRJet.clear();
    charge.clear();
    genMuonMatched.clear();
    genPdgId.clear();
    genDR.clear();
    genMotherPdgId.clear();
    simExtType.clear();
    simPdgId.clear();
    isGlobalMuon.clear();
    isTrackerMuon.clear();
    isStandaloneMuon.clear();
    numberOfMatches.clear();
    segmentCompatibility.clear();
    dxy.clear();
    dz.clear();
    recHitsSize.clear();
    sharedSegmentsMax.clear();
    sharedPartner.clear();
    hasSharedSegment.clear();
    numberOfMatchedStations.clear();
    numberOfChambers.clear();
    numberOfChambersCSCorDT.clear();
    numberOfMatchedRPCLayers.clear();
    numberOfStandAloneMuonHits.clear();
    numberOfStandAloneMuonMatchedStations.clear();
    standAloneMuonNormChi2.clear();
    timeNdof.clear();
    timeAtIpInOut.clear();
    rpcTimeNdof.clear();
    rpcTimeAtIpInOut.clear();
    inverseBeta.clear();
    inverseBetaErr.clear();
  }
};

struct TriggerPathSpec {
  std::string pathPrefix;
  std::string branchPrefix;
  std::string objectCountName;
};

const std::array<TriggerPathSpec, 4> kTriggerPathSpecs = {{
    {"HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm_v",
     "hltDoubleL2Mu10NoVtx2Cha",
     "nHltDoubleL2Mu10NoVtx2ChaObject"},
    {"HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm_v",
     "hltDoubleL2Mu10NoVtx2ChaCosmicSeed",
     "nHltDoubleL2Mu10NoVtx2ChaCosmicSeedObject"},
    {"HLT_DoubleMu4_LowMass_Displaced_v",
     "hltDoubleMu4LowMassDisplaced",
     "nHltDoubleMu4LowMassDisplacedObject"},
    {"HLT_DoubleMu4_3_LowMass_v",
     "hltDoubleMu43LowMass",
     "nHltDoubleMu43LowMassObject"},
}};

struct TriggerBranches {
  int found = 0;
  int pass = 0;
  int index = -1;
  int nObject = 0;
  std::vector<int> objectIndex;
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<int> firstTypeId;

  void clear() {
    found = 0;
    pass = 0;
    index = -1;
    nObject = 0;
    objectIndex.clear();
    pt.clear();
    eta.clear();
    phi.clear();
    firstTypeId.clear();
  }
};

bool pathMatchesTarget(const std::string& pathName, const std::string& targetPath) {
  if (pathName == targetPath) return true;
  if (targetPath.size() < 2) return false;
  if (targetPath.compare(targetPath.size() - 2, 2, "_v") != 0) return false;
  return pathName.compare(0, targetPath.size(), targetPath) == 0;
}

void bookMuonBranches(TTree* tree, const std::string& prefix, const std::string& countName, MuonBranches& branches) {
  tree->Branch(countName.c_str(), &branches.nMuon, (countName + "/I").c_str());
  tree->Branch((prefix + "_pt").c_str(), &branches.pt);
  tree->Branch((prefix + "_eta").c_str(), &branches.eta);
  tree->Branch((prefix + "_phi").c_str(), &branches.phi);
  tree->Branch((prefix + "_dRJet").c_str(), &branches.dRJet);
  tree->Branch((prefix + "_charge").c_str(), &branches.charge);
  tree->Branch((prefix + "_genMuonMatched").c_str(), &branches.genMuonMatched);
  tree->Branch((prefix + "_genPdgId").c_str(), &branches.genPdgId);
  tree->Branch((prefix + "_genDR").c_str(), &branches.genDR);
  tree->Branch((prefix + "_genMotherPdgId").c_str(), &branches.genMotherPdgId);
  tree->Branch((prefix + "_simExtType").c_str(), &branches.simExtType);
  tree->Branch((prefix + "_simPdgId").c_str(), &branches.simPdgId);
  tree->Branch((prefix + "_isGlobalMuon").c_str(), &branches.isGlobalMuon);
  tree->Branch((prefix + "_isTrackerMuon").c_str(), &branches.isTrackerMuon);
  tree->Branch((prefix + "_isStandaloneMuon").c_str(), &branches.isStandaloneMuon);
  tree->Branch((prefix + "_numberOfMatches").c_str(), &branches.numberOfMatches);
  tree->Branch((prefix + "_segmentCompatibility").c_str(), &branches.segmentCompatibility);
  tree->Branch((prefix + "_dxy").c_str(), &branches.dxy);
  tree->Branch((prefix + "_dz").c_str(), &branches.dz);
  tree->Branch((prefix + "_recHitsSize").c_str(), &branches.recHitsSize);
  tree->Branch((prefix + "_sharedSegmentsMax").c_str(), &branches.sharedSegmentsMax);
  tree->Branch((prefix + "_sharedPartner").c_str(), &branches.sharedPartner);
  tree->Branch((prefix + "_hasSharedSegment").c_str(), &branches.hasSharedSegment);
  tree->Branch((prefix + "_numberOfMatchedStations").c_str(), &branches.numberOfMatchedStations);
  tree->Branch((prefix + "_numberOfChambers").c_str(), &branches.numberOfChambers);
  tree->Branch((prefix + "_numberOfChambersCSCorDT").c_str(), &branches.numberOfChambersCSCorDT);
  tree->Branch((prefix + "_numberOfMatchedRPCLayers").c_str(), &branches.numberOfMatchedRPCLayers);
  tree->Branch((prefix + "_numberOfStandAloneMuonHits").c_str(), &branches.numberOfStandAloneMuonHits);
  tree->Branch((prefix + "_numberOfStandAloneMuonMatchedStations").c_str(),
               &branches.numberOfStandAloneMuonMatchedStations);
  tree->Branch((prefix + "_standAloneMuonNormChi2").c_str(), &branches.standAloneMuonNormChi2);
  tree->Branch((prefix + "_timeNdof").c_str(), &branches.timeNdof);
  tree->Branch((prefix + "_timeAtIpInOut").c_str(), &branches.timeAtIpInOut);
  tree->Branch((prefix + "_rpcTimeNdof").c_str(), &branches.rpcTimeNdof);
  tree->Branch((prefix + "_rpcTimeAtIpInOut").c_str(), &branches.rpcTimeAtIpInOut);
  tree->Branch((prefix + "_inverseBeta").c_str(), &branches.inverseBeta);
  tree->Branch((prefix + "_inverseBetaErr").c_str(), &branches.inverseBetaErr);
}

void bookTriggerBranches(TTree* tree, const TriggerPathSpec& spec, TriggerBranches& branches) {
  tree->Branch((spec.branchPrefix + "_found").c_str(), &branches.found, (spec.branchPrefix + "_found/I").c_str());
  tree->Branch((spec.branchPrefix + "_pass").c_str(), &branches.pass, (spec.branchPrefix + "_pass/I").c_str());
  tree->Branch((spec.branchPrefix + "_index").c_str(), &branches.index, (spec.branchPrefix + "_index/I").c_str());
  tree->Branch(spec.objectCountName.c_str(), &branches.nObject, (spec.objectCountName + "/I").c_str());
  tree->Branch((spec.branchPrefix + "Object_index").c_str(), &branches.objectIndex);
  tree->Branch((spec.branchPrefix + "Object_pt").c_str(), &branches.pt);
  tree->Branch((spec.branchPrefix + "Object_eta").c_str(), &branches.eta);
  tree->Branch((spec.branchPrefix + "Object_phi").c_str(), &branches.phi);
  tree->Branch((spec.branchPrefix + "Object_firstTypeId").c_str(), &branches.firstTypeId);
}

void fillMuonBranches(const std::vector<pat::Muon>& muons,
                      const pat::Jet& leadingJet,
                      const edm::Handle<reco::GenParticleCollection>& genParticles,
                      MuonBranches& branches) {
  std::vector<const pat::Muon*> selectedMuons;
  selectedMuons.reserve(muons.size());

  for (const auto& recoMuon : muons) {
    const float dRJet = static_cast<float>(
        reco::deltaR(recoMuon.eta(), recoMuon.phi(), leadingJet.eta(), leadingJet.phi()));
    if (dRJet > 0.4f) continue;

    selectedMuons.push_back(&recoMuon);
    branches.pt.push_back(recoMuon.pt());
    branches.eta.push_back(recoMuon.eta());
    branches.phi.push_back(recoMuon.phi());
    branches.dRJet.push_back(dRJet);
    branches.charge.push_back(recoMuon.charge());

    int genMuonMatched = 0;
    int genPdgId = 0;
    int genMotherPdgId = 0;
    float genDR = -999.f;
    float bestGenMuonDR = 0.1f;
    float bestAnyGenDR = 999.f;
    int closestGenPdgId = 0;
    if (genParticles.isValid()) {
      for (const auto& gen : *genParticles) {
        const float dRGen = static_cast<float>(
            reco::deltaR(recoMuon.eta(), recoMuon.phi(), gen.eta(), gen.phi()));

        if (dRGen < bestAnyGenDR) {
          bestAnyGenDR = dRGen;
          closestGenPdgId = gen.pdgId();
        }

        if (std::abs(gen.pdgId()) != 13) continue;
        if (gen.status() != 1) continue;
        if (dRGen >= bestGenMuonDR) continue;

        genMuonMatched = 1;
        bestGenMuonDR = dRGen;
        genPdgId = gen.pdgId();
        genDR = dRGen;
        genMotherPdgId = gen.numberOfMothers() > 0 && gen.mother(0) ? gen.mother(0)->pdgId() : 0;
      }

      if (!genMuonMatched && closestGenPdgId != 0) {
        genPdgId = closestGenPdgId;
        genDR = bestAnyGenDR;
      }
    }
    branches.genMuonMatched.push_back(genMuonMatched);
    branches.genPdgId.push_back(genPdgId);
    branches.genDR.push_back(genDR);
    branches.genMotherPdgId.push_back(genMotherPdgId);

    branches.simExtType.push_back(static_cast<int>(recoMuon.simExtType()));
    branches.simPdgId.push_back(recoMuon.simPdgId());
    branches.isGlobalMuon.push_back(static_cast<int>(recoMuon.isGlobalMuon()));
    branches.isTrackerMuon.push_back(static_cast<int>(recoMuon.isTrackerMuon()));
    branches.isStandaloneMuon.push_back(static_cast<int>(recoMuon.isStandAloneMuon()));
    branches.numberOfMatches.push_back(recoMuon.numberOfMatches());
    branches.segmentCompatibility.push_back(recoMuon.segmentCompatibility());

    float dxy = -999.f;
    float dz = -999.f;
    int recHitsSize = -1;
    const auto bestTrack = recoMuon.bestTrack();
    if (bestTrack) {
      dxy = bestTrack->dxy();
      dz = bestTrack->dz();
      ignoreMissingProduct([&]() { recHitsSize = bestTrack->recHitsSize(); });
    }
    branches.dxy.push_back(dxy);
    branches.dz.push_back(dz);
    branches.recHitsSize.push_back(recHitsSize);

    branches.numberOfMatchedStations.push_back(recoMuon.numberOfMatchedStations());
    branches.numberOfChambers.push_back(recoMuon.numberOfChambers());
    branches.numberOfChambersCSCorDT.push_back(recoMuon.numberOfChambersCSCorDT());
    branches.numberOfMatchedRPCLayers.push_back(recoMuon.numberOfMatchedRPCLayers());

    float standAloneMuonNormChi2 = -999.f;
    int numberOfStandAloneMuonHits = -1;
    int numberOfStandAloneMuonMatchedStations = -1;
    if (recoMuon.standAloneMuon().isNonnull()) {
      ignoreMissingProduct([&]() {
        standAloneMuonNormChi2 = recoMuon.standAloneMuon()->normalizedChi2();
        const auto& hitPattern = recoMuon.standAloneMuon()->hitPattern();
        numberOfStandAloneMuonHits = hitPattern.numberOfValidMuonHits();
        numberOfStandAloneMuonMatchedStations = hitPattern.muonStationsWithValidHits();
      });
    }
    branches.standAloneMuonNormChi2.push_back(standAloneMuonNormChi2);
    branches.numberOfStandAloneMuonHits.push_back(numberOfStandAloneMuonHits);
    branches.numberOfStandAloneMuonMatchedStations.push_back(numberOfStandAloneMuonMatchedStations);

    int timeNdof = -1;
    float timeAtIpInOut = -999.f;
    if (recoMuon.isTimeValid()) {
      timeNdof = recoMuon.time().nDof;
      timeAtIpInOut = recoMuon.time().timeAtIpInOut;
    }
    branches.timeNdof.push_back(timeNdof);
    branches.timeAtIpInOut.push_back(timeAtIpInOut);

    int rpcTimeNdof = -1;
    float rpcTimeAtIpInOut = -999.f;
    if (recoMuon.rpcTime().nDof > 0) {
      rpcTimeNdof = recoMuon.rpcTime().nDof;
      rpcTimeAtIpInOut = recoMuon.rpcTime().timeAtIpInOut;
    }
    branches.rpcTimeNdof.push_back(rpcTimeNdof);
    branches.rpcTimeAtIpInOut.push_back(rpcTimeAtIpInOut);

    branches.inverseBeta.push_back(recoMuon.inverseBeta());
    branches.inverseBetaErr.push_back(recoMuon.inverseBetaErr());
  }

  branches.nMuon = static_cast<int>(selectedMuons.size());
  branches.sharedSegmentsMax.assign(selectedMuons.size(), 0);
  branches.sharedPartner.assign(selectedMuons.size(), -1);
  branches.hasSharedSegment.assign(selectedMuons.size(), 0);

  for (size_t i = 0; i < selectedMuons.size(); ++i) {
    for (size_t j = i + 1; j < selectedMuons.size(); ++j) {
      const int nShared = muon::sharedSegments(*selectedMuons[i], *selectedMuons[j]);
      if (nShared > branches.sharedSegmentsMax[i]) {
        branches.sharedSegmentsMax[i] = nShared;
        branches.sharedPartner[i] = static_cast<int>(j);
      }
      if (nShared > branches.sharedSegmentsMax[j]) {
        branches.sharedSegmentsMax[j] = nShared;
        branches.sharedPartner[j] = static_cast<int>(i);
      }
      if (nShared > 0) {
        branches.hasSharedSegment[i] = 1;
        branches.hasSharedSegment[j] = 1;
      }
    }
  }
}
}  // namespace

class bkgMuon : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit bkgMuon(const edm::ParameterSet&);
  ~bkgMuon() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void clearBranches();

  edm::EDGetTokenT<std::vector<pat::Jet>> jetToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> displacedMuonToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken_;

  TTree* tree_;

  int run_;
  int lumi_;
  unsigned long long event_;

  int hasLeadingJet_;
  float jet_pt_;
  float jet_eta_;
  float jet_phi_;
  int jet_partonFlavour_;

  std::array<TriggerBranches, kTriggerPathSpecs.size()> triggerBranches_;
  MuonBranches muonBranches_;
  MuonBranches disMuonBranches_;
};

bkgMuon::bkgMuon(const edm::ParameterSet& iConfig)
    : jetToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      muonToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      displacedMuonToken_(
          consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("displacedMuons"))),
      genParticleToken_(
          consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
      triggerResultsToken_(
          consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(
          iConfig.getParameter<edm::InputTag>("triggerObjects"))),
      tree_(nullptr),
      run_(0),
      lumi_(0),
      event_(0),
      hasLeadingJet_(0),
      jet_pt_(-999.f),
      jet_eta_(-999.f),
      jet_phi_(-999.f),
      jet_partonFlavour_(0) {
  usesResource("TFileService");
}

void bkgMuon::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "Leading jet and muons inside dR < 0.4");

  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");
  tree_->Branch("event", &event_, "event/l");

  tree_->Branch("hasLeadingJet", &hasLeadingJet_, "hasLeadingJet/I");
  tree_->Branch("jet_pt", &jet_pt_, "jet_pt/F");
  tree_->Branch("jet_eta", &jet_eta_, "jet_eta/F");
  tree_->Branch("jet_phi", &jet_phi_, "jet_phi/F");
  tree_->Branch("jet_partonFlavour", &jet_partonFlavour_, "jet_partonFlavour/I");

  for (size_t i = 0; i < kTriggerPathSpecs.size(); ++i) {
    bookTriggerBranches(tree_, kTriggerPathSpecs[i], triggerBranches_[i]);
  }

  bookMuonBranches(tree_, "muon", "nMuon", muonBranches_);
  bookMuonBranches(tree_, "disMuon", "nDisMuon", disMuonBranches_);
}

void bkgMuon::clearBranches() {
  hasLeadingJet_ = 0;
  jet_pt_ = -999.f;
  jet_eta_ = -999.f;
  jet_phi_ = -999.f;
  jet_partonFlavour_ = 0;

  for (auto& triggerBranches : triggerBranches_) {
    triggerBranches.clear();
  }

  muonBranches_.clear();
  disMuonBranches_.clear();
}

void bkgMuon::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  clearBranches();

  run_ = iEvent.id().run();
  lumi_ = iEvent.id().luminosityBlock();
  event_ = iEvent.id().event();

  edm::Handle<std::vector<pat::Jet>> jets;
  edm::Handle<std::vector<pat::Muon>> muons;
  edm::Handle<std::vector<pat::Muon>> displacedMuons;
  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<edm::TriggerResults> triggerResults;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(jetToken_, jets);
  iEvent.getByToken(muonToken_, muons);
  iEvent.getByToken(displacedMuonToken_, displacedMuons);
  iEvent.getByToken(genParticleToken_, genParticles);
  iEvent.getByToken(triggerResultsToken_, triggerResults);
  iEvent.getByToken(triggerObjectsToken_, triggerObjects);

  if (triggerResults.isValid()) {
    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
    std::array<std::string, kTriggerPathSpecs.size()> resolvedPaths;
    for (unsigned int i = 0; i < triggerNames.size(); ++i) {
      const std::string& pathName = triggerNames.triggerName(i);
      for (size_t pathIndex = 0; pathIndex < kTriggerPathSpecs.size(); ++pathIndex) {
        if (triggerBranches_[pathIndex].found) continue;
        if (!pathMatchesTarget(pathName, kTriggerPathSpecs[pathIndex].pathPrefix)) continue;

        triggerBranches_[pathIndex].found = 1;
        triggerBranches_[pathIndex].index = static_cast<int>(i);
        triggerBranches_[pathIndex].pass = triggerResults->accept(i) ? 1 : 0;
        resolvedPaths[pathIndex] = pathName;
      }
    }

    if (triggerObjects.isValid()) {
      for (size_t objectIndex = 0; objectIndex < triggerObjects->size(); ++objectIndex) {
        const auto& triggerObject = triggerObjects->at(objectIndex);
        auto unpackedObject = triggerObject;
        unpackedObject.unpackPathNames(triggerNames);
        const auto typeIds = unpackedObject.triggerObjectTypes();

        for (size_t pathIndex = 0; pathIndex < kTriggerPathSpecs.size(); ++pathIndex) {
          if (!triggerBranches_[pathIndex].found) continue;
          if (!unpackedObject.hasPathName(resolvedPaths[pathIndex], true, false)) continue;

          auto& branches = triggerBranches_[pathIndex];
          branches.objectIndex.push_back(static_cast<int>(objectIndex));
          branches.pt.push_back(unpackedObject.pt());
          branches.eta.push_back(unpackedObject.eta());
          branches.phi.push_back(unpackedObject.phi());
          branches.firstTypeId.push_back(typeIds.empty() ? 0 : typeIds.front());
          branches.nObject = static_cast<int>(branches.objectIndex.size());
        }
      }
    }
  }

  if (!jets.isValid()) {
    tree_->Fill();
    return;
  }

  const pat::Jet* leadingJet = nullptr;
  for (const auto& jet : *jets) {
    if (jet.pt() <= 30.f) continue;
    if (std::abs(jet.eta()) >= 2.4f) continue;
    if (!leadingJet || jet.pt() > leadingJet->pt()) {
      leadingJet = &jet;
    }
  }

  if (!leadingJet) {
    tree_->Fill();
    return;
  }

  hasLeadingJet_ = 1;
  jet_pt_ = leadingJet->pt();
  jet_eta_ = leadingJet->eta();
  jet_phi_ = leadingJet->phi();
  jet_partonFlavour_ = leadingJet->partonFlavour();

  if (muons.isValid()) {
    fillMuonBranches(*muons, *leadingJet, genParticles, muonBranches_);
  }
  if (displacedMuons.isValid()) {
    fillMuonBranches(*displacedMuons, *leadingJet, genParticles, disMuonBranches_);
  }

  tree_->Fill();
}

void bkgMuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("displacedMuons", edm::InputTag("slimmedDisplacedMuons"));
  desc.add<edm::InputTag>("genParticles", edm::InputTag("prunedGenParticles"));
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  desc.add<edm::InputTag>("triggerObjects", edm::InputTag("slimmedPatTrigger"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(bkgMuon);
