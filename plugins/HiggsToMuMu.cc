#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TMath.h"
#include "TLorentzVector.h"

class HiggsToMuMu : public edm::EDAnalyzer {
public:
  explicit HiggsToMuMu(const edm::ParameterSet&);
  ~HiggsToMuMu() {}
  void analyze(const edm::Event&, const edm::EventSetup&);
  TString replace_all(const TString& a, const TString& b, const TString& c);
private:
  struct tree_t {
   
    unsigned run;
    unsigned lumi;
    unsigned event;
    int nvertices;

    float mu1_pt, mu1_eta, mu1_phi;
    float mu2_pt, mu2_eta, mu2_phi;

    float mu1_dpv, mu2_dpv;
    float mu1_dz, mu2_dz;
    float mu1_dr, mu2_dr;

    float dimu_mass;
    float dimu_pt;
    float dimu_eta;
    float dimu_phi;
    float dimu_dR;
    float dimu_dPhi;
    float dimu_dEta;

    float met_pt;
    float met_phi;
    int nJets;
    int nJets_lt_2p4;
    int nJets_gt_2p4; 
///Two higest pt jets
    float jet1_pt, jet1_energy, jet1_eta, jet1_phi, jet1_btag;
    float jet2_pt, jet2_energy, jet2_eta, jet2_phi, jet2_btag;
    float jet3_pt, jet3_energy, jet3_eta, jet3_phi, jet3_btag;
    float jet4_pt, jet4_energy, jet4_eta, jet4_phi, jet4_btag;
    float jet1_dpv, jet2_dpv;
    float jet1_dz, jet2_dz;
    float jet1_dr, jet2_dr;
    float j1j2_mass;

///Two higest mass dijet pairs
    float jj1_mass;
    float jj1_dEta;
    float jj2_mass;
    float jj2_dEta;

    int nJets_CSV;
  };

  tree_t t;
  TTree* tree;

  const edm::InputTag dimu_src;
  const edm::InputTag beamspot_src;
  const edm::InputTag met_src;
  const edm::InputTag jet_src;
  const edm::InputTag vertices_src;
  std::vector<edm::InputTag> filterTags;
};

TString HiggsToMuMu::replace_all(const TString& a, const TString& b, const TString& c) {
  TString ret = a;
  ret.ReplaceAll(b, c);
  return ret;
}

HiggsToMuMu::HiggsToMuMu(const edm::ParameterSet& cfg)
  : dimu_src(cfg.getParameter<edm::InputTag>("dimu_src")),
    met_src(cfg.getParameter<edm::InputTag>("met_src")),
    jet_src(cfg.getParameter<edm::InputTag>("jet_src")),
    vertices_src(cfg.getParameter<edm::InputTag>("vertices_src"))
{
 
  consumes<pat::CompositeCandidateCollection>(dimu_src);
  consumes<std::vector<pat::MET>>(met_src);
  consumes<std::vector<pat::Jet>>(jet_src);
  consumes<reco::VertexCollection>(vertices_src);
 
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  tree->Branch("run", &t.run, "run/i");
  tree->Branch("lumi", &t.lumi, "lumi/i");
  tree->Branch("event", &t.event, "event/i");

  tree->Branch("nvertices", &t.nvertices, "nvertices/I");

  tree->Branch("mu1_pt", &t.mu1_pt, "mu1_pt/F");
  tree->Branch("mu1_eta", &t.mu1_eta, "mu1_eta/F");
  tree->Branch("mu1_phi", &t.mu1_phi, "mu1_phi/F");
  tree->Branch("mu2_pt", &t.mu2_pt, "mu2_pt/F");
  tree->Branch("mu2_eta", &t.mu2_eta, "mu2_eta/F");
  tree->Branch("mu2_phi", &t.mu2_phi, "mu2_phi/F");

  tree->Branch("mu1_dpv", &t.mu1_dpv, "mu1_dpv/F");
  tree->Branch("mu2_dpv", &t.mu2_dpv, "mu2_dpv/F");
  tree->Branch("mu1_dz", &t.mu1_dz, "mu1_dz/F");
  tree->Branch("mu2_dz", &t.mu2_dz, "mu2_dz/F");
  tree->Branch("mu1_dr", &t.mu1_dr, "mu1_dr/F");
  tree->Branch("mu2_dr", &t.mu2_dr, "mu2_dr/F");

  tree->Branch("dimu_mass", &t.dimu_mass, "dimu_mass/F");
  tree->Branch("dimu_pt", &t.dimu_pt, "dimu_pt/F");
  tree->Branch("dimu_eta", &t.dimu_eta, "dimu_eta/F");
  tree->Branch("dimu_phi", &t.dimu_phi, "dimu_phi/F");
  tree->Branch("dimu_dR", &t.dimu_dR, "dimu_dR/F");
  tree->Branch("dimu_dPhi", &t.dimu_dPhi, "dimu_dPhi/F");
  tree->Branch("dimu_dEta", &t.dimu_dEta, "dimu_dEta/F");

  tree->Branch("met_pt", &t.met_pt, "met_pt/F");
  tree->Branch("met_phi", &t.met_phi, "met_phi/F");
  tree->Branch("nJets", &t.nJets, "nJets/I");
  tree->Branch("nJets_lt_2p4", &t.nJets_lt_2p4, "nJets_lt_2p4/I");
  tree->Branch("nJets_gt_2p4", &t.nJets_gt_2p4, "nJets_gt_2p4/I");

  tree->Branch("jet1_pt", &t.jet1_pt, "jet1_pt/F");
  tree->Branch("jet1_energy", &t.jet1_energy, "jet1_energy/F");
  tree->Branch("jet1_eta", &t.jet1_eta, "jet1_eta/F");  
  tree->Branch("jet1_phi", &t.jet1_phi, "jet1_phi/F");
  tree->Branch("jet1_btag", &t.jet1_btag, "jet1_btag/F");   

  tree->Branch("jet2_pt", &t.jet2_pt, "jet2_pt/F");
  tree->Branch("jet2_energy", &t.jet2_energy, "jet2_energy/F");
  tree->Branch("jet2_eta", &t.jet2_eta, "jet2_eta/F");
  tree->Branch("jet2_phi", &t.jet2_phi, "jet2_phi/F"); 
  tree->Branch("jet2_btag", &t.jet2_btag, "jet2_btag/F");

  tree->Branch("jet3_pt", &t.jet3_pt, "jet3_pt/F");
  tree->Branch("jet3_energy", &t.jet3_energy, "jet3_energy/F");
  tree->Branch("jet3_eta", &t.jet3_eta, "jet3_eta/F");
  tree->Branch("jet3_phi", &t.jet3_phi, "jet3_phi/F"); 
  tree->Branch("jet3_btag", &t.jet3_btag, "jet3_btag/F");

  tree->Branch("jet4_pt", &t.jet4_pt, "jet4_pt/F");
  tree->Branch("jet4_energy", &t.jet4_energy, "jet4_energy/F");
  tree->Branch("jet4_eta", &t.jet4_eta, "jet4_eta/F");
  tree->Branch("jet4_phi", &t.jet4_phi, "jet4_phi/F"); 
  tree->Branch("jet4_btag", &t.jet4_btag, "jet4_btag/F");

  tree->Branch("j1j2_mass", &t.j1j2_mass, "j1j2_mass/F");

  tree->Branch("jet1_dpv", &t.jet1_dpv, "jet1_dpv/F");
  tree->Branch("jet2_dpv", &t.jet2_dpv, "jet2_dpv/F");
  tree->Branch("jet1_dz", &t.jet1_dz, "jet1_dz/F");
  tree->Branch("jet2_dz", &t.jet2_dz, "jet2_dz/F");
  tree->Branch("jet1_dr", &t.jet1_dr, "jet1_dr/F");
  tree->Branch("jet2_dr", &t.jet2_dr, "jet2_dr/F");

  tree->Branch("jj1_mass", &t.jj1_mass, "jj1_mass/F");
  tree->Branch("jj1_dEta", &t.jj1_dEta, "jj1_dEta/F");
  tree->Branch("jj2_mass", &t.jj2_mass, "jj2_mass/F");
  tree->Branch("jj2_dEta", &t.jj2_dEta, "jj2_dEta/F");

  tree->Branch("nJets_CSV", &t.nJets_CSV, "nJets_CSV/I");
}


void HiggsToMuMu::analyze(const edm::Event& event, const edm::EventSetup&) {
  memset(&t, 0, sizeof(tree_t));

  //
  // Event Information
  //
  t.run = event.id().run();
  t.lumi = event.luminosityBlock();
  t.event = event.id().event();

  // Get Vertex information
  edm::Handle<reco::VertexCollection> pvs;
  event.getByLabel(vertices_src, pvs);
  math::XYZPoint primary_vertex = math::XYZPoint(0.,0.,0.);
  primary_vertex = (*pvs)[0].position();
  t.nvertices = 0;
  BOOST_FOREACH(const reco::Vertex& vtx, *pvs)
    if (vtx.ndof() > 4 && fabs(vtx.z()) <= 24 && fabs(vtx.position().rho()) <= 2)
      t.nvertices += 1;
  


  //
  // Get dilepton collection
  //
  edm::Handle<pat::CompositeCandidateCollection> dils;
  event.getByLabel(dimu_src, dils);

  bool dimu_flag = dils->size() > 0;
  pat::CompositeCandidate dil = (*dils)[0];
  
  // The dils come pre-sorted so that the first in the list is the one to use

  t.mu1_pt = dil.daughter(0)->pt();
  t.mu1_eta = dil.daughter(0)->eta();
  t.mu1_phi = dil.daughter(0)->phi();
  t.mu2_pt = dil.daughter(1)->pt();
  t.mu2_eta = dil.daughter(1)->eta();
  t.mu2_phi = dil.daughter(1)->phi();

  double xx1 = primary_vertex.x() - dil.daughter(0)->vx();
  double yy1 = primary_vertex.y() - dil.daughter(0)->vy();
  double zz1 = primary_vertex.z() - dil.daughter(0)->vz();
  t.mu1_dpv = sqrt(xx1*xx1+yy1*yy1+zz1*zz1);
  t.mu1_dz = fabs(zz1);
  t.mu1_dr = sqrt(xx1*xx1+yy1*yy1);

  double xx2 = primary_vertex.x() - dil.daughter(1)->vx();
  double yy2 = primary_vertex.y() - dil.daughter(1)->vy();
  double zz2 = primary_vertex.z() - dil.daughter(1)->vz();
  t.mu2_dpv = sqrt(xx2*xx2+yy2*yy2+zz2*zz2);
  t.mu2_dz = fabs(zz2);
  t.mu2_dr = sqrt(xx2*xx2+yy2*yy2);

  t.dimu_mass = dil.mass();
  t.dimu_pt = dil.pt();
  t.dimu_eta = dil.eta();
  t.dimu_phi = dil.phi();
  t.dimu_dR = deltaR(*dil.daughter(0), *dil.daughter(1));
  t.dimu_dPhi = fabs(deltaPhi(*dil.daughter(0), *dil.daughter(1)));
  t.dimu_dEta = fabs(dil.daughter(0)->eta() - dil.daughter(1)->eta());

  // set opp_sign and diff_flavor
  const bool opp_sign = dil.daughter(0)->charge() + dil.daughter(1)->charge() == 0;
  if((!opp_sign)||(dil.numberOfDaughters() != 2)||(fabs(dil.daughter(0)->pdgId()) != 13)||(fabs(dil.daughter(1)->pdgId()) != 13)){dimu_flag = false;}
  

  // 
  // Put additional event info here
  // MET, Jets, etc.
  //

  edm::Handle< std::vector< pat::MET > > mets;
  event.getByLabel(met_src, mets);
  t.met_pt = mets->front().pt();
  t.met_phi = mets->front().phi();

  edm::Handle< std::vector< pat::Jet > > jets;
  event.getByLabel(jet_src, jets);


  int iterator = 0;
  int nJets = 0;
  int nJets_lt_2p4 = 0;
  int nJets_gt_2p4 = 0;
  t.jet1_pt = -999;
  t.jet1_energy = -999;
  t.jet1_eta = -999;
  t.jet1_phi = -999;
  t.jet1_btag = -999;
  t.jet1_dpv = -999;

  t.jet2_pt = -999;
  t.jet2_energy = -999;
  t.jet2_eta = -999;
  t.jet2_phi = -999;
  t.jet2_btag = -999;
  t.jet2_dpv = -999;

  t.jet3_pt = -999;
  t.jet3_energy = -999;
  t.jet3_eta = -999;
  t.jet3_phi = -999;
  t.jet3_btag = -999;

  t.jet4_pt = -999;
  t.jet4_energy = -999;
  t.jet4_eta = -999;
  t.jet4_phi = -999;
  t.jet4_btag = -999;

  t.j1j2_mass = -999;

  float jj1_mass = -999;
  float jj1_dEta = -999;
  float jj2_mass = -999;
  float jj2_dEta = -999;

  int nJets_CSV = 0;

  for (std::vector<pat::Jet>::const_iterator itJet = jets->begin(); itJet != jets->end(); itJet++) {
   bool flag1 = ((itJet->neutralHadronEnergyFraction()<0.99) && (itJet->neutralEmEnergyFraction()<0.99) &&
    ((itJet->chargedMultiplicity()+itJet->neutralMultiplicity())>1) && (((abs(itJet->eta())<=2.4) && (itJet->chargedHadronEnergyFraction()>0)
    && (itJet->chargedMultiplicity()>0) && (itJet->chargedEmEnergyFraction()<0.99)) || (abs(itJet->eta())>2.4)));

   bool flag2 =  ((itJet->neutralHadronEnergyFraction()<0.98) && (itJet->neutralEmEnergyFraction()>0.01) && (itJet->neutralMultiplicity()>2));
   bool flag3 =  ((itJet->neutralEmEnergyFraction()<0.90) && (itJet->neutralMultiplicity()>10));
   bool jet1_flag = (((fabs(itJet->eta())<=2.7)&&(flag1))||
                ((fabs(itJet->eta())>2.7)&&(fabs(itJet->eta())<=3.0)&&(flag2))||
                ((fabs(itJet->eta())>3.0)&&(fabs(itJet->eta())<4.7)&&(flag3)));
 
      if (jet1_flag){

        TLorentzVector jet0, jet1;

      if (iterator == 0) {
        t.jet1_pt = itJet->pt();
        t.jet1_energy = itJet->energy();
        t.jet1_eta = itJet->eta();
        t.jet1_phi = itJet->phi();
        t.jet1_btag = itJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        double xx1 = primary_vertex.x() - itJet->vx();
        double yy1 = primary_vertex.y() - itJet->vy();
        double zz1 = primary_vertex.z() - itJet->vz();
        t.jet1_dpv = sqrt(xx1*xx1+yy1*yy1+zz1*zz1);
        t.jet1_dz = fabs(zz1);
        t.jet1_dr = sqrt(xx1*xx1+yy1*yy1);
        jet0.SetPtEtaPhiE(itJet->pt(), itJet->eta(), itJet->phi(), itJet->energy());
        // jet0 = itJet->p4();
      }
      if (iterator == 1) {
        t.jet2_pt = itJet->pt();
        t.jet2_energy = itJet->energy();
        t.jet2_eta = itJet->eta();
        t.jet2_phi = itJet->phi();
        t.jet2_btag = itJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        double xx2 = primary_vertex.x() - itJet->vx();
        double yy2 = primary_vertex.y() - itJet->vy();
        double zz2 = primary_vertex.z() - itJet->vz();
        t.jet2_dpv = sqrt(xx2*xx2+yy2*yy2+zz2*zz2);
        t.jet2_dz = fabs(zz2);
        t.jet2_dr = sqrt(xx2*xx2+yy2*yy2);
        // jet1 = itJet->p4();
        jet1.SetPtEtaPhiE(itJet->pt(), itJet->eta(), itJet->phi(), itJet->energy());
        t.j1j2_mass = (jet0+jet1).M();
      }

      if (iterator == 2){
        t.jet3_pt = itJet->pt();
        t.jet3_energy = itJet->energy();
        t.jet3_eta = itJet->eta();
        t.jet3_phi = itJet->phi();
        t.jet3_btag = itJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      }

      if (iterator == 3) {
        t.jet4_pt = itJet->pt();
        t.jet4_energy = itJet->energy();
        t.jet4_eta = itJet->eta();
        t.jet4_phi = itJet->phi();
        t.jet4_btag = itJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      }



      int nSecondJet = 0; //counter for second jet in dijet pair (only those passing criteria are counted)
      // as both of the loops over jets and the selections of the jets are similar, the counters nJets and nSecondJet are in agreement
      for(std::vector<pat::Jet>::const_iterator secondJet = itJet; secondJet != jets->end(); secondJet++){
           bool flag1_2 = ((secondJet->neutralHadronEnergyFraction()<0.99) && (secondJet->neutralEmEnergyFraction()<0.99) &&
               ((secondJet->chargedMultiplicity()+secondJet->neutralMultiplicity())>1) && (((abs(secondJet->eta())<=2.4) && (secondJet->chargedHadronEnergyFraction()>0)
                && (secondJet->chargedMultiplicity()>0) && (secondJet->chargedEmEnergyFraction()<0.99)) || (abs(secondJet->eta())>2.4)));

           bool flag2_2 =  ((secondJet->neutralHadronEnergyFraction()<0.98) && (secondJet->neutralEmEnergyFraction()>0.01) && (secondJet->neutralMultiplicity()>2));
           bool flag3_2 =  ((secondJet->neutralEmEnergyFraction()<0.90) && (secondJet->neutralMultiplicity()>10));
           bool jet2_flag = (((fabs(secondJet->eta())<=2.7)&&(flag1_2))||
                ((fabs(secondJet->eta())>2.7)&&(fabs(secondJet->eta())<=3.0)&&(flag2_2))||
                ((fabs(secondJet->eta())>3.0)&&(fabs(secondJet->eta())<4.7)&&(flag3_2)));

           if (jet2_flag){
        bool sameJet= ((fabs(itJet->pt()-secondJet->pt())<0.01)&&(fabs(itJet->eta()-secondJet->eta())<0.01)&&(fabs(itJet->phi()-secondJet->phi())<0.01));
        float dijet_mass = (itJet->p4()+secondJet->p4()).M();
        if(!sameJet){
            if (dijet_mass > jj1_mass){
              jj2_mass = jj1_mass;
              jj2_dEta = jj1_dEta;

              jj1_mass = dijet_mass;
              jj1_dEta = fabs(itJet->eta() - secondJet->eta());

        
            } else if (dijet_mass>jj2_mass) {
              jj2_mass = dijet_mass;
              jj2_dEta = fabs(itJet->eta() - secondJet->eta());
            }   
          }
        nSecondJet++;
        }
      }
      iterator = iterator+1;
      if(itJet->pt()>30){
        if (itJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.8484){
          nJets_CSV++;
        }
        if (fabs(itJet->eta())<2.4){
          nJets_lt_2p4++;
        } else {
          nJets_gt_2p4++;
        }
        nJets++; 
      }
 
   }
  }
    //
  t.nJets = nJets;
  t.nJets_lt_2p4 = nJets_lt_2p4;
  t.nJets_gt_2p4 = nJets_gt_2p4;
  t.nJets_CSV = nJets_CSV;
  t.jj1_mass = jj1_mass;
  t.jj1_dEta = jj1_dEta;
  t.jj2_mass = jj2_mass;
  t.jj2_dEta = jj2_dEta;

  if(dimu_flag){
    tree->Fill();
  }
  
} // end HiggsToMuMu::analyze

DEFINE_FWK_MODULE(HiggsToMuMu);
