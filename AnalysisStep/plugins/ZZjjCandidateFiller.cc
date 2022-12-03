/** \class ZZjjCandidateFiller
 *
 *
 *  \author N. Amapane - Torino
 *  \author C. Botta - Torino
 *  \author G. Ortona - LLR
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Jet.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
//#include <ZZAnalysis/AnalysisStep/interface/ZZMassErrors.h>
//#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/CompositeCandMassResolution.h>
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include <DataFormats/GeometryVector/interface/Point3DBase.h>
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h>
#include <RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h>
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>

#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <ZZAnalysis/AnalysisStep/interface/METCorrectionHandler.h>
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include <ZZAnalysis/AnalysisStep/interface/Comparators.h>
#include <ZZAnalysis/AnalysisStep/interface/utils.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/JetCleaner.h>
#include <MelaAnalytics/GenericMEComputer/interface/GMECHelperFunctions.h>
#include <KinZfitter/KinZfitter/interface/KinZfitter.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>

using namespace zzanalysis;
using namespace BranchHelpers;

//bool doVtxFit = false;

class ZZjjCandidateFiller : public edm::EDProducer {
public:
  /// Constructor
  explicit ZZjjCandidateFiller(const edm::ParameterSet&);

  /// Destructor
  virtual ~ZZjjCandidateFiller();

private:
  typedef map<const reco::Candidate*, const pat::PFParticle*> FSRToLepMap;

  virtual void beginJob(){};
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  void getPairMass(const reco::Candidate* lp, const reco::Candidate* lm, FSRToLepMap& photons, float& mass, int& ID);
  void buildMELA();
  void computeMELABranches();
  void updateMELAClusters_Common();
  void updateMELAClusters_J1JEC();
  void updateMELAClusters_J2JEC();
  void updateMELAClusters_LepWH();
  void updateMELAClusters_LepZH();
  void pushMELABranches(pat::CompositeCandidate& myCand);
  void clearMELA();

  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > candidateToken;
  const CutSet<pat::CompositeCandidate> preBestCandSelection;
  const CutSet<pat::CompositeCandidate> cuts;
  int sampleType;
  int setup;
  float superMelaMass;
  Mela* mela;
  std::vector<std::string> recoMElist;
  std::vector<MELAOptionParser*> me_copyopts;
  std::vector<MELAHypothesis*> me_units;
  std::vector<MELAHypothesis*> me_aliased_units;
  std::vector<MELAComputation*> me_computers;
  std::vector<MELACluster*> me_clusters;
  std::vector<MELABranch*> me_branches;

  bool embedDaughterFloats;
  bool isMerged;
  //bool ZRolesByMass;
  reco::CompositeCandidate::role_collection rolesZ1Z2;
  reco::CompositeCandidate::role_collection rolesZ2Z1;
  bool isMC;
  bool recomputeIsoForFSR;
  //bool doKinFit;
  // float muon_iso_cut, electron_iso_cut;
  TH2F* corrSigmaMu;
  TH2F* corrSigmaEle;
  Comparators::ComparatorTypes bestCandType;
  KinZfitter *kinZfitter;
  edm::EDGetTokenT<double> rhoForMuToken;
  edm::EDGetTokenT<double> rhoForEleToken;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  edm::EDGetTokenT<pat::METCollection> metToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > softLeptonToken;
  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > ZCandToken;
  std::string candidateLabel;  
  int _num_of_JEC_variations;
};


ZZjjCandidateFiller::ZZjjCandidateFiller(const edm::ParameterSet& iConfig) :
  candidateToken(consumes<edm::View<reco::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("src"))),
  preBestCandSelection(iConfig.getParameter<edm::ParameterSet>("bestCandAmong")),
  cuts(iConfig.getParameter<edm::ParameterSet>("flags")),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  superMelaMass(iConfig.getParameter<double>("superMelaMass")),
//  mela(0),
  recoMElist(iConfig.getParameter<std::vector<std::string>>("recoProbabilities")),
  embedDaughterFloats(iConfig.getUntrackedParameter<bool>("embedDaughterFloats", true)),
  //ZRolesByMass(iConfig.getParameter<bool>("ZRolesByMass")),
  isMerged(iConfig.getParameter<bool>("isMerged")),
  isMC(iConfig.getParameter<bool>("isMC")),
 // doKinFit(iConfig.getParameter<bool>("doKinFit")),
  recomputeIsoForFSR(iConfig.getParameter<bool>("recomputeIsoForFSR")),
  corrSigmaMu(0),
  corrSigmaEle(0),
  kinZfitter(0)
{
  //mela = new Mela(SetupToSqrts(setup), superMelaMass, TVar::SILENT); //FIXME //aloke
  //mela->setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  produces<pat::CompositeCandidateCollection>();
  rhoForMuToken = consumes<double>(LeptonIsoHelper::getMuRhoTag(sampleType, setup));
  rhoForEleToken = consumes<double>(LeptonIsoHelper::getEleRhoTag(sampleType, setup));

  buildMELA(); 

  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("cleanJets"));
  metToken = consumes<pat::METCollection>(edm::InputTag("slimmedMETs"));
  softLeptonToken = consumes<edm::View<reco::Candidate> >(edm::InputTag("softLeptons"));
  ZCandToken = consumes<edm::View<reco::CompositeCandidate> >(edm::InputTag("ZCand"));

  rolesZ1Z2 = {"Z1", "Z2"};
  rolesZ2Z1 = {"Z2", "Z1"};

  _num_of_JEC_variations = 5; // Define number of total JEC variations 4 (JESUp, JESDn, JERUp, JERDn) + 1 for Nominal

  if (setup < 2015) {// FIXME:  EbE corrections to be updated for Run II
    // Run I ebe corrections; obsolete
    edm::FileInPath fip("ZZAnalysis/AnalysisStep/data/ebeOverallCorrections.Legacy2013.v0.root");
    std::string ebePath=fip.fullPath();

    // EbE corrections
    TFile* fCorrSigma = new TFile(ebePath.data()); // FIXME: is leaked
    std::string sigmaCorrType = (isMC?"mc":"reco");
    std::string sigmaCorrYear = "";
    if (setup==2011) sigmaCorrYear = "42x";
    else if (setup==2012) sigmaCorrYear = "53x";

    corrSigmaMu=  (TH2F*)fCorrSigma->Get(("mu_"+sigmaCorrType+sigmaCorrYear).data());
    corrSigmaEle= (TH2F*)fCorrSigma->Get(("el_"+sigmaCorrType+sigmaCorrYear).data());
  }

  string cmp=iConfig.getParameter<string>("bestCandComparator");
  if      (cmp=="byBestZ1bestZ2") bestCandType=Comparators::byBestZ1bestZ2;
  else if (cmp=="byBestZqq")      bestCandType=Comparators::byBestZqq;
  else if (cmp=="byBestKD")       bestCandType=Comparators::byBestKD;
  else if (cmp=="byBestKD_VH")    bestCandType=Comparators::byBestKD_VH;
  else if (cmp=="byBestPsig")    bestCandType=Comparators::byBestPsig;
  else if (cmp=="byMHWindow")    bestCandType=Comparators::byMHWindow;
  else abort();

  //-- kinematic refitter
  kinZfitter = new KinZfitter(!isMC);
  // No longer used, but keept for future needs
//   muon_iso_cut = iConfig.getParameter<double>("muon_iso_cut");
//   electron_iso_cut = iConfig.getParameter<double>("electron_iso_cut");
 candidateLabel = iConfig.getParameter<edm::InputTag>("src").label();
}

ZZjjCandidateFiller::~ZZjjCandidateFiller(){
  delete kinZfitter;

  clearMELA();
}


void ZZjjCandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  using namespace std;
  using namespace reco;

  auto result = std::make_unique<pat::CompositeCandidateCollection>();

  //const float ZmassValue = PDGHelpers::Zmass;


  double rhoForMu, rhoForEle;
  {
    edm::Handle<double> rhoHandle;
    iEvent.getByToken(rhoForMuToken, rhoHandle);
    rhoForMu = *rhoHandle;
    iEvent.getByToken(rhoForEleToken, rhoHandle);
    rhoForEle = *rhoHandle;
  }


  //--- JEC uncertanties for fat jets


  ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK8PFchs",JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc(JetCorPar);

  JME::JetResolution resolution_pt, resolution_phi;

  resolution_pt = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
  resolution_phi = JME::JetResolution::get(iSetup, "AK4PFchs_phi");



  // Get LLLL candidates
  Handle<edm::View<CompositeCandidate> > LLLLCands;
  iEvent.getByToken(candidateToken, LLLLCands);

  // Get jets
  Handle<edm::View<pat::Jet> > CleanJets;
  iEvent.getByToken(jetToken, CleanJets);

// FIXME aloke FIXME
  // Get MET
  Handle<pat::METCollection> metHandle;
  //Handle<vector<reco::MET> > metHandle;
  iEvent.getByToken(metToken, metHandle);
  //float PFMET = 0.;
  //float PFMETPhi = 0.;
  math::XYZTLorentzVector pfmet;
  if(metHandle.isValid()){
   // PFMET = metHandle->front().pt();
   // PFMETPhi = metHandle->front().phi();
   pfmet=metHandle->front().p4();
   //PFMET=pfmet.pt();
   //PFMET=pfmet.phi();
  }
  // FIXME: May need to correct MET in the MC and use the corrected MET to calculate the VH MEs

  // Get leptons (in order to store extra leptons)
  Handle<View<reco::Candidate> > softleptoncoll;
  iEvent.getByToken(softLeptonToken, softleptoncoll);
  vector<reco::CandidatePtr> goodisoleptonPtrs;
  for( View<reco::Candidate>::const_iterator lep = softleptoncoll->begin(); lep != softleptoncoll->end(); ++ lep ){
    if((bool)userdatahelpers::getUserFloat(&*lep,"isGood")
       //       && (bool)userdatahelpers::getUserFloat(&*lep,"isIsoFSRUncorr") // with old FSR strategy
       && (bool)userdatahelpers::getUserFloat(&*lep,"passCombRelIsoPFFSRCorr") // with new FSR strategy
       ){
      const reco::CandidatePtr lepPtr(softleptoncoll,lep-softleptoncoll->begin());
      goodisoleptonPtrs.push_back(lepPtr);
    }
  }
//  cout<<"good lep size = "<<goodisoleptonPtrs.size()<<endl; 
  //cout <<"softleptoncoll->size() "<<softleptoncoll->size()<<endl; //test aloke FIXME

  // Get Z Candidates
  Handle<View<CompositeCandidate> > ZCands;
  iEvent.getByToken(ZCandToken, ZCands);

  // Get processID
//   edm::Handle<GenEventInfoProduct> gen;
//   iEvent.getByLabel( "generator", gen );
//   int processID = gen->signalProcessID();

  // to calculate mass resolution
  CompositeCandMassResolution errorBuilder;
  errorBuilder.init(iSetup);

  vector<int> bestCandIdx(preBestCandSelection.size(),-1);
  vector<float> maxPtSum(preBestCandSelection.size(),-1);
  vector< vector<int> > preSelCands(preBestCandSelection.size());
  //----------------------------------------------------------------------
  //--- Loop over input candidates
//  cout <<"input candidate size "<<LLLLCands->size()<<endl;
  for( View<CompositeCandidate>::const_iterator cand = LLLLCands->begin(); cand != LLLLCands->end(); ++ cand ) {
 
   int icand = distance(LLLLCands->begin(),cand);

    pat::CompositeCandidate myCand(*cand);

    if (embedDaughterFloats){
      userdatahelpers::embedDaughterData(myCand);
    }

    //--- Set id of the Z1 and "Z1"/"Z2" labels. This allows to call e.g. aHiggs->daughter("Z1").
    // if ZRolesByMass is true, 'Z1' and iZ1 refer to the  Z closest to mZ; otherwise Z1 = daughter(0). The latter is used for control regions.
    const reco::CompositeCandidate::role_collection* ZRoles = &rolesZ1Z2;
//    int iZ1 = 0;
//    int iZ2 = 1;
/*    if (ZRolesByMass) {
      if(std::abs(myCand.daughter(0)->mass()-ZmassValue)>=std::abs(myCand.daughter(1)->mass()-ZmassValue)){
        swap(iZ1,iZ2);
        ZRoles = &rolesZ2Z1;
      }
    }*/
    myCand.setRoles(*ZRoles);
    myCand.applyRoles();

    //--- Z pointers
    const reco::Candidate* Z1= myCand.daughter(0);
    const reco::Candidate* Z2= myCand.daughter(1);
    vector<const reco::Candidate*> Zs = {Z1, Z2}; // in the original order

    //--- Lepton pointers in the original order
    const reco::Candidate* Z1J1= Z1->daughter(0);
    const reco::Candidate* Z1J2= Z1->daughter(1);
    const reco::Candidate* Z2L1= Z2->daughter(0);
    const reco::Candidate* Z2L2= Z2->daughter(1);
    vector<const reco::Candidate*> ZZLeps = {Z1J1,Z1J2,Z2L1,Z2L2}; // array, in the original order

    // Create corresponding array of fourmomenta; will add FSR (below)
    vector<math::XYZTLorentzVector> pij(4);
    std::transform(ZZLeps.begin(), ZZLeps.end(),pij.begin(), [](const reco::Candidate* c){return c->p4();});

    //--- Collect FSR photons and map them to the corresponding leptons
    FSRToLepMap FSRMap;
   // for (unsigned iZ=0; iZ<2; ++iZ) {
      //for (unsigned ifsr=2; ifsr<Zs[1]->numberOfDaughters(); ++ifsr) { 
      for (unsigned ifsr=2; ifsr<Z2->numberOfDaughters(); ++ifsr) {  //FIXME
    const pat::PFParticle* fsr = static_cast<const pat::PFParticle*>(Zs[1]->daughter(ifsr));
    int ilep = 2+fsr->userFloat("leptIdx");
    FSRMap[ZZLeps[ilep]]= fsr;
    pij[ilep]+=fsr->p4();
      }
    //}

    //--- Lepton four-vectors in the original order; with FSR added
    math::XYZTLorentzVector p11 = pij[0];
    math::XYZTLorentzVector p12 = pij[1];
    math::XYZTLorentzVector p21 = pij[2];
    math::XYZTLorentzVector p22 = pij[3];
    //int id11 = Z1L1->pdgId();
    //int id12 = Z1L2->pdgId();
    int id21 = Z2L1->pdgId();
    int id22 = Z2L2->pdgId();
    int candChannel = id21*id22;
    float rho=0.;
 //   if((id11 == 22 && id12 == 22) || (id21 == 22 && id22 == 22)) LogError("Z with 2 tle") << "Found a Z candidate made up of 2 trackless electrons";

//    if(id11 == 22) id11 = -1 * id12;
//    if(id12 == 22) id12 = -1 * id11;
//   if(id21 == 22) id21 = -1 * id22;
//   if(id22 == 22) id22 = -1 * id21;


    // Compute worst-lepton isolation
  //  for (int zIdx=0; zIdx<2; ++zIdx) {
      float worstMuIso=0;
      float worstEleIso=0;
      for (int dauIdx=0; dauIdx<2; ++dauIdx) {
    const reco::Candidate* z = myCand.daughter(1);
    const reco::Candidate* d = z->daughter(dauIdx);
    rho = ((d->isMuon()) ? rhoForMu : rhoForEle);
    float combRelIsoPFCorr = 0;
    if (recomputeIsoForFSR) {  //FIXME: will recompute iso for individual leptons in the new scheme
    float fsrCorr = 0; // The correction to PFPhotonIso
    for (FSRToLepMap::const_iterator ifsr=FSRMap.begin(); ifsr!=FSRMap.end(); ++ifsr) {
          double dR = ROOT::Math::VectorUtil::DeltaR(ifsr->second->p4(), d->momentum());
    // Check if the photon is in the lepton's iso cone and not vetoed
    if (dR<0.3 && ((d->isMuon() && dR > 0.01) ||
            (d->isElectron() && (fabs((static_cast<const pat::Electron*>(d->masterClone().get()))->superCluster()->eta()) < 1.479 || dR > 0.08)))) {
            fsrCorr += ifsr->second->pt();
          }
        }
        combRelIsoPFCorr =  LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, d, fsrCorr);
        string base;
        stringstream str;
        str << "d1." << "d" << dauIdx << ".";
        str >> base;
       // cout<<base+"combRelIsoPFFSRCorr"<<endl;
        myCand.addUserFloat(base+"combRelIsoPFFSRCorr", combRelIsoPFCorr);
        myCand.addUserFloat(base+"passCombRelIsoPFFSRCorr", combRelIsoPFCorr < LeptonIsoHelper::isoCut(d)); // FIXME: not the most elegant solution; hard coded right now to see how things evolve about lepton isolation requirements. 
       }

    else {
         combRelIsoPFCorr = userdatahelpers::getUserFloat(d,"combRelIsoPFFSRCorr");
      }
    if (d->isMuon()) worstMuIso  = max(worstMuIso,  combRelIsoPFCorr);
    else             worstEleIso = max(worstEleIso, combRelIsoPFCorr);
      }
      string base = "d1.";
      myCand.addUserFloat(base+"worstMuIso",worstMuIso);
      myCand.addUserFloat(base+"worstEleIso",worstEleIso);
   // }


    //----------------------------------------------------------------------
    //--- Alternative lepton pairings: "smart cut" and QCD suppression and

    //--- Sign-ordered leptons and leptopn four-vectors (without FSR), to be used to compute mZa, mZb, mZalpha, mZbeta
    const reco::Candidate* Z1Lp(Z1J1);
    const reco::Candidate* Z1Lm(Z1J2);
    const reco::Candidate* Z2Lp(Z2L1);
    const reco::Candidate* Z2Lm(Z2L2);

    // Sort leptons for OS Z candidates; no sorting for the same-sign collections used for CRs
    //if (Z1Lp->charge() < 0 && Z1Lp->charge()*Z1Lm->charge()<0) {
   //   swap(Z1Lp,Z1Lm);
   // }
    if (Z2Lp->charge() < 0 && Z2Lp->charge()*Z2Lm->charge()<0) {
      swap(Z2Lp,Z2Lm);
    }

    math::XYZTLorentzVector p1p(Z1Lp->p4());
    math::XYZTLorentzVector p1m(Z1Lm->p4());
    math::XYZTLorentzVector p2p(Z2Lp->p4());
    math::XYZTLorentzVector p2m(Z2Lm->p4());
    // Build the other SF/OS combination
    /*float mZ1= Z1->mass();
    float mZa, mZb;
    int ZaID, ZbID;
    getPairMass(Z1Lp,Z2Lm,FSRMap,mZa,ZaID);
    getPairMass(Z1Lm,Z2Lp,FSRMap,mZb,ZbID);

    // For same-sign CRs, the Z2 leptons are same sign, so we need to check also the other combination.
    float mZalpha, mZbeta;
    int ZalphaID, ZbetaID;
    getPairMass(Z1Lp,Z2Lp,FSRMap,mZalpha,ZalphaID);
    getPairMass(Z1Lm,Z2Lm,FSRMap,mZbeta,ZbetaID);

    // Sort (mZa,mZb and) (mZalpha,mZbeta) so that a and alpha are the ones closest to mZ
    if (std::abs(mZa-ZmassValue)>=std::abs(mZb-ZmassValue)) {
      swap(mZa,mZb);
      swap(ZaID,ZbID);
    }
    if (std::abs(mZalpha-ZmassValue)>=std::abs(mZbeta-ZmassValue)) {
      swap(mZalpha,mZbeta);
      swap(ZalphaID,ZbetaID);
    }*/

    // "smart cut" mll logic: veto the candidate if by swapping leptons we find a better Z1 and the Z2 is below 12 GeV.
    // To handle same-sign CRs, we have to check both alternate pairings, and consider those that have a SF/OS Z1.
    bool passSmartMLL = true;
    /*if (((ZaID==-121||ZaID==-169) && std::abs(mZa-ZmassValue)<std::abs(mZ1-ZmassValue) && mZb<12) ||
        ((ZalphaID==-121||ZalphaID==-169) && std::abs(mZalpha-ZmassValue)<std::abs(mZ1-ZmassValue) && mZbeta<12)) passSmartMLL = false;
*/

    //--- QCD suppression cut
    vector<const reco::Candidate*> lep;
    lep.push_back(Z1Lm);
    lep.push_back(Z1Lp);
    lep.push_back(Z2Lm);
    lep.push_back(Z2Lp);

  /*  float mll6 = 9999;
    float mll4 = 9999;
    for (int i=0;i<4;++i) {
      for (int j=i+1;j<4;++j) {
        float mll = (lep[i]->p4()+lep[j]->p4()).mass();
        mll6 = min(mll, mll6);
        if (lep[i]->charge()*lep[j]->charge()<0) { //OS
          mll4 = min (mll,mll4);
        }
      }
    }*/


    //--- worst SIP value
    vector<double> SIPS ={ myCand.userFloat("d1.d0.SIP"), myCand.userFloat("d1.d1.SIP") };
    sort(SIPS.begin(),SIPS.end());
    //cout<<"SIPS.size()  "<<SIPS.size()<<endl;// test aloke FIXME
    float SIP4 = SIPS[1];

    //--- Sorted pTs
    vector<pair<double, int>> ptS;
    //ptS.push_back(make_pair(Z1Lm->pt(), Z1Lm->pdgId()));
    //ptS.push_back(make_pair(Z1Lp->pt(), Z1Lp->pdgId()));
    ptS.push_back(make_pair(Z2Lm->pt(), Z2Lm->pdgId()));
    ptS.push_back(make_pair(Z2Lp->pt(),Z2Lp->pdgId()));
    sort(ptS.begin(),ptS.end());
    //cout<<" ptS.size()  "<<ptS.size()<<endl;//test aloke FIXME

    //--- Mass and Lepton uncertainties
    /*std::vector<double> errs;
    float massError = errorBuilder.getMassResolutionWithComponents(myCand, errs);
    int offset =0;
    float sigma[2][3] = {{0,0,0}, {0,0,0}};

    myCand.addUserFloat("massError",      massError);
    myCand.addUserFloat("massError11",    errs[0]);
    sigma[0][0] = errs[0];
    myCand.addUserFloat("massError12",    errs[1]);
    sigma[0][1] = errs[1];
    if (myCand.daughter(0)->numberOfDaughters()==3){
      myCand.addUserFloat("massError13",    errs[2]);
      sigma[0][2] = errs[2];
      offset = 1;
    }
    myCand.addUserFloat("massError21",    errs[2+offset]);
    sigma[1][0] = errs[2+offset];
    myCand.addUserFloat("massError22",    errs[3+offset]);
    sigma[1][1] = errs[3+offset];
    if (myCand.daughter(1)->numberOfDaughters()==3){
      myCand.addUserFloat("massError23",    errs[4+offset]);
      sigma[1][2]=errs[4+offset];
    }

    float massErrorCorr=0;
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        const reco::Candidate* l=cand->daughter(i)->daughter(j);
        const TH2F* h;
        if (l->isMuon()) h = corrSigmaMu;
        else             h = corrSigmaEle;
    float ebecorr=1.;
    if (h!=0) {
      int ptBin  = min(max(1,h->GetXaxis()->FindBin(l->pt())), h->GetNbinsX());
      int etaBin = min(max(1,h->GetYaxis()->FindBin(fabs(l->eta()))), h->GetNbinsY());
      ebecorr = h->GetBinContent(ptBin, etaBin);
    }
        massErrorCorr+= (sigma[i][j]*ebecorr)*(sigma[i][j]*ebecorr);
      }
    }
    massErrorCorr += (sigma[0][2])*(sigma[0][2]);
    massErrorCorr += (sigma[1][2])*(sigma[1][2]);
    massErrorCorr = sqrt(massErrorCorr);
    myCand.addUserFloat("massErrorCorr",      massErrorCorr);


    //--- store good isolated leptons that are not involved in the current ZZ candidate
    int nExtraLep = 0;
    SimpleParticleCollection_t associatedLeptons;
    for (vector<reco::CandidatePtr>::const_iterator lepPtr = goodisoleptonPtrs.begin(); lepPtr != goodisoleptonPtrs.end(); ++lepPtr){
      const reco::Candidate* lep = lepPtr->get();
      if (
        reco::deltaR(lep->p4(), Z1L1->p4()) > 0.02 &&
        reco::deltaR(lep->p4(), Z1L2->p4()) > 0.02 &&
        reco::deltaR(lep->p4(), Z2L1->p4()) > 0.02 &&
        reco::deltaR(lep->p4(), Z2L2->p4()) > 0.02
        ){
        nExtraLep++;
        myCand.addUserCand("ExtraLep"+to_string(nExtraLep), *lepPtr);

        SimpleParticle_t theLepton(
        lep->pdgId(),
        TLorentzVector(lep->p4().x(), lep->p4().y(), lep->p4().z(), lep->p4().t())
        );
        bool inserted=false;
        for (SimpleParticleCollection_t::iterator ielo=associatedLeptons.begin(); ielo<associatedLeptons.end(); ielo++){
          if (lep->pt()>(*ielo).second.Pt()){
            inserted=true;
            associatedLeptons.insert(ielo, theLepton);
            break;
          }
        }
        if (!inserted) associatedLeptons.push_back(theLepton);
      }
    }
    myCand.addUserFloat("nExtraLep",nExtraLep);

    // Leptonically decaying WH
    if (nExtraLep>=1){
      // Take leading-pT lepton to compute fake neutrino
      int nuid = -associatedLeptons.at(0).first + (associatedLeptons.at(0).first>0 ? -1 : +1);

      // Take neutrino momentum from the MET, using a W mass constraint to solve for the z component
      float a = associatedLeptons.at(0).second.X();
      float b = associatedLeptons.at(0).second.Y();
      float c = associatedLeptons.at(0).second.Z();
      float f = associatedLeptons.at(0).second.T();
      TLorentzVector myLep(a, b, c, f);
      float x = PFMET*cos(PFMETPhi);
      float y = PFMET*sin(PFMETPhi);
      float m = PDGHelpers::Wmass;
      float delta = pow(c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y), 2) - 4*(-4*c*c + 4*f*f)*(-a*a*a*a - 2*a*a*b*b - b*b*b*b - 2*a*a*c*c - 2*b*b*c*c - c*c*c*c + 2*a*a*f*f + 2*b*b*f*f + 2*c*c*f*f - f*f*f*f - 2*a*a*m*m - 2*b*b*m*m - 2*c*c*m*m + 2*f*f*m*m - m*m*m*m - 4*a*a*a*x - 4*a*b*b*x - 4*a*c*c*x + 4*a*f*f*x - 4*a*m*m*x - 4*a*a*x*x + 4*f*f*x*x - 4*a*a*b*y - 4*b*b*b*y - 4*b*c*c*y + 4*b*f*f*y - 4*b*m*m*y - 8*a*b*x*y - 4*b*b*y*y + 4*f*f*y*y);

      if (delta>=0.){
        float z1 = (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - sqrt(delta)) / (2*(-4*c*c + 4*f*f));
        float z2 = (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) + sqrt(delta)) / (2*(-4*c*c + 4*f*f));
        TLorentzVector myNu0(x, y, z1, sqrt(x*x+y*y+z1*z1));
        TLorentzVector myNu1(x, y, z2, sqrt(x*x+y*y+z2*z2));
        associatedLeptons.push_back(SimpleParticle_t(nuid, myNu0));
        associatedLeptons.push_back(SimpleParticle_t(nuid, myNu1));
      }
      else{
        TLorentzVector myNu(x, y, 0, TMath::Sqrt(x*x+y*y));
        associatedLeptons.push_back(SimpleParticle_t(nuid, myNu));
      }
    }


    //--- store Z candidates whose leptons are not involved in the current ZZ candidate
    int nExtraZ = 0;
    vector<const CompositeCandidate*> extraZs;
    for( View<CompositeCandidate>::const_iterator zcand = ZCands->begin(); zcand != ZCands->end(); ++ zcand ) {
      if( reco::deltaR( zcand->daughter(0)->p4(), Z1L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(0)->p4(), Z1L2->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(0)->p4(), Z2L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(0)->p4(), Z2L2->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z1L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z1L2->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z2L1->p4() ) > 0.02 &&
      reco::deltaR( zcand->daughter(1)->p4(), Z2L2->p4() ) > 0.02    ){
    const reco::CandidatePtr myZCand(ZCands,zcand-ZCands->begin());
    if((bool)userdatahelpers::getUserFloat(&*myZCand,"GoodIsoLeptons")){
      nExtraZ++;
      extraZs.push_back(&*zcand);
      myCand.addUserCand("assocZ"+to_string(nExtraZ),myZCand);
    }
      }
    }
    myCand.addUserFloat("nExtraZ",nExtraZ);

*/
    /**********************/
    /**********************/
    /***** BEGIN MELA *****/
    /**********************/
    /**********************/
/*
    // Lepton TLorentzVectors, including FSR
    SimpleParticleCollection_t daughters;
//    daughters.push_back(SimpleParticle_t(47, TLorentzVector(p11.x(), p11.y(), p11.z(), p11.t())));  //47 in  2l2q branch
 //   daughters.push_back(SimpleParticle_t(-47, TLorentzVector(p12.x(), p12.y(), p12.z(), p12.t())));

    daughters.push_back(SimpleParticle_t(0, TLorentzVector(p11.x(), p11.y(), p11.z(), p11.t()))); //0 in 2l2q_80x branch
    daughters.push_back(SimpleParticle_t(0, TLorentzVector(p12.x(), p12.y(), p12.z(), p12.t())));
    daughters.push_back(SimpleParticle_t(id21, TLorentzVector(p21.x(), p21.y(), p21.z(), p21.t())));
    daughters.push_back(SimpleParticle_t(id22, TLorentzVector(p22.x(), p22.y(), p22.z(), p22.t())));

    //--- Compute angles, better done here
    float costheta1=0, costheta2=0, phi=0, costhetastar=0, phistar1=0;
    TUtil::computeAngles(
      costhetastar, costheta1, costheta2, phi, phistar1,
      daughters.at(0).second, daughters.at(0).first,
      daughters.at(1).second, daughters.at(1).first,
      daughters.at(2).second, daughters.at(2).first,
      daughters.at(3).second, daughters.at(3).first
    );
    //--- compute higgs azimuthal angles, xi
    TLorentzVector Z14vec = daughters.at(0).second + daughters.at(1).second;
    TLorentzVector higgs = Z14vec + daughters.at(2).second + daughters.at(3).second;
    TVector3 Xaxis(1, 0, 0);
    float xi = higgs.Phi();
    // boost Z1 into rest frame of higgs
    // xistar is the angle between Z decay plane and x-axis
    Z14vec.Boost(-higgs.BoostVector());
    float xistar = Z14vec.Phi();
    // detaJJ, Mjj and Fisher. These are per-event variables in the SR, but not necessarily in the CR as we clean jets also
    // for loose-but-not-tight leptons.
    float DiJetMass  = -99;
    float DiJetDEta  = -99;
    float DiJetFisher  = -99;
    float ZZjjPt     = -99;

    unsigned int nCandidates=0; // Should equal jecnum after the loop below
    
    // Loop over following JEC variations:
    // JES
    // JER
    for (int jecnum = 0; jecnum < _num_of_JEC_variations; jecnum++){
      SimpleParticleCollection_t associated;

      vector<const pat::Jet*> cleanedJetsPt30Jec;
      vector<float> jec_ratio;
      for (edm::View<pat::Jet>::const_iterator jet = CleanJets->begin(); jet != CleanJets->end(); ++jet){
        
        // Nominal jet
        float ratio = 1.;
        float newPt = jet->pt();
        
        // calculate all JEC uncertainty up/down
        
        //JES Up uncertainty
        if (jecnum == 1 )      {
           ratio = 1. + jet->userFloat("jes_unc");
           newPt = jet->pt() * ratio;
        }
        
        //JES Down uncertainty
        else if (jecnum == 2 ) {
           ratio = 1. - jet->userFloat("jes_unc");
           newPt = jet->pt() * ratio;
        }
         
        //JER Up uncertainty
        else if (jecnum == 3 ) {
           ratio = jet->userFloat("pt_jerup") / jet->pt();
           newPt = jet->userFloat("pt_jerup");
        }
         
        //JER Down uncertainty
        else if (jecnum == 4 ) {
           ratio = jet->userFloat("pt_jerdn") / jet->pt();
           newPt = jet->userFloat("pt_jerdn");
        }
         
        // To add new uncertainties you have to:
        // 1) Update _num_of_JEC_variations
        // 2) Add their ration calculation here
        // 3) Update updateMELAClusters_J1JEC and updateMELAClusters_J2JEC functions with names for new variations
        // 4) Update RecoProbabilities.py and RecoProbabilities_minimal.py with those same names
         
        // apply pt>30GeV cut
        if (newPt<=30.) continue;
        // additional jets cleaning for loose leptons belonging to this candidate (for CRs only;
        // does nothing for the SR as jets are already cleaned with all tight isolated leptons )
        if (!jetCleaner::isGood(myCand, *jet)) continue;
         // remove jets belonging to the candidate 
         bool belongs = false;
         for (int dauIdx=0; dauIdx<2; ++dauIdx) {
         const reco::Candidate* z = myCand.daughter(1); 
         const reco::Candidate* d = z->daughter(dauIdx);
         double dR = ROOT::Math::VectorUtil::DeltaR(jet->p4(),d->momentum()); 
         if (dR < 0.05) belongs = true; 
         }
         if (belongs) continue; 
        // store jets and up/down ratio
        cleanedJetsPt30Jec.push_back(&*jet);
        jec_ratio.push_back(ratio);
      }
      if (jecnum==0 && cleanedJetsPt30Jec.size()>1){
        const pat::Jet& jet1 = *(cleanedJetsPt30Jec.at(0));
        const pat::Jet& jet2 = *(cleanedJetsPt30Jec.at(1));
        DiJetDEta = jet1.eta()-jet2.eta();
        DiJetMass = (jet1.p4()+jet2.p4()).M();
        DiJetFisher = fisher(DiJetMass, DiJetDEta);
      }
      
      vector<const pat::Jet*> cleanedJetsPt30;
      //cout<<" CleanJets->size()  "<<CleanJets->size()<<endl; //test aloke FIXME
      for (edm::View<pat::Jet>::const_iterator jet = CleanJets->begin(); jet != CleanJets->end(); ++jet){
        if (jet->pt() > 30.) cleanedJetsPt30.push_back(&*jet);
      }
      
      if(cleanedJetsPt30.size() > 1)
      {
        const pat::Jet& jet1 = *(cleanedJetsPt30.at(0));
        const pat::Jet& jet2 = *(cleanedJetsPt30.at(1));
        ZZjjPt = (Z1Lm->p4()+Z1Lp->p4()+Z2Lm->p4()+Z2Lp->p4()+jet1.p4()+jet2.p4()).pt();
      }
      
      for (unsigned int ijet = 0; ijet<cleanedJetsPt30Jec.size(); ijet++){
        TLorentzVector jet(
          cleanedJetsPt30Jec[ijet]->p4().x()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().y()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().z()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().t()*jec_ratio.at(ijet)
          );
        associated.push_back(SimpleParticle_t(0, jet));
      }
  //    for (unsigned int ilep=0; ilep<associatedLeptons.size(); ilep++) associated.push_back(associatedLeptons.at(ilep));
    //  mela->setInputEvent(&daughters, &associated, 0, 0);
    //  for (unsigned int ijet = 0; ijet<cleanedJetsPt30Jec.size(); ijet++){
      //  TLorentzVector jet(
      //    cleanedJetsPt30Jec[ijet]->p4().x()*jec_ratio.at(ijet),
      //    cleanedJetsPt30Jec[ijet]->p4().y()*jec_ratio.at(ijet),
      //    cleanedJetsPt30Jec[ijet]->p4().z()*jec_ratio.at(ijet),
     //     cleanedJetsPt30Jec[ijet]->p4().t()*jec_ratio.at(ijet)
     //     );
    //    SimpleParticleCollection_t stableTopDaughters; // Just a collection with one jet as the top
     //   stableTopDaughters.push_back(SimpleParticle_t(0, jet));
    //    mela->appendTopCandidate(&stableTopDaughters);
   //   }
 
     mela->setInputEvent(&daughters, &associated, 0, 0); nCandidates++;
    }
    computeMELABranches();
    // IMPORTANT: Reset input events at the end all calculations!
    mela->resetInputEvent();
*/
    /********************/
    /********************/
    /***** END MELA *****/
    /********************/
    /********************/

    //----------------------------------------------------------------------
    //--- Embed angular information and probabilities to build discriminants
    float costheta1=0, costheta2=0, phi=0, costhetastar=0, phistar1=0;
    float xi=0, xistar=0;

    // detaJJ, Mjj and Fisher. These are per-event variables in the SR, but not necessarily in the CR as we clean jets also 
    // for loose-but-not-tight leptons.
    
    float DiJetMass  = -99;
    float DiJetDEta  = -99;
    float DiJetFisher  = -99;
    float ZZjjPt     = -99;

    unsigned int nCandidates=0; // Should equal 3 after the loop below
    for (int jecnum = 0; jecnum < 3; jecnum++){

      // Lepton TLorentzVectors, including FSR 
      SimpleParticleCollection_t daughters;
      
       // multiplier: +1 means JEC up, -1 means JEC down
      double jecnum_multiplier = 0;
      if (jecnum==1) jecnum_multiplier = 1.;
      else if (jecnum==2) jecnum_multiplier = -1.;

       // JEC corrected masses
      float newZ1Mass = 0.;
      float newZZMass = 0.;
      float jratio=0., jratio1=0., jratio2=0.;
      const reco::Candidate* z = myCand.daughter(0);
      if (isMerged) {
	const pat::Jet* thejet = dynamic_cast <const pat::Jet*> (z->masterClone().get());
	jecUnc.setJetEta(thejet->eta());
	jecUnc.setJetPt(thejet->pt());
	float jec_unc = jecUnc.getUncertainty(true);
        jratio = 1. + jecnum_multiplier * jec_unc;
	jecUnc.setJetEta(Z1J1->eta());
	jecUnc.setJetPt(Z1J1->pt());
	float jec_unc1 = jecUnc.getUncertainty(true);
	jratio1 = 1. + jecnum_multiplier * jec_unc1;
	jecUnc.setJetEta(Z1J2->eta());
	jecUnc.setJetPt(Z1J2->pt());
	float jec_unc2 = jecUnc.getUncertainty(true);
        jratio2 = 1. + jecnum_multiplier * jec_unc2;
	newZ1Mass = jratio*myCand.userFloat("d0.ak8PFJetsCHSCorrPrunedMass");
	newZZMass = (jratio*z->p4() + myCand.daughter(1)->p4()).mass();
	
      } else {	  
	const reco::Candidate* d1 = z->daughter(0);
	const pat::Jet* thejet1 = dynamic_cast <const pat::Jet*> (d1->masterClone().get());
	float jec_unc1 = thejet1->userFloat("jes_unc");
	jratio1 = 1. + jecnum_multiplier * jec_unc1;
	const reco::Candidate* d2 = z->daughter(1);
	const pat::Jet* thejet2 = dynamic_cast <const pat::Jet*> (d2->masterClone().get());
	float jec_unc2 = thejet2->userFloat("jes_unc");
	jratio2 = 1. + jecnum_multiplier * jec_unc2;
	newZ1Mass = (jratio1*d1->p4() + jratio2*d2->p4()).mass();
	newZZMass = (jratio1*d1->p4() + jratio2*d2->p4() + myCand.daughter(1)->p4()).mass();
      }
      if (jecnum == 1) {
	myCand.addUserFloat("Z1Mass_JecUp", newZ1Mass);
	myCand.addUserFloat("ZZMass_JecUp", newZZMass);
      } else if (jecnum == 2) {
	myCand.addUserFloat("Z1Mass_JecDown", newZ1Mass);
	myCand.addUserFloat("ZZMass_JecDown", newZZMass);
      } 

      daughters.push_back(SimpleParticle_t(0, TLorentzVector(p11.x()*jratio1, p11.y()*jratio1, p11.z()*jratio1, p11.t()*jratio1)));
      daughters.push_back(SimpleParticle_t(0, TLorentzVector(p12.x()*jratio2, p12.y()*jratio2, p12.z()*jratio2, p12.t()*jratio2)));
      daughters.push_back(SimpleParticle_t(id21, TLorentzVector(p21.x(), p21.y(), p21.z(), p21.t())));
      daughters.push_back(SimpleParticle_t(id22, TLorentzVector(p22.x(), p22.y(), p22.z(), p22.t())));
      
      //--- Compute angles, better done here
      if (jecnum == 0) {
     		TUtil::computeAngles(
      		costhetastar, costheta1, costheta2, phi, phistar1,
      		daughters.at(0).second, daughters.at(0).first,
      		daughters.at(1).second, daughters.at(1).first,
      		daughters.at(2).second, daughters.at(2).first,
      		daughters.at(3).second, daughters.at(3).first
    		);
	//--- compute higgs azimuthal angles, xi
	TLorentzVector Z14vec = daughters.at(0).second + daughters.at(1).second;
	TLorentzVector higgs = Z14vec + daughters.at(2).second + daughters.at(3).second;
	TVector3 Xaxis(1, 0, 0);
	xi = higgs.Phi();
	// boost Z1 into rest frame of higgs
	// xistar is the angle between Z decay plane and x-axis 
	Z14vec.Boost(-higgs.BoostVector());
 	xistar = Z14vec.Phi();
      }
      
      SimpleParticleCollection_t associated;
      
      vector<const pat::Jet*> cleanedJetsPt30Jec;
      vector<float> jec_ratio;
    
      // MELA variables
      for (edm::View<pat::Jet>::const_iterator jet = CleanJets->begin(); jet != CleanJets->end(); ++jet){
        // calculate JEC uncertainty up/down
        float jec_unc = jet->userFloat("jes_unc");
        float ratio = 1. + jecnum_multiplier * jec_unc;
        float newPt = jet->pt() * ratio;
        // apply pt>30GeV cut
        if (newPt<=30.0) continue;
        // additional jets cleaning 
        if (!jetCleaner::isGood(myCand, *jet)) continue;
        // remove jets belonging to the candidate
        bool belongs = false;
        if (isMerged) {
          for (int dauIdx=0; dauIdx<2; ++dauIdx) {       
            const reco::Candidate* z = myCand.daughter(1);
            const reco::Candidate* d = z->daughter(dauIdx);
            double dR = ROOT::Math::VectorUtil::DeltaR(jet->p4(), d->momentum());
            if (dR < 0.4) belongs = true;
          }
          const reco::Candidate* z2 = myCand.daughter(0);
          double dR2 = ROOT::Math::VectorUtil::DeltaR(jet->p4(), z2->momentum());
          // cout << "signal merged jet pt = " << z2->pt() << " eta = " << z2->eta() << std::endl;
          if (dR2 < 0.8) belongs = true;
        } else {
          for (int theCand=0; theCand<2; ++theCand) {
            for (int dauIdx=0; dauIdx<2; ++dauIdx) {
              const reco::Candidate* z = myCand.daughter(theCand);
              const reco::Candidate* d = z->daughter(dauIdx);
              double dR = ROOT::Math::VectorUtil::DeltaR(jet->p4(), d->momentum());
              // if (theCand == 0) cout << "resolved merged jet pt = " << d->pt() << " eta = " << d->eta() << std::endl;
              if (dR < 0.4) belongs = true;
            }
          }
        }
        if (belongs) continue;
        // store jets and up/down ratio
        // cout << "other jet pt = " << jet->pt() << " eta = " << jet->eta() << std::endl;
        cleanedJetsPt30Jec.push_back(&*jet);
        jec_ratio.push_back(ratio);
      }
      if (jecnum==0 && cleanedJetsPt30Jec.size()>1){
        const pat::Jet& jet1 = *(cleanedJetsPt30Jec.at(0));
        const pat::Jet& jet2 = *(cleanedJetsPt30Jec.at(1));
        DiJetDEta = jet1.eta()-jet2.eta();
        DiJetMass = (jet1.p4()+jet2.p4()).M();
        DiJetFisher = fisher(DiJetMass, DiJetDEta);
        ZZjjPt = (Z1Lm->p4()+Z1Lp->p4()+Z2Lm->p4()+Z2Lp->p4()+jet1.p4()+jet2.p4()).pt();
      }
      for (unsigned int ijet = 0; ijet<cleanedJetsPt30Jec.size(); ijet++){
        TLorentzVector jet(
          cleanedJetsPt30Jec[ijet]->p4().x()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().y()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().z()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().t()*jec_ratio.at(ijet)
          );
        associated.push_back(SimpleParticle_t(0, jet));
      }
  //   for (unsigned int ilep=0; ilep<associatedLeptons.size(); ilep++) associated.push_back(associatedLeptons.at(ilep));
      mela->setInputEvent(&daughters, &associated, 0, 0); 
/*     
      for (unsigned int ijet = 0; ijet<cleanedJetsPt30Jec.size(); ijet++){
        TLorentzVector jet(
          cleanedJetsPt30Jec[ijet]->p4().x()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().y()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().z()*jec_ratio.at(ijet),
          cleanedJetsPt30Jec[ijet]->p4().t()*jec_ratio.at(ijet)
          );

          SimpleParticleCollection_t stableTopDaughters; // Just a collection with one jet as the top
          stableTopDaughters.push_back(SimpleParticle_t(0, jet));
          mela->appendTopCandidate(&stableTopDaughters);
      }
*/
     nCandidates++;
    }

 computeMELABranches();
 // IMPORTANT: Reset input events at the end all calculations!
 mela->resetInputEvent();
 /********************/
 /********************/
 /***** END MELA *****/
 /********************/
 /********************/



//added by aloke from 2l2q_80X branch
/*
    mela->setCurrentCandidateFromIndex(0);

    // float p0plus_VAJHU=0;
    // mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    // mela->computeP(p0plus_VAJHU, true);
    float p0minus_VAJHU=0;
    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0minus_VAJHU, true);
    float p0hplus_VAJHU=0;
    mela->setProcess(TVar::H0hplus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0hplus_VAJHU, true);
    // float p2bplus_VAJHU=0;
    // mela->setProcess(TVar::H2_g5, TVar::JHUGen, TVar::ZZGG);
    // mela->computeP(p2bplus_VAJHU, true);
    float p2_VAJHU=0;
    mela->setProcess(TVar::H2_g1g5, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p2_VAJHU, true);

    // float pqqZJJ_VAMCFM=0;
    // mela->setProcess(TVar::bkgZJets, TVar::MCFM, TVar::JJQCD);
    // mela->computeP(pqqZJJ_VAMCFM, true);
    float bkg_VAMCFM=0;
    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela->computeP(bkg_VAMCFM, true);
    float ggzz_VAMCFM=0;
    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    mela->computeP(ggzz_VAMCFM, true);
    // float p0plus_VAMCFM=0;
    // mela->setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    // mela->computeP(p0plus_VAMCFM, true);
    float ggzz_p0plus_VAMCFM=0;
    mela->setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela->computeP(ggzz_p0plus_VAMCFM, true);
    float Dgg10_VAMCFM=0;
    mela->computeD_gg(TVar::MCFM, TVar::D_gg10, Dgg10_VAMCFM);

    float pvbf_VAJHU_highestPTJets=-1;
    float phjj_VAJHU_highestPTJets=-1;
    float pvbf_VAJHU_highestPTJets_up=-1;
    float phjj_VAJHU_highestPTJets_up=-1;
    float pvbf_VAJHU_highestPTJets_dn=-1;
    float phjj_VAJHU_highestPTJets_dn=-1;
    float pvbs_VAMCFM_highestPTJets=-1;
    float pzzjj_VAMCFM_highestPTJets=-1;
    float pvbs_VAMCFM_highestPTJets_up=-1;
    float pzzjj_VAMCFM_highestPTJets_up=-1;
    float pvbs_VAMCFM_highestPTJets_dn=-1;
    float pzzjj_VAMCFM_highestPTJets_dn=-1;
    float pqqZJJ_VAMCFM=-1;
    float p0plus_VAJHU=-1;
    float pqqZJJ_VAMCFM_up=-1;
    float p0plus_VAJHU_up=-1;
    float pqqZJJ_VAMCFM_dn=-1;
    float p0plus_VAJHU_dn=-1;
    float p0plus_VAMCFM=-1;
    float p2bplus_VAJHU=-1;
    float p0plus_VAMCFM_up=-1;
    float p2bplus_VAJHU_up=-1;
    float p0plus_VAMCFM_dn=-1;
    float p2bplus_VAJHU_dn=-1;
    
    // Do these loops at the end to avoid switching particles off first and then on again
    for (unsigned int jecnum=0; jecnum<nCandidates; jecnum++){
      mela->setCurrentCandidateFromIndex(jecnum);
      MELACandidate* melaCand = mela->getCurrentCandidate();

      if (melaCand!=0){
        unsigned int nGoodJets=melaCand->getNAssociatedJets();
        bool hasAtLeastOneJet = (nGoodJets>0);
        //bool hasAtLeastTwoJets = (nGoodJets>1);
	float pqqZJJ_VAMCFM_temp=-1;
	float p0plus_VAJHU_temp=-1;
	float p0plus_VAMCFM_temp=-1;
	float p2bplus_VAJHU_temp=-1;
        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
	mela->computeP(p0plus_VAJHU_temp, true);
        mela->setProcess(TVar::H2_g5, TVar::JHUGen, TVar::ZZGG);
	mela->computeP(p2bplus_VAJHU_temp, true);
        mela->setProcess(TVar::bkgZJets, TVar::MCFM, TVar::JJQCD);
	mela->computeP(pqqZJJ_VAMCFM_temp, true);
        mela->setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
	mela->computeP(p0plus_VAMCFM_temp, true);

        if (hasAtLeastOneJet){
          float phjj_VAJHU_highestPTJets_temp = -1;
          float pvbf_VAJHU_highestPTJets_temp = -1;
          float pzzjj_VAMCFM_highestPTJets_temp = -1;
          float pvbs_VAMCFM_highestPTJets_temp = -1;
      
          for (unsigned int firstjet = 0; firstjet < nGoodJets; firstjet++){
            for (unsigned int secondjet = 1; secondjet < nGoodJets; secondjet++){
              if (secondjet<=firstjet) continue;
              for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++){
                bool flag=false;
                if (disableJet==firstjet || disableJet==secondjet) flag=true;
                melaCand->getAssociatedJet(disableJet)->setSelected(flag);
              }

              float pvbf_temp = -1;
              mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
              mela->computeProdP(pvbf_temp, true);
              float phjj_temp = -1;
              mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
              mela->computeProdP(phjj_temp, true);
              float pvbs_temp = -1;
              mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::JJEW);
              mela->computeProdP(pvbs_temp, true);
              float pzzjj_temp = -1;
              mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::JJQCD);
              mela->computeProdP(pzzjj_temp, true);
            
              if (firstjet == 0 && secondjet == 1){
                phjj_VAJHU_highestPTJets_temp = phjj_temp;
                pvbf_VAJHU_highestPTJets_temp = pvbf_temp;
                pzzjj_VAMCFM_highestPTJets_temp = pzzjj_temp;
                pvbs_VAMCFM_highestPTJets_temp = pvbs_temp;
              
              }
            }
          
          }

          for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(true); // Turn everything back on

          if (jecnum == 0){
            phjj_VAJHU_highestPTJets = phjj_VAJHU_highestPTJets_temp;
            pvbf_VAJHU_highestPTJets = pvbf_VAJHU_highestPTJets_temp;
            pzzjj_VAMCFM_highestPTJets = pzzjj_VAMCFM_highestPTJets_temp;
            pvbs_VAMCFM_highestPTJets = pvbs_VAMCFM_highestPTJets_temp;
          }
          else if (jecnum == 1){
            phjj_VAJHU_highestPTJets_up = phjj_VAJHU_highestPTJets_temp;
            pvbf_VAJHU_highestPTJets_up = pvbf_VAJHU_highestPTJets_temp;
            pzzjj_VAMCFM_highestPTJets_up = pzzjj_VAMCFM_highestPTJets_temp;
            pvbs_VAMCFM_highestPTJets_up = pvbs_VAMCFM_highestPTJets_temp;
          }
          else if (jecnum == 2){
            phjj_VAJHU_highestPTJets_dn = phjj_VAJHU_highestPTJets_temp;
            pvbf_VAJHU_highestPTJets_dn = pvbf_VAJHU_highestPTJets_temp;
            pzzjj_VAMCFM_highestPTJets_dn = pzzjj_VAMCFM_highestPTJets_temp;
            pvbs_VAMCFM_highestPTJets_dn = pvbs_VAMCFM_highestPTJets_temp;
          }
        } // if hasAtLeastOneJet
        if (jecnum == 0){
	  p0plus_VAJHU = p0plus_VAJHU_temp;
	  p0plus_VAMCFM = p0plus_VAMCFM_temp;
	  p2bplus_VAJHU = p2bplus_VAJHU_temp;
	  pqqZJJ_VAMCFM = pqqZJJ_VAMCFM_temp;
        }
        else if (jecnum == 1){
	  p0plus_VAJHU_up = p0plus_VAJHU_temp;
	  p0plus_VAMCFM_up = p0plus_VAMCFM_temp;
	  p2bplus_VAJHU_up = p2bplus_VAJHU_temp;
	  pqqZJJ_VAMCFM_up = pqqZJJ_VAMCFM_temp;
        }
        else if (jecnum == 2){
	  p0plus_VAJHU_dn = p0plus_VAJHU_temp;
	  p0plus_VAMCFM_dn = p0plus_VAMCFM_temp;
	  p2bplus_VAJHU_dn = p2bplus_VAJHU_temp;
	  pqqZJJ_VAMCFM_dn = pqqZJJ_VAMCFM_temp; 
        }

      } // End if melaCand!=0
    } // for jecnum = 0 to 2

*/
// End of aloke MELA adding 



    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    //--- kinematic refitting using Z mass constraint

    float ZZMassRefit = -1;
    float Z1MassRefit = -1;
    //float ZZMassRefitErr = -1;
    //float ZZMassUnrefitErr = -1;
    if(!isMerged){

      vector<TLorentzVector> selectedLeptons;
      vector<TLorentzVector> selectedJets;
//std::map<unsigned int, TLorentzVector> selectedFsrMap;

 for(unsigned ilep=0; ilep<4; ilep++){

        reco::Candidate* oneLep = (reco::Candidate*)ZZLeps[ilep];
        if (oneLep->hasMasterClone()) oneLep = (reco::Candidate*)oneLep->masterClone().get();
        TLorentzVector p4;
        p4.SetPxPyPzE(oneLep->px(),oneLep->py(),oneLep->pz(),oneLep->energy());

        if(FSRMap.find(ZZLeps[ilep])!=FSRMap.end()){
          pat::PFParticle fsr = *(FSRMap[ZZLeps[ilep]]);
          TLorentzVector p4fsr;
          p4fsr.SetPxPyPzE(fsr.px(),fsr.py(),fsr.pz(),fsr.energy());
          p4 += p4fsr;
        }
        if (ilep<2) selectedJets.push_back(p4);
        else selectedLeptons.push_back(p4);
      }

  // cout <<"running fine 1 "<<endl; //test aloke FIXME



/*
      vector<reco::Candidate *> selectedLeptons;
      vector<reco::Candidate *> selectedJets;

      for(unsigned ilep=0; ilep<4; ilep++){

    selectedLeptons.push_back((reco::Candidate*)(ZZLeps[ilep]->masterClone().get()));

    if(FSRMap.find(ZZLeps[ilep])!=FSRMap.end()){
      pat::PFParticle fsr = *(FSRMap[ZZLeps[ilep]]);
      TLorentzVector p4;
      p4.SetPxPyPzE(fsr.px(),fsr.py(),fsr.pz(),fsr.energy());
      selectedFsrMap[ilep] = p4;
    }
   if (ilep<2) selectedJets.push_back(p4);
        else selectedLeptons.push_back(p4);

   }
*/
/*
      kinZfitter->Setup(selectedLeptons, selectedFsrMap);
      kinZfitter->KinRefitZ();
*/
      kinZfitter->Setup2L2Q(selectedLeptons,selectedJets,resolution_pt,resolution_phi,rho);
      kinZfitter->KinRefitZlepZhad();


      /*ZZMassRefit = kinZfitter->GetRefitM4l();
      ZZMassRefitErr = kinZfitter->GetRefitM4lErrFullCov();
      ZZMassUnrefitErr = kinZfitter->GetM4lErr();
	*/
    // To get refit mZZ
     ZZMassRefit = kinZfitter->GetRefitMZZ2L2Q();
 // To get refit hadronic mZ (mjj)
     Z1MassRefit = kinZfitter->GetRefitMZhad();
     //ZZMassUnrefitErr = -1.;



  //  std::cout <<"running fully fine 2 "<<endl; //test aloke FIXME

      // four 4-vectors after refitting order by Z1_1,Z1_2,Z2_1,Z2_2
      //vector<TLorentzVector> p4 = kinZfitter->GetRefitP4s();

    }


    //----------------------------------------------------------------------
    //--- 4l vertex fits (experimental)

    //CandConstraintFit::fit(&myCand, iSetup);
    /*if (doVtxFit && abs(id11)==13 && abs(id12)==13 && abs(id21)==13 && abs(id22)==13) {

      edm::ESHandle<TransientTrackBuilder> theTTBuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

      //Creating a KinematicParticleFactory
      KinematicParticleFactoryFromTransientTrack factory;

      const ParticleMass muon_mass = 0.1056583;
      const ParticleMass electron_mass = 0.0005;
      float muon_sigma = 0.0000000001;
      float electron_sigma = 0.0000000001;

      //initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered
      float chi;
      float ndof;

      vector<RefCountedKinematicParticle> Particles;
      //vector<TransientTrack> t_tks;

      for (unsigned k = 0; k < myCand.numberOfDaughters(); ++k ) {
        const reco::Candidate* Z = myCand.daughter(k);
        for (unsigned l = 0; l < Z->numberOfDaughters(); ++l ) {
          chi = 0.; ndof = 0.;

          const reco::Candidate* lepton= Z->daughter(l);

          if (lepton->isGlobalMuon() || lepton->isTrackerMuon()){
            TransientTrack tt = theTTBuilder->build(lepton->get<TrackRef>());
            Particles.push_back(factory.particle (tt,muon_mass,chi,ndof,muon_sigma));
            //t_tks.push_back(tt);
          }
          else if (lepton->isElectron()){
            TransientTrack tt = theTTBuilder->build(lepton->get<GsfTrackRef>());
            Particles.push_back(factory.particle (tt,electron_mass,chi,ndof,electron_sigma));
            //t_tks.push_back(tt);
          }
        }
      }

      //cout << "Number of particle for constrain fitter= " << Particles.size()<< endl;

      if (Particles.size()>=4){
        KinematicParticleVertexFitter fitter;
        RefCountedKinematicTree myTree = fitter.fit(Particles);

        if ( !myTree->isEmpty()) {
          //accessing the tree components
          myTree->movePointerToTheTop();

          RefCountedKinematicParticle allLeptonsCand     = myTree->currentParticle();
          RefCountedKinematicVertex allLeptonsVertex     = myTree->currentDecayVertex();

          // if(dbg) cout << "m(" << myCand->numberOfDaughters() << "l): " << allLeptonsCand->currentState().mass() << " +- "
          //                   << sqrt(allLeptonsCand->currentState().kinematicParametersError().matrix()(6,6) ) << endl;

          reco::Vertex constrainedVertex(reco::Vertex::Point(allLeptonsVertex->position()),
                                         allLeptonsVertex->error().matrix(),
                                         allLeptonsVertex->chiSquared(),
                                         allLeptonsVertex->degreesOfFreedom(),0);

          //      if(dbg) cout << "kinematicFit vertex, ndof, chi2, prob: "
          //                   << allLeptonsVertex->position() << " , "
          //                   << allLeptonsVertex->degreesOfFreedom() << " , "
          //                   << allLeptonsVertex->chiSquared()   << " , "
          //                   << TMath::Prob(allLeptonsVertex->chiSquared(),allLeptonsVertex->degreesOfFreedom()) << endl;

          //myCand->addUserData("ConstrainedCandVtx",constrainedVertex);
          myCand.addUserFloat("CFitM",allLeptonsCand->currentState().mass());
          myCand.addUserFloat("CFitSigmaM",allLeptonsCand->currentState().kinematicParametersError().matrix()(6,6));
          myCand.addUserFloat("CFitNdof",allLeptonsVertex->degreesOfFreedom());
          myCand.addUserFloat("CFitChi2",allLeptonsVertex->chiSquared());

        } else {
          cout << " ERROR CandConstraintFit: KinematicParticleVertexFitter failed " << endl;
        }
      }
    }*/

    //----------------------------------------------------------------------
    //--- Embed variables
 //   std::cout <<"running fully fine 3 "<<endl; //test aloke FIXME
    myCand.addUserFloat("candChannel",    candChannel);
    myCand.addUserFloat("SIP4",           SIP4);
    //myCand.addUserFloat("pt1",            ptS.at(3).first); // leading-pT  wrong
   // myCand.addUserFloat("pt2",            ptS.at(2).first); // sub-leading pT wrong
    myCand.addUserFloat("pt1",            ptS.at(1).first); // leading-pT 
    myCand.addUserFloat("pt2",            ptS.at(0).first); // sub-leading pT
    //myCand.addUserFloat("pdgId2",         ptS.at(2).second); // sub-leading pT
   // myCand.addUserFloat("mZa",            mZa);
   // myCand.addUserFloat("mZb",            mZb);
   // myCand.addUserFloat("ZaID",           ZaID);
   // myCand.addUserFloat("ZbID",           ZbID);
   // myCand.addUserFloat("mZalpha",        mZalpha);
   // myCand.addUserFloat("mZbeta",         mZbeta);
   // myCand.addUserFloat("ZalphaID",       ZalphaID);
   // myCand.addUserFloat("ZbetaID",        ZbetaID);
    myCand.addUserFloat("mLL4",           Z2->mass()); // smallest mass of any AF/OS pair
    myCand.addUserFloat("mLL6",           Z2->mass());   // smallest mass of any AF/AS pair
    myCand.addUserFloat("passSmartMLL",   passSmartMLL);
    myCand.addUserFloat("costheta1",      costheta1);
    myCand.addUserFloat("costheta2",      costheta2);
    myCand.addUserFloat("phi",            phi);
    myCand.addUserFloat("costhetastar",   costhetastar);
    myCand.addUserFloat("phistar1",       phistar1);
    myCand.addUserFloat("xistar",         xistar);  //azimuthal angle of higgs in rest frame of higgs
    myCand.addUserFloat("xi",             xi);      //azimuthal angle of higgs in lab frame

    myCand.addUserFloat("m4l",            (Z1Lm->p4()+Z1Lp->p4()+Z2Lm->p4()+Z2Lp->p4()).mass()); // mass without FSR
    if(!isMerged) {
      myCand.addUserFloat("ZZMassRefit"   , ZZMassRefit);
      myCand.addUserFloat("Z1MassRefit"   , Z1MassRefit);
     // myCand.addUserFloat("ZZMassRefitErr", ZZMassRefitErr);
     // myCand.addUserFloat("ZZMassUnrefitErr", ZZMassUnrefitErr);
    }

    // Jet quantities
    myCand.addUserFloat("DiJetMass", DiJetMass);
    myCand.addUserFloat("DiJetDEta", DiJetDEta);
    myCand.addUserFloat("DiJetFisher", DiJetFisher);
    
    myCand.addUserFloat("ZZjjPt", ZZjjPt );

  //  std::cout <<"running fully fine 4 "<<endl; //test aloke FIXME
    // MELA branches
    pushMELABranches(myCand); //FIXME
  //Mela v2 added by aloke from 2l2q_80x
/*
    myCand.addUserFloat("p0plus_VAJHU", p0plus_VAJHU);
    myCand.addUserFloat("p0minus_VAJHU", p0minus_VAJHU);
    myCand.addUserFloat("p0hplus_VAJHU", p0hplus_VAJHU);
    myCand.addUserFloat("p2bplus_VAJHU", p2bplus_VAJHU);
    myCand.addUserFloat("p2_VAJHU", p2_VAJHU);
    myCand.addUserFloat("phjj_VAJHU_highestPTJets", phjj_VAJHU_highestPTJets);
    myCand.addUserFloat("pvbf_VAJHU_highestPTJets", pvbf_VAJHU_highestPTJets);
    myCand.addUserFloat("phjj_VAJHU_highestPTJets_up", phjj_VAJHU_highestPTJets_up);
    myCand.addUserFloat("pzzjj_VAMCFM_highestPTJets", pzzjj_VAMCFM_highestPTJets);
    myCand.addUserFloat("pvbs_VAMCFM_highestPTJets", pvbs_VAMCFM_highestPTJets);
    myCand.addUserFloat("pzzjj_VAMCFM_highestPTJets_up", pzzjj_VAMCFM_highestPTJets_up);
    myCand.addUserFloat("p0plus_VAJHU_up", p0plus_VAJHU_up);
    myCand.addUserFloat("p0plus_VAMCFM_up", p0plus_VAMCFM_up);
    myCand.addUserFloat("p2bplus_VAJHU_up", p2bplus_VAJHU_up);
    myCand.addUserFloat("pqqZJJ_VAMCFM_up", pqqZJJ_VAMCFM_up);
    myCand.addUserFloat("pvbf_VAJHU_highestPTJets_up", pvbf_VAJHU_highestPTJets_up);
    myCand.addUserFloat("phjj_VAJHU_highestPTJets_dn", phjj_VAJHU_highestPTJets_dn);
    myCand.addUserFloat("pvbs_VAMCFM_highestPTJets_up", pvbs_VAMCFM_highestPTJets_up);
    myCand.addUserFloat("pzzjj_VAMCFM_highestPTJets_dn", pzzjj_VAMCFM_highestPTJets_dn);
    myCand.addUserFloat("p0plus_VAJHU_dn", p0plus_VAJHU_dn);
    myCand.addUserFloat("p0plus_VAMCFM_dn", p0plus_VAMCFM_dn);
    myCand.addUserFloat("p2bplus_VAJHU_dn", p2bplus_VAJHU_dn);
    myCand.addUserFloat("pqqZJJ_VAMCFM_dn", pqqZJJ_VAMCFM_dn);
    myCand.addUserFloat("pvbf_VAJHU_highestPTJets_dn", pvbf_VAJHU_highestPTJets_dn);
    myCand.addUserFloat("pvbs_VAMCFM_highestPTJets_dn", pvbs_VAMCFM_highestPTJets_dn);

    myCand.addUserFloat("pqqZJJ_VAMCFM", pqqZJJ_VAMCFM);
    myCand.addUserFloat("bkg_VAMCFM", bkg_VAMCFM);
    myCand.addUserFloat("p0plus_VAMCFM", p0plus_VAMCFM);
    myCand.addUserFloat("ggzz_VAMCFM", ggzz_VAMCFM);
    myCand.addUserFloat("ggzz_p0plus_VAMCFM", ggzz_p0plus_VAMCFM);
    myCand.addUserFloat("Dgg10_VAMCFM", Dgg10_VAMCFM);
*/


    //--- MC matching. To be revised, cf. MuFiller, EleFiller
//     if (isMC) {
//       int refID = 25; // FIXME: handle ZZ (sigId = 23)
//       bool MC_isRight = (myCand.userFloat("d0.d0.MCParentCode")==refID &&
//                       myCand.userFloat("d0.d1.MCParentCode")==refID &&
//                       myCand.userFloat("d1.d0.MCParentCode")==refID &&
//                       myCand.userFloat("d1.d1.MCParentCode")==refID);
//       bool MC_isRightPair = false; //FIXME to be

//       myCand.addUserFloat("MC_isRight",     MC_isRight);
//       myCand.addUserFloat("MC_isRightPair", MC_isRightPair);
//     }



    //----------------------------------------------------------------------
    //--- Check if candedate passes the "bestCandAmong" selections (2011 PRL logic)
 //  std::cout <<"running fully fine 5 "<<endl; //test aloke FIXME
    int iCRname=0;
    for(CutSet<pat::CompositeCandidate>::const_iterator bca = preBestCandSelection.begin(); bca != preBestCandSelection.end(); ++bca){
      int preBestCandResult= int((*(bca->second))(myCand));

      if (preBestCandResult){
        // Fill preSelCands matrix
        preSelCands[iCRname].push_back(icand);
      }
      iCRname++;
    }
    result->push_back(myCand);
   //cout<<"preBestCandSelection.size()  "<<preBestCandSelection.size()<<endl;
 //  std::cout <<"running fully fine 5 "<<endl; //test aloke FIXME
  } // End of loop over input candidates


  //--- For each of the bestCandAmong preselections, find the best candidate and store its index (bestCandIdx)
  Comparators::BestCandComparator myComp(*result, bestCandType);
  for (int iCRname=0; iCRname<(int)preSelCands.size(); ++iCRname) {
    if (preSelCands[iCRname].size() > 0) {
      bestCandIdx[iCRname] = *std::min_element( preSelCands[iCRname].begin(), preSelCands[iCRname].end(), myComp);
    }
  }

  //--- Embed best candidate flag (must be done in a separate loop)
  for (int i = 0; i< (int)result->size(); ++i) {
    pat::CompositeCandidate& myCand = (*result)[i];

    int iCRname=0;
    for(CutSet<pat::CompositeCandidate>::const_iterator bca = preBestCandSelection.begin(); bca != preBestCandSelection.end(); ++bca){
      bool isBestCand = (i==bestCandIdx[iCRname]);
      myCand.addUserFloat(bca->first,isBestCand);
      iCRname++;
    }

    //--- Embed flags (ie cuts specified in the "flags" pset).
    //    We do this here so that isBestCand is available within the cuts.
    for(CutSet<pat::CompositeCandidate>::const_iterator cut = cuts.begin(); cut != cuts.end(); ++cut) {
      myCand.addUserFloat(cut->first,int((*(cut->second))(myCand)));
    }
  }

  iEvent.put(std::move(result));

}


void
ZZjjCandidateFiller::getPairMass(const reco::Candidate* lp, const reco::Candidate* lm, ZZjjCandidateFiller::FSRToLepMap& photons, float& mass, int& ID){
  math::XYZTLorentzVector llp4 = lp->p4()+lm->p4();
  auto lpp = photons.find(lp);
  auto lmp = photons.find(lm);
  if (lpp!=photons.end()) llp4+=lpp->second->p4();
  if (lmp!=photons.end()) llp4+=lmp->second->p4();
  mass=llp4.mass();
  ID=lp->pdgId()*lm->pdgId();
}

void ZZjjCandidateFiller::buildMELA(){
  //mela = new Mela(SetupToSqrts(setup), superMelaMass, TVar::ERROR);
  mela = new Mela(SetupToSqrts(setup), superMelaMass, TVar::SILENT);
  mela->setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  for (unsigned int it=0; it<recoMElist.size(); it++){
    MELAOptionParser* me_opt;
    // First find out if the option has a copy specification
    // These copy options will be evaulated in a separate loop
    if (recoMElist.at(it).find("Copy")!=string::npos){
      me_opt = new MELAOptionParser(recoMElist.at(it));
      me_copyopts.push_back(me_opt);
      continue;
    }

    // Create a hypothesis for each option
    MELAHypothesis* me_hypo = new MELAHypothesis(mela, recoMElist.at(it));
    me_units.push_back(me_hypo);

    me_opt = me_hypo->getOption();
    if (me_opt->isAliased()) me_aliased_units.push_back(me_hypo);

    // Create a computation for each hypothesis
    MELAComputation* me_computer = new MELAComputation(me_hypo);
    me_computers.push_back(me_computer);

    // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
    GMECHelperFunctions::addToMELACluster(me_computer, me_clusters);

    // Create the necessary branches for each computation
    // Notice that no tree is passed, so no TBranches are created.
    if (me_opt->doBranch()){
      string basename = me_opt->getName();
      if (me_opt->isGen()) basename = string("Gen_") + basename;
      MELABranch* tmpbranch;
      if (me_opt->hasPAux()){
        tmpbranch = new MELABranch(
          (TTree*)0, TString((string("pAux_") + basename).c_str()),
          me_computer->getVal(MELAHypothesis::UsePAux), me_computer
          );
        me_branches.push_back(tmpbranch);
      }
      if (me_opt->hasPConst()){
        tmpbranch = new MELABranch(
          (TTree*)0, TString((string("pConst_") + basename).c_str()),
          me_computer->getVal(MELAHypothesis::UsePConstant), me_computer
          );
        me_branches.push_back(tmpbranch);
      }
      tmpbranch = new MELABranch(
        (TTree*)0, TString((string("p_") + basename).c_str()),
        me_computer->getVal(MELAHypothesis::UseME), me_computer
        );
      me_branches.push_back(tmpbranch);
    }
  }
  // Resolve copy options
  for (unsigned int it=0; it<me_copyopts.size(); it++){
    MELAOptionParser* me_opt = me_copyopts.at(it);
    MELAHypothesis* original_hypo=0;
    MELAOptionParser* original_opt=0;
    // Find the original options
    for (unsigned int ih=0; ih<me_aliased_units.size(); ih++){
      if (me_opt->testCopyAlias(me_aliased_units.at(ih)->getOption()->getAlias())){
        original_hypo = me_aliased_units.at(ih);
        original_opt = original_hypo->getOption();
        break;
      }
    }
    if (original_opt==0) continue;
    else me_opt->pickOriginalOptions(original_opt);
    // Create a new computation for the copy options
    MELAComputation* me_computer = new MELAComputation(original_hypo);
    me_computer->setOption(me_opt);
    me_computers.push_back(me_computer);

    // The rest is the same story...
    // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
    GMECHelperFunctions::addToMELACluster(me_computer, me_clusters);

    // Create the necessary branches for each computation
    // Notice that no tree is passed, so no TBranches are created.
    if (me_opt->doBranch()){
      string basename = me_opt->getName();
      if (me_opt->isGen()) basename = string("Gen_") + basename;
      MELABranch* tmpbranch;
      if (me_opt->hasPAux()){
        tmpbranch = new MELABranch(
          (TTree*)0, TString((string("pAux_") + basename).c_str()),
          me_computer->getVal(MELAHypothesis::UsePAux), me_computer
          );
        me_branches.push_back(tmpbranch);
      }
      if (me_opt->hasPConst()){
        tmpbranch = new MELABranch(
          (TTree*)0, TString((string("pConst_") + basename).c_str()),
          me_computer->getVal(MELAHypothesis::UsePConstant), me_computer
          );
        me_branches.push_back(tmpbranch);
      }
      tmpbranch = new MELABranch(
        (TTree*)0, TString((string("p_") + basename).c_str()),
        me_computer->getVal(MELAHypothesis::UseME), me_computer
        );
      me_branches.push_back(tmpbranch);
    }
  }
  // Loop over the computations to add any contingencies to aliased hypotheses
  for (unsigned int it=0; it<me_computers.size(); it++) me_computers.at(it)->addContingencies(me_aliased_units);

  if (DEBUG_MB){
    for (unsigned int ib=0; ib<me_branches.size(); ib++) me_branches.at(ib)->Print();
    for (unsigned int icl=0; icl<me_clusters.size(); icl++) cout << "Reco ME cluster " << me_clusters.at(icl)->getName() << " is present in " << me_clusters.size() << " clusters with #Computations = " << me_clusters.at(icl)->getComputations()->size() << endl;
  }
}
void ZZjjCandidateFiller::computeMELABranches(){
  updateMELAClusters_Common(); // "Common"
  updateMELAClusters_J1JEC(); // "J1JECNominal/Up/Dn"
  updateMELAClusters_J2JEC(); // "J2JECNominal/Up/Dn"
  updateMELAClusters_LepWH(); // "LepWH"
  updateMELAClusters_LepZH(); // "LepZH"
}
// Common ME computations with index=0
void ZZjjCandidateFiller::updateMELAClusters_Common(){
  mela->setCurrentCandidateFromIndex(0);
  MELACandidate* melaCand = mela->getCurrentCandidate();
  if (melaCand==0) return;

  for (unsigned int ic=0; ic<me_clusters.size(); ic++){
    MELACluster* theCluster = me_clusters.at(ic);
    if (theCluster->getName()=="Common"){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }
}
// Common ME computations for leptonic WH: Loops over possible fake neutrinos
void ZZjjCandidateFiller::updateMELAClusters_LepWH(){
  mela->setCurrentCandidateFromIndex(0);
  MELACandidate* melaCand = mela->getCurrentCandidate();
  if (melaCand==0) return;

  int nNeutrinos = melaCand->getNAssociatedNeutrinos();
  for (int inu=0; inu<nNeutrinos; inu++){
    // Notice: Looping over Ws does not make much sense unless you have more than one lepton since the fake neutrino is already calculated from the available lepton with W mass constraint.
    // Such a loop over Ws only makes sense if there are more than one lepton in the event, but in that case, it still does not make sense to cross-match neutrinos and leptons.
    for (int disableNu=0; disableNu<nNeutrinos; disableNu++) melaCand->getAssociatedNeutrino(disableNu)->setSelected(disableNu==inu); // Disable all neutrinos other than index==inu
    for (unsigned int icl=0; icl<me_clusters.size(); icl++){ // Loop over clusters to update them
      MELACluster* theCluster = me_clusters.at(icl);
      if (theCluster->getName()=="LepWH"){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then force an update the cluster
        theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
      }
    } // End loop over clusters
  } // End loop over possible neutrinos
  // Re-enable all neutrinos
  for (int disableNu=0; disableNu<nNeutrinos; disableNu++) melaCand->getAssociatedNeutrino(disableNu)->setSelected(true);
}
// Common ME computations for leptonic ZH: Picks best Z3
void ZZjjCandidateFiller::updateMELAClusters_LepZH(){
  mela->setCurrentCandidateFromIndex(0);
  MELACandidate* melaCand = mela->getCurrentCandidate();
  if (melaCand==0) return;

  unsigned int nSortedVs = melaCand->getNSortedVs(); // Be careful, sortedV==0,1 are (guaranteed to be) the ZZ daughters! One needs to start any loop from index=2.
  const unsigned int iSortedVstart=2;
  double dZmass=100000; // Hopefully will the greatest achievable mass for a long time...
  int chosenZ=-1;
  // Choose the Z by mass closest to mZ (~equivalent to ordering by best SM ME but would be equally valid for BSM MEs as well)
  for (unsigned int iV=iSortedVstart; iV<nSortedVs; iV++){
    MELAParticle* associatedV = melaCand->getSortedV(iV);
    if (!PDGHelpers::isAZBoson(associatedV->id)) continue;
    if (!PDGHelpers::isALepton(associatedV->getDaughter(0)->id)) continue;
    if (fabs(associatedV->m()-PDGHelpers::Zmass)<dZmass){ dZmass=associatedV->m()-PDGHelpers::Zmass; chosenZ=(int)iV; }
  }
  if (chosenZ>=0){
    // Disable every associated Z boson (and not its daughters!) unless it is the chosen one
    for (unsigned int disableV=iSortedVstart; disableV<nSortedVs; disableV++){
      bool flag=(((int)disableV)==chosenZ);
      MELAParticle* einV = melaCand->getSortedV(disableV);
      if (PDGHelpers::isAZBoson(einV->id)) einV->setSelected(flag);
    }

    for (unsigned int icl=0; icl<me_clusters.size(); icl++){ // Loop over clusters to update them
      MELACluster* theCluster = me_clusters.at(icl);
      if (theCluster->getName()=="LepZH"){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then force an update the cluster
        theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
      }
    } // End loop over clusters

  } // End if chosenZ>=0
  // Re-enable every associated Z boson and its daughters unless it is the chosen one
  for (unsigned int disableV=iSortedVstart; disableV<nSortedVs; disableV++){
    bool flag=true;
    MELAParticle* einV = melaCand->getSortedV(disableV);
    if (PDGHelpers::isAZBoson(einV->id)) einV->setSelected(flag);
  }
}
// Common ME computations for JECNominal, Up and Down variations, case where ME requires 2 jets
void ZZjjCandidateFiller::updateMELAClusters_J2JEC(){
  for (int jecnum=0; jecnum < _num_of_JEC_variations; jecnum++){
    mela->setCurrentCandidateFromIndex(jecnum);
    MELACandidate* melaCand = mela->getCurrentCandidate();
    if (melaCand==0) continue;

    unsigned int nGoodJets=melaCand->getNAssociatedJets();
    for (unsigned int firstjet = 0; firstjet < nGoodJets; firstjet++){ // Loop over first jet
      for (unsigned int secondjet = firstjet+1; secondjet < nGoodJets; secondjet++){ // Loop over second jet

        // Disable jets and tops
        for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected((disableJet==firstjet || disableJet==secondjet)); // Disable the other jets
        unsigned int nDisabledStableTops=0;
        for (int itop=0; itop<melaCand->getNAssociatedTops(); itop++){
          MELATopCandidate_t* einTop = melaCand->getAssociatedTop(itop);
          if (einTop->getNDaughters()==3) einTop->setSelected(false); // All unstable tops are disabled in the loop for jets (where "jet"=="stable top") since we are looping over jecnum
          else{
            einTop->setSelected((nDisabledStableTops==firstjet || nDisabledStableTops==secondjet)); // Disable the other stable tops
            nDisabledStableTops++;
          }
        }

        for (unsigned int icl=0; icl<me_clusters.size(); icl++){ // Loop over clusters to update them
          MELACluster* theCluster = me_clusters.at(icl);
          if (
            (theCluster->getName()=="J2JECNominal" && jecnum==0) ||
            (theCluster->getName()=="J2JESUp" && jecnum==1) ||
            (theCluster->getName()=="J2JESDn" && jecnum==2) ||
            (theCluster->getName()=="J2JERUp" && jecnum==3) ||
            (theCluster->getName()=="J2JERDn" && jecnum==4)
            ){
            // Re-compute all related hypotheses first...
            theCluster->computeAll();
            // ...then force an update the cluster
            theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
          }
        } // End loop over clusters

      } // End loop over second jet
    } // End loop over first jet
    // Turn associated jets/tops back on
    for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(true); // Turn all jets back on
    for (int itop=0; itop<melaCand->getNAssociatedTops(); itop++) melaCand->getAssociatedTop(itop)->setSelected(true); // Turn all tops back on
  } // End jecnum loop
}
// Common ME computations for JECNominal, Up and Down variations, case where ME requires 1 jet
void ZZjjCandidateFiller::updateMELAClusters_J1JEC(){
  // First determine if any of the candidates has only one jet
  bool doSkip=true;
  for (int jecnum=0; jecnum < _num_of_JEC_variations; jecnum++){
    mela->setCurrentCandidateFromIndex(jecnum);
    MELACandidate* melaCand = mela->getCurrentCandidate();
    if (melaCand==0) continue;

    unsigned int nGoodJets=melaCand->getNAssociatedJets();
    doSkip = doSkip && (nGoodJets!=1);
  }
  if (doSkip) return; // If none of the candidates have exactly 1 jet, skip the computations
  for (int jecnum=0; jecnum < _num_of_JEC_variations; jecnum++){
    mela->setCurrentCandidateFromIndex(jecnum);
    MELACandidate* melaCand = mela->getCurrentCandidate();
    if (melaCand==0) continue;

    unsigned int nGoodJets=min(1, melaCand->getNAssociatedJets());
    for (unsigned int firstjet = 0; firstjet < nGoodJets; firstjet++){ // Loop over first jet

      // Disable jets
      for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(disableJet==firstjet); // Disable the other jets

      for (unsigned int icl=0; icl<me_clusters.size(); icl++){ // Loop over clusters to update them
        MELACluster* theCluster = me_clusters.at(icl);
        if (
          (theCluster->getName()=="J1JECNominal" && jecnum==0) ||
          (theCluster->getName()=="J1JESUp" && jecnum==1) ||
          (theCluster->getName()=="J1JESDn" && jecnum==2) ||
          (theCluster->getName()=="J1JERUp" && jecnum==3) ||
          (theCluster->getName()=="J1JERDn" && jecnum==4)
          ){
          // Re-compute all related hypotheses first...
          theCluster->computeAll();
          // ...then force an update the cluster
          theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
        }
      } // End loop over clusters

    } // End loop over first jet
    // Turn associated jets back on
    for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(true); // Turn all jets back on
  } // End jecnum loop
}
void ZZjjCandidateFiller::pushMELABranches(pat::CompositeCandidate& myCand){
  for (unsigned int ib=0; ib<me_branches.size(); ib++){
    // Pull...
    me_branches.at(ib)->setVal();
    // ...push...
    myCand.addUserFloat(string(me_branches.at(ib)->bname.Data()), (float)me_branches.at(ib)->getVal());
  }
  // ...then reset
  for (unsigned int ic=0; ic<me_clusters.size(); ic++) me_clusters.at(ic)->reset();
}
void ZZjjCandidateFiller::clearMELA(){
  for (unsigned int it=0; it<me_branches.size(); it++) delete me_branches.at(it);
  for (unsigned int it=0; it<me_clusters.size(); it++) delete me_clusters.at(it);
  for (unsigned int it=0; it<me_computers.size(); it++) delete me_computers.at(it);
  for (unsigned int it=0; it<me_copyopts.size(); it++) delete me_copyopts.at(it);
  //for (unsigned int it=0; it<me_aliased_units.size(); it++) delete me_aliased_units.at(it); // DO NOT DELETE THIS, WILL BE DELETED WITH me_units!
  for (unsigned int it=0; it<me_units.size(); it++) delete me_units.at(it);
  delete mela;
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZZjjCandidateFiller);

