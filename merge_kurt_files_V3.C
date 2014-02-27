//  Raghav Kunnawalkam Elayavalli
//  created: 21th Feb 2014

//  sample macro to read in forest files from a file list which is given to us through a condor script. the output file will 
//  be sent to hadoop. 

#include <iostream>
#include <fstream>
#include <TRandom.h>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"

using namespace std;

typedef std::vector<trigger::TriggerObject> trigO;

TStopwatch timer;

float triggerMatch(int trgObjSize, float* trigPhi, float* trigEta, float *trigPt, float jtphi, float jteta, float jtpt){
  float triggerPt=0;
  float closestMatch=999999.;
  
  for(int iObj=0; iObj<trgObjSize; iObj++){
    if(abs(trigPhi[iObj]-jtphi)<0.2 &&  abs(trigEta[iObj]-jteta)<0.2 && abs(trigPt[iObj]-jtpt)/jtpt<1. && abs(trigPt[iObj]-TMath::Floor(trigPt[iObj]))>0.0001){
      if(sqrt(pow((trigPhi[iObj]-jtphi),2)+pow((trigEta[iObj]-jteta),2))<closestMatch){
	closestMatch = sqrt(pow((trigPhi[iObj]-jtphi),2)+pow((trigEta[iObj]-jteta),2));
	triggerPt = trigPt[iObj];
      }
    }
  }
  return triggerPt;
}

void merge_kurt_files_V3(const int startfile=0, const int endfile=1){
  
  TH1::SetDefaultSumw2();

  timer.Start();
  //cout<<"Macro start"<<endl;

  /*
  const int N = 5;
  //Create chain
  TChain* ch[N];

  string dir[N] = {
    "hltanalysis",
    "skimanalysis",
    //"hcalNoise",
    "akPu3PFJetAnalyzer",
    //"akPu5PFJetAnalyzer",
    //"multiPhotonAnalyzer",
    "ppTrack",
    //"pfcandAnalyzer",
    //"anaMET",
    //"muonTree",
    "hiEvtAnalyzer",
  };
    
  string trees[N] = {
    "HltTree",
    "HltTree",
    //"hbhenoise",
    "t",
    //"t",
    //"photon",
    "trackTree",
    //"pfTree",
    //"metTree",
    //"HLTMuTree",
    "HiTree"
  };
  */
  double         triggerPt;

  //output file:
  TFile f(Form("trig_merge_crosscheck_purdueforests_%d.root",endfile),"RECREATE");

  //create the trees and set the branch address
  //jet tree
  int nrefe3;
  float pt3[1000];
  //float old_pt3[1000];
  float raw3[1000];
  float eta3[1000];
  float eta3_CM[1000];
  float phi3[1000];
  float chMax3[1000];
  float trkMax3[1000];
  float chSum3[1000];
  float phSum3[1000];
  float neSum3[1000];
  float trkSum3[1000];
  float phMax3[1000];
  float neMax3[1000];
  /*
  TTree *jetR3Tree = new TTree("jetR3","akPu3PF");
  jetR3Tree->Branch("nrefe",&nrefe3,"nrefe/I");
  jetR3Tree->Branch("pt",&pt3,"pt[nrefe]/F");
  //jetR3Tree->Branch("old_pt",&old_pt3,"old_pt[nrefe]/F");
  jetR3Tree->Branch("raw",&raw3,"raw[nrefe]/F");
  jetR3Tree->Branch("eta",&eta3,"eta[nrefe]/F");
  jetR3Tree->Branch("eta_CM",&eta3_CM,"eta_CM[nrefe]/F");
  jetR3Tree->Branch("phi",&phi3,"phi[nrefe]/F");
  jetR3Tree->Branch("chMax",&chMax3,"chMax[nrefe]/F");
  jetR3Tree->Branch("trkMax",&trkMax3,"trkMax[nrefe]/F");
  jetR3Tree->Branch("phMax",&phMax3,"phMax[nrefe]/F");
  jetR3Tree->Branch("neMax",&neMax3,"neMax[nrefe]/F");
  jetR3Tree->Branch("chSum",&chSum3,"chSum[nrefe]/F");
  jetR3Tree->Branch("phSum",&phSum3,"phSum[nrefe]/F");
  jetR3Tree->Branch("neSum",&neSum3,"neSum[nrefe]/F");
  jetR3Tree->Branch("trkSum",&trkSum3,"trkSum[nrefe]/F");
  */

  //event tree
  int evt;
  int run;
  int lumi;
  int hiBin;
  float vx;
  float vy;
  float vz;
  int hiNtracks;
  float hiHFminus;
  float hiHFplus;
  float hiHFplusEta4;
  float hiHFminusEta4;
  int pPAcollisionEventSelectionPA;
  int pHBHENoiseFilter;
  int pprimaryvertexFilter;
  int pVertexFilterCutGplus;
  /*
  TTree *evtTree = new TTree("evt","evt");
  evtTree->Branch("evt",&evt,"evt/I");
  evtTree->Branch("run",&run,"run/I");
  evtTree->Branch("lumi",&lumi,"lumi/I");
  evtTree->Branch("hiBin",&hiBin,"hiBin/I");
  evtTree->Branch("vx",&vx,"vx/F");
  evtTree->Branch("vy",&vy,"vy/F");
  evtTree->Branch("vz",&vz,"vz/F");
  evtTree->Branch("hiNtracks",&hiNtracks,"hiNtracks/I");
  evtTree->Branch("hiHFminus",&hiHFplus,"hiHFplus/F");
  evtTree->Branch("hiHFplus",&hiHFminus,"hiHFminus/F");
  evtTree->Branch("hiHFplusEta4",&hiHFplusEta4,"hiHFplusEta4/F");
  evtTree->Branch("hiHFminusEta4",&hiHFminusEta4,"hiHFminusEta4/F");
  */
  /* //for the leading,subleading and sub sub leading jets information:  
  float ptlead;
  float etalead;
  float philead;
  float rawlead;
  float ptsublead;
  float etasublead;
  float phisublead;
  float rawsublead;
  float ptsubsublead;
  float etasubsublead;
  float phisubsublead;
  float rawsubsublead;
  evtTree->Branch("ptlead",&ptlead,"ptlead/F");
  evtTree->Branch("rawlead",&rawlead,"rawlead/F");
  evtTree->Branch("etalead",&etalead,"etalead/F");
  evtTree->Branch("philead",&philead,"philead/F");
  evtTree->Branch("ptsublead",&ptsublead,"ptsublead/F");
  evtTree->Branch("rawsublead",&rawsublead,"rawsublead/F");
  evtTree->Branch("etasublead",&etasublead,"etasublead/F");
  evtTree->Branch("phisublead",&phisublead,"phisublead/F");
  evtTree->Branch("ptsubsublead",&ptsubsublead,"ptsubsublead/F");
  evtTree->Branch("rawsubsublead",&rawsubsublead,"rawsubsublead/F");
  evtTree->Branch("etasubsublead",&etasubsublead,"etasubsublead/F");
  evtTree->Branch("phisubsublead",&phisubsublead,"phisubsublead/F");
*/

  //track Tree
  int nTrack;
  float trkPt[10000];
  float trkEta[10000];
  float trkPhi[10000];
  float highPurity[10000];
  float trkDz1[10000];
  float trkDzError1[10000];
  float trkDxy1[10000];
  float trkDxyError1[10000];
  float trkPtError1[10000];
  /*
  TTree *trkTree = new TTree("trk","Track");
  trkTree->Branch("nTrack",&nTrack,"nTrack/I");
  trkTree->Branch("trkPt",&trkPt,"trkPt[nTrack]/F");
  trkTree->Branch("trkEta",&trkEta,"trkEta[nTrack]/F");
  trkTree->Branch("trkPhi",&trkPhi,"trkPhi[nTrack]/F");
  trkTree->Branch("highPurity",&highPurity,"highPurity[nTrack]/F");
  trkTree->Branch("trkDz1",&trkDz1,"trkDz1[nTrack]/F");
  trkTree->Branch("trkDzError1",&trkDzError1,"trkDzError1[nTrack]/F");
  trkTree->Branch("trkDxy1",&trkDxy1,"trkDxy1[nTrack]/F");
  trkTree->Branch("trkDxyError1",&trkDxyError1,"trkDxyError1[nTrack]/F");
  trkTree->Branch("trkPtError1",&trkPtError1,"trkPtError1[nTrack]/F");
  */
  //trigger tree
  int L1_MB;
  int L1_MB_p;
  int jetMB;
  int jet20;
  int jet40;
  int jet60;
  int jet80;
  int jet100;
  int jet120;
  int jetMB_p;
  int jet20_p;
  int jet40_p;
  int jet60_p;
  int jet80_p;
  int jet100_p;
  int jet120_p;
  
  int trigObjSize[5];

  float trigObjId[5][1000];
  float trigObjPt[5][1000];
  float trigObjEta[5][1000];
  float trigObjPhi[5][1000];
  float trigObjMass[5][1000];
  /*
  float trigObj_20_mass[1000];
  float trigObj_40_id[1000];
  float trigObj_40_pt[1000];
  float trigObj_40_eta[1000];
  float trigObj_40_phi[1000];
  float trigObj_40_mass[1000];
  float trigObj_60_id[1000];
  float trigObj_60_pt[1000];
  float trigObj_60_eta[1000];
  float trigObj_60_phi[1000];
  float trigObj_60_mass[1000];
  float trigObj_80_id[1000];
  float trigObj_80_pt[1000];
  float trigObj_80_eta[1000];
  float trigObj_80_phi[1000];
  float trigObj_80_mass[1000];
  float trigObj_100_id[1000];
  float trigObj_100_pt[1000];
  float trigObj_100_eta[1000];
  float trigObj_100_phi[1000];
  float trigObj_100_mass[1000];
  float trigObj_120_id[1000];
  float trigObj_120_pt[1000];
  float trigObj_120_eta[1000];
  float trigObj_120_phi[1000];
  float trigObj_120_mass[1000]; 
  */
  /*
  TTree *hltTree = new TTree("hlt","HLT");
  hltTree->Branch("jetMB",&jetMB,"jetMB/I");
  hltTree->Branch("jet20",&jet20,"jet20/I");
  hltTree->Branch("jet40",&jet40,"jet40/I");
  hltTree->Branch("jet60",&jet60,"jet60/I");
  hltTree->Branch("jet80",&jet80,"jet80/I");
  hltTree->Branch("jet100",&jet100,"jet100/I");
  hltTree->Branch("jet120",&jet120,"jet120/I");
  hltTree->Branch("jetMB_p",&jetMB_p,"jetMB_p/I");
  hltTree->Branch("jet20_p",&jet20_p,"jet20_p/I");
  hltTree->Branch("jet40_p",&jet40_p,"jet40_p/I");
  hltTree->Branch("jet60_p",&jet60_p,"jet60_p/I");
  hltTree->Branch("jet80_p",&jet80_p,"jet80_p/I");
  hltTree->Branch("jet100_p",&jet100_p,"jet100_p/I");
  hltTree->Branch("jet120_p",&jet120_p,"jet120_p/I");
  hltTree->Branch("trigObj_20_size",&trigObj_20_size,"trigObj_20_size/I");
  hltTree->Branch("trigObj_20_id",&trigObj_20_id,"trigObj_20_id[trigObj_20_size]/F");
  hltTree->Branch("trigObj_20_pt",&trigObj_20_pt,"trigObj_20_pt[trigObj_20_size]/F");
  hltTree->Branch("trigObj_20_eta",&trigObj_20_eta,"trigObj_20_eta[trigObj_20_size]/F");
  hltTree->Branch("trigObj_20_phi",&trigObj_20_phi,"trigObj_20_phi[trigObj_20_size]/F");
  hltTree->Branch("trigObj_20_mass",&trigObj_20_mass,"trigObj_20_mass[trigObj_20_size]/F");
  hltTree->Branch("trigObj_40_size",&trigObj_40_size,"trigObj_40_size/I");
  hltTree->Branch("trigObj_40_id",&trigObj_40_id,"trigObj_40_id[trigObj_40_size]/F");
  hltTree->Branch("trigObj_40_pt",&trigObj_40_pt,"trigObj_40_pt[trigObj_40_size]/F");
  hltTree->Branch("trigObj_40_eta",&trigObj_40_eta,"trigObj_40_eta[trigObj_40_size]/F");
  hltTree->Branch("trigObj_40_phi",&trigObj_40_phi,"trigObj_40_phi[trigObj_40_size]/F");
  hltTree->Branch("trigObj_40_mass",&trigObj_40_mass,"trigObj_40_mass[trigObj_40_size]/F");
  hltTree->Branch("trigObj_60_size",&trigObj_60_size,"trigObj_60_size/I");
  hltTree->Branch("trigObj_60_id",&trigObj_60_id,"trigObj_60_id[trigObj_60_size]/F");
  hltTree->Branch("trigObj_60_pt",&trigObj_60_pt,"trigObj_60_pt[trigObj_60_size]/F");
  hltTree->Branch("trigObj_60_eta",&trigObj_60_eta,"trigObj_60_eta[trigObj_60_size]/F");
  hltTree->Branch("trigObj_60_phi",&trigObj_60_phi,"trigObj_60_phi[trigObj_60_size]/F");
  hltTree->Branch("trigObj_60_mass",&trigObj_60_mass,"trigObj_60_mass[trigObj_60_size]/F");
  hltTree->Branch("trigObj_80_size",&trigObj_80_size,"trigObj_80_size/I");
  hltTree->Branch("trigObj_80_id",&trigObj_80_id,"trigObj_80_id[trigObj_80_size]/F");
  hltTree->Branch("trigObj_80_pt",&trigObj_80_pt,"trigObj_80_pt[trigObj_80_size]/F");
  hltTree->Branch("trigObj_80_eta",&trigObj_80_eta,"trigObj_80_eta[trigObj_80_size]/F");
  hltTree->Branch("trigObj_80_phi",&trigObj_80_phi,"trigObj_80_phi[trigObj_80_size]/F");
  hltTree->Branch("trigObj_80_mass",&trigObj_80_mass,"trigObj_80_mass[trigObj_80_size]/F");
  hltTree->Branch("trigObj_100_size",&trigObj_100_size,"trigObj_100_size/I");
  hltTree->Branch("trigObj_100_id",&trigObj_100_id,"trigObj_100_id[trigObj_100_size]/F");
  hltTree->Branch("trigObj_100_pt",&trigObj_100_pt,"trigObj_100_pt[trigObj_100_size]/F");
  hltTree->Branch("trigObj_100_eta",&trigObj_100_eta,"trigObj_100_eta[trigObj_100_size]/F");
  hltTree->Branch("trigObj_100_phi",&trigObj_100_phi,"trigObj_100_phi[trigObj_100_size]/F");
  hltTree->Branch("trigObj_100_mass",&trigObj_100_mass,"trigObj_100_mass[trigObj_100_size]/F");
  hltTree->Branch("trigObj_120_size",&trigObj_120_size,"trigObj_120_size/I");
  hltTree->Branch("trigObj_120_id",&trigObj_120_id,"trigObj_120_id[trigObj_120_size]/F");
  hltTree->Branch("trigObj_120_pt",&trigObj_120_pt,"trigObj_120_pt[trigObj_120_size]/F");
  hltTree->Branch("trigObj_120_eta",&trigObj_120_eta,"trigObj_120_eta[trigObj_120_size]/F");
  hltTree->Branch("trigObj_120_phi",&trigObj_120_phi,"trigObj_120_phi[trigObj_120_size]/F");
  hltTree->Branch("trigObj_120_mass",&trigObj_120_mass,"trigObj_120_mass[trigObj_120_size]/F");  
  */
  trigO *HLT_PAJet_NoJetID_v1_trigObject[6];
  for(int i=0; i<5; i++){
    HLT_PAJet_NoJetID_v1_trigObject[i] = new trigO;
  } 
    
  //cout<<"finished setting up the branches for the new trees"<<endl;
  
  std::string infile;
  
  //change the following file list to the one required. 
  //infile = "kurt_small_filelist.txt"; 
  infile = "PAHighPtPurdueForest_4th.txt";
  
  std::ifstream instr(infile.c_str(), std::ifstream::in);
  
  std::string filename;
  int nFiles=endfile - startfile;
  
  cout<<"Running on "<<nFiles<<" forest files"<<endl;
  
  //for(int i = 0;i<N;i++){
  //  ch[i] = new TChain(string(dir[i]+"/"+trees[i]).data());
  //}
  
  //just to read the files till the start number. 
  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;  
  for(int ifile=0;ifile<startfile;ifile++){
    instr >> filename;
  }

  TH1F* hpPb_Jet100 = new TH1F("hpPb_Jet100","",1000,0,1000);
  TH1F* hpPb_Jet80 = new TH1F("hpPb_Jet80","",1000,0,1000);
  TH1F* hpPb_JetMB = new TH1F("hpPb_JetMB","",1000,0,1000);
  TH1F* hpPb_Comb = new TH1F("hpPb_Comb","",1000,0,1000);

  TH1F* hpPb_Trk40_60 = new TH1F("hpPb_Trk40_60","",1000,0,1000);
  TH1F* hpPb_Trk60_75 = new TH1F("hpPb_Trk60_75","",1000,0,1000);
  TH1F* hpPb_Trk75_95 = new TH1F("hpPb_Trk75_95","",1000,0,1000);
  TH1F* hpPb_Trk95_120 = new TH1F("hpPb_Trk95_120","",1000,0,1000);
  TH1F* hpPb_Trk120 = new TH1F("hpPb_Trk120","",1000,0,1000);
  TH1F* hpPb_TrkComb = new TH1F("hpPb_TrkComb","",1000,0,1000);

  TH1F* hpPb_Kurt100 = new TH1F("hpPb_Kurt100","",1000,0,1000);
  TH1F* hpPb_Kurt80_100 = new TH1F("hpPb_Kurt80_100","",1000,0,1000);
  TH1F* hpPb_Kurt60_80 = new TH1F("hpPb_Kurt60_80","",1000,0,1000);
  TH1F* hpPb_Kurt40_60 = new TH1F("hpPb_Kurt40_60","",1000,0,1000);
  TH1F* hpPb_Kurt20_40 = new TH1F("hpPb_Kurt20_40","",1000,0,1000);
  TH1F* hpPb_KurtComb = new TH1F("hpPb_KurtComb","",1000,0,1000);

  Float_t N_MB_pPb = 2.6026e10; //taken from the merged_MinBiasCentrality_Histo.root
  
  //now we are taking only the files from the given start number to the end number. 
  for(int ifile=startfile;ifile<endfile;ifile++){
    instr >> filename;
    
    cout<<"File: "<<filename<<std::endl;
    
    //cout<<" i = "<<i<<endl;
    
    //ch[i]->Add(filename.c_str());
    //cout << "Tree loaded  " << string(dir[i]+"/"+trees[i]).data() << endl;
    //cout << "Entries : " << ch[i]->GetEntries() << endl;
    
    TFile *fin = TFile::Open(filename.c_str());
    
    TTree* jetTree = (TTree*)fin->Get("akPu3PFJetAnalyzer/t");
    TTree* skimTree_in = (TTree*)fin->Get("skimanalysis/HltTree");
    TTree* trackTree_in = (TTree*)fin->Get("ppTrack/trackTree");
    TTree* evtTree_in = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
    TTree* hltTree_in = (TTree*)fin->Get("hltanalysis/HltTree");
    
    jetTree->AddFriend(skimTree_in);
    jetTree->AddFriend(trackTree_in);
    jetTree->AddFriend(evtTree_in);
    jetTree->AddFriend(hltTree_in);

    cout<<"start setting branch address"<<endl;

    //set the branch addresses:  - one of the most boring parts of the code: 
    jetTree->SetBranchAddress("evt",&evt);
    jetTree->SetBranchAddress("run",&run);
    jetTree->SetBranchAddress("lumi",&lumi);
    jetTree->SetBranchAddress("hiBin",&hiBin);
    jetTree->SetBranchAddress("vz",&vz);
    jetTree->SetBranchAddress("vx",&vx);
    jetTree->SetBranchAddress("vy",&vy);
    jetTree->SetBranchAddress("hiNtracks",&hiNtracks);
    jetTree->SetBranchAddress("hiHFminus",&hiHFminus);
    jetTree->SetBranchAddress("hiHFplus",&hiHFplus);
    jetTree->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
    jetTree->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
    jetTree->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
    jetTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
    jetTree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
    jetTree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);

    jetTree->SetBranchAddress("nref",&nrefe3);
    jetTree->SetBranchAddress("jtpt",&pt3);
    jetTree->SetBranchAddress("jteta",&eta3);
    jetTree->SetBranchAddress("jtphi",&phi3);
    jetTree->SetBranchAddress("rawpt",&raw3);
    jetTree->SetBranchAddress("chargedMax",&chMax3);
    jetTree->SetBranchAddress("chargedSum",&chSum3);
    jetTree->SetBranchAddress("trackMax",&trkMax3);
    jetTree->SetBranchAddress("trackSum",&trkSum3);
    jetTree->SetBranchAddress("photonMax",&phMax3);
    jetTree->SetBranchAddress("photonSum",&phSum3);
    jetTree->SetBranchAddress("neutralMax",&neMax3);
    jetTree->SetBranchAddress("neutralSum",&neSum3);

    jetTree->SetBranchAddress("nTrk",&nTrack);
    jetTree->SetBranchAddress("trkPt",&trkPt);
    jetTree->SetBranchAddress("trkEta",&trkEta);
    jetTree->SetBranchAddress("trkPhi",&trkPhi);
    jetTree->SetBranchAddress("highPurity",&highPurity);
    jetTree->SetBranchAddress("trkDz1",&trkDz1);
    jetTree->SetBranchAddress("trkDzError1",&trkDzError1);
    jetTree->SetBranchAddress("trkDxy1",&trkDxy1);
    jetTree->SetBranchAddress("trkDxyError1",&trkDxyError1);

    jetTree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB);
    jetTree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p);
    jetTree->SetBranchAddress("L1_ZeroBias",&L1_MB);
    jetTree->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p);
    jetTree->SetBranchAddress("HLT_PAJet20_NoJetID_v1",&jet20);
    jetTree->SetBranchAddress("HLT_PAJet20_NoJetID_v1_Prescl",&jet20_p);
    jetTree->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jet40);
    jetTree->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jet40_p);
    jetTree->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jet60);
    jetTree->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jet60_p);
    jetTree->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jet80);
    jetTree->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jet80_p);
    jetTree->SetBranchAddress("HLT_PAJet100_NoJetID_v1",&jet100);
    jetTree->SetBranchAddress("HLT_PAJet100_NoJetID_v1_Prescl",&jet100_p);
    jetTree->SetBranchAddress("HLT_PAJet120_NoJetID_v1",&jet120);
    jetTree->SetBranchAddress("HLT_PAJet120_NoJetID_v1_Prescl",&jet120_p);
    jetTree->SetBranchAddress("HLT_PAJet20_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[0]);
    jetTree->SetBranchAddress("HLT_PAJet40_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[1]);
    jetTree->SetBranchAddress("HLT_PAJet60_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[2]);
    jetTree->SetBranchAddress("HLT_PAJet80_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[3]);
    jetTree->SetBranchAddress("HLT_PAJet100_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[4]);
    //jetTree->SetBranchAddress("HLT_PAJet120_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[5]);

    //now that we have all the branch addresses set! 
    //lets start the event loop process. 
    cout<<"finished setting branch address"<<endl;


    //number convention is going to change here to better understand the merging. 
    // 0 - MB
    // 1 - 20
    // 2 - 40
    // 3 - 60
    // 4 - 80
    // 5 - 100
    // 6 - 120
    
    //before we go into the event and jet loop lets get the histogram for the 12-003 merging which we can just use the project. 
    TH1F *hppb0 = new TH1F("hppb0","",1000,0,1000);
    TH1F *hppb4 = new TH1F("hppb4","",1000,0,1000);
    TH1F *hppb5 = new TH1F("hppb5","",1000,0,1000);
    TH1F *hppbComb = new TH1F("hppbComb","",1000,0,1000);
    
    //histogams for track merging 12-017
    TH1F* hppb_trk40_60 = new TH1F("hppb_trk40_60","",1000,0,1000);
    TH1F* hppb_trk60_75 = new TH1F("hppb_trk60_75","",1000,0,1000);
    TH1F* hppb_trk75_95 = new TH1F("hppb_trk75_95","",1000,0,1000);
    TH1F* hppb_trk95_120 = new TH1F("hppb_trk95_120","",1000,0,1000);
    TH1F* hppb_trk120 = new TH1F("hppb_trk120","",1000,0,1000);
    TH1F* hppb_trkComb = new TH1F("hppb_trkComb","",1000,0,1000);
    
    //histograms for kurt's merging 
    TH1F* hppb_kurt100 = new TH1F("hppb_kurt100","",1000,0,1000);
    TH1F* hppb_kurt80_100 = new TH1F("hppb_kurt80_100","",1000,0,1000);
    TH1F* hppb_kurt60_80 = new TH1F("hppb_kurt60_80","",1000,0,1000);
    TH1F* hppb_kurt40_60 = new TH1F("hppb_kurt40_60","",1000,0,1000);
    TH1F* hppb_kurt20_40 = new TH1F("hppb_kurt20_40","",1000,0,1000);
    TH1F* hppb_kurtComb = new TH1F("hppb_kurtComb","",1000,0,1000);
  
    //this is going to be problematic since we dont have eta_CM here. 
    /*
    TCut eventcut = "fabs(vz)<15&&pHBHENoiseFilter&&pprimaryvertexFilter&&pPAcollisionEventSelectionPA&&pVertexFilterCutGplus";
    TCut ppb0 = "abs(eta_CM)<1&&jetMB&&!jet80&&!jet100&&raw>20";
    TCut ppb1 = "abs(eta_CM)<1&&jet20&&!jet80&&!jet100&&raw>20";
    TCut ppb2 = "abs(eta_CM)<1&&jet40&&!jet80&&!jet100&&raw>20";
    TCut ppb3 = "abs(eta_CM)<1&&jet60&&!jet80&&!jet100&&raw>20";
    TCut ppb4 = "abs(eta_CM)<1&&jet80&&!jet100&&raw>20";//try with or without 60 first and then try without 60
    TCut ppb5 = "abs(eta_CM)<1&&jet100&&raw>20";
    
    jetTree->Project("hppb0","jtpt","jetMB_p*L1_MB_p"*ppb0);
    jetTree->Project("hppb4","jtpt",ppb4);
    jetTree->Project("hppb5","jtpt",ppb5);

    hpPb_Jet100->Add(hppb5);
    hpPb_Jet80->Add(hppb4);
    hpPb_JetMb->Add(hppb0);
    hppbComb->Add(hppb0);
    hppbComb->Add(hppb4);
    hppbComb->Add(hppb5);
    hpPb_Comb->Add(hppbComb);
    */
      
    Long64_t nentries = jetTree->GetEntries();
    cout<<"nentries = "<<nentries<<endl;

    for (int i = 0; i < nentries; i++){ //start of event loop. 
    //for(int i = 0;i<100;i++){ //event loop test 
    
      jetTree->GetEntry(i);
      cout<<vz<<endl;
      if(i%1000==0)cout<<"event = "<<i<<"; run = "<<run<<endl;

         
      if(!pHBHENoiseFilter || !pprimaryvertexFilter || !pPAcollisionEventSelectionPA) continue;
      if(!pVertexFilterCutGplus) continue;
      //if(vz>15. || vz<-15.) continue;
      if(fabs(vz)<15)cout<<"hi inside vz<15"<<endl;
      //cout<<"hi"<<endl;
      for(int j=0; j<5; j++){ trigObjSize[j] = HLT_PAJet_NoJetID_v1_trigObject[j]->size();}
      
      //Fill the trigger Pt/Eta/Phi from the TriggerObjects
      for(int ii=0; ii<5; ii++){
	for(unsigned int iObj=0; iObj<trigObjSize[ii]; iObj++){
	  trigObjPt[ii][iObj] = HLT_PAJet_NoJetID_v1_trigObject[ii]->at(iObj).pt();
	  trigObjEta[ii][iObj] = HLT_PAJet_NoJetID_v1_trigObject[ii]->at(iObj).eta();
	  trigObjPhi[ii][iObj] = HLT_PAJet_NoJetID_v1_trigObject[ii]->at(iObj).phi();
	}
      }
      
      /*
      //cout<<"A"<<endl;
      
      trigObj_20_size = HLT_PAJet_NoJetID_v1_trigObject[0]->size();
      
      //cout<<" = "<<trigObj_20_size<<endl;
      for (int j = 0; j < trigObj_20_size; j++){
        trigObj_20_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).id();
        trigObj_20_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).pt();
        trigObj_20_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).eta();
        trigObj_20_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).phi();
        trigObj_20_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).mass();
      }

      //cout<<"B"<<endl;

      trigObj_40_size = HLT_PAJet_NoJetID_v1_trigObject[1]->size();
      //cout<<" = "<<trigObj_40_size<<endl;
      for (int j = 0; j < trigObj_40_size; j++){
	trigObj_40_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).id();
	trigObj_40_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).pt();
	trigObj_40_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).eta();
	trigObj_40_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).phi();
	trigObj_40_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).mass();
      }

      //cout<<"C"<<endl;

      trigObj_60_size = HLT_PAJet_NoJetID_v1_trigObject[2]->size();
      //cout<<"trigObj_60_size = "<<trigObj_60_size<<endl;
      for (int j = 0; j < trigObj_60_size; j++){
	trigObj_60_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).id();
	trigObj_60_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).pt();
	trigObj_60_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).eta();
	trigObj_60_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).phi();
	trigObj_60_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).mass();
      }

      //cout<<"D"<<endl;

      trigObj_80_size = HLT_PAJet_NoJetID_v1_trigObject[3]->size();
      //cout<<" = "<<trigObj_80_size<<endl;
      for (int j = 0; j < trigObj_80_size; j++){
	trigObj_80_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).id();
	trigObj_80_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).pt();
	trigObj_80_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).eta();
	trigObj_80_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).phi();
	trigObj_80_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).mass();
      }

      //cout<<"E"<<endl;

      trigObj_100_size = HLT_PAJet_NoJetID_v1_trigObject[4]->size();
      for (int j = 0; j < trigObj_100_size; j++){
	trigObj_100_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).id();
	trigObj_100_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).pt();
	trigObj_100_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).eta();
	trigObj_100_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).phi();
	trigObj_100_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).mass();
      }

      //cout<<"F"<<endl;

      trigObj_120_size = HLT_PAJet_NoJetID_v1_trigObject[5]->size();
      for (int j = 0; j < trigObj_120_size; j++){
	trigObj_120_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).id();
	trigObj_120_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).pt();
	trigObj_120_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).eta();
	trigObj_120_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).phi();
	trigObj_120_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).mass();
      }
      
       */
 
      
      //cout<<"event = "<<evt<<endl;
      //cout<<"run = "<<run<<endl;
      
      float etashift = 0.;
      if(run>210497 && run<211300) etashift = 0.465;
      if(run>211300 && run<211800) etashift = -0.465;

      //here add all the different trigger combination methods 
      double treePrescl[6] = {jetMB_p,jet20_p,jet40_p,jet60_p,jet80_p,jet100_p};
      bool trgDec[5] = {(bool)jet20, (bool)jet40, (bool)jet60, (bool)jet80, (bool)jet100};
      
      if(i&1000==0)cout<<"leading jet step"<<endl;
      //lets get the leading jet information here
      int leadJet = -1;
      float temppt = -9;
      for(int j = 0;j<nrefe3;j++){
	if(raw3[j]<20) continue;
	if(fabs(eta3[j]+etashift)>1) continue;
	float jetpt = pt3[j];
	if(jetpt>temppt){
	  temppt = jetpt;
	  leadJet = j;
	}
      
      }

      if(i%1000==0)cout<<"start of jet loop"<<endl;
      //cout<<"number of jets = "<<nrefe3<<endl;
      for(int j = 0;j<nrefe3;j++){//start of jet loop

      	//raw pt cut - keep that for the analysis level.  
	if(raw3[j]<20) continue;

	//first do the 12-003 merging 
	if(fabs(eta3[j]+etashift)<1){
	  //cout<<"inside 12-003 method"<<endl;
	  if(jet100) hppb5->Fill(pt3[j],jet100_p);
	  if(jet80 && !jet100) hppb4->Fill(pt3[j],1);
	  if(jetMB && !jet80 && !jet100) hppb0->Fill(pt3[j],L1_MB_p*jetMB_p);
	}

	//now do 12-017 
	if(fabs(eta3[j]+etashift)<1){
	  //cout<<"inside 12-017 method"<<endl;
	  if(jet20 && pt3[j]>40 && pt3[j]<60) hppb_trk40_60->Fill(pt3[j],jet20_p);
	  if(jet40 && pt3[j]>60 && pt3[j]<75) hppb_trk60_75->Fill(pt3[j],jet40_p);
	  if(jet60 && pt3[j]>75 && pt3[j]<95) hppb_trk75_95->Fill(pt3[j],jet60_p);
	  if(jet80 && pt3[j]>95 && pt3[j]<120) hppb_trk95_120->Fill(pt3[j],jet80_p);
	  if(jet100 && pt3[j]>120) hppb_trk120->Fill(pt3[j],jet100_p);
	}

	//now do it for kurt's method
	if(fabs(eta3[j]+etashift)<1){
	  for(int ii=0; ii<5; ii++){
	    if(trgDec[ii]){
	      triggerPt = triggerMatch(trigObjSize[ii], trigObjPhi[ii], trigObjEta[ii], trigObjPt[ii], phi3[j], eta3[j], pt3[j]);
	      //if you find a trigger match that has the right pt window for that particular trigger, break out of the for loop!
	      if((triggerPt>(ii+1)*20 && triggerPt<(ii+2)*20) || (triggerPt>100 && ii==4)){
		break;
	      }
	    }
	  }
	  if(jet100 && triggerPt>=100) hppb_kurt100->Fill(pt3[j],jet100_p);
	  if(jet80 && triggerPt>=80 && triggerPt<100) hppb_kurt80_100->Fill(pt3[j],jet80_p);
	  if(jet60 && triggerPt>=60 && triggerPt<80) hppb_kurt60_80->Fill(pt3[j],jet60_p);
	  if(jet40 && triggerPt>=40 && triggerPt<60) hppb_kurt40_60->Fill(pt3[j],jet40_p);
	  if(jet20 && triggerPt>=20 && triggerPt<40) hppb_kurt20_40->Fill(pt3[j],jet20_p);
	}
	


	
      }//end jet loop 

      
      if(i%1000==0)cout<<"finished jet loop"<<endl;
      
      
      //add the histograms from the 12-003 method. 
      hpPb_Jet100->Add(hppb5);
      hpPb_Jet80->Add(hppb4);
      hpPb_JetMB->Add(hppb0);
      hppbComb->Add(hppb0);
      hppbComb->Add(hppb4);
      hppbComb->Add(hppb5);
      hpPb_Comb->Add(hppbComb);
      if(i%1000==0)hppbComb->Print("base");

      //add the histograms from the 12-017 method
      hpPb_Trk40_60->Add(hppb_trk40_60);
      hpPb_Trk60_75->Add(hppb_trk60_75);
      hpPb_Trk75_95->Add(hppb_trk75_95);
      hpPb_Trk95_120->Add(hppb_trk95_120);
      hpPb_Trk120->Add(hppb_trk120);
      hppb_trkComb->Add(hppb_trk40_60);
      hppb_trkComb->Add(hppb_trk60_75);
      hppb_trkComb->Add(hppb_trk75_95);
      hppb_trkComb->Add(hppb_trk95_120);
      hppb_trkComb->Add(hppb_trk120);
      hpPb_TrkComb->Add(hppb_trkComb);
      if(i%1000==0)hppb_trkComb->Print("base");

      //add the histogram from kurt's method 
      hpPb_Kurt100->Add(hppb_kurt100);
      hpPb_Kurt80_100->Add(hppb_kurt80_100);
      hpPb_Kurt60_80->Add(hppb_kurt60_80);
      hpPb_Kurt40_60->Add(hppb_kurt40_60);
      hpPb_Kurt20_40->Add(hppb_kurt20_40);
      hppb_kurtComb->Add(hppb_kurt100);
      hppb_kurtComb->Add(hppb_kurt80_100);
      hppb_kurtComb->Add(hppb_kurt60_80);
      hppb_kurtComb->Add(hppb_kurt40_60);
      hppb_kurtComb->Add(hppb_kurt20_40);
      hpPb_KurtComb->Add(hppb_kurtComb);
      if(i%1000==0)hppb_kurtComb->Print("base");

      
      //jetR3Tree->Fill();
      //evtTree->Fill();
      //hltTree->Fill();
      //trkTree->Fill();
    
    }//new event loop 
    
    fin->Close();
    
  }//end of file loop
  /*

  //now we have loaded all the trees into the important tchains. 

  ch[2]->AddFriend(ch[1]);
  ch[2]->AddFriend(ch[0]);
  ch[2]->AddFriend(ch[3]);
  ch[2]->AddFriend(ch[4]);

  cout<<"start setting branch address"<<endl;

  //set the branch addresses:  - one of the most boring parts of the code: 
  ch[2]->SetBranchAddress("evt",&evt);
  ch[2]->SetBranchAddress("run",&run);
  ch[2]->SetBranchAddress("lumi",&lumi);
  ch[2]->SetBranchAddress("hiBin",&hiBin);
  ch[2]->SetBranchAddress("vz",&vz);
  ch[2]->SetBranchAddress("vx",&vx);
  ch[2]->SetBranchAddress("vy",&vy);
  ch[2]->SetBranchAddress("hiNtracks",&hiNtracks);
  ch[2]->SetBranchAddress("hiHFminus",&hiHFminus);
  ch[2]->SetBranchAddress("hiHFplus",&hiHFplus);
  ch[2]->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
  ch[2]->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
  ch[2]->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
  ch[2]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
  ch[2]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
  ch[2]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);

  ch[2]->SetBranchAddress("nref",&nrefe3);
  ch[2]->SetBranchAddress("jtpt",&pt3);
  ch[2]->SetBranchAddress("jteta",&eta3);
  ch[2]->SetBranchAddress("jtphi",&phi3);
  ch[2]->SetBranchAddress("rawpt",&raw3);
  ch[2]->SetBranchAddress("chargedMax",&chMax3);
  ch[2]->SetBranchAddress("chargedSum",&chSum3);
  ch[2]->SetBranchAddress("trackMax",&trkMax3);
  ch[2]->SetBranchAddress("trackSum",&trkSum3);
  ch[2]->SetBranchAddress("photonMax",&phMax3);
  ch[2]->SetBranchAddress("photonSum",&phSum3);
  ch[2]->SetBranchAddress("neutralMax",&neMax3);
  ch[2]->SetBranchAddress("neutralSum",&neSum3);

  ch[2]->SetBranchAddress("nTrk",&nTrack);
  ch[2]->SetBranchAddress("trkPt",&trkPt);
  ch[2]->SetBranchAddress("trkEta",&trkEta);
  ch[2]->SetBranchAddress("trkPhi",&trkPhi);
  ch[2]->SetBranchAddress("highPurity",&highPurity);
  ch[2]->SetBranchAddress("trkDz1",&trkDz1);
  ch[2]->SetBranchAddress("trkDzError1",&trkDzError1);
  ch[2]->SetBranchAddress("trkDxy1",&trkDxy1);
  ch[2]->SetBranchAddress("trkDxyError1",&trkDxyError1);

  ch[2]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB);
  ch[2]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p);
  ch[2]->SetBranchAddress("HLT_PAJet20_NoJetID_v1",&jet20);
  ch[2]->SetBranchAddress("HLT_PAJet20_NoJetID_v1_Prescl",&jet20_p);
  ch[2]->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jet40);
  ch[2]->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jet40_p);
  ch[2]->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jet60);
  ch[2]->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jet60_p);
  ch[2]->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jet80);
  ch[2]->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jet80_p);
  ch[2]->SetBranchAddress("HLT_PAJet100_NoJetID_v1",&jet100);
  ch[2]->SetBranchAddress("HLT_PAJet100_NoJetID_v1_Prescl",&jet100_p);
  ch[2]->SetBranchAddress("HLT_PAJet120_NoJetID_v1",&jet120);
  ch[2]->SetBranchAddress("HLT_PAJet120_NoJetID_v1_Prescl",&jet120_p);
  ch[2]->SetBranchAddress("HLT_PAJet20_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[0]);
  ch[2]->SetBranchAddress("HLT_PAJet40_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[1]);
  ch[2]->SetBranchAddress("HLT_PAJet60_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[2]);
  ch[2]->SetBranchAddress("HLT_PAJet80_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[3]);
  ch[2]->SetBranchAddress("HLT_PAJet100_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[4]);
  ch[2]->SetBranchAddress("HLT_PAJet120_NoJetID_v1_trigObject",&HLT_PAJet_NoJetID_v1_trigObject[5]);

  //now that we have all the branch addresses set! 
  //lets start the event loop process. 
  cout<<"finished setting branch address"<<endl;

  Long64_t nentries = ch[2]->GetEntries();

  for (int i = 0; i < nentries; i++){ //start of event loop. 
    //cout<<"event = "<<i<<endl;
    
    ch[2]->GetEntry(i);

    //cout<<"A"<<endl;

    trigObj_20_size = HLT_PAJet_NoJetID_v1_trigObject[0]->size();

    //cout<<" = "<<trigObj_20_size<<endl;
    for (int j = 0; j < trigObj_20_size; j++){
      trigObj_20_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).id();
      trigObj_20_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).pt();
      trigObj_20_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).eta();
      trigObj_20_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).phi();
      trigObj_20_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).mass();
    }

    //cout<<"B"<<endl;

    trigObj_40_size = HLT_PAJet_NoJetID_v1_trigObject[1]->size();
    //cout<<" = "<<trigObj_40_size<<endl;
    for (int j = 0; j < trigObj_40_size; j++){
      trigObj_40_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).id();
      trigObj_40_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).pt();
      trigObj_40_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).eta();
      trigObj_40_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).phi();
      trigObj_40_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).mass();
    }

    //cout<<"C"<<endl;

    trigObj_60_size = HLT_PAJet_NoJetID_v1_trigObject[2]->size();
    //cout<<"trigObj_60_size = "<<trigObj_60_size<<endl;
    for (int j = 0; j < trigObj_60_size; j++){
      trigObj_60_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).id();
      trigObj_60_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).pt();
      trigObj_60_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).eta();
      trigObj_60_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).phi();
      trigObj_60_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).mass();
    }

    //cout<<"D"<<endl;

    trigObj_80_size = HLT_PAJet_NoJetID_v1_trigObject[3]->size();
    //cout<<" = "<<trigObj_80_size<<endl;
    for (int j = 0; j < trigObj_80_size; j++){
      trigObj_80_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).id();
      trigObj_80_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).pt();
      trigObj_80_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).eta();
      trigObj_80_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).phi();
      trigObj_80_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).mass();
    }

    //cout<<"E"<<endl;

    trigObj_100_size = HLT_PAJet_NoJetID_v1_trigObject[4]->size();
    for (int j = 0; j < trigObj_100_size; j++){
      trigObj_100_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).id();
      trigObj_100_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).pt();
      trigObj_100_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).eta();
      trigObj_100_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).phi();
      trigObj_100_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).mass();
    }

    //cout<<"F"<<endl;

    trigObj_120_size = HLT_PAJet_NoJetID_v1_trigObject[5]->size();
    for (int j = 0; j < trigObj_120_size; j++){
      trigObj_120_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).id();
      trigObj_120_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).pt();
      trigObj_120_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).eta();
      trigObj_120_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).phi();
      trigObj_120_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).mass();
    }

    //cout<<"E"<<endl;

    if(!pHBHENoiseFilter || !pprimaryvertexFilter || !pPAcollisionEventSelectionPA) continue;
    if(!pVertexFilterCutGplus) continue;

    //cout<<"event = "<<evt<<endl;
    //cout<<"run = "<<run<<endl;

    for(int j = 0;j<nrefe3;j++){//start of jet loop

      if(run>211256){
	eta3_CM[j] = eta3[j] - 0.465;
	
      }else {
	eta3_CM[i] = eta3[j] + 0.465;
      }

      //apply the JEC here? 

      //raw pt cut - keep that for the analysis level.  
      //if(raw3<20) continue;
      

    }
    // fill the ntuples here: 
    jetR3Tree->Fill();
    evtTree->Fill();
    hltTree->Fill();
    trkTree->Fill();

  }//end of event loop

  */

  //cout<<"finished the event loop"<<endl;

  //output file definition 
  
  f.cd();

  //jetR3Tree->Write();
  //evtTree->Write();
  //hltTree->Write();
  //trkTree->Write();
  
  f.Write();
  
  //hpPb_Comb->Write();
  //hpPb_Jet100->Write();
  //hpPb_Jet80->Write();
  //hpPb_JetMB->Write();

  f.Close();

  timer.Stop();
  float rtime  = timer.RealTime();
  float ctime  = timer.CpuTime();
  
  std::cout<<"\t"<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;
  
}
