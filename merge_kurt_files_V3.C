
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

void merge_kurt_files_V3(const int startfile=0, const int endfile=1){
  timer.Start();
  cout<<"Macro start"<<endl;

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
  
  //double         triggerPt;

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


  //event tree
  int evt;
  int run;
  int lumi;
  int hiBin;
  float vx;
  float vy;
  float vz;
  float hiNtracks;
  float hiHFminus;
  float hiHFplus;
  float hiHFplusEta4;
  float hiHFminusEta4;
  int pPAcollisionEventSelectionPA;
  int pHBHENoiseFilter;
  int pprimaryvertexFilter;
  int pVertexFilterCutGplus;

  TTree *evtTree = new TTree("evt","evt");
  evtTree->Branch("evt",&evt,"evt/I");
  evtTree->Branch("run",&run,"run/I");
  evtTree->Branch("lumi",&lumi,"lumi/I");
  evtTree->Branch("hiBin",&hiBin,"hiBin/I");
  evtTree->Branch("vx",&vx,"vx/F");
  evtTree->Branch("vy",&vy,"vy/F");
  evtTree->Branch("vz",&vz,"vz/F");
  evtTree->Branch("hiNtracks",&hiNtracks,"hiNtracks/F");
  evtTree->Branch("hiHFminus",&hiHFplus,"hiHFplus/F");
  evtTree->Branch("hiHFplus",&hiHFminus,"hiHFminus/F");
  evtTree->Branch("hiHFplusEta4",&hiHFplusEta4,"hiHFplusEta4/F");
  evtTree->Branch("hiHFminusEta4",&hiHFminusEta4,"hiHFminusEta4/F");
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
  float trkPt[1000];
  float trkEta[1000];
  float trkPhi[1000];
  float highPurity[1000];
  float trkDz1[1000];
  float trkDzError1[1000];
  float trkDxy1[1000];
  float trkDxyError1[1000];
  float trkPtError1[1000];
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

  //trigger tree
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
  int trigObj_20_size;
  int trigObj_40_size;
  int trigObj_60_size;
  int trigObj_80_size;
  int trigObj_100_size;
  int trigObj_120_size;
  float trigObj_20_id[1000];
  float trigObj_20_pt[1000];
  float trigObj_20_eta[1000];
  float trigObj_20_phi[1000];
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
  hltTree->Branch("trigObj_40_id",&trigObj_40_id,"trigObj_40_id[trgObj_40_size]/F");
  hltTree->Branch("trigObj_40_pt",&trigObj_40_pt,"trigObj_40_pt[trgObj_40_size]/F");
  hltTree->Branch("trigObj_40_eta",&trigObj_40_eta,"trigObj_40_eta[trgObj_40_size]/F");
  hltTree->Branch("trigObj_40_phi",&trigObj_40_phi,"trigObj_40_phi[trgObj_40_size]/F");
  hltTree->Branch("trigObj_40_mass",&trigObj_40_mass,"trigObj_40_mass[trgObj_40_size]/F");
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

  trigO *HLT_PAJet_NoJetID_v1_trigObject[6];
  for(int i=0; i<6; i++){
    HLT_PAJet_NoJetID_v1_trigObject[i] = new trigO;
  } 
    
  cout<<"finished setting up the branches for the new trees"<<endl;

  std::string infile;
  
  //change the following file list to the one required. 
  //infile = "kurt_small_filelist.txt"; 
  infile = "PAHighPtPurdueForest_4th.txt";

  std::ifstream instr(infile.c_str(), std::ifstream::in);

  std::string filename;
  int nFiles=endfile - startfile;

  cout<<"Running on "<<nFiles<<" forest files"<<endl;

  for(int i = 0;i<N;i++){
    ch[i] = new TChain(string(dir[i]+"/"+trees[i]).data());
  }
  
  //just to read the files till the start number. 
  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;  
  for(int ifile=0;ifile<startfile;ifile++){
    instr >> filename;
  }
  
  //now we are taking only the files from the given start number to the end number. 
  for(int ifile=startfile;ifile<endfile;ifile++){
    instr >> filename;
    
    cout<<"File: "<<filename<<std::endl;
    
    for (int i = 0; i < N; i++) {
      //cout<<" i = "<<i<<endl;
      
      ch[i]->Add(filename.c_str());
      cout << "Tree loaded  " << string(dir[i]+"/"+trees[i]).data() << endl;
      cout << "Entries : " << ch[i]->GetEntries() << endl;
      
    }
    
  }

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
    cout<<"event = "<<i<<endl;
    
    ch[2]->GetEntry(i);

    cout<<"A"<<endl;

    trigObj_20_size = HLT_PAJet_NoJetID_v1_trigObject[0]->size();

    cout<<" = "<<trigObj_20_size<<endl;
    for (int j = 0; j < trigObj_20_size; j++){
      trigObj_20_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).id();
      trigObj_20_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).pt();
      trigObj_20_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).eta();
      trigObj_20_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).phi();
      trigObj_20_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[0]->at(j).mass();
    }

    cout<<"B"<<endl;

    trigObj_40_size = HLT_PAJet_NoJetID_v1_trigObject[1]->size();
    cout<<" = "<<trigObj_40_size<<endl;
    for (int j = 0; j < trigObj_40_size; j++){
      trigObj_40_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).id();
      trigObj_40_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).pt();
      trigObj_40_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).eta();
      trigObj_40_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).phi();
      trigObj_40_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[1]->at(j).mass();
    }

    cout<<"C"<<endl;

    trigObj_60_size = HLT_PAJet_NoJetID_v1_trigObject[2]->size();
    cout<<"trigObj_60_size = "<<trigObj_60_size<<endl;
    for (int j = 0; j < trigObj_60_size; j++){
      trigObj_60_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).id();
      trigObj_60_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).pt();
      trigObj_60_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).eta();
      trigObj_60_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).phi();
      trigObj_60_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[2]->at(j).mass();
    }

    cout<<"D"<<endl;

    trigObj_80_size = HLT_PAJet_NoJetID_v1_trigObject[3]->size();
    cout<<" = "<<trigObj_80_size<<endl;
    for (int j = 0; j < trigObj_80_size; j++){
      trigObj_80_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).id();
      trigObj_80_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).pt();
      trigObj_80_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).eta();
      trigObj_80_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).phi();
      trigObj_80_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[3]->at(j).mass();
    }

    cout<<"E"<<endl;

    trigObj_100_size = HLT_PAJet_NoJetID_v1_trigObject[4]->size();
    for (int j = 0; j < trigObj_100_size; j++){
      trigObj_100_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).id();
      trigObj_100_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).pt();
      trigObj_100_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).eta();
      trigObj_100_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).phi();
      trigObj_100_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[4]->at(j).mass();
    }

    cout<<"F"<<endl;

    trigObj_120_size = HLT_PAJet_NoJetID_v1_trigObject[5]->size();
    for (int j = 0; j < trigObj_120_size; j++){
      trigObj_120_id[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).id();
      trigObj_120_pt[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).pt();
      trigObj_120_eta[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).eta();
      trigObj_120_phi[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).phi();
      trigObj_120_mass[j] =  HLT_PAJet_NoJetID_v1_trigObject[5]->at(j).mass();
    }

    cout<<"E"<<endl;

    bool evSel = false;
    evSel = fabs(vz)<15 && pHBHENoiseFilter && pPAcollisionEventSelectionPA && pprimaryvertexFilter && pVertexFilterCutGplus;
    if(!evSel) continue;


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

  cout<<"finished the event loop"<<endl;

  //output file definition 
  TFile f(Form("ntuple_HighPt_prod_JSON_nofilter_pPb_v8JEC_%d",endfile),"RECREATE");
  jetR3Tree->Write();
  evtTree->Write();
  hltTree->Write();
  trkTree->Write();
  f.Close();

  timer.Stop();
  float rtime  = timer.RealTime();
  float ctime  = timer.CpuTime();
  
  std::cout<<"\t"<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;
  
}
