
/*
  Raghav Kunnawalkam Elayavalli
  created: 17th Feb 2014

  sample macro to read in forest files from a file list which is given to us through a condor script. the output file will 
  be sent to hadoop. 
*/



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

using namespace std;

void merge_kurt_files(const int startfile=0, const int endfile=1){

  const int N = 4;
  

  //Create chain
  TChain* ch[N];

  string dir[N] = {
    "hltanalysis",
    "skimanalysis",
    //"hcalNoise",
    "akPu3PFJetAnalyzer",
    //"akPu5PFJetAnalyzer",
    //"multiPhotonAnalyzer",
    //"ppTrack",
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
    //"trackTree",
    //"pfTree",
    //"metTree",
    //"HLTMuTree",
    "HiTree"
  };
    

  
  std::string infile;
  
  //change the following file list to the one required. 
  infile = "kurt_small_filelist.txt"; 

  std::ifstream instr(infile.c_str(), std::ifstream::in);

  std::string filename;
  int nFiles=0;

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
      cout<<" i = "<<i<<endl;
      ch[i] = new TChain(string(dir[i]+"/"+trees[i]).data());
      ch[i]->Add(filename.c_str());
      cout << "Tree loaded  " << string(dir[i]+"/"+trees[i]).data() << endl;
      cout << "Entries : " << ch[i]->GetEntries() << endl;
      
    }
    
  }


  //friend the trees now you can use them as usual, ch[2]->JetTree
  ch[2]->AddFriend(ch[0]);
  ch[2]->AddFriend(ch[1]);
  ch[2]->AddFriend(ch[3]);

  ch[2]->Project("h1","jtpt","abs(vz)<15&&pPAcollisionEventSelectionPA&&pHBHENoiseFilter&&HLT_PAJet100_NoJetID_v1");

  TFile f(Form("kurt_sample_output_%d.root",endfile),"RECREATE");
  h1->Write();
  f.Close();
  
  //}

  // Start working on the chains as trees. 

  /*
  // Output chains to files
  TFile *file = new TFile(outfile, "RECREATE");
  file->cd();
  for (int i = 0; i < N; i++) {
    cout << string(dir[i]+"/"+trees[i]).data() << endl;
    if (i == 0) {
       file->mkdir(dir[i].data())->cd();
    }
    else {
       if (TString(dir[i].data()) != TString(dir[i-1].data()))
         file->mkdir(dir[i].data())->cd();
       else
         file->cd(dir[i].data());
    }
    ch[i]->Merge(file,0,"keep");
  }
  //file->Write();

  file->Close();

  */

}
