// Raghav Kunnawalkam Elayavalli
// created 4/22/2014

// macro to read in the second part of the pPb data files.to add to the RpA analysis. do it for all different trigger 
// combination methods. 


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


TH1F *functionHist(TF1 *f, TH1F* h,char *fHistname)
{
	TH1F *hF = (TH1F*)h->Clone(fHistname);
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double var = f->Integral(h->GetBinLowEdge(i),h->GetBinLowEdge(i+1))/h->GetBinWidth(i);
		hF->SetBinContent(i,var);
		hF->SetBinError(i,0);
	}
	return hF;
}


// divide by bin width
void divideBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val/=h->GetBinWidth(i);
		valErr/=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}



void drawText(const char *text, float xp, float yp, int size){
	TLatex *tex = new TLatex(xp,yp,text);
	tex->SetTextFont(63);
	tex->SetTextSize(size);
	tex->SetTextColor(kBlack);
	tex->SetLineWidth(1);
	//tex->SetTextFont(42);
	tex->SetNDC();
	tex->Draw();
}


void putCMSPrel(double x, double y, double size){
	TLatex *tex=0;
	tex = new TLatex(x,y,"CMS Preliminary");
	tex->SetTextSize(size);
	tex->SetLineWidth(2);
	tex->SetNDC();
	tex->Draw();
}


TLegend *myLegend(double x1,double y1,double x2, double y2)
{
	TLegend *leg = new TLegend(x1,y1,x2,y2);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	return leg; 
	
}

// Remove bins with error > central value
void cleanup(TH1F *h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val1 = h->GetBinContent(i);
		double valErr1 = h->GetBinError(i);
		if (valErr1>=val1) {
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}
	}   
	
}


void formatCanvas(TCanvas *c){
  c->Divide(1,2,0.01,0.01);
  c->cd(1);
  //c->GetPad(1)->SetLogy();
  c->GetPad(1)->SetPad(0.,0.425,1.,1.);
  c->GetPad(2)->SetPad(0.,0.0,1.,0.425);
  c->GetPad(2)->SetBottomMargin(0.3);
  c->GetPad(2)->SetGridy(1);
}


//static const int nbins_yaxian = 29;
//static const double boundaries_yaxian[nbins_yaxian+1] = {3,4,5,7,9,12,15,18,22,27,33,39,47,55,64,74,84,97,114,133,153,174,196,220,245,272,300,429,692,1000};

static const int nbins_yaxian = 29;
static const double boundaries_yaxian[nbins_yaxian+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

//static const int nbins_rec = 50;
//statis const int boundaries_rec[]

// rebin the spectra
TH1F *rebin(TH1F *h, char *histName)
{
  TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),50,0,200);
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      double val=h->GetBinContent(i);
      double valErr=h->GetBinError(i);
      int binNum = hRebin->FindBin(h->GetBinCenter(i));
      double val1 = hRebin->GetBinContent(binNum);
      double valErr1 = hRebin->GetBinError(binNum);
      hRebin->SetBinContent(binNum,val+val1);
      hRebin->SetBinError(binNum,sqrt(valErr1*valErr1+valErr*valErr));
    }
  cleanup(hRebin);
  hRebin->SetName(histName);
  return hRebin;
}


// rebin the spectra
TH1F *rebin_yaxian(TH1F *h, char *histName)
{
  TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),nbins_yaxian,boundaries_yaxian);
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      double val=h->GetBinContent(i);
      double valErr=h->GetBinError(i);
      int binNum = hRebin->FindBin(h->GetBinCenter(i));
      double val1 = hRebin->GetBinContent(binNum);
      double valErr1 = hRebin->GetBinError(binNum);
      hRebin->SetBinContent(binNum,val+val1);
      hRebin->SetBinError(binNum,sqrt(valErr1*valErr1+valErr*valErr));
    }
  cleanup(hRebin);
  hRebin->SetName(histName);
  return hRebin;
}


using namespace std;

TStopwatch timer;


void readpPbdata(){

  TH1::SetDefaultSumw2();

  TFile *fin = TFile::Open("/net/hidsk0001/d00/scratch/maoyx/pPb/12-017/Unfold/Qiao/PbPdata_ppReco_akPu3PF_AlljetTrigKurtTrCombFile0_6815_JetPt0noIPupperCut.root");

  TTree* jetTree = (TTree*)fin->Get("nt");

  TFile fout("Pbp_spectrahistos_v1.root","RECREATE");

  // we have to do all this business here inside the event loop

  TH1F* hPbp_Jet100 = new TH1F("hPbp_Jet100","",1000,0,1000);
  TH1F* hPbp_Jet80 = new TH1F("hPbp_Jet80","",1000,0,1000);
  TH1F* hPbp_Jet60 = new TH1F("hPbp_Jet60","",1000,0,1000);
  TH1F* hPbp_Jet40 = new TH1F("hPbp_Jet40","",1000,0,1000);
  TH1F* hPbp_Jet20 = new TH1F("hPbp_Jet20","",1000,0,1000);
  TH1F* hPbp_Jet20_v2 = new TH1F("hPbp_Jet20_v2","",1000,0,1000);
  TH1F* hPbp_Jet60_v2 = new TH1F("hPbp_Jet60_v2","",1000,0,1000);
  //TH1F* hpPb_JetMB = new TH1F("hpPb_JetMB","",1000,0,1000);
  //TH1F* hpPb_Comb_MB = new TH1F("hpPb_Comb_MB","",1000,0,1000);
  TH1F* hPbp_Comb_20 = new TH1F("hPbp_Comb_20","",1000,0,1000);
  TH1F* hPbp_Comb_40 = new TH1F("hPbp_Comb_40","",1000,0,1000);
  TH1F* hPbp_Comb_60 = new TH1F("hPbp_Comb_60","",1000,0,1000);
  TH1F* hPbp_Comb_20_60_100 = new TH1F("hPbp_Comb_20_60_100","",1000,0,1000);

  TH1F* hPbp_Trk40_60 = new TH1F("hPbp_Trk40_60","",1000,0,1000);
  TH1F* hPbp_Trk60_75 = new TH1F("hPbp_Trk60_75","",1000,0,1000);
  TH1F* hPbp_Trk75_95 = new TH1F("hPbp_Trk75_95","",1000,0,1000);
  TH1F* hPbp_Trk95_120 = new TH1F("hPbp_Trk95_120","",1000,0,1000);
  TH1F* hPbp_Trk120 = new TH1F("hPbp_Trk120","",1000,0,1000);
  TH1F* hPbp_TrkComb = new TH1F("hPbp_TrkComb","",1000,0,1000);
  
  /*
  TH1F* hpPb_Trk40_60_v2 = new TH1F("hpPb_Trk40_60_v2","",1000,0,1000);
  TH1F* hpPb_Trk60_75_v2 = new TH1F("hpPb_Trk60_75_v2","",1000,0,1000);
  TH1F* hpPb_Trk75_95_v2 = new TH1F("hpPb_Trk75_95_v2","",1000,0,1000);
  TH1F* hpPb_Trk95_120_v2 = new TH1F("hpPb_Trk95_120_v2","",1000,0,1000);
  TH1F* hpPb_Trk120_v2 = new TH1F("hpPb_Trk120_v2","",1000,0,1000);
  TH1F* hpPb_TrkComb_v2 = new TH1F("hpPb_TrkComb_v2","",1000,0,1000);
  */
  /*
  TH1F* hpPb_Kurt100 = new TH1F("hpPb_Kurt100","",1000,0,1000);
  TH1F* hpPb_Kurt80_100 = new TH1F("hpPb_Kurt80_100","",1000,0,1000);
  TH1F* hpPb_Kurt60_80 = new TH1F("hpPb_Kurt60_80","",1000,0,1000);
  TH1F* hpPb_Kurt40_60 = new TH1F("hpPb_Kurt40_60","",1000,0,1000);
  TH1F* hpPb_Kurt20_40 = new TH1F("hpPb_Kurt20_40","",1000,0,1000);
  TH1F* hpPb_KurtComb = new TH1F("hpPb_KurtComb","",1000,0,1000);
  */
  
  TH1F* hPbp_Kurt100 = new TH1F("hPbp_Kurt100","",1000,0,1000);
  TH1F* hPbp_Kurt80_100 = new TH1F("hPbp_Kurt80_100","",1000,0,1000);
  TH1F* hPbp_Kurt60_80 = new TH1F("hPbp_Kurt60_80","",1000,0,1000);
  TH1F* hPbp_Kurt40_60 = new TH1F("hPbp_Kurt40_60","",1000,0,1000);
  TH1F* hPbp_Kurt20_40 = new TH1F("hPbp_Kurt20_40","",1000,0,1000);
  TH1F* hPbp_KurtComb = new TH1F("hPbp_KurtComb","",1000,0,1000);

  TH1F* hPbp_KurtComb_test2 = new TH1F("hPbp_KurtComb_test2","",1000,0,1000);


  int nrefe;
  Double_t pt[1000];
  Double_t raw[1000];
  Double_t eta[1000];
  Double_t eta_CM[1000];
  Double_t phi[1000];
  Double_t chMax[1000];
  Double_t trkMax[1000];
  Double_t chSum[1000];
  Double_t phSum[1000];
  Double_t neSum[1000];
  Double_t trkSum[1000];
  Double_t phMax[1000];
  Double_t neMax[1000];

  int evt;
  int run;
  Double_t vz;

  int pPAcollisionEventSelectionPA;
  int pHBHENoiseFilter;
  int pprimaryvertexFilter;
  int pVertexFilterCutGplus;
  
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

  Double_t weight;

  jetTree->SetBranchAddress("evt",&evt);
  jetTree->SetBranchAddress("run",&run);
  jetTree->SetBranchAddress("vz",&vz);

  jetTree->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
  jetTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
  jetTree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
  jetTree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
  
  jetTree->SetBranchAddress("nref",&nrefe);
  jetTree->SetBranchAddress("jtpt",&pt);
  jetTree->SetBranchAddress("jteta",&eta);
  jetTree->SetBranchAddress("jtphi",&phi);
  jetTree->SetBranchAddress("rawpt",&raw);
  jetTree->SetBranchAddress("chargedMax",&chMax);
  jetTree->SetBranchAddress("chargedSum",&chSum);
  //jetTree->SetBranchAddress("trackMax",&trkMax);
  //jetTree->SetBranchAddress("trackSum",&trkSum);
  jetTree->SetBranchAddress("photonMax",&phMax);
  jetTree->SetBranchAddress("photonSum",&phSum);
  jetTree->SetBranchAddress("neutralMax",&neMax);
  jetTree->SetBranchAddress("neutralSum",&neSum);

  jetTree->SetBranchAddress("weight",&weight);

  //jetTree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB);
  //jetTree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p);
  jetTree->SetBranchAddress("HLT_PAJet20_noJetID_v1",&jet20);
  jetTree->SetBranchAddress("HLT_PAJet20_noJetID_v1_Prescl",&jet20_p);
  jetTree->SetBranchAddress("HLT_PAJet40_noJetID_v1",&jet40);
  jetTree->SetBranchAddress("HLT_PAJet40_noJetID_v1_Prescl",&jet40_p);
  jetTree->SetBranchAddress("HLT_PAJet60_noJetID_v1",&jet60);
  jetTree->SetBranchAddress("HLT_PAJet60_noJetID_v1_Prescl",&jet60_p);
  jetTree->SetBranchAddress("HLT_PAJet80_noJetID_v1",&jet80);
  jetTree->SetBranchAddress("HLT_PAJet80_noJetID_v1_Prescl",&jet80_p);
  jetTree->SetBranchAddress("HLT_PAJet100_noJetID_v1",&jet100);
  jetTree->SetBranchAddress("HLT_PAJet100_noJetID_v1_Prescl",&jet100_p);
  //jetTree->SetBranchAddress("HLT_PAJet120_noJetID_v1",&jet120);
  //jetTree->SetBranchAddress("HLT_PAJet120_noJetID_v1_Prescl",&jet120_p);
  

  Long64_t nentries = jetTree->GetEntries();
  cout<<"nentries = "<<nentries<<endl;

  for(int i = 0;i<nentries;i++){
  //for(int i = 0;i<100;i++){

     jetTree->GetEntry(i);

     if(i%1000==0)cout<<"event = "<<i<<"; run = "<<run<<endl;

     float etashift = 0.;
     if(run>210497 && run<211300) etashift = 0.465;
     if(run>211300 && run<211800) etashift = -0.465;
     if(etashift==0)break;
     
     if(!pHBHENoiseFilter || !pprimaryvertexFilter || !pPAcollisionEventSelectionPA) continue;
     //if(!pVertexFilterCutGplus) continue;
     //if(vz>15. || vz<-15.) continue;
     if(fabs(vz)>15) continue;
     
     if(!jetMB && !jet20 && !jet40 && !jet60 && !jet80 && !jet100) continue;

     int leadJet = 0;
     float temppt = -9;
     for(int j = 0;j<nrefe;j++){
       if(raw[j]<30) continue;
       if(fabs(eta[j]+etashift)>1) continue;
       float jetpt = pt[j];
       if(jetpt>temppt){
	 temppt = jetpt;
	 leadJet = j;
       }
     }
     //cout<<"leadjet = "<<leadJet;
     
     for(int j = 0;j<nrefe;j++){
       
       if(fabs(eta[j]+etashift)<1){
	 //cout<<"inside 12-003 method"<<endl;
         if(jet100) hPbp_Jet100->Fill(pt[j],1);
         if(jet80 && !jet100) hPbp_Jet80->Fill(pt[j],1);
         if(jet20 && !jet80 && !jet100) hPbp_Jet20->Fill(pt[j],jet20_p);
         if(jet40 && !jet80 && !jet100) hPbp_Jet40->Fill(pt[j],jet40_p);
         if(jet60 && !jet80 && !jet100) hPbp_Jet60->Fill(pt[j],jet60_p);
         if(jet20 && !jet60 && !jet100) hPbp_Jet20_v2->Fill(pt[j],jet20_p);
         if(jet60 && !jet100) hPbp_Jet60_v2->Fill(pt[j],jet60_p);
       }
       
       //now do 12-017 
       if(fabs(eta[j]+etashift)<1){
	 //cout<<"inside 12-017 method"<<endl;
         if(jet20 && pt[leadJet]>=40 && pt[leadJet]<60) hPbp_Trk40_60->Fill(pt[j],jet20_p);
         if(jet40 && pt[leadJet]>=60 && pt[leadJet]<75) hPbp_Trk60_75->Fill(pt[j],jet40_p);
         if(jet60 && pt[leadJet]>=75 && pt[leadJet]<95) hPbp_Trk75_95->Fill(pt[j],jet60_p);
         if(jet80 && pt[leadJet]>=95 && pt[leadJet]<120) hPbp_Trk95_120->Fill(pt[j],jet80_p);
         if(jet100 && pt[leadJet]>=120) hPbp_Trk120->Fill(pt[j],jet100_p);
       }
       
       //now do it for kurt's method
       if(fabs(eta[j]+etashift)<1){
	 if(weight == jet100_p) hPbp_Kurt100->Fill(pt[j],weight);
         if(weight == jet80_p) hPbp_Kurt80_100->Fill(pt[j],weight);
         if(weight == jet60_p) hPbp_Kurt60_80->Fill(pt[j],weight);
         if(weight == jet40_p) hPbp_Kurt40_60->Fill(pt[j],weight);
         if(weight == jet20_p) hPbp_Kurt20_40->Fill(pt[j],weight);
	 hPbp_KurtComb_test2->Fill(pt[j],weight);
       }     
       
     }// end jet loop
  
  
  }// end event loop


  // hpPb_Comb_MB->Add(hpPb_JetMB);
  //hpPb_Comb_MB->Add(hpPb_Jet80);
  //hpPb_Comb_MB->Add(hpPb_Jet100);
  //hpPb_Comb_MB->Print();

  hPbp_Comb_20->Add(hPbp_Jet20);
  hPbp_Comb_20->Add(hPbp_Jet80);
  hPbp_Comb_20->Add(hPbp_Jet100);
  hPbp_Comb_20->Print();

  hPbp_Comb_40->Add(hPbp_Jet40);
  hPbp_Comb_40->Add(hPbp_Jet80);
  hPbp_Comb_40->Add(hPbp_Jet100);
  hPbp_Comb_40->Print();

  hPbp_Comb_60->Add(hPbp_Jet60);
  hPbp_Comb_60->Add(hPbp_Jet80);
  hPbp_Comb_60->Add(hPbp_Jet100);
  hPbp_Comb_60->Print();

  hPbp_Comb_20_60_100->Add(hPbp_Jet20_v2);
  hPbp_Comb_20_60_100->Add(hPbp_Jet60_v2);
  hPbp_Comb_20_60_100->Add(hPbp_Jet100);
  hPbp_Comb_20_60_100->Print();
  
  //add the histograms from the 12-017 method
  hPbp_TrkComb->Add(hPbp_Trk40_60);
  hPbp_TrkComb->Add(hPbp_Trk60_75);
  hPbp_TrkComb->Add(hPbp_Trk75_95);
  hPbp_TrkComb->Add(hPbp_Trk95_120);
  hPbp_TrkComb->Add(hPbp_Trk120);
  hPbp_TrkComb->Print();
  
  //add the histogram from kurt's method 
  hPbp_KurtComb->Add(hPbp_Kurt100);
  hPbp_KurtComb->Add(hPbp_Kurt80_100);
  hPbp_KurtComb->Add(hPbp_Kurt60_80);
  hPbp_KurtComb->Add(hPbp_Kurt40_60);
  hPbp_KurtComb->Add(hPbp_Kurt20_40);
  hPbp_KurtComb->Print();

  hPbp_KurtComb_test2->Print();

  fout.Write();
  fout.Close();
  
  timer.Stop();
  float rtime  = timer.RealTime();
  float ctime  = timer.CpuTime();
  
  std::cout<<"\t"<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;

}
