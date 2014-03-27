// Raghav Kunnawalkam Elayavalli
// created: March 2nd 2014

// macro to perform simple cross check on different trigger methods for the RpA analysis 


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
  c->GetPad(1)->SetLogy();
  c->GetPad(1)->SetPad(0.,0.425,1.,1.);
  c->GetPad(2)->SetPad(0.,0.0,1.,0.45);
  c->GetPad(2)->SetBottomMargin(0.3);
  c->GetPad(2)->SetGridy(1);
}


//static const int nbins_yaxian = 29;
//static const double boundaries_yaxian[nbins_yaxian+1] = {3,4,5,7,9,12,15,18,22,27,33,39,47,55,64,74,84,97,114,133,153,174,196,220,245,272,300,429,692,1000};

static const int nbins_yaxian = 29;
static const double boundaries_yaxian[nbins_yaxian+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

//static const int nbins_rec = 50;
//statis const int boundaries_rec[]

static const int nbins_trigger = 26;
static const double boundaries_trigger[nbins_trigger+1] = {0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114, 120, 126, 132, 138,144,150,156};

static const int nbins_trigger20 = 15;
static const double boundaries_trigger20[nbins_trigger20+1] = {0, 6, 12, 18, 24, 30, 36, 42, 54, 66, 72, 84, 102, 120, 138,156};

//static const int nbins_trigger20 = 8;
//static const double boundaries_trigger20[nbins_trigger20+1] = {0, 6, 12, 18, 24, 30, 36, 42,48};

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

//TStopwatch timer;

void RpA_crosscheck(){
  
  TH1::SetDefaultSumw2();

  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);

  //timer.start();

  TFile *fin = TFile::Open("pAforest_merge_output_akPu3PF_v2.root");
  //TFile *fin_2 = TFile::Open("pAforest_trig_nocut.root");
  //TFile *fin_rerun = TFile::Open("trig_merge_crosscheck_purdueforests_2368.root");
  TFile * fYaxian = TFile::Open("PPbJetTrigPYTHIAak3PFJetSpectraRpAHFsumEta4Bin1.root");
  //TFile *fMB = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/2013/data/ntuple_2013_JEC_applied_ppbJetMB_v2.root");

  //TTree *jet = (TTree*)fMB->Get("jetR3");
  //TTree *evt = (TTree*)fMB->Get("evt");

  //jet->AddFriend(evt);

  //TH1F* hMB = new TH1F("hMB","",1000,0,1000);
  //TCut ppb0 = "abs(eta_CM)<1&&jetMB&&!jet80&&!jet100&&raw>20&&run>210658";
  //jet->Project("hMB","pt","jetMB_p"*ppb0);

  TH1F* hRaghav_test_20 = (TH1F*)fin->Get("hpPb_Comb_20");
  TH1F* hKurt_test = (TH1F*)fin->Get("hpPb_KurtComb");
  TH1F* hYaxian_test = (TH1F*)fin->Get("hpPb_TrkComb");
  TH1F* hpPb_Comb_MB = (TH1F*)fin->Get("hpPb_Comb_MB");

  //TH1F* hRaghav_rerun = (TH1F*)fin_rerun->Get("hpPb_Comb_MB");
  //TH1F* hKurt_rerun = (TH1F*)fin_rerun->Get("hpPb_KurtComb");
  //TH1F* hYaxian_rerun = (TH1F*)fin_rerun->Get("hpPb_TrkComb");

  //TH1F* hYaxian_spectra = (TH1F*)fYaxian->Get("");

  //hRaghav_test->Add(hRaghav_rerun);
  //hYaxian_test->Add(hYaxian_rerun);
  //hKurt_test->Add(hKurt_rerun);

  hRaghav_test_20->Print();
  hKurt_test->Print();
  hYaxian_test->Print();

  TH1F* hpPb_Jet100 = (TH1F*)fin->Get("hpPb_100");
  TH1F* hpPb_Jet80 = (TH1F*)fin->Get("hpPb_80");
  TH1F* hpPb_Jet60 = (TH1F*)fin->Get("hpPb_60");
  TH1F* hpPb_Jet40 = (TH1F*)fin->Get("hpPb_40");
  TH1F* hpPb_Jet20 = (TH1F*)fin->Get("hpPb_20");
  TH1F* hpPb_JetMB = (TH1F*)fin->Get("hpPb_MB");

  TH1F* hpPb_comb = (TH1F*)fin->Get("hpPb_KurtComb_nocut");

  TH1F* hpPb_raghav_MB = (TH1F*)fin->Get("hpPb_JetMB");
  TH1F* hpPb_raghav_20 = (TH1F*)fin->Get("hpPb_Jet20");
  TH1F* hpPb_raghav_40 = (TH1F*)fin->Get("hpPb_Jet40");
  TH1F* hpPb_raghav_60 = (TH1F*)fin->Get("hpPb_Jet60");
  TH1F* hpPb_raghav_80 = (TH1F*)fin->Get("hpPb_Jet80");
  TH1F* hpPb_raghav_100 = (TH1F*)fin->Get("hpPb_Jet100");

  TH1F* hpPb_raghav_comb_60 = (TH1F*)fin->Get("hpPb_Comb_60");
  TH1F* hpPb_raghav_comb_20 = (TH1F*)hRaghav_test_20->Clone("hpPb_raghav_comb_20");
  TH1F* hpPb_raghav_comb_40 = (TH1F*)fin->Get("hpPb_Comb_40");
  TH1F* hpPb_raghav_comb_MB = (TH1F*)fin->Get("hpPb_Comb_MB");

  TH1F* hpPb_kurt_comb = (TH1F*)hKurt_test->Clone("hpPb_kurt_comb");
  TH1F* hpPb_kurt_20_40 = (TH1F*)fin->Get("hpPb_Kurt20_40");
  TH1F* hpPb_kurt_40_60 = (TH1F*)fin->Get("hpPb_Kurt40_60");
  TH1F* hpPb_kurt_60_80 = (TH1F*)fin->Get("hpPb_Kurt60_80");
  TH1F* hpPb_kurt_80_100 = (TH1F*)fin->Get("hpPb_Kurt80_100");
  TH1F* hpPb_kurt_100 = (TH1F*)fin->Get("hpPb_Kurt100");
  /*
  hpPb_raghav_MB = (TH1F*)hpPb_raghav_MB->Rebin(nbins_yaxian,"hpPb_raghav_MB",boundaries_yaxian);
  hpPb_raghav_20 = (TH1F*)hpPb_raghav_20->Rebin(nbins_yaxian,"hpPb_raghav_20",boundaries_yaxian);
  hpPb_raghav_40 = (TH1F*)hpPb_raghav_40->Rebin(nbins_yaxian,"hpPb_raghav_40",boundaries_yaxian);
  hpPb_raghav_60 = (TH1F*)hpPb_raghav_60->Rebin(nbins_yaxian,"hpPb_raghav_60",boundaries_yaxian);
  hpPb_raghav_80 = (TH1F*)hpPb_raghav_80->Rebin(nbins_yaxian,"hpPb_raghav_80",boundaries_yaxian);
  hpPb_raghav_100 = (TH1F*)hpPb_raghav_100->Rebin(nbins_yaxian,"hpPb_raghav_100",boundaries_yaxian);

  hpPb_raghav_comb_60 = (TH1F*)hpPb_raghav_comb_60->Rebin(nbins_yaxian,"hpPb_raghav_comb_60",boundaries_yaxian);
  hpPb_raghav_comb_40 = (TH1F*)hpPb_raghav_comb_40->Rebin(nbins_yaxian,"hpPb_raghav_comb_40",boundaries_yaxian);
  hpPb_raghav_comb_20 = (TH1F*)hpPb_raghav_comb_20->Rebin(nbins_yaxian,"hpPb_raghav_comb_20",boundaries_yaxian);
  hpPb_raghav_comb_MB = (TH1F*)hpPb_raghav_comb_MB->Rebin(nbins_yaxian,"hpPb_raghav_comb_MB",boundaries_yaxian);

  
  hpPb_kurt_comb = (TH1F*)hpPb_kurt_comb->Rebin(nbins_yaxian,"hpPb_kurt_comb",boundaries_yaxian);
  hpPb_kurt_20_40 = (TH1F*)hpPb_kurt_20_40->Rebin(nbins_yaxian,"hpPb_kurt_20_40",boundaries_yaxian);
  hpPb_kurt_40_60 = (TH1F*)hpPb_kurt_40_60->Rebin(nbins_yaxian,"hpPb_kurt_40_60",boundaries_yaxian);
  hpPb_kurt_60_80 = (TH1F*)hpPb_kurt_60_80->Rebin(nbins_yaxian,"hpPb_kurt_comb",boundaries_yaxian);
  hpPb_kurt_80_100 = (TH1F*)hpPb_kurt_80_100->Rebin(nbins_yaxian,"hpPb_kurt_comb",boundaries_yaxian);
  hpPb_kurt_100 = (TH1F*)hpPb_kurt_100->Rebin(nbins_yaxian,"hpPb_kurt_comb",boundaries_yaxian);
  */
  divideBinWidth(hpPb_raghav_MB);
  divideBinWidth(hpPb_raghav_20);
  divideBinWidth(hpPb_raghav_40);
  divideBinWidth(hpPb_raghav_60);
  divideBinWidth(hpPb_raghav_80);
  divideBinWidth(hpPb_raghav_100);

  divideBinWidth(hpPb_raghav_comb_60);
  divideBinWidth(hpPb_raghav_comb_40);
  divideBinWidth(hpPb_raghav_comb_20);
  divideBinWidth(hpPb_raghav_comb_MB);

  divideBinWidth(hpPb_kurt_comb);
  divideBinWidth(hpPb_kurt_20_40);
  divideBinWidth(hpPb_kurt_40_60);
  divideBinWidth(hpPb_kurt_60_80);
  divideBinWidth(hpPb_kurt_80_100);
  divideBinWidth(hpPb_kurt_100);

  TH1F* hRaghav = rebin_yaxian(hRaghav_test_20,"hRaghav");
  TH1F* hYaxian = rebin_yaxian(hYaxian_test,"hYaxian");
  TH1F* hKurt = rebin_yaxian(hKurt_test,"hKurt");
  /*
  hpPb_Jet100 = (TH1F*)Rebin(nbins_yaxian,"hpPb_Jet100",boundaries_yaxian);
  hpPb_Jet80 = (TH1F*)Rebin(nbins_yaxian,"hpPb_Jet80",boundaries_yaxian);
  hpPb_Jet60 = (TH1F*)Rebin(nbins_yaxian,"hpPb_Jet60",boundaries_yaxian);
  hpPb_Jet40 = (TH1F*)Rebin(nbins_yaxian,"hpPb_Jet40",boundaries_yaxian);
  hpPb_Jet20 = (TH1F*)Rebin(nbins_yaxian,"hpPb_Jet20",boundaries_yaxian);
  hpPb_JetMB = (TH1F*)Rebin(nbins_yaxian,"hpPb_JetMB",boundaries_yaxian);
  */

  hpPb_Jet100->Print("base");
  hpPb_Jet80->Print("base");
  hpPb_Jet60->Print("base");
  hpPb_Jet40->Print("base");
  hpPb_Jet20->Print("base");
  hpPb_JetMB->Print("base");
  hpPb_comb->Print("base");

  divideBinWidth(hpPb_Jet100);
  divideBinWidth(hpPb_Jet80);
  divideBinWidth(hpPb_Jet60);
  divideBinWidth(hpPb_Jet40);
  divideBinWidth(hpPb_Jet20);
  divideBinWidth(hpPb_JetMB);
  divideBinWidth(hpPb_comb);

  /*
  TH1F *hComb = new TH1F("hComb","",1000,0,1000);
  hComb->Print("base");
  hComb->Add(hpPb_Jet100);
  hComb->Add(hpPb_Jet80);
  hComb->Add(hpPb_Jet60);
  hComb->Add(hpPb_Jet40);
  hComb->Add(hpPb_Jet20);
  hComb->Print("base");
  */

  TH1F* hTurnon100 = (TH1F*)hpPb_Jet100->Clone("hTurnon100");
  TH1F* hTurnon80 = (TH1F*)hpPb_Jet80->Clone("hTurnon80");
  TH1F* hTurnon60 = (TH1F*)hpPb_Jet60->Clone("hTurnon60");
  TH1F* hTurnon40 = (TH1F*)hpPb_Jet40->Clone("hTurnon40");
  TH1F* hTurnon20 = (TH1F*)hpPb_Jet20->Clone("hTurnon20");

  hTurnon100->Divide(hpPb_comb);
  hTurnon80->Divide(hpPb_comb);
  hTurnon60->Divide(hpPb_comb);
  hTurnon40->Divide(hpPb_comb);
  hTurnon20->Divide(hpPb_comb);
  
  //hKurt->Add(hKurt_rerun);

  TCanvas *c1 = new TCanvas("c1","",800,600);
  //c1->Divide(2,1);
  //c1->cd(1)->SetLogy();
  //c1->cd(1)->SetGridy();
  formatCanvas(c1);

  divideBinWidth(hKurt);
  divideBinWidth(hRaghav);
  divideBinWidth(hYaxian);

  hRaghav->SetMarkerStyle(20);
  hRaghav->SetMarkerColor(kCyan);
  hKurt->SetMarkerStyle(20);
  hKurt->SetMarkerColor(kBlack);
  hYaxian->SetMarkerColor(kRed);
  hYaxian->SetMarkerStyle(20);

  hRaghav->SetTitle("");
  hRaghav->SetXTitle("p_{T} GeV/c");
  hRaghav->SetYTitle("arbirary counts");
 
  hRaghav->Draw();
  hKurt->Draw("same");
  //hYaxian->Draw("same");

  TLegend *title = myLegend(0.3,0.65,0.75,0.9);
  title->AddEntry(hRaghav,"12-003 (HLT_100, HLT_80 and HLT_20)","pl");
  title->AddEntry(hKurt,"14-007 (latest event by event method)","pl");
  //title->AddEntry(hYaxian,"12-017 (modified method using jet by jet weighting - not correct)","pl");
  title->SetTextSize(0.04);
  title->Draw();

  //drawText("pPb 2013 Data",0.3,0.65,20);
  //drawText("Anti-k_{T} Pu PF Jets R = 0.3, |#eta_{CM}<1|, |vz|<15",0.3,0.56,20);

  putCMSPrel(0.1,0.9,0.06);
  drawText("#sqrt{s_{NN}}=5.02 TeV",0.6,0.91,16);
  drawText("akPu3PF, |#eta_{CM}|<1, |vz|<15, 0-100%",0.55,0.55,16);

  c1->cd(2);
  //c1->cd(2)->SetGridy();
  TH1F* hRaghav_ratio = rebin_yaxian(hKurt,"hRaghav_ratio");
  TH1F* hYaxian_ratio = rebin_yaxian(hYaxian,"hYaxian_ratio");

  hRaghav_ratio->Print();
  hYaxian_ratio->Print();
  hRaghav_ratio->Divide(hRaghav);
  hYaxian_ratio->Divide(hKurt);

  hRaghav_ratio->SetMarkerColor(kCyan);
  hRaghav_ratio->SetMarkerStyle(20);
  hYaxian_ratio->SetMarkerColor(kRed);
  hYaxian_ratio->SetMarkerStyle(21);

  hRaghav_ratio->SetXTitle("p_{T} GeV/c");
  hRaghav_ratio->SetYTitle(" ");
  hRaghav_ratio->SetTitle(" ");
  hRaghav_ratio->SetAxisRange(20,300,"X");
  hRaghav_ratio->SetAxisRange(0.8,1.2,"Y");
  hRaghav_ratio->Draw();
  //hYaxian_ratio->Draw("same");

  TLegend *title2 = myLegend(0.54,0.25,0.75,0.45);
  title2->AddEntry(hRaghav_ratio,"14-007/12-003","pl");
  //title2->AddEntry(hYaxian_ratio,"12-017/14-007","pl");
  title2->SetTextSize(0.04);
  title2->Draw();

  c1->SaveAs("RpA_trigger_Crosscheck.pdf","RECREATE");

  //plot the trigger turn on curves: 
  TCanvas *c2 = new TCanvas("c2","",800,600);
  //TH1F* hTurnon20_values = (TH1F*)hTurnon20->Rebin(nbins_yaxian,"hTurnon20_values",boundaries_yaxian);
  //divideBinWidth(hTurnon20_values);
  //for(int i = 0;i<nbins_yaxian;i++){
    //if(boundaries_yaxian[i]<=40){
  //    cout<<boundaries_yaxian[i]<<" :  "<<hTurnon20_values->GetBinContent(i)<<endl;
      //}else
      //cout<<boundaries_yaxian[i]<<" : "<<"1"<<endl;
  //}

  hTurnon100 = (TH1F*)hTurnon100->Rebin(nbins_trigger,"hTurnon100",boundaries_trigger);
  hTurnon80 = (TH1F*)hTurnon80->Rebin(nbins_trigger,"hTurnon80",boundaries_trigger);
  hTurnon60 = (TH1F*)hTurnon60->Rebin(nbins_trigger,"hTurnon60",boundaries_trigger);
  hTurnon40 = (TH1F*)hTurnon40->Rebin(nbins_trigger,"hTurnon40",boundaries_trigger);
  hTurnon20 = (TH1F*)hTurnon20->Rebin(nbins_trigger20,"hTurnon20",boundaries_trigger20);

  divideBinWidth(hTurnon100);
  divideBinWidth(hTurnon80);
  divideBinWidth(hTurnon60);
  divideBinWidth(hTurnon40);
  divideBinWidth(hTurnon20);

  for(int i = 0;i<hTurnon20->GetNbinsX();i++){
    if(hTurnon20->GetBinContent(i)>1.0) hTurnon20->SetBinContent(i,1);
  }
  
  for(int i = 0;i<hTurnon40->GetNbinsX();i++){
    if(hTurnon40->GetBinContent(i)>1.0) hTurnon40->SetBinContent(i,1);
    if(hTurnon60->GetBinContent(i)>1.0) hTurnon60->SetBinContent(i,1);
    if(hTurnon80->GetBinContent(i)>1.0) hTurnon80->SetBinContent(i,1);
    if(hTurnon100->GetBinContent(i)>1.0) hTurnon100->SetBinContent(i,1);
  }

  hTurnon100->SetMarkerStyle(20);
  hTurnon100->SetMarkerColor(kCyan);
  hTurnon80->SetMarkerStyle(20);
  hTurnon80->SetMarkerColor(kRed);
  hTurnon60->SetMarkerStyle(20);
  hTurnon60->SetMarkerColor(kYellow);
  hTurnon40->SetMarkerStyle(20);
  hTurnon40->SetMarkerColor(kGreen);
  hTurnon20->SetMarkerStyle(20);
  hTurnon20->SetMarkerColor(kBlue);

  hTurnon100->SetAxisRange(0,150,"X");
  hTurnon100->SetAxisRange(0,1.3,"Y");
  hTurnon100->SetXTitle("offline Jet p_{T} GeV/c");
  hTurnon100->SetYTitle("Trigger Efficiency");
  hTurnon100->GetXaxis()->SetTitleOffset(1.3);
  hTurnon100->Draw("");
  hTurnon80->Draw("same");
  hTurnon60->Draw("same");
  hTurnon40->Draw("same");
  hTurnon20->Draw("same");				

  TLegend *title3 = myLegend(0.68,0.15,0.92,0.45);
  title3->AddEntry(hTurnon20,"HLT_PAJet20","pl");
  title3->AddEntry(hTurnon40,"HLT_PAJet40","pl");
  title3->AddEntry(hTurnon60,"HLT_PAJet60","pl");
  title3->AddEntry(hTurnon80,"HLT_PAJet80","pl");
  title3->AddEntry(hTurnon100,"HLT_PAJet100","pl");
  title3->SetTextSize(0.03);
  title3->Draw();
  
  putCMSPrel(0.1,0.92,0.06);
  drawText("pPb #int dt = 20.7 nb^{-1}, #sqrt{s_{NN}}=5.02 TeV",0.5,0.93,16);

  c2->SaveAs("RpA_trigger_turnon.pdf","RECREATE");

  // plot the hlt trigger merging for 12-003 and 14-007
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->SetLogy();
  TF1 *fPowerLaw = new TF1("fPowerLaw","[0]/pow(x,[1])");
  hpPb_kurt_comb->Fit("fPowerLaw","","",35,350);
  hpPb_kurt_comb->Fit("fPowerLaw","","",35,350);
  hpPb_kurt_comb->Fit("fPowerLaw","","",35,350);
  hpPb_kurt_comb->Fit("fPowerLaw","","",35,350);
  hpPb_kurt_comb->Fit("fPowerLaw","","",35,350);

  hpPb_kurt_comb->SetMarkerStyle(20);
  hpPb_kurt_comb->SetMarkerColor(kBlack);
  hpPb_kurt_20_40->SetMarkerStyle(20);
  hpPb_kurt_20_40->SetMarkerColor(kBlue);
  hpPb_kurt_40_60->SetMarkerStyle(20);
  hpPb_kurt_40_60->SetMarkerColor(kGreen);
  hpPb_kurt_60_80->SetMarkerStyle(20);
  hpPb_kurt_60_80->SetMarkerColor(kYellow);
  hpPb_kurt_80_100->SetMarkerStyle(20);
  hpPb_kurt_80_100->SetMarkerColor(kRed);
  hpPb_kurt_100->SetMarkerStyle(20);
  hpPb_kurt_100->SetMarkerColor(kCyan);

  hpPb_kurt_comb->SetXTitle("p_{T} GeV/c");
  hpPb_kurt_comb->SetYTitle("Arbitrary counts");
  hpPb_kurt_comb->SetAxisRange(20,350,"X");
  hpPb_kurt_comb->Draw();
  hpPb_kurt_20_40->Draw("same");
  hpPb_kurt_40_60->Draw("same");
  hpPb_kurt_60_80->Draw("same");
  hpPb_kurt_80_100->Draw("same");
  hpPb_kurt_100->Draw("same");
  hpPb_kurt_comb->Draw("same");

  TLegend* title4 = myLegend(0.5,0.6,0.9,0.9);
  title4->AddEntry(hpPb_kurt_comb,"14-007 method combined spectra","pl");
  title4->AddEntry(hpPb_kurt_20_40,"HLT_PAJet20, 20<= Max trigger pt < 40","pl");
  title4->AddEntry(hpPb_kurt_40_60,"HLT_PAJet40, 40<= Max trigger pt < 60","pl");
  title4->AddEntry(hpPb_kurt_60_80,"HLT_PAJet60, 60<= Max trigger pt < 80","pl");
  title4->AddEntry(hpPb_kurt_80_100,"HLT_PAJet80, 80<= Max trigger pt < 100","pl");
  title4->AddEntry(hpPb_kurt_100,"HLT_PAJet100, 100<= Max trigger pt","pl");
  title4->Draw();

  putCMSPrel(0.1,0.9,0.06);
  drawText("#sqrt{s_{NN}}=5.02 TeV",0.5,0.91,16);
  drawText("pPb #int dt = 15.784 nb^{-1}, |#eta_{CM}|<1, 0-100%",0.5,0.5,16);
  //drawText();

  c3->SaveAs("RpA_14007_method_combination.pdf","RECREATE");

  TCanvas *c4 = new TCanvas("c4","",800,600);
  c4->SetLogy();
  //TF1* fPowerLaw2 = new TF1("fPowerLaw2","[0]/pow(x,[1])");
  hpPb_raghav_comb_MB->SetMarkerStyle(20);
  hpPb_raghav_comb_MB->SetMarkerColor(kBlue);
  hpPb_raghav_comb_20->SetMarkerStyle(20);
  hpPb_raghav_comb_20->SetMarkerColor(kBlue);
  hpPb_raghav_comb_40->SetMarkerStyle(20);
  hpPb_raghav_comb_40->SetMarkerColor(kGreen);
  hpPb_raghav_comb_60->SetMarkerStyle(20);
  hpPb_raghav_comb_60->SetMarkerColor(kYellow);
  hpPb_kurt_comb->SetMarkerStyle(20);
  hpPb_kurt_comb->SetMarkerColor(kBlack);

  hpPb_kurt_comb->SetAxisRange(30,150,"X");
  hpPb_kurt_comb->SetYTitle("arbitrary counts");
  hpPb_kurt_comb->SetXTitle("p_{T} GeV/c");

  hpPb_kurt_comb->Draw();
  //hpPb_raghav_comb_MB->Draw("same");
  hpPb_raghav_comb_20->Draw("same");
  hpPb_raghav_comb_40->Draw("same");
  hpPb_raghav_comb_60->Draw("same");
  
  TLegend* title5 = myLegend(0.5,0.6,0.9,0.9);
  title5->AddEntry(hpPb_kurt_comb,"14-007 method","pl");
  //title5->AddEntry(hpPb_raghav_comb_MB,"12-003 HLT100 + HLT80 + HLT MB","pl");
  title5->AddEntry(hpPb_raghav_comb_20,"12-003 HLT100 + HLT80 + HLT 20","pl");
  title5->AddEntry(hpPb_raghav_comb_40,"12-003 HLT100 + HLT80 + HLT 40","pl");
  title5->AddEntry(hpPb_raghav_comb_60,"12-003 HLT100 + HLT80 + HLT 60","pl");
  title5->Draw();

  putCMSPrel(0.1,0.9,0.06);
  drawText("#sqrt{s_{NN}}=5.02 TeV",0.5,0.91,16);
  drawText("pPb #int dt = 15.784 nb^{-1}, |#eta_{CM}|<1, 0-100%",0.5,0.5,16);
  c4->SaveAs("RpA_12003_comparingwith_14007.pdf","RECREATE");

  TCanvas *c5 = new TCanvas("c5","",800,600);
  c5->SetLogy();
  hpPb_raghav_comb_20->SetMarkerStyle(20);
  hpPb_raghav_comb_20->SetMarkerColor(kBlack);
  hpPb_raghav_20->SetMarkerStyle(20);
  hpPb_raghav_20->SetMarkerColor(kBlue);
  hpPb_raghav_80->SetMarkerStyle(20);
  hpPb_raghav_80->SetMarkerColor(kRed);
  hpPb_raghav_100->SetMarkerStyle(20);
  hpPb_raghav_100->SetMarkerColor(kCyan);

  hpPb_raghav_comb_20->SetAxisRange(30,150,"X");
  hpPb_raghav_comb_20->SetYTitle("arbitrary weighted counts");
  hpPb_raghav_comb_20->SetXTitle("p_{T} GeV/c");

  hpPb_raghav_comb_20->Draw();
  hpPb_raghav_20->Draw("same");
  hpPb_raghav_80->Draw("same");
  hpPb_raghav_100->Draw("same");
  hpPb_raghav_comb_20->Draw("same");

  TLegend* title6 = myLegend(0.5,0.6,0.9,0.9);
  title6->AddEntry(hpPb_kurt_comb,"14-007 method","pl");
  //title5->AddEntry(hpPb_raghav_comb_MB,"12-003 HLT100 + HLT80 + HLT MB","pl");
  title6->AddEntry(hpPb_raghav_20,"HLT_PAJet20 && !HLT_PAJet80 && !HLT_PAJet100","pl");
  title6->AddEntry(hpPb_raghav_80,"HLT_PAJet80 && !HLT_PAJet100","pl");
  title6->AddEntry(hpPb_raghav_100,"HLT_PAJet100","pl");
  title6->Draw();

  putCMSPrel(0.1,0.9,0.06);
  drawText("#sqrt{s_{NN}}=5.02 TeV",0.5,0.91,16);
  drawText("pPb #int dt = 15.784 nb^{-1}, |#eta_{CM}|<1, 0-100%",0.5,0.5,16);
  c5->SaveAs("RpA_12003_method_combination.pdf","RECREATE");

  TFile fout("RpA_trigcombination_check.root","RECREATE");
  hRaghav->Write();
  hKurt->Write();
  hYaxian->Write();
  hRaghav_ratio->Write();
  hYaxian_ratio->Write();
  fout.Write();
  fout.Close();

  
 


}
