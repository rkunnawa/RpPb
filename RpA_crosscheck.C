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
  c->GetPad(1)->SetPad(0.,0.45,1.,1.);
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

static const int nbins_trigger = 25;
static const double boundaries_trigger[nbins_trigger+1] = {0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114, 120, 126, 132, 138,144,150};

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
  //gStyle->SetErrorX(0);

  //timer.start();

  TFile *fin = TFile::Open("pPb_spectrahistos_v4.root");
  TFile *fin2 = TFile::Open("Pbp_spectrahistos_v2.root");
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
  TH1F* hRaghav_test_20_60_100 = (TH1F*)fin->Get("hpPb_Comb_40");
  //TH1F* hRaghav_test_20_60_100 = hRaghav_test_20->Clone("hpPb_Comb_20_60_100");
  TH1F* hKurt_test = (TH1F*)fin->Get("hpPb_KurtComb");
  TH1F* hYaxian_test = (TH1F*)fin->Get("hpPb_TrkComb");

  TH1F* hRaghav_test_20_B = (TH1F*)fin2->Get("hpPb_Comb_20");
  TH1F* hRaghav_test_20_60_100_B = (TH1F*)fin2->Get("hpPb_Comb_20_60_100");
  //TH1F* hRaghav_test_20_60_100_B = hRaghav_test_20_60_100->Clone("hPbp_Comb_20_60_100");
  TH1F* hKurt_test_B = (TH1F*)fin2->Get("hpPb_KurtComb");
  TH1F* hYaxian_test_B = (TH1F*)fin2->Get("hpPb_TrkComb");

  TH1F* hpPb_Comb_MB = (TH1F*)fin->Get("hpPb_Comb_MB");

  //cout<<"hRaghav entries = "<< hRaghav_test_20->GetEntries();
  //cout<<"hRaghav 20_60_100 = entriea = "<<hRaghav_test_20_60_100->GetEntries();
  //cout<<"hKurt entries = "<<hKurt_test->GetEntries();
  //cout<<"hYaxian entries = "<<hYaxian_test->GetEntries();

  //cout<<"hRaghav integral = "<< hRaghav_test_20->Integral();
  //cout<<"hKurt integral = "<<hKurt_test->Integral();
  //cout<<"hYaxian integral = "<<hYaxian_test->Integral();

  //TH1F* hRaghav_rerun = (TH1F*)fin_rerun->Get("hpPb_Comb_MB");
  //TH1F* hKurt_rerun = (TH1F*)fin_rerun->Get("hpPb_KurtComb");
  //TH1F* hYaxian_rerun = (TH1F*)fin_rerun->Get("hpPb_TrkComb");

  //TH1F* hYaxian_spectra = (TH1F*)fYaxian->Get("");

  //hRaghav_test->Add(hRaghav_rerun);
  //hYaxian_test->Add(hYaxian_rerun);
  //hKurt_test->Add(hKurt_rerun);

  hRaghav_test_20->Print();
  hRaghav_test_20_60_100->Print();
  hKurt_test->Print();
  hYaxian_test->Print();

  TH1F* hpPb_Jet100 = (TH1F*)fin->Get("hpPb_100");
  TH1F* hpPb_Jet80 = (TH1F*)fin->Get("hpPb_80");
  TH1F* hpPb_Jet60 = (TH1F*)fin->Get("hpPb_60");
  TH1F* hpPb_Jet40 = (TH1F*)fin->Get("hpPb_40");
  TH1F* hpPb_Jet20 = (TH1F*)fin->Get("hpPb_20");
  TH1F* hpPb_JetMB = (TH1F*)fin->Get("hpPb_MB");

  TH1F* hPbp_Jet100 = (TH1F*)fin2->Get("hpPb_100");
  TH1F* hPbp_Jet80 = (TH1F*)fin2->Get("hpPb_80");
  TH1F* hPbp_Jet60 = (TH1F*)fin2->Get("hpPb_60");
  TH1F* hPbp_Jet40 = (TH1F*)fin2->Get("hpPb_40");
  TH1F* hPbp_Jet20 = (TH1F*)fin2->Get("hpPb_20");
  TH1F* hPbp_JetMB = (TH1F*)fin2->Get("hpPb_MB");

  TH1F* hpPb_comb = (TH1F*)fin->Get("hpPb_KurtComb_nocut");
  TH1F* hPbp_comb = (TH1F*)fin2->Get("hpPb_KurtComb_nocut");

  TH1F* hpPb_raghav_MB = (TH1F*)fin->Get("hpPb_JetMB");
  TH1F* hpPb_raghav_20 = (TH1F*)fin->Get("hpPb_Jet20");
  TH1F* hpPb_raghav_40 = (TH1F*)fin->Get("hpPb_Jet40");
  TH1F* hpPb_raghav_60 = (TH1F*)fin->Get("hpPb_Jet60");
  TH1F* hpPb_raghav_80 = (TH1F*)fin->Get("hpPb_Jet80");
  TH1F* hpPb_raghav_100 = (TH1F*)fin->Get("hpPb_Jet100");

  TH1F* hpPb_raghav_comb_60 = (TH1F*)fin->Get("hpPb_Comb_60");
  TH1F* hpPb_raghav_comb_20 = (TH1F*)hRaghav_test_20->Clone("hpPb_raghav_comb_20");
  TH1F* hpPb_raghav_comb_20_60_100 = (TH1F*)hRaghav_test_20_60_100->Clone("hpPb_raghav_comb_20_60_100");
  TH1F* hpPb_raghav_comb_40 = (TH1F*)fin->Get("hpPb_Comb_40");
  TH1F* hpPb_raghav_comb_MB = (TH1F*)fin->Get("hpPb_Comb_MB");

  TH1F* hpPb_kurt_comb = (TH1F*)hKurt_test->Clone("hpPb_kurt_comb");
  TH1F* hpPb_kurt_20_40 = (TH1F*)fin->Get("hpPb_Kurt20_40");
  TH1F* hpPb_kurt_40_60 = (TH1F*)fin->Get("hpPb_Kurt40_60");
  TH1F* hpPb_kurt_60_80 = (TH1F*)fin->Get("hpPb_Kurt60_80");
  TH1F* hpPb_kurt_80_100 = (TH1F*)fin->Get("hpPb_Kurt80_100");
  TH1F* hpPb_kurt_100 = (TH1F*)fin->Get("hpPb_Kurt100");
  /*
  TH1F* hpPb_Jet100 = (TH1F*)fin->Get("hpPb_100");
  TH1F* hpPb_Jet80 = (TH1F*)fin->Get("hpPb_80");
  TH1F* hpPb_Jet60 = (TH1F*)fin->Get("hpPb_60");
  TH1F* hpPb_Jet40 = (TH1F*)fin->Get("hpPb_40");
  TH1F* hpPb_Jet20 = (TH1F*)fin->Get("hpPb_20");
  TH1F* hpPb_JetMB = (TH1F*)fin->Get("hpPb_MB");

  TH1F* hpPb_comb = (TH1F*)fin->Get("hpPb_KurtComb_nocut");
  */
  //TH1F* hpPb_raghav_20 = (TH1F*)fin->Get("hpPb_Jet20");
  //TH1F* hpPb_raghav_40 = (TH1F*)fin->Get("hpPb_Jet40");
  //TH1F* hpPb_raghav_60 = (TH1F*)fin->Get("hpPb_Jet60");
  //TH1F* hpPb_raghav_80 = (TH1F*)fin->Get("hpPb_Jet80");
  //TH1F* hpPb_raghav_100 = (TH1F*)fin->Get("hpPb_Jet100");

  cout<<"got all the histograms from file pPb"<<endl;

  TH1F* hPbp_raghav_comb_60 = (TH1F*)fin2->Get("hpPb_Comb_60");
  TH1F* hPbp_raghav_comb_20 = (TH1F*)hRaghav_test_20->Clone("hPbp_raghav_comb_20");
  TH1F* hPbp_raghav_comb_20_60_100 = (TH1F*)hRaghav_test_20_60_100->Clone("hPbp_raghav_comb_20_60_100");
  TH1F* hPbp_raghav_comb_40 = (TH1F*)fin2->Get("hpPb_Comb_40");
  TH1F* hPbp_raghav_comb_MB = (TH1F*)fin2->Get("hpPb_Comb_MB");

  TH1F* hPbp_raghav_20 = (TH1F*)fin2->Get("hpPb_Jet20");
  TH1F* hPbp_raghav_80 = (TH1F*)fin2->Get("hpPb_Jet80");
  TH1F* hPbp_raghav_100 = (TH1F*)fin2->Get("hpPb_Jet100");

  TH1F* hPbp_kurt_comb = (TH1F*)hKurt_test->Clone("hPbp_kurt_comb");
  TH1F* hPbp_kurt_20_40 = (TH1F*)fin2->Get("hpPb_Kurt20_40");
  TH1F* hPbp_kurt_40_60 = (TH1F*)fin2->Get("hpPb_Kurt40_60");
  TH1F* hPbp_kurt_60_80 = (TH1F*)fin2->Get("hpPb_Kurt60_80");
  TH1F* hPbp_kurt_80_100 = (TH1F*)fin2->Get("hpPb_Kurt80_100");
  TH1F* hPbp_kurt_100 = (TH1F*)fin2->Get("hpPb_Kurt100");

  cout<<"got all the histograms from file Pbp"<<endl;


  //hpPb_raghav_20 = (TH1F*)hpPb_raghav_20->Rebin(nbins_yaxian,"hpPb_raghav_20",boundaries_yaxian);
  //hpPb_raghav_60 = (TH1F*)hpPb_raghav_60->Rebin(nbins_yaxian,"hpPb_raghav_60",boundaries_yaxian);
  //hpPb_raghav_100 = (TH1F*)hpPb_raghav_100->Rebin(nbins_yaxian,"hpPb_raghav_100",boundaries_yaxian);

  //TH1F* hpPb_raghav_20_60_100 = new TH1F("hpPb_raghav_20_60_100","",1000,0,1000);
  //hpPb_raghav_20_60_100->Add(hpPb_raghav_20);
  //hpPb_raghav_20_60_100->Add(hpPb_raghav_60);
  //hpPb_raghav_20_60_100->Add(hpPb_raghav_100);

  //hpPb_raghav_comb_20_60_100 = (TH1F*)hpPb_raghav_comb_20_60_100->Rebin(nbins_yaxian,"hpPb_raghav_comb_20_60_100",boundaries_yaxian);

  //divideBinWidth(hpPb_raghav_comb_20_60_100);

  //hPbp_raghav_comb_20_60_100 = (TH1F*)hPbp_raghav_comb_20_60_100->Rebin(nbins_yaxian,"hPbp_raghav_comb_20_60_100",boundaries_yaxian);

  //divideBinWidth(hPbp_raghav_comb_20_60_100);


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

  cout<<"passed the divide bin width stuff"<<endl;


  TH1F* hRaghav = rebin_yaxian(hRaghav_test_20,"hRaghav");
  divideBinWidth(hRaghav);
  //TH1F* hRaghav = (TH1F*)hpPb_raghav_comb_20_60_100->Clone("hRaghav");

  TH1F* hRaghav_B = rebin_yaxian(hRaghav_test_20_B,"hRaghav_B");
  divideBinWidth(hRaghav_B);
  TH1F* hYaxian = rebin_yaxian(hYaxian_test,"hYaxian");
  divideBinWidth(hYaxian);
  TH1F* hYaxian_B = rebin_yaxian(hYaxian_test_B,"hYaxian_B"); 
  divideBinWidth(hYaxian_B);
  TH1F* hKurt = rebin_yaxian(hKurt_test,"hKurt");
  divideBinWidth(hKurt);
  TH1F* hKurt_B = rebin_yaxian(hKurt_test_B,"hKurt_B");
  divideBinWidth(hKurt_B);

  cout<<"passed getting histograms for the three trigger merging methods"<<endl;

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

  hTurnon100->Add(hPbp_Jet100);
  hTurnon80->Add(hPbp_Jet80);
  hTurnon60->Add(hPbp_Jet60);
  hTurnon40->Add(hPbp_Jet40);
  hTurnon20->Add(hPbp_Jet20);

  hpPb_comb->Add(hPbp_comb);

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

  hKurt->Scale(1./2);
  hRaghav->Scale(1./2);
  hYaxian->Scale(1./2);

  hKurt->Scale(1./20.97e6);
  hRaghav->Scale(1./20.97e6);
  hYaxian->Scale(1./20.97e6);

  hRaghav_B->Scale(1./2);
  hRaghav_B->Scale(1./14.1e6);

  hYaxian_B->Scale(1./2);
  hYaxian_B->Scale(1./14.1e6);

  hKurt_B->Scale(1./2);
  hKurt_B->Scale(1./14.1e6);
  
  
  hRaghav->Add(hRaghav_B);
  hYaxian->Add(hYaxian_B);
  hKurt->Add(hKurt_B);
  

  hRaghav->SetMarkerStyle(20);
  hRaghav->SetMarkerColor(kCyan);
  hKurt->SetMarkerStyle(20);
  hKurt->SetMarkerColor(kBlack);
  hYaxian->SetMarkerColor(kRed);
  hYaxian->SetMarkerStyle(20);

  hRaghav->SetTitle("");
  hRaghav->SetXTitle("Jet p_{T} (GeV/c)");
  hRaghav->SetYTitle("#frac{d^{2}#sigma}{dp_{T} d#eta} (mb/GeV/c)");
 
  hRaghav->Draw();
  hKurt->Draw("same");
  hYaxian->Draw("same");

  TLegend *title = myLegend(0.3,0.65,0.75,0.9);
  title->AddEntry(hRaghav,"12-003 (HLT_100, HLT_80 and HLT_20)","pl");
  //title->AddEntry(hRaghav,"12-003 (HLT_100, HLT_80 and HLT_20)","pl");
  title->AddEntry(hKurt,"14-007","pl");
  title->AddEntry(hYaxian,"12-017","pl");
  title->SetTextSize(0.04);
  title->Draw();

  //drawText("pPb 2013 Data",0.3,0.65,20);
  //drawText("Anti-k_{T} Pu PF Jets R = 0.3, |#eta_{CM}<1|, |vz|<15",0.3,0.56,20);

  putCMSPrel(0.1,0.9,0.06);
  drawText("pPb #int dt = 35 nb^{-1}, #sqrt{s_{NN}}=5.02 TeV",0.6,0.91,16);
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

  hRaghav_ratio->SetXTitle("Jet p_{T} (GeV/c)");
  hRaghav_ratio->SetYTitle("Ratio");
  hRaghav_ratio->SetTitle(" ");
  //hRaghav_ratio->SetAxisRange(20,300,"X");
  hRaghav_ratio->SetAxisRange(0.8,1.2,"Y");
  hRaghav_ratio->Draw();
  hYaxian_ratio->Draw("same");

  TLegend *title2 = myLegend(0.14,0.65,0.45,0.85);
  title2->AddEntry(hRaghav_ratio,"12-003/14-007","pl");
  title2->AddEntry(hYaxian_ratio,"12-017/14-007","pl");
  title2->SetTextSize(0.04);
  title2->Draw();

  c1->SaveAs("RpA_trigger_Crosscheck_full.pdf","RECREATE");

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
  hTurnon100->SetXTitle("offline Jet p_{T} (GeV/c)");
  hTurnon100->SetYTitle("Trigger Efficiency");
  hTurnon100->GetXaxis()->SetTitleOffset(1.3);
  hTurnon100->Draw("");
  hTurnon80->Draw("same");
  hTurnon60->Draw("same");
  hTurnon40->Draw("same");
  hTurnon20->Draw("same");
  hTurnon40->Draw("same");
  hTurnon60->Draw("same");
  hTurnon80->Draw("same");
  hTurnon100->Draw("same");
  
  TLine *line = new TLine(0,1,150,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();

  TLegend *title3 = myLegend(0.68,0.15,0.92,0.45);
  title3->AddEntry(hTurnon20,"HLT_PAJet20","pl");
  title3->AddEntry(hTurnon40,"HLT_PAJet40","pl");
  title3->AddEntry(hTurnon60,"HLT_PAJet60","pl");
  title3->AddEntry(hTurnon80,"HLT_PAJet80","pl");
  title3->AddEntry(hTurnon100,"HLT_PAJet100","pl");
  title3->SetTextSize(0.03);
  title3->Draw();
  
  putCMSPrel(0.1,0.92,0.06);
  drawText("pPb #int dt = 35 nb^{-1}, #sqrt{s_{NN}}=5.02 TeV",0.5,0.93,16);

  c2->SaveAs("RpA_trigger_turnon.pdf","RECREATE");

  // plot the hlt trigger merging for 12-003 and 14-007
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->SetLogy();
  //TF1 *fPowerLaw = new TF1("fPowerLaw","[0]/pow(x,[1])");
  //hpPb_kurt_comb->Fit("fPowerLaw","","",35,350);
  //hpPb_kurt_comb->Fit("fPowerLaw","","",35,350);
  //hpPb_kurt_comb->Fit("fPowerLaw","","",35,350);
  //hpPb_kurt_comb->Fit("fPowerLaw","","",35,350);
  //hpPb_kurt_comb->Fit("fPowerLaw","","",35,350);

  hpPb_kurt_comb->Scale(1./2);
  hpPb_kurt_comb->Scale(1./20.97e6);

  hpPb_kurt_20_40->Scale(1./2);
  hpPb_kurt_20_40->Scale(1./20.97e6);

  hpPb_kurt_40_60->Scale(1./2);
  hpPb_kurt_40_60->Scale(1./20.97e6);

  hpPb_kurt_60_80->Scale(1./2);
  hpPb_kurt_60_80->Scale(1./20.97e6);

  hpPb_kurt_80_100->Scale(1./2);
  hpPb_kurt_80_100->Scale(1./20.97e6);

  hpPb_kurt_100->Scale(1./2);
  hpPb_kurt_100->Scale(1./20.97e6);

  //hPbp_kurt_comb->Scale(1./2);
  //hPbp_kurt_comb->Scale(1./14.1e6);

  hPbp_kurt_20_40->Scale(1./2);
  hPbp_kurt_20_40->Scale(1./14.1e6);

  hPbp_kurt_40_60->Scale(1./2);
  hPbp_kurt_40_60->Scale(1./14.1e6);

  hPbp_kurt_60_80->Scale(1./2);
  hPbp_kurt_60_80->Scale(1./14.1e6);

  hPbp_kurt_80_100->Scale(1./2);
  hPbp_kurt_80_100->Scale(1./14.1e6);

  hPbp_kurt_100->Scale(1./2);
  hPbp_kurt_100->Scale(1./14.1e6);

  TH1F* hPbp_my_Kurt_Comb = (TH1F*)hPbp_kurt_100->Clone("hPbp_my_Kurt_Comb");
  hPbp_my_Kurt_Comb->Add(hPbp_kurt_80_100);
  hPbp_my_Kurt_Comb->Add(hPbp_kurt_60_80);
  hPbp_my_Kurt_Comb->Add(hPbp_kurt_40_60);
  hPbp_my_Kurt_Comb->Add(hPbp_kurt_20_40);


  
  hpPb_kurt_comb->Add(hPbp_my_Kurt_Comb);
  hpPb_kurt_20_40->Add(hPbp_kurt_20_40);
  hpPb_kurt_40_60->Add(hPbp_kurt_40_60);
  hpPb_kurt_60_80->Add(hPbp_kurt_60_80);
  hpPb_kurt_80_100->Add(hPbp_kurt_80_100);
  hpPb_kurt_100->Add(hPbp_kurt_100);
  

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
  hpPb_kurt_comb->SetYTitle("#frac{d^{2}#sigma}{dp_{T} d#eta} (mb/GeV/c)");
  hpPb_kurt_comb->SetAxisRange(20,350,"X");

  
  hpPb_kurt_comb->Draw();
  hpPb_kurt_20_40->Draw("same");
  hpPb_kurt_40_60->Draw("same");
  hpPb_kurt_60_80->Draw("same");
  hpPb_kurt_80_100->Draw("same");
  hpPb_kurt_100->Draw("same");
  hpPb_kurt_comb->Draw("same");
  
  /*
  hPbp_my_Kurt_Comb->SetMarkerColor(kRed);
  hPbp_my_Kurt_Comb->Draw();
  hPbp_kurt_20_40->Draw("same");
  hPbp_kurt_40_60->Draw("same");
  hPbp_kurt_60_80->Draw("same");
  hPbp_kurt_80_100->Draw("same");
  hPbp_kurt_100->Draw("same");
  hPbp_my_Kurt_Comb->Draw("same");
  */


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
  drawText("pPb #int dt = 35 nb^{-1}, |#eta_{CM}|<1, 0-100%",0.5,0.5,16);
  //drawText();

  c3->SaveAs("RpA_14007_method_combination_full.pdf","RECREATE");


  /*
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
  hRaghav->SetMarkerStyle(20);
  hRaghav->SetMarkerColor(kRed);
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
  hpPb_kurt_comb->Draw("same");
  //hRaghav->Draw("same");
  
  TLegend* title5 = myLegend(0.5,0.6,0.9,0.9);
  title5->AddEntry(hpPb_kurt_comb,"14-007 method","pl");
  //title5->AddEntry(hpPb_raghav_comb_MB,"12-003 HLT100 + HLT80 + HLT MB","pl");
  title5->AddEntry(hpPb_raghav_comb_20,"12-003 HLT100 + HLT80 + HLT 20","pl");
  title5->AddEntry(hpPb_raghav_comb_40,"12-003 HLT100 + HLT80 + HLT 40","pl");
  title5->AddEntry(hpPb_raghav_comb_60,"12-003 HLT100 + HLT80 + HLT 60","pl");
  //title5->AddEntry(hRaghav,"12-003 HLT100 + HLT60 + HLT 20","pl");
  title5->Draw();

  putCMSPrel(0.1,0.9,0.06);
  drawText("#sqrt{s_{NN}}=5.02 TeV",0.5,0.91,16);
  drawText("pPb #int dt = 20.7 nb^{-1}, |#eta_{CM}|<1, 0-100%",0.5,0.5,16);
  c4->SaveAs("RpA_12003_comparingwith_14007.pdf","RECREATE");
  */

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

  hpPb_raghav_comb_20->SetAxisRange(30,300,"X");
  hpPb_raghav_comb_20->SetYTitle("#frac{d^{2}#sigma}{dp_{T} d#eta} (mb/GeV/c)");
  hpPb_raghav_comb_20->SetXTitle("Jet p_{T} (GeV/c)");

  hpPb_raghav_comb_20->Scale(1./2);
  hpPb_raghav_comb_20->Scale(1./20.97e6);
  
  hpPb_raghav_20->Scale(1./2);
  hpPb_raghav_20->Scale(1./20.97e6);

  hpPb_raghav_80->Scale(1./2);
  hpPb_raghav_80->Scale(1./20.97e6);

  hpPb_raghav_100->Scale(1./2);
  hpPb_raghav_100->Scale(1./20.97e6);

  //hPbp_raghav_comb_20->Scale(1./2);
  //hPbp_raghav_comb_20->Scale(1./14.1e6);
  
  hPbp_raghav_20->Scale(1./2);
  hPbp_raghav_20->Scale(1./14.1e6);

  hPbp_raghav_80->Scale(1./2);
  hPbp_raghav_80->Scale(1./14.1e6);

  hPbp_raghav_100->Scale(1./2);
  hPbp_raghav_100->Scale(1./14.1e6);

  TH1F* hPbp_my_raghav_comb = (TH1F*)hPbp_raghav_100->Clone("hPbp_my_raghav_comb");
  hPbp_my_raghav_comb->Add(hPbp_raghav_80);
  hPbp_my_raghav_comb->Add(hPbp_raghav_20);

  
  hpPb_raghav_comb_20->Add(hPbp_my_raghav_comb);
  hpPb_raghav_20->Add(hPbp_raghav_20);
  hpPb_raghav_80->Add(hPbp_raghav_80);
  hpPb_raghav_100->Add(hPbp_raghav_100);
  
  
  hpPb_raghav_comb_20->Draw();
  hpPb_raghav_20->Draw("same");
  hpPb_raghav_80->Draw("same");
  hpPb_raghav_100->Draw("same");
  hpPb_raghav_comb_20->Draw("same");
  
  /*
  hPbp_my_raghav_comb->SetMarkerColor(kRed);
  hPbp_my_raghav_comb->Draw();
  hPbp_raghav_20->Draw("same");
  hPbp_raghav_80->Draw("same");
  hPbp_raghav_100->Draw("same");
  hPbp_my_raghav_comb->Draw("same");
  */

  TLegend* title6 = myLegend(0.5,0.6,0.9,0.9);
  //title6->AddEntry(hpPb_kurt_comb,"14-007 method","pl");
  title6->AddEntry(hpPb_raghav_comb_20,"12-003 method trigger combination","pl");
  title6->AddEntry(hpPb_raghav_20,"HLT_PAJet20 && !HLT_PAJet80 && !HLT_PAJet100","pl");
  title6->AddEntry(hpPb_raghav_80,"HLT_PAJet80 && !HLT_PAJet100","pl");
  title6->AddEntry(hpPb_raghav_100,"HLT_PAJet100","pl");
  title6->Draw();

  putCMSPrel(0.1,0.9,0.06);
  drawText("#sqrt{s_{NN}}=5.02 TeV",0.5,0.91,16);
  drawText("pPb #int dt = 35 nb^{-1}, |#eta_{CM}|<1, 0-100%",0.5,0.5,16);
  c5->SaveAs("RpA_12003_method_combination_full.pdf","RECREATE");


  //TFile* fin2 = TFile::Open("pPb_eta_rap_vals.root");
  
  TH2F* heta_rap_ratio = new TH2F("heta_rap_ratio","",100,0,500,100,0.9,1.03);
  TH2F* heta_rap_ratio_2 = new TH2F("heta_rap_ratio_2","",100,0,500,100,0.93,1.03);
  TH1F* heta_rap_diff = new TH2F("heta_rap_diff","",100,0,500,100,0,0.03);
  cout<<"A"<<endl;
  TTree* jet = (TTree*)fin->Get("t1");
  jet->Print();
  //jet->Draw("(jeteta/jetrap):jetpt>>heta_rap_ratio","abs(jeteta)<1&&jetpt>=30&&jetpt<500","colz");
  jet->Project("heta_rap_ratio","(jeteta/jetrap):jetpt","abs(jeteta)<1&&jetpt>=30&&jetpt<500");
  cout<<"B"<<endl;
  heta_rap_ratio->Print("base");
  cout<<"C"<<endl;
  //jet->Draw("(jeteta/jetrap):jetpt>>heta_rap_ratio_2","abs(jeteta)<0.5&&jetpt>=30&&jetpt<500","colz");
  jet->Project("heta_rap_ratio_2","(jeteta/jetrap):jetpt","abs(jeteta)<0.5&&jetpt>=30&&jetpt<500");
  heta_rap_ratio_2->Print("base");
  //jet->Draw("((jeteta-jetrap)/jeteta):jetpt>>heta_rap_diff","abs(jeteta)<1&&jetpt>=30&&jetpt<500","colz");
  jet->Project("heta_rap_diff","((jeteta-jetrap)/jeteta):jetpt","abs(jeteta)<1&&jetpt>=30&&jetpt<500");
  heta_rap_diff->Print("base");

  //TH2F* heta_rap_ratio2,*heta_rap_ratio2_2,*heta_rap_diff2;
  TH2D* heta_rap_ratio2 = new TH2D("heta_rap_ratio2","",100,0,500,100,0.9,1.03);
  TH2D* heta_rap_ratio2_2 = new TH2D("heta_rap_ratio2_2","",100,0,500,100,0.93,1.03);
  TH1D* heta_rap_diff2 = new TH2D("heta_rap_diff2","",100,0,500,100,0,0.03);
  TTree* jet2 = (TTree*)fin2->Get("t1");
  jet2->Print();
  jet2->Draw("(jeteta/jetrap):jetpt>>heta_rap_ratio2","abs(jeteta)<1&&jetpt>=30&&jetpt<500","colz");
  heta_rap_ratio2->Print("base");
  jet2->Draw("(jeteta/jetrap):jetpt>>heta_rap_ratio2_2","abs(jeteta)<0.5&&jetpt>=30&&jetpt<500","colz");
  heta_rap_ratio2_2->Print("base");
  jet2->Draw("((jeteta-jetrap)/jeteta):jetpt>>heta_rap_diff2","abs(jeteta)<1&&jetpt>=30&&jetpt<500","colz");
  heta_rap_diff2->Print("base");

  heta_rap_ratio->Add(heta_rap_ratio2);
  heta_rap_ratio_2->Add(heta_rap_ratio2_2);
  heta_rap_diff->Add(heta_rap_diff2);

  
  TH2F* heta_vs_rap = (TH2F*)fin->Get("heta_vs_rap");

  heta_vs_rap->Print("base");
  
  TH1F* hpPb_eta_05 = new TH1F("hpPb_eta_05","",nbins_yaxian,boundaries_yaxian); 
  TH1F* hpPb_rap_05 = new TH1F("hpPb_rap_05","",nbins_yaxian,boundaries_yaxian);

  TH1F* hPbp_eta_05 = new TH1F("hPbp_eta_05","",nbins_yaxian,boundaries_yaxian); 
  TH1F* hPbp_rap_05 = new TH1F("hPbp_rap_05","",nbins_yaxian,boundaries_yaxian);

  jet->Draw("jetpt>>hpPb_eta_05","fabs(jeteta)<0.5&&jetpt>=30","");
  jet->Draw("jetpt>>hpPb_rap_05","fabs(jetrap)<0.5&&jetpt>=30","");

  jet2->Draw("jetpt>>hPbp_eta_05","fabs(jeteta)<0.5&&jetpt>=30","");
  jet2->Draw("jetpt>>hPbp_rap_05","fabs(jetrap)<0.5&&jetpt>=30","");

  hpPb_eta_05->Print("base");
  hpPb_rap_05->Print("base");

  divideBinWidth(hpPb_eta_05); 
  divideBinWidth(hpPb_rap_05);

  hpPb_eta_05->Scale(1./1); // delta eta 
  hpPb_eta_05->Scale(1./20.97); //lumi value to take it to cross section. 
  
  hpPb_rap_05->Scale(1./1);
  hpPb_rap_05->Scale(1./20.97);

  hPbp_eta_05->Print("base");
  hPbp_rap_05->Print("base");

  divideBinWidth(hPbp_eta_05); 
  divideBinWidth(hPbp_rap_05);

  hPbp_eta_05->Scale(1./1); // delta eta 
  hPbp_eta_05->Scale(1./14.1); //lumi value to take it to cross section. 
  
  hPbp_rap_05->Scale(1./1);
  hPbp_rap_05->Scale(1./14.1);

  hpPb_eta_05->Add(hPbp_eta_05);
  hpPb_rap_05->Add(hPbp_rap_05);
 
  
  TCanvas *c10 = new TCanvas("c10","",800,600);

  formatCanvas(c10);


  c10->cd(1);


  hpPb_eta_05->SetMarkerStyle(24);
  hpPb_eta_05->SetMarkerColor(kBlue);
  hpPb_eta_05->SetYTitle("#frac{d #sigma}{dp_{T} d #eta} (nb)");
  hpPb_eta_05->SetXTitle("Jet p_{T} (GeV/c)");
  hpPb_eta_05->SetAxisRange(20,600,"X");
  hpPb_eta_05->Draw("");

  hpPb_rap_05->SetMarkerStyle(20);
  hpPb_rap_05->SetMarkerColor(kRed);
  hpPb_rap_05->Draw("same");

  putCMSPrel(0.1,0.91,0.06);
  drawText("#sqrt{s_{NN}}=5.02 (TeV)",0.6,0.91,16);
  drawText("pPb #int dt = 35 nb^{-1} 0-100%",0.6,0.8,16);

  TLegend* title7 = myLegend(0.6,0.5,0.8,0.75);
  title7->AddEntry(hpPb_eta_05,"|#eta_{CM}|<0.5","pl");
  title7->AddEntry(hpPb_rap_05,"|y_{CM}|<0.5","pl");
  title7->SetTextSize(0.06);
  title7->Draw();

  c10->cd(2);

  TH1F* hpPb_eta_rap_ratio = (TH1F*)hpPb_eta_05->Clone("hpPb_eta_rap_ratio");
  hpPb_eta_rap_ratio->Divide(hpPb_rap_05);
  hpPb_eta_rap_ratio->SetMarkerStyle(20);
  hpPb_eta_rap_ratio->SetMarkerColor(kBlack);
  hpPb_eta_rap_ratio->SetXTitle("Jet p_{T} (GeV/c)");
  hpPb_eta_rap_ratio->SetYTitle("#eta selection over y selection");
  hpPb_eta_rap_ratio->SetAxisRange(0.95,1.05,"Y");
  hpPb_eta_rap_ratio->SetAxisRange(20,600,"X");
  hpPb_eta_rap_ratio->Draw();

  c10->SaveAs("RpA_pt_eta_vs_rap_comparison.pdf","RECREATE");


  TCanvas *c6 = new TCanvas("c6","",800,600);

  c6->SetLogz();
  heta_vs_rap->SetXTitle("#eta Pseudorapidity");
  heta_vs_rap->SetYTitle("y Rapidity");
  heta_vs_rap->Draw("colz");
  putCMSPrel(0.1,0.9,0.06);
  drawText("#sqrt{s_{NN}}=5.02 TeV",0.5,0.91,16);
  drawText("pPb #int dt = 35 nb^{-1}, |#eta_{CM}|<5, 0-100%",0.5,0.85,16);
  c6->SaveAs("RpA_eta_vs_rapidity_comp.pdf","RECREATE");

  TCanvas *c7 = new TCanvas("c7","",800,600);
  c7->SetLogz();
  heta_rap_ratio->SetXTitle("Jet p_{T} (GeV/c)");
  heta_rap_ratio->SetYTitle("#frac{#eta}{y}");
  heta_rap_ratio->SetTitle(" ");
  heta_rap_ratio->GetYaxis()->SetTitleOffset(1.4);
  heta_rap_ratio->Draw("colz");
  heta_rap_ratio->ProfileX()->DrawClone("same");
  putCMSPrel(0.1,0.9,0.06);
  drawText("#sqrt{s_{NN}}=5.02 TeV",0.6,0.85,16);
  drawText("pPb #int dt = 35 nb^{-1}, |#eta_{CM}|<1, 0-100%",0.5,0.92,16);
  c7->SaveAs("RpA_eta_vs_rapidity_ratio.pdf","RECREATE");

  TCanvas *c12 = new TCanvas("c12","",800,600);
  heta_rap_ratio->ProfileX()->SetAxisRange(1.0,1.008,"Y");
  heta_rap_ratio->ProfileX()->DrawClone();
  putCMSPrel(0.1,0.9,0.06);
  drawText("#sqrt{s_{NN}}=5.02 TeV",0.6,0.85,16);
  drawText("Profile of #frac{#eta}{y}",0.5,0.75,16);
  drawText("pPb #int dt = 35 nb^{-1}, |#eta_{CM}|<1, 0-100%",0.5,0.92,16);
  c12->SaveAs("RpA_eta_rap_ratio_profile.pdf","RECREATE");

  TCanvas *c9 = new TCanvas("c9","",800,600);
  c9->SetLogz();
  heta_rap_ratio_2->SetXTitle("Jet p_{T} (GeV/c)");
  heta_rap_ratio_2->SetYTitle("#frac{#eta}{y}");
  heta_rap_ratio_2->SetTitle(" ");
  heta_rap_ratio_2->GetYaxis()->SetTitleOffset(1.4);
  heta_rap_ratio_2->Draw("colz");
  heta_rap_ratio_2->ProfileX()->DrawClone("same");
  putCMSPrel(0.1,0.9,0.06);
  drawText("#sqrt{s_{NN}}=5.02 TeV",0.6,0.85,16);
  drawText("pPb #int dt = 35 nb^{-1}, |#eta_{CM}|<0.5, 0-100%",0.5,0.92,16);
  c9->SaveAs("RpA_eta_vs_rapidity_ratio_2.pdf","RECREATE");


  TCanvas* c8 = new TCanvas("c8","",800,600);
  c8->SetLogz();
  heta_rap_diff->SetXTitle("Jet p_{T} (GeV/c)");
  heta_rap_diff->SetYTitle("#frac{#eta - y}{#eta}");
  heta_rap_diff->SetTitle(" ");
  heta_rap_diff->GetYaxis()->SetTitleOffset(1.4);
  heta_rap_diff->Draw("colz");
  heta_rap_diff->ProfileX()->DrawClone("same");
  putCMSPrel(0.1,0.915,0.06);
  drawText("#sqrt{s_{NN}}=5.02 TeV",0.5,0.85,16);
  drawText("pPb #int dt = 35 nb^{-1}, |#eta_{CM}|<1, 0-100%",0.5,0.915,16);
  c8->SaveAs("RpA_eta_vs_rapidity_diff.pdf","RECREATE");

  
  
  
  
  TFile fout("RpA_trigcombination_check.root","RECREATE");
  hRaghav->Write();
  hKurt->Write();
  hYaxian->Write();
  hRaghav_ratio->Write();
  hYaxian_ratio->Write();
  fout.Write();
  fout.Close();

  
 


}
