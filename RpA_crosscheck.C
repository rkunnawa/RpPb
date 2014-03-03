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


static const int nbins_yaxian = 29;
static const double boundaries_yaxian[nbins_yaxian+1] = {3,4,5,7,9,12,15,18,22,27,33,39,47,55,64,74,84,97,114,133,153,174,196,220,245,272,300,429,692,1000};

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

//TStopwatch timer;

void RpA_crosscheck(){
  
  TH1::SetDefaultSumw2();

  //timer.start();


  TFile *fin = TFile::Open("RpA_trig_merge_crosscheck_purdue_forests_merge.root");

  TH1F* hRaghav_test = (TH1F*)fin->Get("hpPb_Comb");
  TH1F* hKurt = (TH1F*)fin->Get("hpPb_KurtComb");
  TH1F* hYaxian_test = (TH1F*)fin->Get("hpPb_TrkComb");

  hRaghav_test->Print();
  hKurt->Print();
  hYaxian_test->Print();

  TH1F* hRaghav = rebin(hRaghav_test,"hRaghav");
  TH1F* hYaxian = rebin(hYaxian_test,"hYaxian");


  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Divide(2,1);
  c1->cd(1)->SetLogy();

  divideBinWidth(hKurt);
  divideBinWidth(hRaghav);
  divideBinWidth(hYaxian);

  hRaghav->SetMarkerStyle(25);
  hRaghav->SetMarkerColor(kGreen);
  hKurt->SetMarkerStyle(24);
  hKurt->SetMarkerColor(kBlack);
  hYaxian->SetMarkerColor(kRed);
  hYaxian->SetMarkerStyle(26);

  hRaghav->SetTitle("");
  hRaghav->SetXTitle("p_{T} GeV/c");
  hRaghav->SetYTitle("counts");

  hRaghav->Draw();
  hKurt->Draw("same");
  hYaxian->Draw("same");

  TLegend *title = myLegend(0.34,0.65,0.65,0.9);
  title->AddEntry(hRaghav,"12-003","pl");
  title->AddEntry(hKurt,"14-007","pl");
  title->AddEntry(hYaxian,"12-017","pl");
  title->SetTextSize(0.04);
  title->Draw();

  //drawText("pPb 2013 Data",0.3,0.65,20);
  //drawText("Anti-k_{T} Pu PF Jets R = 0.3, |#eta_{CM}<1|, |vz|<15",0.3,0.56,20);

  c1->cd(2);
  TH1F* hRaghav_ratio = rebin(hRaghav,"hRaghav_ratio");
  TH1F* hYaxian_ratio = rebin(hYaxian,"hYaxian_ratio");

  hRaghav_ratio->Print();
  hYaxian_ratio->Print();
  hRaghav_ratio->Divide(hKurt);
  hYaxian_ratio->Divide(hKurt);

  hRaghav_ratio->SetMarkerColor(25);
  hRaghav_ratio->SetMarkerStyle(kGreen);
  hYaxian_ratio->SetMarkerColor(kRed);
  hYaxian_ratio->SetMarkerStyle(26);

  hRaghav_ratio->SetXTitle("p_{T} GeV/c");
  hRaghav_ratio->SetTitle("");

  hRaghav_ratio->Draw();
  hYaxian_ratio->Draw("same");

  TLegend *title2 = myLegend(0.34,0.65,0.65,0.9);
  title2->AddEntry(hRaghav_ratio,"12-003/14-007","pl");
  title2->AddEntry(hYaxian_ratio,"12-017/14-007","pl");
  title2->SetTextSize(0.04);
  title2->Draw();

  c1->SaveAs("RpA_trigger_Crosscheck.pdf","RECREATE");
  
 


}
