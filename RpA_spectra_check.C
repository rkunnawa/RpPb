// Raghav Kunnawalkam Elayavalli
// created 4th March 2014 

// macro to cross check different spectras 


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

void RpA_spectra_check(){

  TH1::SetDefaultSumw2();

  // RpA = 1/A * 1/L * [d^2 N_pA / (dpt deta)] / [d^2 pp / (dpt deta) ]
  // where A = 208
  //       L = 15.612 e-9 barns
  //       N_pA -> pA spectra normalized with pt bin width and eta width = 2. 
  //       pp -> pp spectra 

  //lets get the required histograms

  TFile *finMC = TFile::Open("RpA_pp_MC_NLO_reference_5020GeV.root");
  TFile *finData = TFile::Open("trig_merge_crosscheck_purdueforests_merge.root");
  TFile *fPP = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/CMSSW_6_0_0/src/result-2013-ppb-ak3PF-cent-1/ppb_merge_ak3PF_MB_correctedMC_weighting_eta_CM_1_lowest_pp_mc_Unfo_ak3PF_cent_1.root");
  TFile *fYaxian = TFile::Open("PPbJetTrigPYTHIAak3PFJetSpectraRpAHFsumEta4Bin1.root");

  TFile fout("RpA_spectra_check.root","RECREATE");

  TH1F* hRaghav_test = (TH1F*)finData->Get("hpPb_Comb");
  TH1F* hKurt_test = (TH1F*)finData->Get("hpPb_KurtComb");
  TH1F* hYaxian_test = (TH1F*)finData->Get("hpPb_TrkComb");

  TH1F* hRaghav = rebin_yaxian(hRaghav_test,"hRaghav");
  TH1F* hKurt = rebin_yaxian(hKurt_test,"hKurt");
  TH1F* hYaxian = rebin_yaxian(hYaxian_test,"hYaxian");

  TH1F* hPP = (TH1F*)fPP->Get("hGen_cent1");
  TH1F* hPPrebin = rebin_yaxian(hPP,"hPPrebin");
  
  hRaghav->Print("base");
  TH1F* hRaghav_v2 = (TH1F*)hRaghav->Clone("hRaghav_v2");

  TH1F* hRpA_12003 = (TH1F*)hRaghav->Clone("hRpA_12003");
  TH1F* hRpA_14007 = (TH1F*)hKurt->Clone("hRpA_14007");
  TH1F* hRpA_12017 = (TH1F*)hYaxian->Clone("hRpA_12017");

  TH1F* hPPMC = (TH1F*)finMC->Get("hPPrebin");
  hPPMC->Print("base");
  TH1F* hPP_nnpdf_NLO = (TH1F*)finMC->Get("hPP_nnpdf21_NLO");
  TH1F* hPP_ct10n_NLO = (TH1F*)finMC->Get("hPP_ct10n_NLO");
  TH1F* hPP_hera15all_NLO = (TH1F*)finMC->Get("hPP_hera15all_NLO");
  //TH1F* hPPMCrebin = rebin_yaxian(hPPMC,"hPPMCrebin");

  hRaghav->Scale(1./208);
  hRaghav->Scale(1./15.612e-9);
  hRaghav->Scale(1./2);
  divideBinWidth(hRaghav);

  hRaghav_v2->Scale(1./208);
  hRaghav_v2->Scale(1./15.612e-9);
  hRaghav_v2->Scale(1./2);
  divideBinWidth(hRaghav_v2);

  hRpA_12003->Scale(1./208);
  hRpA_12003->Scale(1./15.612e-9);
  hRpA_12003->Scale(1./2);
  divideBinWidth(hRpA_12003);

  hRpA_12017->Scale(1./208);
  hRpA_12017->Scale(1./15.612e-9);
  hRpA_12017->Scale(1./2);
  divideBinWidth(hRpA_12017);

  hRpA_14007->Scale(1./208);
  hRpA_14007->Scale(1./15.612e-9);
  hRpA_14007->Scale(1./2);
  divideBinWidth(hRpA_14007);

  TH1F* hYaxian_2 = (TH1F*)fYaxian->Get("DataJetInEtaBin-10_10;1");
  hYaxian_2->Scale(1./208);
  hYaxian_2->Scale(1./15.612e-9);
  hYaxian_2->Scale(1./2);
  divideBinWidth(hYaxian_2);
  hYaxian_2->Print("base");

  //hRaghav_v2->Scale(1./15.612e-9);
  //hRaghav_v2->Scale(1./208);

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->SetLogy();
  //hRaghav_v2->Draw();
  //hYaxian_2->SetMarkerColor(kRed);
  hYaxian_2->Draw();
  hRaghav_v2->SetMarkerStyle(22);
  hRaghav_v2->SetMarkerColor(kRed);
  TLegend *title2 = myLegend(0.34,0.65,0.65,0.9);
  title2->AddEntry(hRaghav_v2,"Raghav code - 12-017","pl");
  title2->AddEntry(hYaxian_2,"Yaxian code - 12-017(unfolded)","pl");
  title2->Draw();
  c1->SaveAs("RpA_12003_yaxian_12017_comparison.pdf","RECREATE");


  hPPMC->Scale(1./1e-6);
  hPP_nnpdf_NLO->Scale(1./1e-6);
  hRpA_12003->Divide(hPP_nnpdf_NLO);
  hRpA_12017->Divide(hPP_nnpdf_NLO);
  hRpA_14007->Divide(hPP_nnpdf_NLO);

  hRaghav->Divide(hPPMC);

  TCanvas *c2 = new TCanvas("c2","",800,600);
  hRaghav->SetMarkerStyle(25);
  hRaghav->SetMarkerColor(kBlack);
  hRaghav->SetTitle("RpA using 12-003 (measured) method with different PP references");
  hRaghav->Draw();
  hRpA_12003->SetMarkerStyle(26);
  hRpA_12003->SetMarkerColor(kRed);
  hRpA_12003->Draw("same");

  TLegend *title = myLegend(0.34,0.65,0.65,0.9);
  title->AddEntry(hRaghav,"MC Gen reference","pl");
  title->AddEntry(hRpA_12003,"NLO pp at 5.02TeV","pl");
  title->Draw();
  
  c2->SaveAs("RpA_nlo_vs_pp_gen_12003.pdf","RECREATE");

  
  
  
  

  fout.cd();
  hYaxian_2->Write();
  //hRaghav->Write();
  fout.Write();
  fout.Close();


}
