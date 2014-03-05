// Raghav Kunnawalkam Elayavalli
// created March 4th 2014

// macro to cross check with kurt the validity of my analysis. 

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

static const int nbins_yaxian = 29;
static const double boundaries_yaxian[nbins_yaxian+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};


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


void cross_check_wth_kurt(){

  TH1::SetDefaultSumw2();

  TFile fout("RpA_kurt_cross_check_histos.root","RECREATE");
  
  TFile * fPP = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/CMSSW_6_0_0/src/result-2013-ppb-ak3PF-cent-1/ppb_merge_ak3PF_MB_correctedMC_weighting_eta_CM_1_lowest_pp_mc_Unfo_ak3PF_cent_1.root");
  //TFile *fPP = TFile::Open("ppb_merge_ak3PF_MB_correctedMC_weighting_eta_CM_1_lowest_pp_mc_Unfo_ak3PF_cent_1-5.root");
  TFile *finData = TFile::Open("trig_merge_crosscheck_purdueforests_merge.root");

  TH1F* hRaghav= (TH1F*)finData->Get("hpPb_Comb");
  TH1F* hKurt= (TH1F*)finData->Get("hpPb_KurtComb");
  TH1F* hYaxian= (TH1F*)finData->Get("hpPb_TrkComb");

  TH1F* hPP = (TH1F*)fPP->Get("hGen_cent1");


  hRaghav = (TH1F*)hRaghav->Rebin(nbins_yaxian,"hRaghav",boundaries_yaxian);
  hKurt = (TH1F*)hKurt->Rebin(nbins_yaxian,"hKurt",boundaries_yaxian);

  hPP = (TH1F*)hPP->Rebin(nbins_yaxian,"hPP",boundaries_yaxian);

  //divideBinWidth(hRaghav);
  //divideBinWidth(hPP);

  hRaghav->Scale(1./15.61E6);
  hRaghav->Divide(hPP);
  hRaghav->Scale(1./208);

  
  

  fout.cd();
  hRaghav->Write();
  hKurt->Write();
  hYaxian->Write();
  hPP->Write();
  fout.Write();
  
  fout.Close();




}
