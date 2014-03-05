// Raghav Kunnawalkam Elayavalli
// created: Marth 3rd 2014

// macro to read the pp NLO at 5.02 TeV.



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


void RpA_NLO_cross_section_comparison(){

  timer.Start();

  TH1::SetDefaultSumw2();

  //output file 
  TFile fout("RpA_pp_MC_NLO_reference_5020GeV.root","RECREATE");

  // load the files required:
  TFile * fnnpdf21_0 = TFile::Open("fnl5350eta0_nnpdf21-nlo_aspdf.root"); 
  TFile * fnnpdf21_1 = TFile::Open("fnl5350eta1_nnpdf21-nlo_aspdf.root"); 
  TFile * fnnpdf21_2 = TFile::Open("fnl5350eta2_nnpdf21-nlo_aspdf.root"); 
  TFile * fnnpdf21_3 = TFile::Open("fnl5350eta3_nnpdf21-nlo_aspdf.root"); 
  TFile * fnnpdf21_4 = TFile::Open("fnl5350eta4_nnpdf21-nlo_aspdf.root"); 

  TFile * fct10n_0 = TFile::Open("fnl5350eta0_ct10n-nlo_aspdf.root"); 
  TFile * fct10n_1 = TFile::Open("fnl5350eta1_ct10n-nlo_aspdf.root"); 
  TFile * fct10n_2 = TFile::Open("fnl5350eta2_ct10n-nlo_aspdf.root"); 
  TFile * fct10n_3 = TFile::Open("fnl5350eta3_ct10n-nlo_aspdf.root"); 
  TFile * fct10n_4 = TFile::Open("fnl5350eta4_ct10n-nlo_aspdf.root");

  TFile * fhera15all_0 = TFile::Open("fnl5350eta0_hera15all-nlo_aspdf.root"); 
  TFile * fhera15all_1 = TFile::Open("fnl5350eta1_hera15all-nlo_aspdf.root"); 
  TFile * fhera15all_2 = TFile::Open("fnl5350eta2_hera15all-nlo_aspdf.root"); 
  TFile * fhera15all_3 = TFile::Open("fnl5350eta3_hera15all-nlo_aspdf.root"); 
  TFile * fhera15all_4 = TFile::Open("fnl5350eta4_hera15all-nlo_aspdf.root");

  TFile * fPP = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/CMSSW_6_0_0/src/result-2013-ppb-ak3PF-cent-1/ppb_merge_ak3PF_MB_correctedMC_weighting_eta_CM_1_lowest_pp_mc_Unfo_ak3PF_cent_1.root");

  // we are going to look at nnpdf21 nlo calculatino and the eta bins we are interested in are 0.0 to 0.3, 0.3 to 0.7 and 0.7 to 1. there are in +eta since for pp the final output is symmetric. the histogram from the root file are named accroding to 
  // The R in histogram number xxxxRxx goes from 0.2 jet size to 0.4 jet size for R=1,2,3.
  // So the NLO for anti-kT R=0.3 is in histogram 100200
  
  // naming convension for the file names appendage corresponding to the eta bin 
  // 0 : 0.0 to 0.3
  // 1 : 0.3 to 0.7 
  // 2 : 0.7 to 1.0 
  // 3 : 1.0 to 1.2
  // 4 : 1.2 to 2.2                    
  
  // for our ourput - -1 to 1. no number
  
  // once we get the histograms from the file we have to add them to get to eta range from -1 to 1. 
  // for this we need to multiply/scale the histograms by delta eta for its respective bin. 
  
  Double_t deta_0 = 0.6;
  Double_t deta_1 = 0.8;
  Double_t deta_2 = 0.6;
  Double_t deta_3 = 0.4;
  Double_t deta_4 = 2.0;

  Double_t deta = 2;


  //get the histograms 

  TH1F* hpt_nnpdf21_0 = (TH1F*)fnnpdf21_0->Get("h100200");
  TH1F* hpt_nnpdf21_1 = (TH1F*)fnnpdf21_1->Get("h100200");
  TH1F* hpt_nnpdf21_2 = (TH1F*)fnnpdf21_2->Get("h100200");
  TH1F* hpt_nnpdf21_3 = (TH1F*)fnnpdf21_3->Get("h100200");
  TH1F* hpt_nnpdf21_4 = (TH1F*)fnnpdf21_4->Get("h100200");

  TH1F* hpt_ct10n_0 = (TH1F*)fct10n_0->Get("h100200");
  TH1F* hpt_ct10n_1 = (TH1F*)fct10n_1->Get("h100200");
  TH1F* hpt_ct10n_2 = (TH1F*)fct10n_2->Get("h100200");
  TH1F* hpt_ct10n_3 = (TH1F*)fct10n_3->Get("h100200");
  TH1F* hpt_ct10n_4 = (TH1F*)fct10n_4->Get("h100200");
  
  TH1F* hpt_hera15all_0 = (TH1F*)fhera15all_0->Get("h100200");
  TH1F* hpt_hera15all_1 = (TH1F*)fhera15all_1->Get("h100200");
  TH1F* hpt_hera15all_2 = (TH1F*)fhera15all_2->Get("h100200");
  TH1F* hpt_hera15all_3 = (TH1F*)fhera15all_3->Get("h100200");
  TH1F* hpt_hera15all_4 = (TH1F*)fhera15all_4->Get("h100200");

  TH1F* hPP = (TH1F*)fPP->Get("hGen_cent1");
  TH1F* hPPrebin = rebin_yaxian(hPP,"hPPrebin");

  /*
  hpt_0->Print("base");
  hpt_1->Print("base");
  hpt_2->Print("base");
  hpt_3->Print("base");
  hpt_4->Print("base");
  */
  

  fout.cd();

  TH1F* hPP_nnpdf21_NLO = new TH1F("hPP_nnpdf21_NLO","",nbins_yaxian,boundaries_yaxian);
  TH1F* hPP_ct10n_NLO = new TH1F("hPP_ct10n_NLO","",nbins_yaxian,boundaries_yaxian);
  TH1F* hPP_hera15all_NLO = new TH1F("hPP_hera15all_NLO","",nbins_yaxian,boundaries_yaxian);

  hpt_nnpdf21_0->Scale(deta_0);
  hpt_nnpdf21_1->Scale(deta_1);
  hpt_nnpdf21_2->Scale(deta_2);
  hpt_nnpdf21_3->Scale(deta_3);
  hpt_nnpdf21_4->Scale(deta_4);

  hPP_nnpdf21_NLO->Add(hpt_nnpdf21_0);
  hPP_nnpdf21_NLO->Add(hpt_nnpdf21_1);
  hPP_nnpdf21_NLO->Add(hpt_nnpdf21_2);

  hPP_nnpdf21_NLO->Scale(1./deta);
  hPP_nnpdf21_NLO->Print("base");

  hpt_ct10n_0->Scale(deta_0);
  hpt_ct10n_1->Scale(deta_1);
  hpt_ct10n_2->Scale(deta_2);
  hpt_ct10n_3->Scale(deta_3);
  hpt_ct10n_4->Scale(deta_4);

  hPP_ct10n_NLO->Add(hpt_ct10n_0);
  hPP_ct10n_NLO->Add(hpt_ct10n_1);
  hPP_ct10n_NLO->Add(hpt_ct10n_2);

  hPP_ct10n_NLO->Scale(1./deta);
  hPP_ct10n_NLO->Print("base");

  hpt_hera15all_0->Scale(deta_0);
  hpt_hera15all_1->Scale(deta_1);
  hpt_hera15all_2->Scale(deta_2);
  hpt_hera15all_3->Scale(deta_3);
  hpt_hera15all_4->Scale(deta_4);

  hPP_hera15all_NLO->Add(hpt_hera15all_0);
  hPP_hera15all_NLO->Add(hpt_hera15all_1);
  hPP_hera15all_NLO->Add(hpt_hera15all_2);

  hPP_hera15all_NLO->Scale(1./deta);
  hPP_hera15all_NLO->Print("base");

  hPPrebin->Scale(1./deta);
  hPPrebin->Scale(1e9);//take it to pico barns
  divideBinWidth(hPPrebin);

  
  //double integral = hPP_refe_NLO->Integral();
  //hPP_refe_NLO->Scale(1./integral);

  //divide the two histograms 

  TH1F* hRatio_nnpdf21 = (TH1F*)hPP_nnpdf21_NLO->Clone("hRatio");
  hRatio_nnpdf21->Divide(hPPrebin);

  hPP_nnpdf21_NLO->Write();

  TH1F* hRatio_ct10n = (TH1F*)hPP_ct10n_NLO->Clone("hRatio");
  hRatio_ct10n->Divide(hPPrebin);

  hPP_ct10n_NLO->Write();

  TH1F* hRatio_hera15all = (TH1F*)hPP_hera15all_NLO->Clone("hRatio");
  hRatio_hera15all->Divide(hPPrebin);

  hPP_hera15all_NLO->Write();
  hPPrebin->Write();
  hRatio->Write();

  fout.Write();
  /*
  int binvalue=hPP_ct10n_NLO->FindBin(300);
 
  Float_t val_a = hPP_ct10n_NLO->GetBinContent(binvalue);
  Float_t val_b = hPPrebin->GetBinContent(binvalue);

  Float_t err_a = hPP_ct10n_NLO->GetBinError(binvalue);
  Float_t err_b = hPPrebin->GetBinError(binvalue);


  //get the bin values to check error propagation
  cout<<"value   "<<hPP_ct10n_NLO->GetBinContent(binvalue)<<" pp mc   "<<hPPrebin->GetBinContent(binvalue)<<endl;
  cout<<" ratio = "<<(Float_t)hPP_ct10n_NLO->GetBinContent(binvalue)/hPPrebin->GetBinContent(binvalue)<<endl;

  cout<<"error:  "<<hPP_ct10n_NLO->GetBinError(binvalue)<<"  PP mc   "<<hPPrebin->GetBinError(binvalue)<<endl;
  cout<<"ratio  "<< (Float_t)(val_a/val_b)*(TMath::Sqrt((err_a/val_a)*(err_a/val_a) + (err_b/val_b)*(err_b/val_b)))<<endl;

  cout<<"check : "<<endl;

  cout<<"ratio value = "<<hRatio_ct10n->GetBinContent(binvalue)<<endl;
  cout<<"ratio error = "<<hRatio_ct10n->GetBinError(binvalue)<<endl;
  */

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Divide(2,1);
  c1->cd(1);
  c1->cd(1)->SetLogy();
  hPP_nnpdf21_NLO->SetLineColor(kRed);
  hPP_ct10n_NLO->SetLineColor(kBlue);
  hPP_hera15all_NLO->SetLineColor(kGreen);
  hPPrebin->SetMarkerColor(kBlack);
  hPPrebin->SetMarkerStyle(25);

  hPP_nnpdf21_NLO->SetYTitle("#sigma pico barns");
  hPP_nnpdf21_NLO->SetXTitle("p_{T} GeV/c");

  hPP_nnpdf21_NLO->Draw();
  hPP_ct10n_NLO->Draw("same");
  hPP_hera15all_NLO->Draw("same");
  hPPrebin->Draw("same");

  TLegend * title = myLegend(0.34, 0.65,0.65, 0.9);
  title->AddEntry(hPP_nnpdf21_NLO,"NLO nnpdf21","l");
  title->AddEntry(hPP_ct10n_NLO,"NLO ct10n","l");
  title->AddEntry(hPP_hera15all_NLO,"NLO hera15all","l");
  title->AddEntry(hPPrebin,"MC Gen spectra","pl");
  title->Draw();

  c1->cd(2);
  hRatio_nnpdf21->SetYTitle("NLO / MC Gen");
  hRatio_nnpdf21->SetXTitle("p_{T} GeV/c");
  hRatio_nnpdf21->SetAxisRange(0,2,"Y");
  hRatio_nnpdf21->SetLineColor(kRed);
  hRatio_hera15all->SetLineColor(kGreen);
  hRatio_ct10n->SetLineColor(kBlue);
  hRatio_nnpdf21->Draw();
  hRatio_hera15all->Draw("same");
  hRatio_ct10n->Draw("same");

  c1->SaveAs("pp_5020GeV_NLO_nnpdf21_vs_MC_gen_spectra.pdf","RECREATE");
  
  hPP_nnpdf21_NLO->Write();
  hPP_ct10n_NLO->Write();
  hPP_hera15all_NLO->Write();
  hPPrebin->Write();


  fout.Close();


  timer.Stop();
  float rtime  = timer.RealTime();
  float ctime  = timer.CpuTime();
  
  std::cout<<"\t"<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;
  

}
