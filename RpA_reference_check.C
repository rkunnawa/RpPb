// Raghav Kunnawalkam Elayavalli
// created April 15th 2014

// macro to compare the references for the RpA analysis 


#include <iostream>
#include <stdio.h>
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

//#include "headers/utilities_V0.h"

void formatCanvas(TCanvas *c){
  c->Divide(1,2,0.01,0.01);
  c->cd(1);
  c->GetPad(1)->SetLogy();
  c->GetPad(1)->SetPad(0.,0.425,1.,1.);
  c->GetPad(2)->SetPad(0.,0.0,1.,0.425);
  c->GetPad(2)->SetBottomMargin(0.3);
  c->GetPad(2)->SetGridy(1);
}

static const int nbins_yaxian = 29;
static const double boundaries_yaxian[nbins_yaxian+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638,790,967};

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

void * convertToInvYield(TH1 *hist) {
    for(int i = 1; i<=hist->GetNbinsX(); i++) {
        double content = hist->GetBinContent(i);
        double pt = hist->GetBinCenter(i);
        double error = hist->GetBinError(i);
        
        double new_content = content/(2.*TMath::Pi()*pt);
        double new_error = error/(2.*TMath::Pi()*pt);
        
        hist->SetBinContent(i,new_content);
        hist->SetBinError(i,new_error);
    }
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

void RpA_reference_check(){

  //lets get all the histograms which we have created and assign them here. 
  gStyle->SetOptStat(0);

  TFile *fpp276 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/CMSSW_6_0_0/src/pp_2760GeV_nlo_histos.root");
  TFile *fpp502_R2 = TFile::Open("RpA_pp_MC_eric_NLO_reference_ak2PF_5020GeV.root");
  TFile *fpp502_R3 = TFile::Open("RpA_pp_MC_eric_NLO_reference_ak3PF_5020GeV.root");
  TFile *fpp502_R4 = TFile::Open("RpA_pp_MC_eric_NLO_reference_ak4PF_5020GeV.root");

  TFile *fppPythia_2760 = TFile::Open("AnaGENJetR357_2760GeV_Apr15_Z2Combined.root");
  TFile *fppPythia_5020 = TFile::Open("AnaGENJetR357_5020GeV_Apr15_Z2Combined.root");
  TFile *fppPythia_7000 = TFile::Open("AnaGENJetR357_7000GeV_Apr15_Z2Combined.root");

  TFile* fpp276_05_ak3 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/CMSSW_6_0_0/src/pp_2013_2760_abs_eta_05_mc_ak3PF.root");
  //TFile* fpp276_05_ak4 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/CMSSW_6_0_0/src/pp_2013_2760_abs_eta_05_mc_ak4PF.root");
  TFile* fpp276_05_ak5 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/CMSSW_6_0_0/src/pp_2013_2760_abs_eta_05_mc_ak5PF.root");

  TFile* fpp276_15_20_ak3 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/CMSSW_6_0_0/src/pp_2013_2760_abs_eta_15_20_mc_ak3PF.root");
  //TFile* fpp276_15_20_ak4 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/CMSSW_6_0_0/src/pp_2013_2760_abs_eta_15_20_mc_ak4PF.root");
  TFile* fpp276_15_20_ak5 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/CMSSW_6_0_0/src/pp_2013_2760_abs_eta_15_20_mc_ak5PF.root");

  TFile foutput("RpA_pp_nlo_ratios.root","RECREATE");

  TH1F* hPP_2760_data_05_R_3 = (TH1F*)fpp276_05_ak3->Get("Unfolded_cent1");
  hPP_2760_data_05_R_3->Scale(1./5300e6);
  divideBinWidth(hPP_2760_data_05_R_3);
  
  //TH1F* hPP_2760_data_05_R_4 = (TH1F*)fpp276_05_ak4->Get("Unfolded_cent1");
  //hPP_2760_data_05_R_4->Scale(1./5300e6);
  //divideBinWidth(hPP_2760_data_05_R_4);

  TH1F* hPP_2760_data_05_R_5 = (TH1F*)fpp276_05_ak5->Get("Unfolded_cent1");
  hPP_2760_data_05_R_5->Scale(1./5300e6);
  divideBinWidth(hPP_2760_data_05_R_5);

  TH1F* hPP_2760_data_15_20_R_3 = (TH1F*)fpp276_15_20_ak3->Get("Unfolded_cent1");
  hPP_2760_data_15_20_R_3->Scale(1./5300e6);
  divideBinWidth(hPP_2760_data_15_20_R_3);

  //TH1F* hPP_2760_data_15_20_R_4 = (TH1F*)fpp276_15_20_ak4->Get("Unfolded_cent1");
  //hPP_2760_data_15_20_R_4->Scale(1./5300e6);
  //divideBinWidth(hPP_2760_data_15_20_R_4);
  
  TH1F* hPP_2760_data_15_20_R_5 = (TH1F*)fpp276_15_20_ak5->Get("Unfolded_cent1");
  hPP_2760_data_15_20_R_5->Scale(1./5300e6);
  divideBinWidth(hPP_2760_data_15_20_R_5);

  TH1F* hPP_2760_data_05_R_3_5 = (TH1F*)hPP_2760_data_05_R_3->Clone("hPP_2760_data_05_R_3_5");
  hPP_2760_data_05_R_3_5->Divide(hPP_2760_data_05_R_5);

  //TH1F* hPP_2760_data_05_R_4_5 = (TH1F*)hPP_2760_data_05_R_4->Clone("hPP_2760_data_05_R_4_5");
  //hPP_2760_data_05_R_4_5->Divide(hPP_2760_data_05_R_5);

  TH1F* hPP_2760_data_15_20_R_3_5 = (TH1F*)hPP_2760_data_15_20_R_3->Clone("hPP_2760_data_15_20_R_3_5");
  hPP_2760_data_15_20_R_3_5->Divide(hPP_2760_data_15_20_R_5);

  //TH1F* hPP_2760_data_15_20_R_4_5 = (TH1F*)hPP_2760_data_15_20_R_4->Clone("hPP_2760_data_15_20_R_4_5");
  //hPP_2760_data_15_20_R_4_5->Divide(hPP_2760_data_15_20_R_5);

  
  //load the pp 2760 unfolded data and 
  //first lets do the NLO at 2.76 to NLO at 5.02 

  TH1F* hPP_5020_nnpdf_22_22_R2 = (TH1F*)fpp502_R2->Get("hPP5020_nnpdf_22_22");
  TH1F* hPP_5020_nnpdf_22_22_R3 = (TH1F*)fpp502_R3->Get("hPP5020_nnpdf_22_22");
  TH1F* hPP_5020_nnpdf_03_03_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_nnpdf_03_03");
  TH1F* hPP_5020_nnpdf_03_07_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_nnpdf_03_07");
  TH1F* hPP_5020_nnpdf_07_10_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_nnpdf_07_10");
  TH1F* hPP_5020_nnpdf_10_12_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_nnpdf_10_12");
  TH1F* hPP_5020_nnpdf_12_22_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_nnpdf_12_22");
  TH1F* hPP_5020_nnpdf_10_10_R3 = (TH1F*)fpp502_R3->Get("hPP_nnpdf21_NLO");
  TH1F* hPP_5020_nnpdf_22_22_R4 = (TH1F*)fpp502_R4->Get("hPP5020_nnpdf_22_22");

  TH1F* hPP_5020_ct10n_22_22_R2 = (TH1F*)fpp502_R2->Get("hPP5020_ct10n_22_22");
  TH1F* hPP_5020_ct10n_22_22_R3 = (TH1F*)fpp502_R3->Get("hPP5020_ct10n_22_22");
  TH1F* hPP_5020_ct10n_03_03_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_ct10n_03_03");
  TH1F* hPP_5020_ct10n_03_07_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_ct10n_03_07");
  TH1F* hPP_5020_ct10n_07_10_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_ct10n_07_10");
  TH1F* hPP_5020_ct10n_10_12_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_ct10n_10_12");
  TH1F* hPP_5020_ct10n_12_22_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_ct10n_12_22");
  TH1F* hPP_5020_ct10n_10_10_R3 = (TH1F*)fpp502_R3->Get("hPP_ct10n_NLO");
  TH1F* hPP_5020_ct10n_22_22_R4 = (TH1F*)fpp502_R4->Get("hPP5020_ct10n_22_22");

  TH1F* hPP_5020_hera_22_22_R2 = (TH1F*)fpp502_R2->Get("hPP5020_hera_22_22");
  TH1F* hPP_5020_hera_22_22_R3 = (TH1F*)fpp502_R3->Get("hPP5020_hera_22_22");
  TH1F* hPP_5020_hera_03_03_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_hera_03_03");
  TH1F* hPP_5020_hera_03_07_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_hera_03_07");
  TH1F* hPP_5020_hera_07_10_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_hera_07_10");
  TH1F* hPP_5020_hera_10_12_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_hera_10_12");
  TH1F* hPP_5020_hera_12_22_R3 = (TH1F*)fpp502_R3->Get("hPP_5020_hera_12_22");
  TH1F* hPP_5020_hera_10_10_R3 = (TH1F*)fpp502_R3->Get("hPP_hera15all_NLO");
  TH1F* hPP_5020_hera_22_22_R4 = (TH1F*)fpp502_R4->Get("hPP5020_hera_22_22"); 

  TH1F* hPP_2760_nnpdf_2_2_R2 = (TH1F*)fpp276->Get("hPP_nnpdf_NLO_R2");
  TH1F* hPP_2760_nnpdf_2_2_R3 = (TH1F*)fpp276->Get("hPP_nnpdf_NLO");
  TH1F* hPP_2760_nnpdf_2_2_R4 = (TH1F*)fpp276->Get("hPP_nnpdf_NLO_R4");


  TH1F* hPP_2760_ct10n_2_2_R2 = (TH1F*)fpp276->Get("hPP_ct10n_NLO_R2");
  TH1F* hPP_2760_ct10n_2_2_R3 = (TH1F*)fpp276->Get("hPP_ct10n_NLO");
  TH1F* hPP_2760_ct10n_2_2_R4 = (TH1F*)fpp276->Get("hPP_ct10n_NLO_R4");


  TH1F* hPP_2760_hera_2_2_R2 = (TH1F*)fpp276->Get("hPP_hera_NLO_R2");
  TH1F* hPP_2760_hera_2_2_R3 = (TH1F*)fpp276->Get("hPP_hera_NLO");
  TH1F* hPP_2760_hera_2_2_R4 = (TH1F*)fpp276->Get("hPP_hera_NLO_R4");

  //we have to get histograms which show the ratio between one NLO pdf to another at the same radius - 0.3 

  TH1F* hPP_2760_ct10n_nnpdf_2_2_R3 = (TH1F*)hPP_2760_ct10n_2_2_R3->Clone("hPP_2760_ct10n_nnpdf_2_2_R3");
  hPP_2760_ct10n_nnpdf_2_2_R3->Divide(hPP_2760_nnpdf_2_2_R3);

  TH1F* hPP_2760_ct10n_hera_2_2_R3 = (TH1F*)hPP_2760_ct10n_2_2_R3->Clone("hPP_2760_ct10n_hera_2_2_R3");
  hPP_2760_ct10n_hera_2_2_R3->Divide(hPP_2760_hera_2_2_R3);

  TH1F* hPP_5020_ct10n_nnpdf_22_22_R3 = (TH1F*)hPP_5020_ct10n_22_22_R3->Clone("hPP_5020_ct10n_nnpdf_22_22_R3");
  hPP_5020_ct10n_nnpdf_22_22_R3->Divide(hPP_5020_nnpdf_22_22_R3);

  TH1F* hPP_5020_ct10n_hera_22_22_R3 = (TH1F*)hPP_5020_ct10n_22_22_R3->Clone("hPP_5020_ct10n_hera_22_22_R3");
  hPP_5020_ct10n_hera_22_22_R3->Divide(hPP_5020_hera_22_22_R3);

  TH1F* hPP_5020_ct10n_nnpdf_03_03_R3 = (TH1F*)hPP_5020_ct10n_03_03_R3->Clone("hPP_5020_ct10n_nnpdf_03_03_R3");
  hPP_5020_ct10n_nnpdf_03_03_R3->Divide(hPP_5020_nnpdf_03_03_R3);

  TH1F* hPP_5020_ct10n_nnpdf_03_07_R3 = (TH1F*)hPP_5020_ct10n_03_07_R3->Clone("hPP_5020_ct10n_nnpdf_03_07_R3");
  hPP_5020_ct10n_nnpdf_03_07_R3->Divide(hPP_5020_nnpdf_03_07_R3);

  TH1F* hPP_5020_ct10n_nnpdf_07_10_R3 = (TH1F*)hPP_5020_ct10n_07_10_R3->Clone("hPP_5020_ct10n_nnpdf_07_10_R3");
  hPP_5020_ct10n_nnpdf_07_10_R3->Divide(hPP_5020_nnpdf_07_10_R3);

  TH1F* hPP_5020_ct10n_nnpdf_10_12_R3 = (TH1F*)hPP_5020_ct10n_10_12_R3->Clone("hPP_5020_ct10n_nnpdf_10_12_R3");
  hPP_5020_ct10n_nnpdf_10_12_R3->Divide(hPP_5020_nnpdf_10_12_R3);

  TH1F* hPP_5020_ct10n_nnpdf_12_22_R3 = (TH1F*)hPP_5020_ct10n_12_22_R3->Clone("hPP_5020_ct10n_nnpdf_12_22_R3");
  hPP_5020_ct10n_nnpdf_12_22_R3->Divide(hPP_5020_nnpdf_12_22_R3);

  TH1F* hPP_5020_ct10n_nnpdf_10_10_R3 = (TH1F*)hPP_5020_ct10n_10_10_R3->Clone("hPP_5020_ct10n_nnpdf_10_10_R3");
  hPP_5020_ct10n_nnpdf_10_10_R3->Divide(hPP_5020_nnpdf_10_10_R3);

  TH1F* hPP_5020_ct10n_hera_03_03_R3 = (TH1F*)hPP_5020_ct10n_03_03_R3->Clone("hPP_5020_ct10n_hera_03_03_R3");
  hPP_5020_ct10n_hera_03_03_R3->Divide(hPP_5020_hera_03_03_R3);

  TH1F* hPP_5020_ct10n_hera_03_07_R3 = (TH1F*)hPP_5020_ct10n_03_07_R3->Clone("hPP_5020_ct10n_hera_03_07_R3");
  hPP_5020_ct10n_hera_03_07_R3->Divide(hPP_5020_hera_03_07_R3);

  TH1F* hPP_5020_ct10n_hera_07_10_R3 = (TH1F*)hPP_5020_ct10n_07_10_R3->Clone("hPP_5020_ct10n_hera_07_10_R3");
  hPP_5020_ct10n_hera_07_10_R3->Divide(hPP_5020_hera_07_10_R3);

  TH1F* hPP_5020_ct10n_hera_10_12_R3 = (TH1F*)hPP_5020_ct10n_10_12_R3->Clone("hPP_5020_ct10n_hera_10_12_R3");
  hPP_5020_ct10n_hera_10_12_R3->Divide(hPP_5020_hera_10_12_R3);

  TH1F* hPP_5020_ct10n_hera_12_22_R3 = (TH1F*)hPP_5020_ct10n_12_22_R3->Clone("hPP_5020_ct10n_hera_12_22_R3");
  hPP_5020_ct10n_hera_12_22_R3->Divide(hPP_5020_hera_12_22_R3);

  TH1F* hPP_5020_ct10n_hera_10_10_R3 = (TH1F*)hPP_5020_ct10n_10_10_R3->Clone("hPP_5020_ct10n_hera_10_10_R3");
  hPP_5020_ct10n_hera_10_10_R3->Divide(hPP_5020_hera_10_10_R3);
  
  //now lets make the text files and get the values required. 
  
  //03_03
  ofstream pp_5020_ct10n_nnpdf_03_03_R3;
  pp_5020_ct10n_nnpdf_03_03_R3.open("pp_5020_ct10n_nnpdf_03_03_R3.txt");
  for(int i = 0;i<hPP_5020_ct10n_nnpdf_03_03_R3->GetNbinsX();i++){
    float bincenter = hPP_5020_ct10n_nnpdf_03_03_R3->GetBinCenter(i);
    float val_1 = hPP_5020_ct10n_nnpdf_03_03_R3->GetBinContent(i);
    val_1 = 100*TMath::Abs(1-val_1);
    float val_2 = hPP_5020_ct10n_hera_03_03_R3->GetBinContent(i);
    val_2 = 100*TMath::Abs(1-val_2);
    pp_5020_ct10n_nnpdf_03_03_R3<<bincenter<<" "<<TMath::Max(val_1,val_2)<<endl;
  }
  pp_5020_ct10n_nnpdf_03_03_R3.close();


  //03_07
  ofstream pp_5020_ct10n_nnpdf_03_07_R3;
  pp_5020_ct10n_nnpdf_03_07_R3.open("pp_5020_ct10n_nnpdf_03_07_R3.txt");
  for(int i = 0;i<hPP_5020_ct10n_nnpdf_03_07_R3->GetNbinsX();i++){
    float bincenter = hPP_5020_ct10n_nnpdf_03_07_R3->GetBinCenter(i);
    float val_1 = hPP_5020_ct10n_nnpdf_03_07_R3->GetBinContent(i);
    val_1 = 100*TMath::Abs(1-val_1);
    float val_2 = hPP_5020_ct10n_hera_03_07_R3->GetBinContent(i);
    val_2 = 100*TMath::Abs(1-val_2);
    pp_5020_ct10n_nnpdf_03_07_R3<<bincenter<<" "<<TMath::Max(val_1,val_2)<<endl;
  }
  pp_5020_ct10n_nnpdf_03_07_R3.close();



  //07_10
  ofstream pp_5020_ct10n_nnpdf_07_10_R3;
  pp_5020_ct10n_nnpdf_07_10_R3.open("pp_5020_ct10n_nnpdf_07_10_R3.txt");
  for(int i = 0;i<hPP_5020_ct10n_nnpdf_07_10_R3->GetNbinsX();i++){
    float bincenter = hPP_5020_ct10n_nnpdf_07_10_R3->GetBinCenter(i);
    float val_1 = hPP_5020_ct10n_nnpdf_07_10_R3->GetBinContent(i);
    val_1 = 100*TMath::Abs(1-val_1);
    float val_2 = hPP_5020_ct10n_hera_07_10_R3->GetBinContent(i);
    val_2 = 100*TMath::Abs(1-val_2);
    pp_5020_ct10n_nnpdf_07_10_R3<<bincenter<<" "<<TMath::Max(val_1,val_2)<<endl;
  }
  pp_5020_ct10n_nnpdf_07_10_R3.close();

 

 //10_12
  ofstream pp_5020_ct10n_nnpdf_10_12_R3;
  pp_5020_ct10n_nnpdf_10_12_R3.open("pp_5020_ct10n_nnpdf_10_12_R3.txt");
  for(int i = 0;i<hPP_5020_ct10n_nnpdf_10_12_R3->GetNbinsX();i++){
    float bincenter = hPP_5020_ct10n_nnpdf_10_12_R3->GetBinCenter(i);
    float val_1 = hPP_5020_ct10n_nnpdf_10_12_R3->GetBinContent(i);
    val_1 = 100*TMath::Abs(1-val_1);
    float val_2 = hPP_5020_ct10n_hera_10_12_R3->GetBinContent(i);
    val_2 = 100*TMath::Abs(1-val_2);
    pp_5020_ct10n_nnpdf_10_12_R3<<bincenter<<" "<<TMath::Max(val_1,val_2)<<endl;
  }
  pp_5020_ct10n_nnpdf_10_12_R3.close();

   //12_22
  ofstream pp_5020_ct10n_nnpdf_12_22_R3;
  pp_5020_ct10n_nnpdf_12_22_R3.open("pp_5020_ct10n_nnpdf_12_22_R3.txt");
  for(int i = 0;i<hPP_5020_ct10n_nnpdf_12_22_R3->GetNbinsX();i++){
    float bincenter = hPP_5020_ct10n_nnpdf_12_22_R3->GetBinCenter(i);
    float val_1 = hPP_5020_ct10n_nnpdf_12_22_R3->GetBinContent(i);
    val_1 = 100*TMath::Abs(1-val_1);
    float val_2 = hPP_5020_ct10n_hera_12_22_R3->GetBinContent(i);
    val_2 = 100*TMath::Abs(1-val_2);
    pp_5020_ct10n_nnpdf_12_22_R3<<bincenter<<" "<<TMath::Max(val_1,val_2)<<endl;
  }
  pp_5020_ct10n_nnpdf_12_22_R3.close();


   //22_22
  ofstream pp_5020_ct10n_nnpdf_22_22_R3;
  pp_5020_ct10n_nnpdf_22_22_R3.open("pp_5020_ct10n_nnpdf_22_22_R3.txt");
  for(int i = 0;i<hPP_5020_ct10n_nnpdf_22_22_R3->GetNbinsX();i++){
    float bincenter = hPP_5020_ct10n_nnpdf_22_22_R3->GetBinCenter(i);
    float val_1 = hPP_5020_ct10n_nnpdf_22_22_R3->GetBinContent(i);
    val_1 = 100*TMath::Abs(1-val_1);
    float val_2 = hPP_5020_ct10n_hera_22_22_R3->GetBinContent(i);
    val_2 = 100*TMath::Abs(1-val_2);
    pp_5020_ct10n_nnpdf_22_22_R3<<bincenter<<" "<<TMath::Max(val_1,val_2)<<endl;
  }
  pp_5020_ct10n_nnpdf_22_22_R3.close();


   //10_10
  ofstream pp_5020_ct10n_nnpdf_10_10_R3;
  pp_5020_ct10n_nnpdf_10_10_R3.open("pp_5020_ct10n_nnpdf_10_10_R3.txt");
  for(int i = 0;i<hPP_5020_ct10n_nnpdf_10_10_R3->GetNbinsX();i++){
    float bincenter = hPP_5020_ct10n_nnpdf_10_10_R3->GetBinCenter(i);
    float val_1 = hPP_5020_ct10n_nnpdf_10_10_R3->GetBinContent(i);
    val_1 = 100*TMath::Abs(1-val_1);
    float val_2 = hPP_5020_ct10n_hera_10_10_R3->GetBinContent(i);
    val_2 = 100*TMath::Abs(1-val_2);
    pp_5020_ct10n_nnpdf_10_10_R3<<bincenter<<" "<<TMath::Max(val_1,val_2)<<endl;
  }
  pp_5020_ct10n_nnpdf_10_10_R3.close();



  //2.76 TeV
  TDirectoryFile* ak3GenJet_2760_05 = (TDirectoryFile*)fppPythia_2760->Get("ak3GenJetSpectrum_QCD10001_00_05");
  TH1F* hPP_2760_Pythia_05_R3 = (TH1F*)ak3GenJet_2760_05->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_2760_05_10 = (TDirectoryFile*)fppPythia_2760->Get("ak3GenJetSpectrum_QCD10001_05_10");
  TH1F* hPP_2760_Pythia_05_10_R3 = (TH1F*)ak3GenJet_2760_05_10->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_2760_10_15 = (TDirectoryFile*)fppPythia_2760->Get("ak3GenJetSpectrum_QCD10001_10_15");
  TH1F* hPP_2760_Pythia_10_15_R3 = (TH1F*)ak3GenJet_2760_10_15->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_2760_15_20 = (TDirectoryFile*)fppPythia_2760->Get("ak3GenJetSpectrum_QCD10001_15_20");
  TH1F* hPP_2760_Pythia_15_20_R3 = (TH1F*)ak3GenJet_2760_15_20->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_2760_20_25 = (TDirectoryFile*)fppPythia_2760->Get("ak3GenJetSpectrum_QCD10001_20_25");
  TH1F* hPP_2760_Pythia_20_25_R3 = (TH1F*)ak3GenJet_2760_20_25->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_2760_25_30 = (TDirectoryFile*)fppPythia_2760->Get("ak3GenJetSpectrum_QCD10001_25_30");
  TH1F* hPP_2760_Pythia_25_30_R3 = (TH1F*)ak3GenJet_2760_25_30->Get("JetSpectrum");

  TDirectoryFile* ak4GenJet_2760_05 = (TDirectoryFile*)fppPythia_2760->Get("ak4GenJetSpectrum_QCD10001_00_05");
  TH1F* hPP_2760_Pythia_05_R4 = (TH1F*)ak4GenJet_2760_05->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_2760_05_10 = (TDirectoryFile*)fppPythia_2760->Get("ak4GenJetSpectrum_QCD10001_05_10");
  TH1F* hPP_2760_Pythia_05_10_R4 = (TH1F*)ak4GenJet_2760_05_10->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_2760_10_15 = (TDirectoryFile*)fppPythia_2760->Get("ak4GenJetSpectrum_QCD10001_10_15");
  TH1F* hPP_2760_Pythia_10_15_R4 = (TH1F*)ak4GenJet_2760_10_15->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_2760_15_20 = (TDirectoryFile*)fppPythia_2760->Get("ak4GenJetSpectrum_QCD10001_15_20");
  TH1F* hPP_2760_Pythia_15_20_R4 = (TH1F*)ak4GenJet_2760_15_20->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_2760_20_25 = (TDirectoryFile*)fppPythia_2760->Get("ak4GenJetSpectrum_QCD10001_20_25");
  TH1F* hPP_2760_Pythia_20_25_R4 = (TH1F*)ak4GenJet_2760_20_25->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_2760_25_30 = (TDirectoryFile*)fppPythia_2760->Get("ak4GenJetSpectrum_QCD10001_25_30");
  TH1F* hPP_2760_Pythia_25_30_R4 = (TH1F*)ak4GenJet_2760_25_30->Get("JetSpectrum");

  TDirectoryFile* ak5GenJet_2760_05 = (TDirectoryFile*)fppPythia_2760->Get("ak5GenJetSpectrum_QCD10001_00_05");
  TH1F* hPP_2760_Pythia_05_R5 = (TH1F*)ak5GenJet_2760_05->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_2760_05_10 = (TDirectoryFile*)fppPythia_2760->Get("ak5GenJetSpectrum_QCD10001_05_10");
  TH1F* hPP_2760_Pythia_05_10_R5 = (TH1F*)ak5GenJet_2760_05_10->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_2760_10_15 = (TDirectoryFile*)fppPythia_2760->Get("ak5GenJetSpectrum_QCD10001_10_15");
  TH1F* hPP_2760_Pythia_10_15_R5 = (TH1F*)ak5GenJet_2760_10_15->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_2760_15_20 = (TDirectoryFile*)fppPythia_2760->Get("ak5GenJetSpectrum_QCD10001_15_20");
  TH1F* hPP_2760_Pythia_15_20_R5 = (TH1F*)ak5GenJet_2760_15_20->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_2760_20_25 = (TDirectoryFile*)fppPythia_2760->Get("ak5GenJetSpectrum_QCD10001_20_25");
  TH1F* hPP_2760_Pythia_20_25_R5 = (TH1F*)ak5GenJet_2760_20_25->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_2760_25_30 = (TDirectoryFile*)fppPythia_2760->Get("ak5GenJetSpectrum_QCD10001_25_30");
  TH1F* hPP_2760_Pythia_25_30_R5 = (TH1F*)ak5GenJet_2760_25_30->Get("JetSpectrum");


  //5.02 TeV
  TDirectoryFile* ak3GenJet_5020_05 = (TDirectoryFile*)fppPythia_5020->Get("ak3GenJetSpectrum_QCD10001_00_05");
  TH1F* hPP_5020_Pythia_05_R3 = (TH1F*)ak3GenJet_5020_05->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_5020_05_10 = (TDirectoryFile*)fppPythia_5020->Get("ak3GenJetSpectrum_QCD10001_05_10");
  TH1F* hPP_5020_Pythia_05_10_R3 = (TH1F*)ak3GenJet_5020_05_10->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_5020_10_15 = (TDirectoryFile*)fppPythia_5020->Get("ak3GenJetSpectrum_QCD10001_10_15");
  TH1F* hPP_5020_Pythia_10_15_R3 = (TH1F*)ak3GenJet_5020_10_15->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_5020_15_20 = (TDirectoryFile*)fppPythia_5020->Get("ak3GenJetSpectrum_QCD10001_15_20");
  TH1F* hPP_5020_Pythia_15_20_R3 = (TH1F*)ak3GenJet_5020_15_20->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_5020_20_25 = (TDirectoryFile*)fppPythia_5020->Get("ak3GenJetSpectrum_QCD10001_20_25");
  TH1F* hPP_5020_Pythia_20_25_R3 = (TH1F*)ak3GenJet_5020_20_25->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_5020_25_30 = (TDirectoryFile*)fppPythia_5020->Get("ak3GenJetSpectrum_QCD10001_25_30");
  TH1F* hPP_5020_Pythia_25_30_R3 = (TH1F*)ak3GenJet_5020_25_30->Get("JetSpectrum");

  TDirectoryFile* ak4GenJet_5020_05 = (TDirectoryFile*)fppPythia_5020->Get("ak4GenJetSpectrum_QCD10001_00_05");
  TH1F* hPP_5020_Pythia_05_R4 = (TH1F*)ak4GenJet_5020_05->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_5020_05_10 = (TDirectoryFile*)fppPythia_5020->Get("ak4GenJetSpectrum_QCD10001_05_10");
  TH1F* hPP_5020_Pythia_05_10_R4 = (TH1F*)ak4GenJet_5020_05_10->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_5020_10_15 = (TDirectoryFile*)fppPythia_5020->Get("ak4GenJetSpectrum_QCD10001_10_15");
  TH1F* hPP_5020_Pythia_10_15_R4 = (TH1F*)ak4GenJet_5020_10_15->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_5020_15_20 = (TDirectoryFile*)fppPythia_5020->Get("ak4GenJetSpectrum_QCD10001_15_20");
  TH1F* hPP_5020_Pythia_15_20_R4 = (TH1F*)ak4GenJet_5020_15_20->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_5020_20_25 = (TDirectoryFile*)fppPythia_5020->Get("ak4GenJetSpectrum_QCD10001_20_25");
  TH1F* hPP_5020_Pythia_20_25_R4 = (TH1F*)ak4GenJet_5020_20_25->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_5020_25_30 = (TDirectoryFile*)fppPythia_5020->Get("ak4GenJetSpectrum_QCD10001_25_30");
  TH1F* hPP_5020_Pythia_25_30_R4 = (TH1F*)ak4GenJet_5020_25_30->Get("JetSpectrum");

  TDirectoryFile* ak5GenJet_5020_05 = (TDirectoryFile*)fppPythia_5020->Get("ak5GenJetSpectrum_QCD10001_00_05");
  TH1F* hPP_5020_Pythia_05_R5 = (TH1F*)ak5GenJet_5020_05->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_5020_05_10 = (TDirectoryFile*)fppPythia_5020->Get("ak5GenJetSpectrum_QCD10001_05_10");
  TH1F* hPP_5020_Pythia_05_10_R5 = (TH1F*)ak5GenJet_5020_05_10->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_5020_10_15 = (TDirectoryFile*)fppPythia_5020->Get("ak5GenJetSpectrum_QCD10001_10_15");
  TH1F* hPP_5020_Pythia_10_15_R5 = (TH1F*)ak5GenJet_5020_10_15->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_5020_15_20 = (TDirectoryFile*)fppPythia_5020->Get("ak5GenJetSpectrum_QCD10001_15_20");
  TH1F* hPP_5020_Pythia_15_20_R5 = (TH1F*)ak5GenJet_5020_15_20->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_5020_20_25 = (TDirectoryFile*)fppPythia_5020->Get("ak5GenJetSpectrum_QCD10001_20_25");
  TH1F* hPP_5020_Pythia_20_25_R5 = (TH1F*)ak5GenJet_5020_20_25->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_5020_25_30 = (TDirectoryFile*)fppPythia_5020->Get("ak5GenJetSpectrum_QCD10001_25_30");
  TH1F* hPP_5020_Pythia_25_30_R5 = (TH1F*)ak5GenJet_5020_25_30->Get("JetSpectrum");


  //7 TeV
  TDirectoryFile* ak3GenJet_7000_05 = (TDirectoryFile*)fppPythia_7000->Get("ak3GenJetSpectrum_QCD10001_00_05");
  TH1F* hPP_7000_Pythia_05_R3 = (TH1F*)ak3GenJet_7000_05->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_7000_05_10 = (TDirectoryFile*)fppPythia_7000->Get("ak3GenJetSpectrum_QCD10001_05_10");
  TH1F* hPP_7000_Pythia_05_10_R3 = (TH1F*)ak3GenJet_7000_05_10->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_7000_10_15 = (TDirectoryFile*)fppPythia_7000->Get("ak3GenJetSpectrum_QCD10001_10_15");
  TH1F* hPP_7000_Pythia_10_15_R3 = (TH1F*)ak3GenJet_7000_10_15->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_7000_15_20 = (TDirectoryFile*)fppPythia_7000->Get("ak3GenJetSpectrum_QCD10001_15_20");
  TH1F* hPP_7000_Pythia_15_20_R3 = (TH1F*)ak3GenJet_7000_15_20->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_7000_20_25 = (TDirectoryFile*)fppPythia_7000->Get("ak3GenJetSpectrum_QCD10001_20_25");
  TH1F* hPP_7000_Pythia_20_25_R3 = (TH1F*)ak3GenJet_7000_20_25->Get("JetSpectrum");
  TDirectoryFile* ak3GenJet_7000_25_30 = (TDirectoryFile*)fppPythia_7000->Get("ak3GenJetSpectrum_QCD10001_25_30");
  TH1F* hPP_7000_Pythia_25_30_R3 = (TH1F*)ak3GenJet_7000_25_30->Get("JetSpectrum");

  TDirectoryFile* ak4GenJet_7000_05 = (TDirectoryFile*)fppPythia_7000->Get("ak4GenJetSpectrum_QCD10001_00_05");
  TH1F* hPP_7000_Pythia_05_R4 = (TH1F*)ak4GenJet_7000_05->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_7000_05_10 = (TDirectoryFile*)fppPythia_7000->Get("ak4GenJetSpectrum_QCD10001_05_10");
  TH1F* hPP_7000_Pythia_05_10_R4 = (TH1F*)ak4GenJet_7000_05_10->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_7000_10_15 = (TDirectoryFile*)fppPythia_7000->Get("ak4GenJetSpectrum_QCD10001_10_15");
  TH1F* hPP_7000_Pythia_10_15_R4 = (TH1F*)ak4GenJet_7000_10_15->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_7000_15_20 = (TDirectoryFile*)fppPythia_7000->Get("ak4GenJetSpectrum_QCD10001_15_20");
  TH1F* hPP_7000_Pythia_15_20_R4 = (TH1F*)ak4GenJet_7000_15_20->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_7000_20_25 = (TDirectoryFile*)fppPythia_7000->Get("ak4GenJetSpectrum_QCD10001_20_25");
  TH1F* hPP_7000_Pythia_20_25_R4 = (TH1F*)ak4GenJet_7000_20_25->Get("JetSpectrum");
  TDirectoryFile* ak4GenJet_7000_25_30 = (TDirectoryFile*)fppPythia_7000->Get("ak4GenJetSpectrum_QCD10001_25_30");
  TH1F* hPP_7000_Pythia_25_30_R4 = (TH1F*)ak4GenJet_7000_25_30->Get("JetSpectrum");

  TDirectoryFile* ak5GenJet_7000_05 = (TDirectoryFile*)fppPythia_7000->Get("ak5GenJetSpectrum_QCD10001_00_05");
  TH1F* hPP_7000_Pythia_05_R5 = (TH1F*)ak5GenJet_7000_05->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_7000_05_10 = (TDirectoryFile*)fppPythia_7000->Get("ak5GenJetSpectrum_QCD10001_05_10");
  TH1F* hPP_7000_Pythia_05_10_R5 = (TH1F*)ak5GenJet_7000_05_10->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_7000_10_15 = (TDirectoryFile*)fppPythia_7000->Get("ak5GenJetSpectrum_QCD10001_10_15");
  TH1F* hPP_7000_Pythia_10_15_R5 = (TH1F*)ak5GenJet_7000_10_15->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_7000_15_20 = (TDirectoryFile*)fppPythia_7000->Get("ak5GenJetSpectrum_QCD10001_15_20");
  TH1F* hPP_7000_Pythia_15_20_R5 = (TH1F*)ak5GenJet_7000_15_20->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_7000_20_25 = (TDirectoryFile*)fppPythia_7000->Get("ak5GenJetSpectrum_QCD10001_20_25");
  TH1F* hPP_7000_Pythia_20_25_R5 = (TH1F*)ak5GenJet_7000_20_25->Get("JetSpectrum");
  TDirectoryFile* ak5GenJet_7000_25_30 = (TDirectoryFile*)fppPythia_7000->Get("ak5GenJetSpectrum_QCD10001_25_30");
  TH1F* hPP_7000_Pythia_25_30_R5 = (TH1F*)ak5GenJet_7000_25_30->Get("JetSpectrum");


  //ok we got the histograms from the file. now lets get the histograms needed from these histograms. 
  //these are all in differnet bins, so no way we can add them. for lets just do it with -0.5 to +0.5 for diff radius and diff energy 
  /*
  TH1F* hPP_2760_Pythia_20_20_R3 = (TH1F*)hPP_2760_Pythia_05_R3->Clone("hPP_2760_Pythia_20_20_R3");
  hPP_2760_Pythia_20_20_R3->Add(hPP_2760_Pythia_05_10_R3);
  hPP_2760_Pythia_20_20_R3->Add(hPP_2760_Pythia_10_15_R3);
  hPP_2760_Pythia_20_20_R3->Add(hPP_2760_Pythia_15_20_R3);
  hPP_2760_Pythia_20_20_R3->Divide(1./4);

  TH1F* hPP_2760_Pythia_20_20_R4 = (TH1F*)hPP_2760_Pythia_05_R4->Clone("hPP_2760_Pythia_20_20_R4");
  hPP_2760_Pythia_20_20_R4->Add(hPP_2760_Pythia_05_10_R4);
  hPP_2760_Pythia_20_20_R4->Add(hPP_2760_Pythia_10_15_R4);
  hPP_2760_Pythia_20_20_R4->Add(hPP_2760_Pythia_15_20_R4);
  hPP_2760_Pythia_20_20_R4->Divide(1./4);

  TH1F* hPP_2760_Pythia_20_20_R4 = (TH1F*)hPP_2760_Pythia_05_R4->Clone("hPP_2760_Pythia_20_20_R4");
  hPP_2760_Pythia_20_20_R4->Add(hPP_2760_Pythia_05_10_R4);
  hPP_2760_Pythia_20_20_R4->Add(hPP_2760_Pythia_10_15_R4);
  hPP_2760_Pythia_20_20_R4->Add(hPP_2760_Pythia_15_20_R4);
  hPP_2760_Pythia_20_20_R4->Divide(1./4);
  */

  TH1F* hPP_2760_Pythia_05_R_3_5 = (TH1F*)hPP_2760_Pythia_05_R3->Clone("hPP_2760_Pythia_05_R_3_5");
  hPP_2760_Pythia_05_R_3_5->Divide(hPP_2760_Pythia_05_R5);

  TH1F* hPP_2760_Pythia_05_R_4_5 = (TH1F*)hPP_2760_Pythia_05_R4->Clone("hPP_2760_Pythia_05_R_4_5");
  hPP_2760_Pythia_05_R_4_5->Divide(hPP_2760_Pythia_05_R5);

  TH1F* hPP_2760_Pythia_05_10_R_3_5 = (TH1F*)hPP_2760_Pythia_05_10_R3->Clone("hPP_2760_Pythia_05_10_R_3_5");
  hPP_2760_Pythia_05_10_R_3_5->Divide(hPP_2760_Pythia_05_10_R5);

  TH1F* hPP_2760_Pythia_05_10_R_4_5 = (TH1F*)hPP_2760_Pythia_05_10_R4->Clone("hPP_2760_Pythia_05_10_R_4_5");
  hPP_2760_Pythia_05_10_R_4_5->Divide(hPP_2760_Pythia_05_10_R5);

  TH1F* hPP_2760_Pythia_15_20_R_3_5 = (TH1F*)hPP_2760_Pythia_15_20_R3->Clone("hPP_2760_Pythia_15_20_R_3_5");
  hPP_2760_Pythia_15_20_R_3_5->Divide(hPP_2760_Pythia_15_20_R5);

  TH1F* hPP_2760_Pythia_15_20_R_4_5 = (TH1F*)hPP_2760_Pythia_15_20_R4->Clone("hPP_2760_Pythia_15_20_R_4_5");
  hPP_2760_Pythia_15_20_R_4_5->Divide(hPP_2760_Pythia_15_20_R5);

  TH1F* hPP_5020_Pythia_05_R_3_5 = (TH1F*)hPP_5020_Pythia_05_R3->Clone("hPP_5020_Pythia_05_R_3_5");
  hPP_5020_Pythia_05_R_3_5->Divide(hPP_5020_Pythia_05_R5);

  TH1F* hPP_5020_Pythia_05_R_4_5 = (TH1F*)hPP_5020_Pythia_05_R4->Clone("hPP_5020_Pythia_05_R_4_5");
  hPP_5020_Pythia_05_R_4_5->Divide(hPP_5020_Pythia_05_R5);

  TH1F* hPP_5020_Pythia_05_10_R_3_5 = (TH1F*)hPP_5020_Pythia_05_10_R3->Clone("hPP_5020_Pythia_05_10_R_3_5");
  hPP_5020_Pythia_05_10_R_3_5->Divide(hPP_5020_Pythia_05_10_R5);

  TH1F* hPP_5020_Pythia_05_10_R_4_5 = (TH1F*)hPP_5020_Pythia_05_10_R4->Clone("hPP_5020_Pythia_05_10_R_4_5");
  hPP_5020_Pythia_05_10_R_4_5->Divide(hPP_5020_Pythia_05_10_R5);

  TH1F* hPP_5020_Pythia_15_20_R_3_5 = (TH1F*)hPP_5020_Pythia_15_20_R3->Clone("hPP_5020_Pythia_15_20_R_3_5");
  hPP_5020_Pythia_15_20_R_3_5->Divide(hPP_5020_Pythia_15_20_R5);

  TH1F* hPP_5020_Pythia_15_20_R_4_5 = (TH1F*)hPP_5020_Pythia_15_20_R4->Clone("hPP_5020_Pythia_15_20_R_4_5");
  hPP_5020_Pythia_15_20_R_4_5->Divide(hPP_5020_Pythia_15_20_R5);

  TH1F* hPP_7000_Pythia_05_R_3_5 = (TH1F*)hPP_7000_Pythia_05_R3->Clone("hPP_7000_Pythia_05_R_3_5");
  hPP_7000_Pythia_05_R_3_5->Divide(hPP_7000_Pythia_05_R5);

  TH1F* hPP_7000_Pythia_05_R_4_5 = (TH1F*)hPP_7000_Pythia_05_R4->Clone("hPP_7000_Pythia_05_R_4_5");
  hPP_7000_Pythia_05_R_4_5->Divide(hPP_7000_Pythia_05_R5);

  TH1F* hPP_7000_Pythia_05_10_R_3_5 = (TH1F*)hPP_7000_Pythia_05_10_R3->Clone("hPP_7000_Pythia_05_10_R_3_5");
  hPP_7000_Pythia_05_10_R_3_5->Divide(hPP_7000_Pythia_05_10_R5);

  TH1F* hPP_7000_Pythia_05_10_R_4_5 = (TH1F*)hPP_7000_Pythia_05_10_R4->Clone("hPP_7000_Pythia_05_10_R_4_5");
  hPP_7000_Pythia_05_10_R_4_5->Divide(hPP_7000_Pythia_05_10_R5);

  TH1F* hPP_7000_Pythia_15_20_R_3_5 = (TH1F*)hPP_7000_Pythia_15_20_R3->Clone("hPP_7000_Pythia_15_20_R_3_5");
  hPP_7000_Pythia_15_20_R_3_5->Divide(hPP_7000_Pythia_15_20_R5);

  TH1F* hPP_7000_Pythia_15_20_R_4_5 = (TH1F*)hPP_7000_Pythia_15_20_R4->Clone("hPP_7000_Pythia_15_20_R_4_5");
  hPP_7000_Pythia_15_20_R_4_5->Divide(hPP_7000_Pythia_15_20_R5);

  TH1F* hPP_2760_5020_Pythia_05_R3 = (TH1F*)hPP_2760_Pythia_05_R3->Clone("hPP_2760_5020_Pythia_05_R3");
  hPP_2760_5020_Pythia_05_R3->Divide(hPP_5020_Pythia_05_R3);

  TH1F* hPP_2760_7000_Pythia_05_R3 = (TH1F*)hPP_2760_Pythia_05_R3->Clone("hPP_2760_7000_Pythia_05_R3");
  hPP_2760_7000_Pythia_05_R3->Divide(hPP_7000_Pythia_05_R3);

  TH1F* hPP_2760_5020_Pythia_05_R4 = (TH1F*)hPP_2760_Pythia_05_R4->Clone("hPP_2760_5020_Pythia_05_R4");
  hPP_2760_5020_Pythia_05_R4->Divide(hPP_5020_Pythia_05_R4);

  TH1F* hPP_2760_7000_Pythia_05_R4 = (TH1F*)hPP_2760_Pythia_05_R4->Clone("hPP_2760_7000_Pythia_05_R4");
  hPP_2760_7000_Pythia_05_R4->Divide(hPP_7000_Pythia_05_R4);

  TH1F* hPP_2760_5020_Pythia_05_R5 = (TH1F*)hPP_2760_Pythia_05_R5->Clone("hPP_2760_5020_Pythia_05_R5");
  hPP_2760_5020_Pythia_05_R5->Divide(hPP_5020_Pythia_05_R5);

  TH1F* hPP_2760_7000_Pythia_05_R5 = (TH1F*)hPP_2760_Pythia_05_R5->Clone("hPP_2760_7000_Pythia_05_R5");
  hPP_2760_7000_Pythia_05_R5->Divide(hPP_7000_Pythia_05_R5);
  
  TH1F* hPP_5020_2760_nnpdf_R2 = (TH1F*)hPP_5020_nnpdf_22_22_R2->Clone("hPP_5020_2760_nnpdf_R2");
  hPP_5020_2760_nnpdf_R2->Divide(hPP_2760_nnpdf_2_2_R2);
  TH1F* hPP_5020_2760_nnpdf_R3 = (TH1F*)hPP_5020_nnpdf_22_22_R3->Clone("hPP_5020_2760_nnpdf_R2");
  hPP_5020_2760_nnpdf_R3->Divide(hPP_2760_nnpdf_2_2_R3);
  TH1F* hPP_5020_2760_nnpdf_R4 = (TH1F*)hPP_5020_nnpdf_22_22_R4->Clone("hPP_5020_2760_nnpdf_R2");
  hPP_5020_2760_nnpdf_R4->Divide(hPP_2760_nnpdf_2_2_R4);

  TH1F* hPP_5020_2760_ct10n_R2 = (TH1F*)hPP_5020_ct10n_22_22_R2->Clone("hPP_5020_2760_ct10n_R2");
  hPP_5020_2760_ct10n_R2->Divide(hPP_2760_ct10n_2_2_R2);
  TH1F* hPP_5020_2760_ct10n_R3 = (TH1F*)hPP_5020_ct10n_22_22_R3->Clone("hPP_5020_2760_ct10n_R2");
  hPP_5020_2760_ct10n_R3->Divide(hPP_2760_ct10n_2_2_R3);
  TH1F* hPP_5020_2760_ct10n_R4 = (TH1F*)hPP_5020_ct10n_22_22_R4->Clone("hPP_5020_2760_ct10n_R2");
  hPP_5020_2760_ct10n_R4->Divide(hPP_2760_ct10n_2_2_R4);

  TH1F* hPP_5020_2760_hera_R2 = (TH1F*)hPP_5020_hera_22_22_R2->Clone("hPP_5020_2760_hera_R2");
  hPP_5020_2760_hera_R2->Divide(hPP_2760_hera_2_2_R2);
  TH1F* hPP_5020_2760_hera_R3 = (TH1F*)hPP_5020_hera_22_22_R3->Clone("hPP_5020_2760_hera_R2");
  hPP_5020_2760_hera_R3->Divide(hPP_2760_hera_2_2_R3);
  TH1F* hPP_5020_2760_hera_R4 = (TH1F*)hPP_5020_hera_22_22_R4->Clone("hPP_5020_2760_hera_R2");
  hPP_5020_2760_hera_R4->Divide(hPP_2760_hera_2_2_R4);

  TH1F* hPP_5020_nnpdf_R_2_4 = (TH1F*)hPP_5020_nnpdf_22_22_R2->Clone("hPP_5020_nnpdf_R_2_4");
  hPP_5020_nnpdf_R_2_4->Divide(hPP_5020_nnpdf_22_22_R4);
  TH1F* hPP_5020_nnpdf_R_3_4 = (TH1F*)hPP_5020_nnpdf_22_22_R3->Clone("hPP_5020_nnpdf_R_3_4");
  hPP_5020_nnpdf_R_3_4->Divide(hPP_5020_nnpdf_22_22_R4);

  TH1F* hPP_2760_nnpdf_R_2_4 = (TH1F*)fpp276->Get("hRatio_nnpdf_R_2_4");
  TH1F* hPP_2760_nnpdf_R_3_4 = (TH1F*)fpp276->Get("hRatio_nnpdf_R_3_4");

  TH1F* hPP_data_2760_R_3_4 = (TH1F*)fpp276->Get("hRatio_data_R_3_4");
  TH1F* hPP_data_2760_R_3_5 = (TH1F*)fpp276->Get("hRatio_data_R_3_5");
  TH1F* hPP_data_2760_R_4_5 = (TH1F*)fpp276->Get("hRatio_data_R_4_5");

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetGridy();
  hPP_5020_nnpdf_R_2_4->SetMarkerStyle(20);
  hPP_5020_nnpdf_R_2_4->SetMarkerColor(3);
  hPP_5020_nnpdf_R_2_4->SetYTitle("Ratio of differential cross sections");
  hPP_5020_nnpdf_R_2_4->SetXTitle("Jet p_{T} (GeV/c)");
  hPP_5020_nnpdf_R_2_4->SetTitle(" ");
  hPP_5020_nnpdf_R_2_4->SetAxisRange(22,500,"X");
  hPP_5020_nnpdf_R_2_4->SetAxisRange(0.75,1.1,"Y");
  hPP_5020_nnpdf_R_2_4->Draw("p");
  hPP_5020_nnpdf_R_3_4->SetMarkerStyle(22);
  hPP_5020_nnpdf_R_3_4->SetMarkerColor(4);
  hPP_5020_nnpdf_R_3_4->Draw("same p");

  hPP_2760_nnpdf_R_2_4->SetMarkerStyle(20);
  hPP_2760_nnpdf_R_2_4->SetMarkerColor(6);
  hPP_2760_nnpdf_R_2_4->Draw("same p");
  hPP_2760_nnpdf_R_3_4->SetMarkerStyle(22);
  hPP_2760_nnpdf_R_3_4->SetMarkerColor(7);
  hPP_2760_nnpdf_R_3_4->Draw("same p");

  hPP_data_2760_R_3_4->SetMarkerStyle(22);
  hPP_data_2760_R_3_4->SetMarkerColor(8);
  hPP_data_2760_R_3_4->Draw("same p");

  TLegend *title3 = myLegend(0.13,0.55,0.33,0.85);
  title3->AddEntry(hPP_5020_nnpdf_R_2_4,"NLO 5.02, R=0.2/R=0.4","pl");
  title3->AddEntry(hPP_5020_nnpdf_R_3_4,"NLO 5.02, R=0.3/R=0.4","pl");
  title3->AddEntry(hPP_2760_nnpdf_R_2_4,"NLO 2.76, R=0.2/R=0.4","pl");
  title3->AddEntry(hPP_2760_nnpdf_R_3_4,"NLO 2.76, R=0.3/R=0.4","pl");
  title3->AddEntry(hPP_data_2760_R_3_4,"Data 2.76 R=0.3/R=0.4","pl");
  title3->SetTextSize(0.04);
  title3->Draw();

  putCMSPrel(0.1,0.92,0.04);
  drawText("pp NLO-NNPDF21 , #sqrt{s}=2.76(TeV) and 5.02(TeV)",0.35,0.92,16);
  //drawText(Form("anti k_{T} Jets",2),0.47,0.83,16);

  c2->SaveAs("RpA_reference_NLO_data_diff_energy_radius_cross_section_ratio.pdf","RECREATE");


  TCanvas *c1 = new TCanvas("c1","",800,600);
  formatCanvas(c1);
  c1->cd(1);

  hPP_5020_nnpdf_22_22_R2->SetMarkerStyle(20);
  hPP_5020_nnpdf_22_22_R2->SetMarkerColor(3);
  hPP_5020_nnpdf_22_22_R2->SetTitle(" ");
  hPP_5020_nnpdf_22_22_R2->SetXTitle("Jet p_{T} (GeV/c)");
  hPP_5020_nnpdf_22_22_R2->SetYTitle("#frac{d^{2} #sigma}{d p_{T} d #eta} (pb#frac{GeV}{c})");
  hPP_5020_nnpdf_22_22_R2->Draw("p");
  hPP_5020_nnpdf_22_22_R3->SetMarkerStyle(20);
  hPP_5020_nnpdf_22_22_R3->SetMarkerColor(4);
  hPP_5020_nnpdf_22_22_R3->Draw("same p");
  hPP_5020_nnpdf_22_22_R4->SetMarkerStyle(20);
  hPP_5020_nnpdf_22_22_R4->SetMarkerColor(6);
  hPP_5020_nnpdf_22_22_R4->Draw("same");

  hPP_2760_nnpdf_2_2_R2->SetMarkerStyle(33);
  hPP_2760_nnpdf_2_2_R2->SetMarkerColor(3);
  hPP_2760_nnpdf_2_2_R2->Draw("same p");
  hPP_2760_nnpdf_2_2_R3->SetMarkerStyle(33);
  hPP_2760_nnpdf_2_2_R3->SetMarkerColor(4);
  hPP_2760_nnpdf_2_2_R3->Draw("same p");
  hPP_2760_nnpdf_2_2_R4->SetMarkerStyle(33);
  hPP_2760_nnpdf_2_2_R4->SetMarkerColor(6);
  hPP_2760_nnpdf_2_2_R3->Draw("same p");

  putCMSPrel(0.1,0.92,0.06);
  drawText("pp NLO-NNPDF21 , #sqrt{s}=2.76(TeV) and 5.02(TeV)",0.35,0.92,16);
  drawText(Form("anti k_{T} Jets",2),0.47,0.83,16);
  
  TLegend *title = myLegend(0.47,0.50,0.67,0.8);
  title->AddEntry(hPP_5020_nnpdf_22_22_R2,"5.02, R=0.2","pl");
  title->AddEntry(hPP_5020_nnpdf_22_22_R3,"5.02, R=0.3","pl");
  title->AddEntry(hPP_5020_nnpdf_22_22_R4,"5.02, R=0.4","pl");
  title->AddEntry(hPP_2760_nnpdf_2_2_R2,"2.76, R=0.2","pl");
  title->AddEntry(hPP_2760_nnpdf_2_2_R3,"2.76, R=0.3","pl");
  title->AddEntry(hPP_2760_nnpdf_2_2_R4,"2.76, R=0.4","pl");
  title->SetTextSize(0.04);
  title->Draw();

  c1->cd(2);
  c1->cd(2)->SetLogy();

  hPP_5020_2760_nnpdf_R2->SetMarkerStyle(20);
  hPP_5020_2760_nnpdf_R2->SetMarkerColor(3);
  hPP_5020_2760_nnpdf_R2->SetTitle(" ");
  hPP_5020_2760_nnpdf_R2->SetXTitle("Jet p_{T} (GeV/c)");
  hPP_5020_2760_nnpdf_R2->SetYTitle("Ratio");
  hPP_5020_2760_nnpdf_R2->Draw("p");
  hPP_5020_2760_nnpdf_R3->SetMarkerStyle(20);
  hPP_5020_2760_nnpdf_R3->SetMarkerColor(4);
  hPP_5020_2760_nnpdf_R3->Draw("same p");
  hPP_5020_2760_nnpdf_R4->SetMarkerStyle(20);
  hPP_5020_2760_nnpdf_R4->SetMarkerColor(6);
  hPP_5020_2760_nnpdf_R4->Draw("same p");

  TLegend *title2 = myLegend(0.1,0.6,0.3,0.9);
  title2->AddEntry(hPP_5020_2760_nnpdf_R2,"R=0.2, 5.02/2.76","pl");
  title2->AddEntry(hPP_5020_2760_nnpdf_R3,"R=0.3, 5.02/2.76","pl");
  title2->AddEntry(hPP_5020_2760_nnpdf_R4,"R=0.4, 5.02/2.76","pl");
  title2->SetTextSize(0.04);
  title2->Draw();

  c1->SaveAs("RpA_reference_check_NLO_diff_energy_diff_radius.pdf","RECREATE");

  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->SetGridy();
  hPP_5020_Pythia_05_R_3_5->SetMarkerStyle(20);
  hPP_5020_Pythia_05_R_3_5->SetMarkerColor(3);
  hPP_5020_Pythia_05_R_3_5->SetYTitle("Ratio of differential cross sections");
  hPP_5020_Pythia_05_R_3_5->SetXTitle("Jet p_{T} (GeV/c)");
  hPP_5020_Pythia_05_R_3_5->SetTitle(" ");
  hPP_5020_Pythia_05_R_3_5->SetAxisRange(22,500,"X");
  hPP_5020_Pythia_05_R_3_5->SetAxisRange(0,2,"Y");
  hPP_5020_Pythia_05_R_3_5->Draw("p");
  hPP_5020_Pythia_05_R_4_5->SetMarkerStyle(22);
  hPP_5020_Pythia_05_R_4_5->SetMarkerColor(4);
  hPP_5020_Pythia_05_R_4_5->Draw("same p");

  hPP_2760_Pythia_05_R_3_5->SetMarkerStyle(20);
  hPP_2760_Pythia_05_R_3_5->SetMarkerColor(6);
  hPP_2760_Pythia_05_R_3_5->Draw("same p");
  hPP_2760_Pythia_05_R_4_5->SetMarkerStyle(22);
  hPP_2760_Pythia_05_R_4_5->SetMarkerColor(7);
  hPP_2760_Pythia_05_R_4_5->Draw("same p");

  hPP_7000_Pythia_05_R_3_5->SetMarkerStyle(20);
  hPP_7000_Pythia_05_R_3_5->SetMarkerColor(12);
  hPP_7000_Pythia_05_R_3_5->Draw("same p");
  hPP_7000_Pythia_05_R_4_5->SetMarkerStyle(22);
  hPP_7000_Pythia_05_R_4_5->SetMarkerColor(13);
  hPP_7000_Pythia_05_R_4_5->Draw("same p");


  TLegend *title4 = myLegend(0.13,0.55,0.33,0.85);
  title4->AddEntry(hPP_5020_Pythia_05_R_3_5,"Z2 5.02, R=0.3/R=0.5","pl");
  title4->AddEntry(hPP_5020_Pythia_05_R_4_5,"Z2 5.02, R=0.4/R=0.5","pl");
  title4->AddEntry(hPP_2760_Pythia_05_R_3_5,"Z2 2.76, R=0.3/R=0.5","pl");
  title4->AddEntry(hPP_2760_Pythia_05_R_4_5,"Z2 2.76, R=0.4/R=0.5","pl");
  title4->AddEntry(hPP_7000_Pythia_05_R_3_5,"Z2 7, R=0.3/R=0.5","pl");
  title4->AddEntry(hPP_7000_Pythia_05_R_4_5,"Z2 7, R=0.4/R=0.5","pl");

  title4->SetTextSize(0.04);
  title4->Draw();

  putCMSPrel(0.1,0.92,0.04);
  drawText("pp Pythia Z2 , #sqrt{s}=2.7, 5.02 and 7(TeV)",0.35,0.92,16);
  drawText("Pythia |y|<0.5, Anti k_{T} PF Jets",0.47,0.83,16);

  c3->SaveAs("RpA_reference_pythia_diff_energy_radius_cross_section_ratio_05.pdf","RECREATE");





  TCanvas *c4 = new TCanvas("c4","",800,600);
  c4->SetGridy();
  hPP_5020_Pythia_05_R_3_5->SetMarkerStyle(20);
  hPP_5020_Pythia_05_R_3_5->SetMarkerColor(3);
  hPP_5020_Pythia_05_R_3_5->SetYTitle("Ratio of differential cross sections");
  hPP_5020_Pythia_05_R_3_5->SetXTitle("Jet p_{T} (GeV/c)");
  hPP_5020_Pythia_05_R_3_5->SetTitle(" ");
  hPP_5020_Pythia_05_R_3_5->SetAxisRange(22,500,"X");
  hPP_5020_Pythia_05_R_3_5->SetAxisRange(0,2,"Y");
  hPP_5020_Pythia_05_R_3_5->Draw("p");
  hPP_5020_Pythia_15_20_R_3_5->SetMarkerStyle(22);
  hPP_5020_Pythia_15_20_R_3_5->SetMarkerColor(4);
  hPP_5020_Pythia_15_20_R_3_5->Draw("same p");

  hPP_2760_Pythia_05_R_3_5->SetMarkerStyle(20);
  hPP_2760_Pythia_05_R_3_5->SetMarkerColor(6);
  hPP_2760_Pythia_05_R_3_5->Draw("same p");
  hPP_2760_Pythia_15_20_R_3_5->SetMarkerStyle(22);
  hPP_2760_Pythia_15_20_R_3_5->SetMarkerColor(7);
  hPP_2760_Pythia_15_20_R_3_5->Draw("same p");

  hPP_7000_Pythia_05_R_3_5->SetMarkerStyle(20);
  hPP_7000_Pythia_05_R_3_5->SetMarkerColor(12);
  hPP_7000_Pythia_05_R_3_5->Draw("same p");
  hPP_7000_Pythia_15_20_R_3_5->SetMarkerStyle(22);
  hPP_7000_Pythia_15_20_R_3_5->SetMarkerColor(13);
  hPP_7000_Pythia_15_20_R_3_5->Draw("same p");

  //hPP_2760_data_05_R_3_5->SetMarkerStyle(20);
  hPP_data_2760_R_3_5->SetMarkerStyle(20);
  hPP_data_2760_R_3_5->SetMarkerColor(8);
  hPP_data_2760_R_3_5->Draw("same p");

  //hPP_2760_data_15_20_R_3_5->SetMarkerStyle(20);
  //hPP_2760_data_15_20_R_3_5->SetMarkerColor(8);
  //hPP_2760_data_15_20_R_3_5->Draw("same p");

  TLegend *title5 = myLegend(0.13,0.55,0.33,0.85);
  title5->AddEntry(hPP_5020_Pythia_05_R_3_5,"Z2 5.02, |y|<0.5","pl");
  title5->AddEntry(hPP_5020_Pythia_15_20_R_3_5,"Z2 5.02, 1.5<|y|<2.0","pl");
  title5->AddEntry(hPP_2760_Pythia_05_R_3_5,"Z2 2.76, |y|<0.5","pl");
  title5->AddEntry(hPP_2760_Pythia_15_20_R_3_5,"Z2 2.76, 1.5<|y|<2.0","pl");
  title5->AddEntry(hPP_7000_Pythia_05_R_3_5,"Z2 7, |y|<0.5","pl");
  title5->AddEntry(hPP_7000_Pythia_15_20_R_3_5,"Z2 7, 1.5<|y|<2.0","pl");
  title5->AddEntry(hPP_data_2760_R_3_5,"Data 2.76 |eta|<2","pl");
  //title5->AddEntry(hPP_2760_data_15_20_R_3_5,"Data 2.76 1.5<|y|<2.0","pl");
  title5->SetTextSize(0.04);
  title5->Draw();

  putCMSPrel(0.1,0.92,0.04);
  drawText("pp Pythia Z2 , #sqrt{s}=2.7, 5.02 and 7(TeV)",0.35,0.92,16);
  drawText("R=0.3/R=0.5, Anti k_{T} PF Jets",0.47,0.83,16);

  c4->SaveAs("RpA_reference_pythia_data_diff_energy_radius_cross_section_ratio_R_3_5.pdf","RECREATE");




  TCanvas *c5 = new TCanvas("c5","",800,600);
  c5->SetGridy();
  hPP_5020_Pythia_05_R_4_5->SetMarkerStyle(20);
  hPP_5020_Pythia_05_R_4_5->SetMarkerColor(3);
  hPP_5020_Pythia_05_R_4_5->SetYTitle("Ratio of differential cross sections");
  hPP_5020_Pythia_05_R_4_5->SetXTitle("Jet p_{T} (GeV/c)");
  hPP_5020_Pythia_05_R_4_5->SetTitle(" ");
  hPP_5020_Pythia_05_R_4_5->SetAxisRange(22,500,"X");
  hPP_5020_Pythia_05_R_4_5->SetAxisRange(0,2,"Y");
  hPP_5020_Pythia_05_R_4_5->Draw("p");
  hPP_5020_Pythia_15_20_R_4_5->SetMarkerStyle(22);
  hPP_5020_Pythia_15_20_R_4_5->SetMarkerColor(4);
  hPP_5020_Pythia_15_20_R_4_5->Draw("same p");

  hPP_2760_Pythia_05_R_4_5->SetMarkerStyle(20);
  hPP_2760_Pythia_05_R_4_5->SetMarkerColor(6);
  hPP_2760_Pythia_05_R_4_5->Draw("same p");
  hPP_2760_Pythia_15_20_R_4_5->SetMarkerStyle(22);
  hPP_2760_Pythia_15_20_R_4_5->SetMarkerColor(7);
  hPP_2760_Pythia_15_20_R_4_5->Draw("same p");

  hPP_7000_Pythia_05_R_4_5->SetMarkerStyle(20);
  hPP_7000_Pythia_05_R_4_5->SetMarkerColor(12);
  hPP_7000_Pythia_05_R_4_5->Draw("same p");
  hPP_7000_Pythia_15_20_R_4_5->SetMarkerStyle(22);
  hPP_7000_Pythia_15_20_R_4_5->SetMarkerColor(13);
  hPP_7000_Pythia_15_20_R_4_5->Draw("same p");

  hPP_data_2760_R_4_5->SetMarkerStyle(20);
  hPP_data_2760_R_4_5->SetMarkerColor(8);
  hPP_data_2760_R_4_5->Draw("same p");

  //hPP_2760_data_15_20_R_4_->SetMarkerStyle(20);
  //hPP_2760_data_15_20_R_4_5->SetMarkerColor(8);
  //hPP_2760_data_15_20_R_4_5->Draw("same p");

  TLegend *title6 = myLegend(0.13,0.55,0.33,0.85);
  title6->AddEntry(hPP_5020_Pythia_05_R_4_5,"Z2 5.02, |y|<0.5","pl");
  title6->AddEntry(hPP_5020_Pythia_15_20_R_4_5,"Z2 5.02, 1.5<|y|<2.0","pl");
  title6->AddEntry(hPP_2760_Pythia_05_R_4_5,"Z2 2.76, |y|<0.5","pl");
  title6->AddEntry(hPP_2760_Pythia_15_20_R_4_5,"Z2 2.76, 1.5<|y|<2.0","pl");
  title6->AddEntry(hPP_7000_Pythia_05_R_4_5,"Z2 7, |y|<0.5","pl");
  title6->AddEntry(hPP_7000_Pythia_15_20_R_4_5,"Z2 7, 1.5<|y|<2.0","pl");
  title6->AddEntry(hPP_data_2760_R_4_5,"Data 2.76 |eta|<2","pl");
  //title6->AddEntry(hPP_2760_data_15_20_R_4_5,"Data 2.76 1.5<|y|<2.0","pl");
  title6->SetTextSize(0.04);
  title6->Draw();

  putCMSPrel(0.1,0.92,0.04);
  drawText("pp Pythia Z2 , #sqrt{s}=2.7, 5.02 and 7(TeV)",0.35,0.92,16);
  drawText("R=0.4/R=0.5, Anti k_{T} PF Jets",0.47,0.83,16);

  c5->SaveAs("RpA_reference_pythia_data_diff_energy_radius_cross_section_ratio_R_4_5.pdf","RECREATE");

  
  TCanvas *c6 = new TCanvas("c6","",800,600);
  
  c6->SetLogy();

  hPP_2760_5020_Pythia_05_R3->SetMarkerStyle(20);
  hPP_2760_5020_Pythia_05_R3->SetMarkerColor(3);
  hPP_2760_5020_Pythia_05_R3->SetTitle(" ");
  hPP_2760_5020_Pythia_05_R3->SetXTitle("Jet p_{T} (GeV/c)");
  hPP_2760_5020_Pythia_05_R3->SetYTitle("Ratio of Diff cross section");
  hPP_2760_5020_Pythia_05_R3->Draw("p");
  hPP_2760_5020_Pythia_05_R4->SetMarkerStyle(27);
  hPP_2760_5020_Pythia_05_R4->SetMarkerColor(4);
  hPP_2760_5020_Pythia_05_R4->Draw("same p");
  hPP_2760_5020_Pythia_05_R5->SetMarkerStyle(33);
  hPP_2760_5020_Pythia_05_R5->SetMarkerColor(6);
  hPP_2760_5020_Pythia_05_R5->Draw("same p");

  hPP_2760_7000_Pythia_05_R3->SetMarkerStyle(20);
  hPP_2760_7000_Pythia_05_R3->SetMarkerColor(7);
  hPP_2760_7000_Pythia_05_R3->Draw("same p");
  hPP_2760_7000_Pythia_05_R4->SetMarkerStyle(27);
  hPP_2760_7000_Pythia_05_R4->SetMarkerColor(8);
  hPP_2760_7000_Pythia_05_R4->Draw("same p");
  hPP_2760_7000_Pythia_05_R5->SetMarkerStyle(33);
  hPP_2760_7000_Pythia_05_R5->SetMarkerColor(13);
  hPP_2760_7000_Pythia_05_R5->Draw("same p");

  TLegend *title8 = myLegend(0.5,0.6,0.7,0.9);
  title8->AddEntry(hPP_2760_5020_Pythia_05_R3,"R=0.3, 2.76/5.02","pl");
  title8->AddEntry(hPP_2760_5020_Pythia_05_R4,"R=0.4, 2.76/5.02","pl");
  title8->AddEntry(hPP_2760_5020_Pythia_05_R5,"R=0.5, 2.76/5.02","pl");
  title8->AddEntry(hPP_2760_7000_Pythia_05_R3,"R=0.3, 2.76/7","pl");
  title8->AddEntry(hPP_2760_7000_Pythia_05_R4,"R=0.4, 2.76/7","pl");
  title8->AddEntry(hPP_2760_7000_Pythia_05_R5,"R=0.5, 2.76/7","pl");
  title8->SetTextSize(0.04);
  title8->Draw();

  putCMSPrel(0.1,0.92,0.04);
  drawText("pp Pythia Z2, #sqrt{s}=2.7, 5.02 and 7(TeV)",0.35,0.92,16);
  drawText("|y|<0.5, Anti k_{T} PF Jets",0.5,0.55,16);

  c6->SaveAs("RpA_reference_check_Pythia_05_diff_energy_diff_radius.pdf","RECREATE");

  TCanvas* c7 = new TCanvas("c7","",800,600);
  c7->SetGridy();

  
  hPP_2760_ct10n_nnpdf_2_2_R3->SetTitle(" ");
  hPP_2760_ct10n_nnpdf_2_2_R3->SetYTitle("Ratio of Diff cross section");
  hPP_2760_ct10n_nnpdf_2_2_R3->SetXTitle("Jet p_{T} (GeV/c)");
  hPP_2760_ct10n_nnpdf_2_2_R3->SetAxisRange(0.9,1.2,"Y");
  hPP_2760_ct10n_nnpdf_2_2_R3->SetAxisRange(30,500,"X");
  hPP_2760_ct10n_nnpdf_2_2_R3->SetMarkerStyle(20);
  hPP_2760_ct10n_nnpdf_2_2_R3->SetMarkerColor(3);
  hPP_2760_ct10n_nnpdf_2_2_R3->Draw("p");

  hPP_2760_ct10n_hera_2_2_R3->SetMarkerStyle(20);
  hPP_2760_ct10n_hera_2_2_R3->SetMarkerColor(4);
  hPP_2760_ct10n_hera_2_2_R3->Draw("same p");

  hPP_5020_ct10n_nnpdf_22_22_R3->SetMarkerStyle(23);
  hPP_5020_ct10n_nnpdf_22_22_R3->SetMarkerColor(6);
  hPP_5020_ct10n_nnpdf_22_22_R3->Draw("same p");

  hPP_5020_ct10n_hera_22_22_R3->SetMarkerStyle(23);
  hPP_5020_ct10n_hera_22_22_R3->SetMarkerColor(8);
  hPP_5020_ct10n_hera_22_22_R3->Draw("same p");

  TLegend *title7 = myLegend(0.1,0.6,0.3,0.85);
  title7->AddEntry(hPP_2760_ct10n_nnpdf_2_2_R3,"ct10n/nnpdf 2.76TeV","pl");
  title7->AddEntry(hPP_2760_ct10n_hera_2_2_R3,"ct10n/hera 2.76TeV","pl");
  title7->AddEntry(hPP_5020_ct10n_nnpdf_22_22_R3,"ct10n/nnpdf 5.02TeV","pl");
  title7->AddEntry(hPP_5020_ct10n_hera_22_22_R3,"ct10n/hera 5.02TeV","pl");
  title7->SetTextSize(0.04);
  title7->Draw();

  putCMSPrel(0.1,0.92,0.04);
  drawText("pp NLO,|#eta|<2 - 2.76TeV, |#eta|<2 - 5.02TeV, Anti k_{T} PF Jets",0.14,0.86,14);
  //drawText("|#eta|<2 - 2.76TeV, |#eta|<2 - 5.02TeV, Anti k_{T} PF Jets",0.2,0.55,16);

  c7->SaveAs("RpA_reference_NLO_systematics_calculation.pdf","RECREATE");

  
  foutput.Write();
  foutput.Close();

}
